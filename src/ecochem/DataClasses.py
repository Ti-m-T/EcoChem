from dataclasses import dataclass, field
from math import exp
from re import search
import re
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, Descriptors
import pubchempy as pcp
from chempy.chemistry import balance_stoichiometry
from thermo import Chemical as density_finder
import requests
import pubchempy as pcp
from tomlkit import item
from sympy import sympify


@dataclass
class Chemical:  # Class grouping all chemicals

    smiles: str

    CID: int | None = None

    GHS: dict = field(default_factory=dict)

    pictograms: list = field(default_factory=list)

    def __post_init__(self):

        self.mol = Chem.MolFromSmiles(self.smiles)  # Gets molecular object from smiles

        if self.mol:

            self.mw = rdMolDescriptors.CalcExactMolWt(
                self.mol
            )  # Calculates molecular weight from smiles

            self.logp = rdMolDescriptors.CalcCrippenDescriptors(self.mol)[
                0
            ]  # Finds the hydrophobicity of the molecule

            self.mol_f = rdMolDescriptors.CalcMolFormula(
                self.mol
            )  # Finds the molecular formula using the smiles

            self.coeff: int = 1

            self.nb_atom = Chem.AddHs(
                self.mol
            ).GetNumAtoms()  # Finds the number of atoms in the molecule including hydrogens

            self.moles: float = 0.0

            self._mass: float = 0.0

        else:
            self.mw = 0.0

            self.smiles = "Invalid"

            self.CID = "None"

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, value: float):
        self._mass = value

        if self.mw > 0:
            self.moles = value / self.mw

    def get_CID(
        self,
    ):  # Searches for the CID of a molecule using the PubChem data base.

        try:
            CID = pcp.get_cids(
                self.smiles, "smiles"
            )  # tries to use the smiles to get the CID of the molecule from PubChem database.

            if CID:  # if the CID is found it's stored in the class and returned.
                self.CID = CID[0]
                return self.CID

            else:  # if the CID is not found the CID is set to None.
                self.CID = None
                return None

        except (
            Exception
        ) as e:  # If there is an error during the search for the CID we print the error and set the CID to None.

            print(f"CID error : {e}")
            self.CID = None

            return None

    def get_GHS(self):

        if (
            self.CID is None
        ):  # If the CID is not set, search for it using the get_CID function.
            self.get_CID()

        if (
            self.CID is None
        ):  # If the CID is not found or an error occurs then we return not GHS.
            return None

        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/"
            f"rest/pug_view/data/compound/{self.CID}/JSON"
        )  # We use the CID information to acces the PubChem database of the said molecule. We use the requests library to get the data in json format.

        data = requests.get(url).json()

        self.GHS = {}  # We create a dictionary to store the GHS information.

        def extract_ghs(section):

            if isinstance(section, dict):  # We asks if the section is a dictionary.

                if (
                    section.get("TOCHeading") == "GHS Classification"
                ):  # We search for the section of the database that contains the GHS information. If we find it we extract the GHS codes and labels and store them in the class lists.

                    def collect_GHS(
                        section,
                    ):  # Recursive function to collect GHS codes and labels from the section of the database that contains the GHS information.

                        if isinstance(section, dict):  # Case our data is a dictionary.

                            for (
                                values
                            ) in (
                                section.values()
                            ):  # Since section is a dictionary we loop through its values to find the GHS codes and labels.

                                if isinstance(
                                    values, str
                                ):  # We check if the values are strings.

                                    if (
                                        values.startswith("H") and ":" in values
                                    ):  # We check for "Hazard statements" , which start with H and are followed by a number, and are separated from their label by a ":"

                                        raw_code = values.split(":", 1)[0]
                                        code = raw_code.split()[0].strip()
                                        label = values.split(":", 1)[
                                            1
                                        ].strip()  # We extract the code and label from the value string.

                                        if code not in self.GHS:

                                            self.GHS[code] = (
                                                label  # We store the code and label in the GHS dictionary of the class if the code is not already in it.
                                            )

                                else:  # If the values are not strings we call the collect function on them to keep searching for GHS codes and labels in the database.
                                    collect_GHS(values)

                        elif isinstance(
                            section, list
                        ):  # We ask if the section is a list.

                            for (
                                item
                            ) in (
                                section
                            ):  # If it is a list we loop through its items and call the collect function on them to keep searching for GHS codes and labels in the database.
                                collect_GHS(item)

                    collect_GHS(
                        section
                    )  # We define the collect_GHS function and call it on the section of the database that contains the GHS information to extract all GHS codes and labels.

                    return  # We return after finding the GHS information to avoid unnecessary searching in the rest of the database.

                for (
                    values
                ) in (
                    section.values()
                ):  # If the current dictionary is not the GHS section we explore the PubChem database deeper.
                    extract_ghs(values)

            elif isinstance(
                section, list
            ):  # If the section wasn't a dictionary but a list we loop through its items and call the extract_ghs function on them to keep searching for the GHS information in the database.

                for item in section:

                    extract_ghs(item)

        extract_ghs(
            data
        )  # We execute the extract_ghs function on the data we got from the PubChem database of the molecule wanted.

        return {
            "pictograms": list(self.pictograms),
            "codes": list(self.GHS.keys()),
            "labels": list(self.GHS.values()),
        }

    def get_pictograms(self) -> list:

        if self.CID is None:
            self.get_CID()

        if self.CID is None:
            return []

        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/"
            f"rest/pug_view/data/compound/{self.CID}/JSON"
        )

        data = requests.get(url).json()

        self.pictograms = set()

        mapping = {
            "GHS01": "Exploding bomb",
            "GHS02": "Flame",
            "GHS03": "Oxidizer (flame over circle)",
            "GHS04": "Gas cylinder",
            "GHS05": "Corrosion",
            "GHS06": "Skull and crossbones",
            "GHS07": "Exclamation mark",
            "GHS08": "Health hazard",
            "GHS09": "Environment",
        }  # We introduce a mapping between the possible GHS pictogram code and their meaning.

        def search_pictograms(section, in_primary=False):

            if isinstance(section, dict):

                if (
                    section.get("TOCHeading") == "Primary Hazards"
                ):  # Only look for primary hazards
                    in_primary = True

                for key, values in section.items():

                    if in_primary and isinstance(
                        values, str
                    ):  # Search strings inside Primary Hazards
                        for code, name in mapping.items():

                            if code in values:

                                self.pictograms.add(name)

                    else:
                        search_pictograms(values, in_primary)

            elif isinstance(section, list):

                for item in section:

                    search_pictograms(item, in_primary)

        search_pictograms(data)

        return list(
            self.pictograms
        )  # We return the list of GHS pictograms found for the molecule.


@dataclass
class ChemswithMass(Chemical):
    initial_mass: float = 0.0

    def __post_init__(self):
        super().__post_init__()

        if self.initial_mass > 0:
            self.mass = self.initial_mass


@dataclass
class LiquidChemical(
    Chemical
):  # Subclass of chemical grouping all chemicals that are liquids, mainly used for solvents and extractants
    volume: float = 0.0
    evaluated_density: float = field(init=False, default=0.0)
    user_density: float = 0.0

    def __post_init__(self):
        super().__post_init__()
        self.density = 0.0

        try:
            chem_data = density_finder(self.smiles, T=298.15)

            if chem_data.rhol is not None:
                self.density = chem_data.rhol / 1000

        except Exception as e:
            print(e)
            self.density = 0.0

    @property
    def m_liquid(self):
        if self.user_density > 0:
            return self.user_density * self.volume
        # return self.density * self.volume


@dataclass
class Solvent(LiquidChemical):  # Subclass of LiquidChemical meant for solvents

    def __post_init__(self):
        super().__post_init__()


@dataclass
class Extractant(LiquidChemical):  # Subclass of LiquidChemical meant for extractants

    def __post_init__(self):
        super().__post_init__()


@dataclass
class Reaction:
    reactants: list[Chemical] = field(
        default_factory=list
    )  # Class that contains all reactants used in a reaction
    byproducts: list[Chemical] = field(
        default_factory=list
    )  # Class that contains all products used in a reaction
    wanted_product: ChemswithMass = field(
        default_factory=lambda: ChemswithMass(smiles="")
    )
    Catalysts: list[ChemswithMass] = field(default_factory=list)
    solvents: list[Solvent] = field(default_factory=list)
    extractants: list[Extractant] = field(default_factory=list)
    Chosen_Yield: float = 1

    def stoich_of_reaction(self):

        reac = {reactant.mol_f for reactant in self.reactants}
        prod = {byproduct.mol_f for byproduct in self.byproducts}
        prod.add(self.wanted_product.mol_f)
        reactants_coeff, products_coeff = balance_stoichiometry(reac, prod)

        def calc_DOF_a_un(Coef_with_DOF: float | str) -> float:

            rech_var = re.compile(
                r"\bx\d+\b"
            )  # We look for the DOF of form x followed by a number in the coefficient expression using regular expressions.

            def evaluer_exp(Coef_with_DOF: float | str) -> float:

                Coef_as_str = str(Coef_with_DOF)

                Coef_as_str = rech_var.sub("1", Coef_as_str)  # Sets DOF to 1.

                try:
                    return float(
                        sympify(Coef_as_str)
                    )  # Tries to evaluate the expression with the DOF set to 1
                except Exception:

                    try:
                        return float(Coef_as_str)
                    except Exception as e:
                        raise ValueError(
                            f"Cannot evaluate expression: {Coef_with_DOF}"
                        ) from e  # If there is an error during the evaluation it raises an error.

            coef_without_DOF = evaluer_exp(Coef_with_DOF)

            return coef_without_DOF

        for reactant in self.reactants:
            reactant.coeff = calc_DOF_a_un(reactants_coeff.get(reactant.mol_f, 1))

        for product in self.byproducts:
            product.coeff = calc_DOF_a_un(products_coeff.get(product.mol_f, 1))

        self.wanted_product.coeff = calc_DOF_a_un(
            products_coeff.get(self.wanted_product.mol_f, 1)
        )

        return reactants_coeff, products_coeff

    def calcul_eco_atom(self):

        self.stoich_of_reaction()

        return (
            (self.wanted_product.coeff * self.wanted_product.nb_atom)
            * 100
            / sum(reactant.coeff * reactant.nb_atom for reactant in self.reactants)
        )

    def total_mass_catalysts(self):
        # self.stoich_of_reaction()
        total_mass_catalyst: float = 0.0

        for catalyst in self.Catalysts:

            total_mass_catalyst += catalyst.mass

        return total_mass_catalyst

    def total_mass_solvents(self):
        # self.stoich_of_reaction()
        total_mass_solvent: float = 0.0

        for solvent in self.solvents:

            total_mass_solvent += solvent.m_liquid

        return total_mass_solvent

    def total_mass_extractants(self):
        # self.stoich_of_reaction()
        total_mass_extractant: float = 0.0

        for extractant in self.extractants:

            total_mass_extractant += extractant.m_liquid

        return total_mass_extractant

    def total_mass_reactants(self):
        self.stoich_of_reaction()

        total_mass_reactant: float = 0.0

        for reactant in self.reactants:

            reactant.moles = self.wanted_product.moles / (
                (self.wanted_product.coeff / reactant.coeff) * self.Chosen_Yield
            )
            reactant.mass = reactant.moles * reactant.mw
            total_mass_reactant += reactant.mass

        return total_mass_reactant

    def total_mass_byproducts(self):
        self.stoich_of_reaction()

        total_mass_byproduct: float = 0.0

        for byproduct in self.byproducts:

            byproduct.moles = self.wanted_product.moles / (
                (self.wanted_product.coeff / byproduct.coeff)
            )
            byproduct.mass = byproduct.moles * byproduct.mw
            total_mass_byproduct += byproduct.mass

        return total_mass_byproduct

    def mass_reactants_left(self):
        self.total_mass_reactants()
        self.stoich_of_reaction()

        tot_mass_reactant_left: float = 0.0

        for reactant in self.reactants:

            mol_reactant_left = (
                reactant.moles
                - (self.wanted_product.moles * reactant.coeff)
                / self.wanted_product.coeff
            )
            mass_reactant_left = mol_reactant_left * reactant.mw
            tot_mass_reactant_left += mass_reactant_left

        return tot_mass_reactant_left

    def PMI(self):
        self.stoich_of_reaction()

        total_input = (
            self.total_mass_reactants()
            + self.total_mass_solvents()
            + self.total_mass_extractants()
            + self.total_mass_catalysts()
        )
        return total_input / self.wanted_product.mass

    def e_factor(self):
        self.stoich_of_reaction()

        waste = (
            self.total_mass_catalysts()
            + self.mass_reactants_left()
            + self.total_mass_extractants()
            + self.total_mass_solvents()
            + self.total_mass_byproducts()
        )
        return waste / self.wanted_product.mass
