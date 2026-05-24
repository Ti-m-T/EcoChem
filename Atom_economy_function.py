import re
from chempy.chemistry import balance_stoichiometry
from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from sympy import sympify
from experiment import Reaction
from experiment import Chemical
from experiment import Extras
from experiment import LiquidChemical
from experiment import Solvent
from experiment import Extractant

dictionary_input_smiles: dict[str, list[str]] = {
    "reactants": ["CCO", "O=O", "CC=O"],
    "products": ["CC(=O)O", "O", "C=O"],
}

###################################################################################################################

def react_builder(dictionary_input_smiles : dict[str,list[str]]) :
    react_list : str = []
    prod_list : str = []

    for react in dictionary_input_smiles["reactants"] :
        react = Chemical(smiles = react, mass = 1)
        react_list.append(react)

    for prod in dictionary_input_smiles["products"] : 
        prod = Chemical(smiles = prod, mass = 1)
        prod_list.append(prod)

    react = Reaction(react_list , prod_list)
    return react

###################################################################################################################

def calcul_eco_atom_react_M(react) :

    eco_atom_reactants_M: float = 0.0
    
    for reactant in react.reactants :
        eco_atom_reactants_M += (reactant.coeff) * reactant.mw

        # print(reactant.coeff ,reactant.mw)

    return ( eco_atom_reactants_M )

###################################################################################################################

def calcul_eco_atom_react_nb_atom(react) : 

    nb_atom_react : int = 0

    for reactant in react.reactants :
        nb_atom_react += (reactant.coeff) * reactant.nb_atom

        # print(reactant.coeff ,reactant.nb_atom)

    return (nb_atom_react)

###################################################################################################################

def calc_DOF_a_un(exp: float | str ) -> float:
   
    rech_var = re.compile(r"\bx\d+\b")

    def evaluer_exp(exp):
        exp_str = str(exp)

        exp_str = rech_var.sub("1", exp_str)

        # DOF a 1

        try:
            return float(sympify(exp_str))
        except Exception : 
            
            try:
                return float(exp_str)
            except Exception as e:
                raise ValueError(f"Cannot evaluate expression: {exp}") from e

        # Evalue l'expression et regarde s'il a un probleme

    equation_sans_DOF = evaluer_exp(exp)

    return equation_sans_DOF

###################################################################################################################
       
def calculate_eco_atm_M_et_nb_atom(dictionary_input_smiles: dict[str, list[str]],) -> float:

    react = react_builder(dictionary_input_smiles)

    react.stoich_of_reaction()

    eco_atom_reactants_M = calcul_eco_atom_react_M(react)

    eco_atom_mp = (react.products[0].coeff) * react.products[0].mw
    
    eco_atom_M = calc_DOF_a_un((eco_atom_mp * 100)/ eco_atom_reactants_M)

    nb_atom_reactants = calcul_eco_atom_react_nb_atom(react)

    nb_atom_mp = (react.products[0].coeff) * react.products[0].nb_atom

    nb_atom = calc_DOF_a_un((nb_atom_mp * 100)/ nb_atom_reactants)
    
    return (eco_atom_M , nb_atom)

print(calculate_eco_atm_M_et_nb_atom(dictionary_input_smiles))

###################################################################################################################

