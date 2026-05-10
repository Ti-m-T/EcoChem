from dataclasses import dataclass, field
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from chempy.chemistry import balance_stoichiometry
from thermo import chemical as density_finder

@dataclass
class Chemical: # Class grouping all chemicals 
    smiles: str
    mass: float = 0.0

    def __post_init__(self):
        self.mol = Chem.MolFromSmiles(self.smiles) # Gets molecular object from smiles
        if self.mol:
            self.mw = rdMolDescriptors.CalcExactMolWt(self.mol) # Calculates molecular weight from smiles
            self.logp = rdMolDescriptors.SlogP_VSA_(self.mol) # Finds the hydrophobicity of the molecule
            self.mol_f = rdMolDescriptors.CalcMolFormula(self.mol) # Finds the molecular formula using the smiles
            self.coeff : int = 1


            mol_H = Chem.AddHs(self.mol)

            self.nb_atom = mol_H.GetNumAtoms()
        else:
            self.mw = 0.0
            self.smiles = "Invalid"
    

@dataclass
class Extras(Chemical):
    pass


@dataclass    
class LiquidChemical(Chemical): # Subclass of chemical grouping all chemicals that are liquids, mainly used for solvents and extractants
    volume: float = 0.0
    density: float = 0.0

    def __post_init__(self):
        super().__post_init__()
        self.density = density_finder(f"SMILES={self.smiles}").rho

    @property
    def m_solvent(self):
        return self.density*self.volume


@dataclass    
class Solvent(LiquidChemical): # Subclass of LiquidChemical meant for solvents
    pass

@dataclass    
class Extractant(LiquidChemical): # Subclass of LiquidChemical meant for extractants
    pass

@dataclass
class Reaction:
    reactants: list[Chemical] = field(default_factory=list)
    products: list[Chemical] = field(default_factory=list)
    #solvents: list[Solvent] = field(default_factory=list)
    #extractants: list[Extractant] = field(default_factory=list)
    
    def stoich_of_reaction(self):
        reac = {reactant.mol_f for reactant in self.reactants}
        prod = {product.mol_f for product in self.products}
        reactants_coeff,products_coeff = balance_stoichiometry(reac,prod)
        
        for reactant in self.reactants:
            reactant.coeff = reactants_coeff.get(reactant.mol_f,1)
        for product in self.products:
            product.coeff = products_coeff.get(product.mol_f,1)
        
#r1 = Chemical(smiles="C",mass=1)
#r2 = Chemical(smiles="O=O",mass=1)
#p1 = Chemical(smiles="C(=O)=O",mass=1)
#p2 = Chemical(smiles="O",mass=1)
#react = Reaction([r1,r2],[p1,p2])
#react.stoich_of_reaction()
#for r in react.reactants:
#    print(r.coeff) 
#for p in react.products:
#    print(p.coeff)

mol = Chemical("CC=O")
print(mol.nb_atom)