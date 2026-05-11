from dataclasses import dataclass, field
import re
from rdkit import Chem
import streamlit as st
from rdkit.Chem import AllChem, rdMolDescriptors, Descriptors
from stmol import showmol
import py3Dmol
from pathlib import Path
import pandas as pd
import os
from streamlit_ketcher import st_ketcher
from rdkit.Chem import rdFingerprintGenerator
import numpy as np
import mols2grid
import streamlit.components.v1 as components
import plotly.figure_factory as ff
from typing import Tuple, List
import pubchempy as pcp
from chempy.chemistry import balance_stoichiometry
from thermo import chemical as density_finder

@dataclass
class Chemical: # Class grouping all chemicals 
    smiles: str
    mass: float=0.0

    def __post_init__(self):
        self.mol=Chem.MolFromSmiles(self.smiles) # Gets molecular object from smiles
        if self.mol:
            self.mw=rdMolDescriptors.CalcExactMolWt(self.mol) # Calculates molecular weight from smiles
            self.logp=rdMolDescriptors.SlogP_VSA_(self.mol) # Finds the hydrophobicity of the molecule
            self.mol_f=rdMolDescriptors.CalcMolFormula(self.mol) # Finds the molecular formula using the smiles
            self.coeff: int=1
            self.nb_atom = Chem.AddHs(self.mol).GetNumAtoms() # Finds the number of atoms in the molecule including hydrogens
        else:
            self.mw=0.0
            self.smiles="Invalid"
    

@dataclass
class Extras(Chemical):
    pass


@dataclass    
class LiquidChemical(Chemical): # Subclass of chemical grouping all chemicals that are liquids, mainly used for solvents and extractants
    volume: float =0.0
    density: float=0.0

    def __post_init__(self):
        super().__post_init__()
        self.density=density_finder(f"SMILES={self.smiles}").rho

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
        reac={reactant.mol_f for reactant in self.reactants}
        prod={product.mol_f for product in self.products}
        reactants_coeff,products_coeff=balance_stoichiometry(reac,prod)
        for reactant in self.reactants:
            reactant.coeff=reactants_coeff.get(reactant.mol_f,1)
        for product in self.products:
            product.coeff=products_coeff.get(product.mol_f,1)
    
    def calcul_eco_atom_reactants_M(self) :
        #if all(isinstance(reactant.mw,(float,int)) for reactant in reactants):
        return sum(reactant.coeff*reactant.mw for reactant in self.reactants)
      
    
    def calcul_eco_atom_products_M(self) :
        #if all(isinstance(product.mw,(float,int)) for product in products):
        return sum(product.coeff*product.mw for product in self.products)
    
    def calcul_eco_atom_react_nb_atom(self):
        #if all(isinstance(reactant.coeff,(float,int)) for reactant in reactants):
        return sum(reactant.coeff*reactant.nb_atom for reactant in self.reactants)
       

    def calcul_eco_atom_product_nb_atom(self):
       # if all(isinstance(product.coeff,(float,int)) for product in products):
        return sum(product.coeff*product.nb_atom for product in self.products)
        

        
'''r1=Chemical(smiles="C",mass=1)
r2=Chemical(smiles="O=O",mass=1)
p1=Chemical(smiles="C(=O)=O",mass=1)
p2=Chemical(smiles="O",mass=1)
react=Reaction([r1,r2],[p1,p2])'''
reactants=[Chemical(smiles="C",mass=1),Chemical(smiles="O=O",mass=1)]
products=[Chemical(smiles="O",mass=1),Chemical(smiles="C(=O)=O",mass=1)]
reacto=Reaction(reactants,products)
reacto.stoich_of_reaction()

for r in reacto.reactants:
    print(r.mol_f,r.coeff)
for p in reacto.products:
    print(p.mol_f,p.coeff)
print(reacto.calcul_eco_atom_react_nb_atom())