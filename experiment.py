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
    mass: float = field(init=False, default=0.0)

    def __post_init__(self):
        self.mol = Chem.MolFromSmiles(self.smiles) # Gets molecular object from smiles
        if self.mol:
            self.mw=rdMolDescriptors.CalcExactMolWt(self.mol) # Calculates molecular weight from smiles
            self.logp=rdMolDescriptors.CalcCrippenDescriptors(self.mol)[0] # Finds the hydrophobicity of the molecule
            self.mol_f=rdMolDescriptors.CalcMolFormula(self.mol) # Finds the molecular formula using the smiles
            self.coeff: int=1
            self.nb_atom = Chem.AddHs(self.mol).GetNumAtoms() # Finds the number of atoms in the molecule including hydrogens
            self.moles:float = self.mass/self.mw
        else:
            self.mw = 0.0
            self.smiles = "Invalid"
    

@dataclass
class ChemswithMass(Chemical):
    mass: float = 0.0

    def __post_init__(self):
        super().__post_init__()


@dataclass    
class LiquidChemical(Chemical): # Subclass of chemical grouping all chemicals that are liquids, mainly used for solvents and extractants
    volume: float = 0.0
    density: float = field(init=False, default=0.0)

    def __post_init__(self):
        super().__post_init__()
        self.density = density_finder(self.smiles).rhols

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
    byproducts: list[Chemical] = field(default_factory=list)
    #wanted_products: list[ChemswithMass] = field(default_factory=list)
    wanted_product: ChemswithMass = 0
    reagents: list[ChemswithMass] = field(default_factory=list)
    solvents: list[Solvent] = field(default_factory=list)
    extractants: list[Extractant] = field(default_factory=list)
    Chosen_Yield: float = 0.0
    
    def stoich_of_reaction(self):
        reac={reactant.mol_f for reactant in self.reactants}
        prod={other_product.mol_f for other_product in self.other_products}
        prod.add(self.wanted_product.mol_f)
        reactants_coeff,products_coeff=balance_stoichiometry(reac,prod)
        for reactant in self.reactants:
            reactant.coeff=reactants_coeff.get(reactant.mol_f,1)
        for product in self.products:
            product.coeff=products_coeff.get(product.mol_f,1)
        self.wanted_product.coeff=products_coeff.get(self.wanted_product.mol_f,1)
        return reactants_coeff,products_coeff
    
    def calcul_eco_atom(self):
        self.stoich_of_reaction()
        return (self.wanted_product.coeff*self.wanted_product.mw)*100/sum(reactant.coeff*reactant.mw for reactant in self.reactants)

    def total_mass_reagents(self):
         self.stoich_of_reaction()
         for reagent in self.reagents:
            reagent.moles = self.wanted_product.moles/((self.wanted_product.coeff/reagent.coeff)*self.Chosen_Yield)
            reagent.mass = reagent.moles*reagent.mw
            total_mass_reagent+=reagent.mass
         return total_mass_reagent
    
    def total_mass_solvents(self):
        self.stoich_of_reaction()
        total_mass_solvent:float = 0.0
        for solvent in self.solvents:
            total_mass_solvent+=solvent.mass
        return total_mass_solvent
    
    def total_mass_extractants(self):
        self.stoich_of_reaction()
        total_mass_extractant:float = 0.0
        for extractant in self.extractants:
            total_mass_extractant+=extractant.mass
        return total_mass_extractant

    def total_mass_reactants(self):
        self.stoich_of_reaction()
        total_mass_reactant:float = 0.0
        for reactant in self.reactants:
            reactant.moles = self.wanted_product.moles/((self.wanted_product.coeff/reactant.coeff)*self.Chosen_Yield)
            reactant.mass = reactant.moles*reactant.mw
            total_mass_reactant+=reactant.mass
        return total_mass_reactant
    
    def total_mass_byproducts(self):
        self.stoich_of_reaction()
        total_mass_byproduct:float = 0.0
        for byproduct in self.byproducts:
            byproduct.moles = self.wanted_product.moles/((self.wanted_product.coeff/byproduct.coeff))
            byproduct.mass = byproduct.moles*byproduct.mw
            total_mass_byproduct+=byproduct.mass
        return total_mass_byproduct
    
    def mass_reactants_left(self):
        self.total_mass_extractants()
        self.stoich_of_reaction()
        tot_mass_reactant_left:float = 0.0
        for reactant in self.reactants:
            mol_reactant_left=reactant.moles-(self.wanted_product.moles*reactant.coeff)/self.wanted_product.coeff
            mass_reactant_left=mol_reactant_left*reactant.mw
            tot_mass_reactant_left+=mass_reactant_left
        return tot_mass_reactant_left
    
    def mass_reagents_left(self):
        self.total_mass_extractants()
        self.stoich_of_reaction()
        tot_mass_reagent_left:float = 0.0
        for reagent in self.reagents:
            mol_reagent_left=reagent.moles-(self.wanted_product.moles*reagent.coeff)/self.wanted_product.coeff
            mass_reagent_left=mol_reagent_left*reagent.mw
            tot_mass_reagent_left+=mass_reagent_left
        return tot_mass_reagent_left

    def PMI(self):
        self.stoich_of_reaction()
        total_input = self.total_mass_reactants() + self.total_mass_solvents() + self.total_mass_extractants() + self.total_mass_reagents()
        return total_input/self.wanted_product.mass
     
    def e_factor(self):
        self.stoich_of_reaction()
        waste=self.mass_reagents_left()+self.mass_reactants_left()+self.total_mass_extractants()+self.total_mass_solvents()
        return waste/self.wanted_product.mass
        

        

        
'''r1=Chemical(smiles="C",mass=1)
r2=Chemical(smiles="O=O",mass=1)
p1=Chemical(smiles="C(=O)=O",mass=1)
p2=Chemical(smiles="O",mass=1)
react=Reaction([r1,r2],[p1,p2])'''
reactants=[Chemical(smiles="C"),Chemical(smiles="O=O")]
byproducts=[Chemical(smiles="O")]
main_product=Chemical(smiles="C(=O)=O")
Solvent
reacto=Reaction()


