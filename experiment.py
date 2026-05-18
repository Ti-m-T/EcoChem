from dataclasses import dataclass, field
#import re
from rdkit import Chem
#import streamlit as st
from rdkit.Chem import AllChem, rdMolDescriptors, Descriptors
#from stmol import showmol
#import py3Dmol
#from pathlib import Path
#import pandas as pd
#import os
#from streamlit_ketcher import st_ketcher
#from rdkit.Chem import rdFingerprintGenerator
#import numpy as np
#import mols2grid  
#import streamlit.components.v1 as components
#import plotly.figure_factory as ff
#from typing import Tuple, List
import pubchempy as pcp
from chempy.chemistry import balance_stoichiometry
from thermo import chemical as density_finder
import requests
import pubchempy as pcp

@dataclass
class Chemical: # Class grouping all chemicals 

    smiles : str

    CID : int | None = None 

    GHS_codes :list = field (default_factory = list)

    GHS_labels : list = field (default_factory = list)
    

    def __post_init__(self):

        self.mol = Chem.MolFromSmiles(self.smiles) # Gets molecular object from smiles

        if self.mol:

            self.mw = rdMolDescriptors.CalcExactMolWt(self.mol) # Calculates molecular weight from smiles

            self.logp = rdMolDescriptors.CalcCrippenDescriptors(self.mol)[0] # Finds the hydrophobicity of the molecule

            self.mol_f = rdMolDescriptors.CalcMolFormula(self.mol) # Finds the molecular formula using the smiles

            self.coeff: int=1

            self.nb_atom = Chem.AddHs(self.mol).GetNumAtoms() # Finds the number of atoms in the molecule including hydrogens

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
    def mass(self, value:float):
        self._mass=value
        if self.mw>0:
         self.moles=value/self.mw

    def get_CID(self) :

        try :
            
            CID = pcp.get_cids(self.smiles, "smiles")
            
            if CID :

                self.CID = CID[0]

                return self.CID

            else : 
                
                self.CID = None 

                return None
            
        except Exception as e :

            print(f"CID error : {e}")
            self.CID = None 

            return None 
        


    def get_GHS(self) :

        

        if self.CID is None :
        
            self.get_CID()

        if self.CID is None :

            return None
        
        print(self.CID)


        url = (f"https://pubchem.ncbi.nlm.nih.gov/" f"rest/pug_view/data/compound/{self.CID}/JSON")

        data = requests.get(url).json()

        self.GHS_codes = []
        self.GHS_labels = []

        def extract_ghs(section):

            if isinstance(section, dict):

                # ONLY process correct section
                if section.get("TOCHeading") == "GHS Classification":

                    def collect(obj):

                        if isinstance(obj, dict):

                            for v in obj.values():

                                if isinstance(v, str):

                                    if v.startswith("H") and ":" in v:

                                        code, label = v.split(":", 1)

                                        self.GHS_codes.append(code.strip())
                                        self.GHS_labels.append(label.strip())

                                else:
                                    collect(v)

                        elif isinstance(obj, list):

                            for item in obj:
                                collect(item)

                    collect(section)

                    return  # stop after finding GHS section

                # continue searching other branches
                for v in section.values():
                    extract_ghs(v)

            elif isinstance(section, list):

                for item in section:
                    extract_ghs(item)

        extract_ghs(data)

        clean_codes = []
        clean_labels = []

        for c in self.GHS_codes:

        # remove "(100%)" etc
            code = c.split()[0]
            clean_codes.append(code)

            self.GHS_codes = list(dict.fromkeys(clean_codes))
            self.GHS_labels = list(dict.fromkeys(self.GHS_labels))

        return {
        "codes": self.GHS_codes,
        "labels": self.GHS_labels
    }





@dataclass
class ChemswithMass(Chemical):
    initial_mass: float = 0.0

    def __post_init__(self):
        super().__post_init__()
        if self.initial_mass>0:
            self.mass=self.initial_mass


@dataclass    
class LiquidChemical(Chemical): # Subclass of chemical grouping all chemicals that are liquids, mainly used for solvents and extractants
    volume: float = 0.0
    density: float = field(init=False, default=0.0)

    def __post_init__(self):
        super().__post_init__()
        self.density = density_finder(self.smiles).rhols /1000

    @property
    def m_liquid(self):
        return self.density*self.volume


@dataclass    
class Solvent(LiquidChemical): # Subclass of LiquidChemical meant for solvents
    pass

@dataclass    
class Extractant(LiquidChemical): # Subclass of LiquidChemical meant for extractants
    pass

@dataclass
class Reaction:
    reactants: list[Chemical] = field(default_factory=list) # Class that contains all reactants used in a reaction
    byproducts: list[Chemical] = field(default_factory=list) # Class that contains all products used in a reaction
    #wanted_products: list[ChemswithMass] = field(default_factory=list)
    wanted_product: ChemswithMass = field(default_factory=lambda: ChemswithMass(smiles=""))
    Catalysts: list[ChemswithMass] = field(default_factory=list)
    solvents: list[Solvent] = field(default_factory=list)
    extractants: list[Extractant] = field(default_factory=list)
    Chosen_Yield: float = 1
    
    def stoich_of_reaction(self):
        reac = {reactant.mol_f for reactant in self.reactants}
        prod = {byproduct.mol_f for byproduct in self.byproducts}
        prod.add(self.wanted_product.mol_f)
        reactants_coeff,products_coeff = balance_stoichiometry(reac,prod)
        for reactant in self.reactants:
            reactant.coeff = reactants_coeff.get(reactant.mol_f,1)
        for product in self.byproducts:
            product.coeff = products_coeff.get(product.mol_f,1)
        self.wanted_product.coeff=products_coeff.get(self.wanted_product.mol_f,1)
        return reactants_coeff,products_coeff
    
    def calcul_eco_atom(self):
        self.stoich_of_reaction()
        return (self.wanted_product.coeff*self.wanted_product.mw)*100/sum(reactant.coeff*reactant.mw for reactant in self.reactants)

    def total_mass_catalysts(self):
         self.stoich_of_reaction()
         total_mass_catalyst=0
         for catalyst in self.Catalysts:
            total_mass_catalyst+=catalyst.mass
         return total_mass_catalyst
    
    def total_mass_solvents(self):
        self.stoich_of_reaction()
        total_mass_solvent:float = 0.0
        for solvent in self.solvents:
            total_mass_solvent +=solvent.m_liquid
        return total_mass_solvent
    
    def total_mass_extractants(self):
        self.stoich_of_reaction()
        total_mass_extractant:float = 0.0
        for extractant in self.extractants:
            total_mass_extractant += extractant.m_liquid
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
        self.total_mass_reactants()
        self.stoich_of_reaction()
        tot_mass_reactant_left:float = 0.0
        for reactant in self.reactants:
            mol_reactant_left=reactant.moles-(self.wanted_product.moles*reactant.coeff)/self.wanted_product.coeff
            mass_reactant_left=mol_reactant_left*reactant.mw
            tot_mass_reactant_left+=mass_reactant_left
        return tot_mass_reactant_left

    def PMI(self):
        self.stoich_of_reaction()
        total_input = self.total_mass_reactants() + self.total_mass_solvents() + self.total_mass_extractants() + self.total_mass_catalysts()
        return total_input/self.wanted_product.mass
     
    def e_factor(self):
        self.stoich_of_reaction()
        waste=self.total_mass_catalysts()+self.mass_reactants_left()+self.total_mass_extractants()+self.total_mass_solvents()+self.total_mass_byproducts
        return waste/self.wanted_product.mass
        

        

     
r1=Chemical(smiles="C")
r2=Chemical(smiles="O=O")
p1=Chemical(smiles="C(=O)=O")
p2=Chemical(smiles="O")
react=Reaction([r1,r2],[p1],p2)
reactants = [Chemical(smiles="C"),Chemical(smiles="O=O")]
byproducts=[Chemical(smiles="O")]
main_product=Chemical(smiles="C(=O)=O")
Solvent
reacto=Reaction()


for chem in react.reactants:
    chem.get_GHS()
    print(chem.GHS_codes)
    print(chem.GHS_labels)

react.wanted_product.get_GHS()
print(react.wanted_product.GHS_codes)
print(react.wanted_product.GHS_labels)

for byprod in react.byproducts :
    byprod.get_GHS()
    print(byprod.GHS_codes)
    print(byprod.GHS_labels)