# tests/test_ecochem.py
from ecochem.DataClasses import (
    Reaction,
    Chemical,
    ChemswithMass,
    Solvent,
    Extractant,
    LiquidChemical,
)

import pytest
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
from unittest.mock import patch
from sympy import sympify
import re
import thermo


def test_chemical_valid_smiles():
    """Test that a valid SMILES creates a proper Chemical object"""
    chem = Chemical(smiles="C")  
    assert chem.mol is not None
    assert chem.mw > 0
    assert chem.mol_f == "CH4"

def test_chemical_invalid_smiles():
    """Test that an invalid SMILES returns 0.0 for molecular weight and invalid as a SMILES"""
    chem = Chemical(smiles="INVALID")
    assert chem.mw == 0.0
    assert chem.smiles == "Invalid"

def test_chemical_molecular_weight():
    """Test molecular weight calculation for water"""
    water = Chemical(smiles="O")  
    assert round(water.mw, 2) == 18.01

def test_chemical_mass_setter():
    """Test that by giving a mass the corresponding moles are calculated"""
    chem = Chemical(smiles="O") 
    chem.mass = 18.0
    assert round(chem.moles, 1) == 1.0

def test_chemical_logp():
    """Test the logP calculation"""
    chem = Chemical(smiles="C")
    assert isinstance(chem.logp, float)

def test_chemical_nb_atom():
    """Test atom count including hydrogens"""
    methane = Chemical(smiles="C")  
    assert methane.nb_atom == 5

def test_chemswithmass_sets_mass():
    """Test that the initial_mass is given correctly"""
    chem = ChemswithMass(smiles="O", initial_mass=18.0)
    assert chem.mass == 18.0

def test_chemswithmass_calculates_moles():
    """Test that moles are calculated from the initial mass given with the ChemswithMass class"""
    chem = ChemswithMass(smiles="O", initial_mass=18.0)
    assert round(chem.moles, 1) == 1.0

def test_chemswithmass_zero_mass():
    """Test that moles are not calculated when the initial mass is 0"""
    chem = ChemswithMass(smiles="O", initial_mass=0.0)
    assert chem.moles == 0.0

"Example of a Combustion reaction to test the Reaction class"
@pytest.fixture
def combustion_reaction():
    """
    Methane combustion: CH4 + 2O2 → CO2 + 2H2O
    """
    reactants = [
        Chemical(smiles="C"),    
        Chemical(smiles="O=O"),  
    ]
    byproducts = [
        Chemical(smiles="O"),    
    ]
    wanted_product = ChemswithMass(smiles="C(=O)=O", initial_mass=44.0) 
    return Reaction(
        reactants=reactants,
        byproducts=byproducts,
        wanted_product=wanted_product
    )

def test_atom_economy(combustion_reaction):
    """Test atom economy calculation"""
    ae = combustion_reaction.calcul_eco_atom()
    assert 0 < ae <= 100  

def test_pmi(combustion_reaction):
    """Test PMI is greater than 1"""
    pmi = combustion_reaction.PMI()
    assert pmi >= 1.0

def test_e_factor(combustion_reaction):
    """Test E-factor is non-negative"""
    ef = combustion_reaction.e_factor()
    assert ef >= 0.0

def test_total_mass_reactants(combustion_reaction):
    """Test that total reactant mass is positive"""
    mass = combustion_reaction.total_mass_reactants()
    assert mass > 0

def test_total_mass_byproducts(combustion_reaction):
    """Test that total byproduct mass is positive"""
    mass = combustion_reaction.total_mass_byproducts()
    assert mass > 0

def test_stoich_of_reaction(combustion_reaction):
    """Test that stoichiometry returns dicts"""
    reactants_coeff, products_coeff = combustion_reaction.stoich_of_reaction()
    assert isinstance(reactants_coeff, dict)
    assert isinstance(products_coeff, dict)

def test_yield_effect():
    """Test that lower yield increases reactant mass needed"""
    def make_reaction(y):
        return Reaction(
            reactants=[Chemical(smiles="C"), Chemical(smiles="O=O")],
            byproducts=[Chemical(smiles="O")],
            wanted_product=ChemswithMass(smiles="C(=O)=O", initial_mass=44.0),
            Chosen_Yield=y
        )
    high_yield = make_reaction(1.0).total_mass_reactants()
    low_yield  = make_reaction(0.5).total_mass_reactants()
    assert low_yield > high_yield  

def test_get_CID_valid():
    """Test that a valid SMILES returns a CID from PubChem"""
    chem = Chemical(smiles="O")
    cid = chem.get_CID()
    assert cid is not None
    assert isinstance(cid, int)  

def test_get_CID_invalid():
    """Test that a valid SMILES does not returns None"""
    chem = Chemical(smiles="C")  
    cid = chem.get_CID()
    assert cid is None

def test_get_CID_stored():
    """Test that CID is stored in the object after calling get_CID"""
    chem = Chemical(smiles="O")
    chem.get_CID()
    assert hasattr(chem, "CID") 
    assert chem.CID is not None

def test_liquid_chemical_density_valid():
    """Test that density is calculated for a valid liquid molecule"""
    liquid = LiquidChemical(smiles="CCO")  
    assert liquid.density > 0.0
    assert isinstance(liquid.density, float)

def test_liquid_chemical_density_value():
    """Test that ethanol density is approximately correct (~0.789 g/mL)"""
    liquid = LiquidChemical(smiles="CCO")  
    assert 0.7 < liquid.density < 0.9     

def test_liquid_chemical_density_invalid_smiles():
    """Test that density defaults to 0.0 for a molecule with invalid SMILES"""
    liquid = LiquidChemical(smiles="FAKE") 
    assert liquid.density == 0.0

def test_evaluer_exp_plain_float():
    """Test that a plain float is returned correctly"""
    rech_var = re.compile(r'[a-zA-Z]')
    
    def evaluer_exp(Coef_with_DOF):
        Coef_as_str = str(Coef_with_DOF)
        Coef_as_str = rech_var.sub("1", Coef_as_str)
        try:
            return float(sympify(Coef_as_str))
        except Exception:
            try:
                return float(Coef_as_str)
            except Exception as e:
                raise ValueError(f"Cannot evaluate expression: {Coef_with_DOF}") from e

    assert evaluer_exp(2.5) == 2.5

def test_evaluer_exp_string_number():
    """Test that a string number is converted to float correctly"""
    rech_var = re.compile(r'[a-zA-Z]')

    def evaluer_exp(Coef_with_DOF):
        Coef_as_str = str(Coef_with_DOF)
        Coef_as_str = rech_var.sub("1", Coef_as_str)
        try:
            return float(sympify(Coef_as_str))
        except Exception:
            try:
                return float(Coef_as_str)
            except Exception as e:
                raise ValueError(f"Cannot evaluate expression: {Coef_with_DOF}") from e

    assert evaluer_exp("3") == 3.0

def test_evaluer_exp_single_variable():
    """Test that a single variable is replaced by 1"""
    rech_var = re.compile(r'[a-zA-Z]')

    def evaluer_exp(Coef_with_DOF):
        Coef_as_str = str(Coef_with_DOF)
        Coef_as_str = rech_var.sub("1", Coef_as_str)
        try:
            return float(sympify(Coef_as_str))
        except Exception:
            try:
                return float(Coef_as_str)
            except Exception as e:
                raise ValueError(f"Cannot evaluate expression: {Coef_with_DOF}") from e

    assert evaluer_exp("x") == 1.0

def test_evaluer_exp_invalid():
    """Test that an invalid expression raises ValueError"""
    rech_var = re.compile(r'[a-zA-Z]')

    def evaluer_exp(Coef_with_DOF):
        Coef_as_str = str(Coef_with_DOF)
        Coef_as_str = rech_var.sub("1", Coef_as_str)
        try:
            return float(sympify(Coef_as_str))
        except Exception:
            try:
                return float(Coef_as_str)
            except Exception as e:
                raise ValueError(f"Cannot evaluate expression: {Coef_with_DOF}") from e

    with pytest.raises(ValueError):
        evaluer_exp("@@@")
 


