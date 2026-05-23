# tests/test_ecochem.py
import pytest
from experiment import Chemical, ChemswithMass, Reaction 
def test_chemical_valid_smiles():
    """Test that a valid SMILES creates a proper Chemical object"""
    chem = Chemical(smiles="C")  
    assert chem.mol is not None
    assert chem.mw > 0
    assert chem.mol_f == "CH4"

def test_chemical_invalid_smiles():
    """Test that an invalid SMILES is handled gracefully"""
    chem = Chemical(smiles="INVALID")
    assert chem.mw == 0.0
    assert chem.smiles == "Invalid"

def test_chemical_molecular_weight():
    """Test molecular weight calculation for water"""
    water = Chemical(smiles="O")  # water
    assert round(water.mw, 2) == 18.01

def test_chemical_mass_setter():
    """Test that setting mass also calculates moles"""
    chem = Chemical(smiles="O")  # water, mw ~18
    chem.mass = 18.0
    assert round(chem.moles, 1) == 1.0

def test_chemical_logp():
    """Test that logP is calculated"""
    chem = Chemical(smiles="C")
    assert isinstance(chem.logp, float)

def test_chemical_nb_atom():
    """Test atom count including hydrogens"""
    methane = Chemical(smiles="C")  # CH4 = 5 atoms
    assert methane.nb_atom == 5

def test_chemswithmass_sets_mass():
    """Test that initial_mass is set correctly"""
    chem = ChemswithMass(smiles="O", initial_mass=18.0)
    assert chem.mass == 18.0

def test_chemswithmass_calculates_moles():
    """Test that moles are calculated from initial_mass"""
    chem = ChemswithMass(smiles="O", initial_mass=18.0)
    assert round(chem.moles, 1) == 1.0

@pytest.fixture
def combustion_reaction():
    """
    Methane combustion: CH4 + 2O2 → CO2 + 2H2O
    Simple, well known reaction good for testing
    """
    reactants = [
        Chemical(smiles="C"),    # methane
        Chemical(smiles="O=O"),  # oxygen
    ]
    byproducts = [
        Chemical(smiles="O"),    # water
    ]
    wanted_product = ChemswithMass(smiles="C(=O)=O", initial_mass=44.0)  # CO2, ~1 mole
    return Reaction(
        reactants=reactants,
        byproducts=byproducts,
        wanted_product=wanted_product
    )

def test_atom_economy(combustion_reaction):
    """Test atom economy calculation"""
    ae = combustion_reaction.calcul_eco_atom()
    assert 0 < ae <= 100  # must be a percentage

def test_pmi(combustion_reaction):
    """Test PMI is greater than 1 (can never be less)"""
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
    assert low_yield > high_yield  # lower yield = more reactants needed

