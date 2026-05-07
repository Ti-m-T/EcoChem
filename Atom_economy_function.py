from unittest import result
from chempy.chemistry import balance_stoichiometry
from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

reactants: list[str] = ["H2", "O2"]
products: list[str] = ["H2O"]
dictionary_input: dict[list[str]] = {"reactants": reactants, "products": products}
dictionary_input_smiles: dict[str, list[str]] = {
    "reactants": ["[H][H]", "O=O"],
    "products": ["O"],
}


def smiles_to_brute(
    dictionary_input_smiles: dict[str, list[str]],
) -> dict[str, list[str]]:

    dictionary_input_brute: dict[str, list[str]] = {}

    for key, liste_smiles in dictionary_input_smiles.items():
        dictionary_input_brute[key] = []

        for smiles in liste_smiles:
            molecule = Chem.MolFromSmiles(smiles)
            formule_brute = CalcMolFormula(molecule)
            dictionary_input_brute[key].append(formule_brute)

    return dictionary_input_brute


# Fonction qui passe de smiles a formule brute et retourne un dictionnarie de listes de formules brutes

print(smiles_to_brute(dictionary_input_smiles))


def split_input(dictionary_input: dict[str, list[str]]) -> tuple[list[str], list[str]]:
    reactants: list[str] = dictionary_input.get("reactants", [])
    products: list[str] = dictionary_input.get("products", [])
    return (reactants, products)


# Split en deux classes le dictionnaire d'input


def calcul_coef_stoechio(
    dictionary_input: dict[str, list[str]],
) -> tuple[dict[str, int], dict[str, int]]:
    reactants, products = split_input(dictionary_input)
    reac, prod = balance_stoichiometry(reactants, products)
    return reac, prod


# print(calcul_coef_stoechio(dictionary_input))
# Fonction qui calcule les coefficients stoechiométriques d'une réaction chimique donnée sous forme de dictionnaire de listes de réactifs et de produits, retourant les coefficients stochiométriques sous forme de dictionaires.


def coef_stoechio_reactants(dictionary_input: dict[str, list[str]]) -> list[int]:
    list_coef_reactants: list[int] = []
    for key in calcul_coef_stoechio(dictionary_input)[0]:
        list_coef_reactants.append(calcul_coef_stoechio(dictionary_input)[0][key])
    return list_coef_reactants


# Fonction qui calcule les coefficients stoechiométriques des réactifs d'une réaction chimique donnée sous forme de dictionnaire de listes de réactifs et de produits, retourant les coefficients stochiométriques des réactifs sous forme de liste d'entiers.

### print(Reaction(*balance_stoichiometry(['H2', 'O2'], ['H2O'])).string())
# Fonction qui calcule les coefficients stoechiométriques d'une réaction chimique donnée sous forme de listes de réactifs et de produits, retournant l'equation chimique équilibrée sous forme de string.


def calculate_molecular_weight(
    dictionary_input: dict[str, list[str]],
) -> dict[str, list[float]]:
    mol_weights_dict: dict[str, list[float]] = {"reactants": [], "products": []}
    for key in ["reactants", "products"]:
        for individual_smiles in dictionary_input.get(key, []):
            molecule = Chem.MolFromSmiles(individual_smiles)
            mol_weight = Descriptors.MolWt(molecule)
            mol_weights_dict[key].append(mol_weight)
    return mol_weights_dict


# Fonction qui calcule le poids moléculaire de chaque réactif et produit d'une réaction chimique apartir d'un dictionnaire contenant les listes de SMILES des réactifs et des produits, retournant un dictionnaire contenant les poids moléculaires correspondants.
### print(calculate_molecular_weight(dictionary_input_smiles))


def calcul_eco_atom_reactants(dictionary_input: dict[str, list[str]]) -> float:
    list_coef_reactants = coef_stoechio_reactants(dictionary_input)
    mol_weight_dict = calculate_molecular_weight(dictionary_input)
    mol_weight_reactants = mol_weight_dict.get("reactants", [])
    eco_atom_reactants: float = 0.0
    for i in list_coef_reactants:
        eco_atom_reactants += list_coef_reactants[i] * mol_weight_reactants[i]
    return eco_atom_reactants


# need of the smiles / brute conversion

# def calculate_eco_atm(smiles_of_reaction: dict[str,list[str]]) -> tuple[dict[str, float], dict[str, float]]:
#    reac,prod = calcul_coef_stoechio(smiles_of_reaction.get("reactants", []), smiles_of_reaction.get("products", []))
#    mol_weights_dict = calculate_molecular_weight(smiles_of_reaction)
# for i in range(len(smiles_of_reaction.get("reactants",[]))):


# test:list[str]= ["ea", "eb", "ec", "ed", "ee"]
# print(test.index("ec"))
