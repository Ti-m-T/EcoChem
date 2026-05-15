import re
from chempy.chemistry import balance_stoichiometry
from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from sympy import sympify

dictionary_input_smiles: dict[str, list[str]] = {
    "reactants": ["CCO", "O=O", "CC=O"],
    "products": ["CC(=O)O", "O", "C=O"],
}


def split_input(
    dictionary_input_smiles: dict[str, list[str]],
) -> tuple[list[str], list[str]]:
    reactants: list[str] = dictionary_input_smiles.get("reactants", [])
    products: list[str] = dictionary_input_smiles.get("products", [])
    return (reactants, products)


# Split en deux classes le dictionnaire d'input


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


# print(smiles_to_brute(dictionary_input_smiles))
# Fonction qui passe de smiles a formule brute et retourne un dictionnarie de listes de formules brutes


def calcul_coef_stoechio(
    dictionary_input_smiles: dict[str, list[str]],
) -> dict[str,dict[str, int | str]]:
    dict_coef_stoechio: dict[str,dict[str, int | str]] = {"reactants": {}, "products": {}}
    dict_input_brute = smiles_to_brute(dictionary_input_smiles)
    reactant, product = split_input(dict_input_brute)
    coef_stoi_reac, coef_stoi_prod = balance_stoichiometry(reactant, product)

    # Calcul des coefficients stoechiométriques à partir des formules brutes des réactifs et des produits

    liste_smiles_tot: list[str] = []
    for key in dictionary_input_smiles["reactants"]:
        liste_smiles_tot.append(key)

    for key in dictionary_input_smiles["products"]:
        liste_smiles_tot.append(key)

    print(liste_smiles_tot)
    # Creation d'une liste des smiles des réactifs et des produits

    

    liste_brute_tot: list[str] = []
    for key in coef_stoi_reac.keys():
        liste_brute_tot.append(key)

    for key in coef_stoi_prod.keys():
        liste_brute_tot.append(key)

    print(coef_stoi_prod, coef_stoi_reac)

    print(liste_brute_tot)

    # Creation d'une liste de formules brutes associées

    for key, coef_stoi in coef_stoi_reac.items():
        Loc = liste_brute_tot.index(key)
        smiles_associé_au_brute = liste_smiles_tot[Loc]
        dict_coef_stoechio["reactants"][smiles_associé_au_brute] = coef_stoi

    for key, coef_stoi in coef_stoi_prod.items():
        Loc = liste_brute_tot.index(key)
        smiles_associé_au_brute = liste_smiles_tot[Loc]
        dict_coef_stoechio["products"][smiles_associé_au_brute] = coef_stoi

    

    # On relie la formule brute au smile par emplacement dans la liste

    return dict_coef_stoechio
# print(calcul_coef_stoechio(dictionary_input_smiles))
# Fonction qui calcule les coefficients stoechiométriques d'une réaction chimique donnée sous forme de dictionnaire de listes de réactifs et de produits, retourant les coefficients stochiométriques sous forme de dictionaires.
####


def calc_DOF_a_un(coef_dict: dict[str, dict[str, int | str]]) -> dict[str,dict[str,int]]:

    rech_var = re.compile(r"x\d+")
    # Recherche de variables de type x1,x2,....,xn

    def evaluer_exp(exp):
        exp_str = str(exp)

        if not rech_var.search(exp_str):
            try:
                return int(exp_str)
            except:
                return float(exp_str)
        
        # Cas ou il n'y pas de degré de libertée
        
        exp_str = rech_var.sub("1", exp_str)

        return int(sympify(exp_str))

        # Retourne les coef. stoi. sans modifications

    equation_sans_DOF : dict[str,dict[str,int]] = {"reactants": {}, "products": {}}

    for side in ["reactants", "products"]:
        for species, coef in coef_dict[side].items():
            equation_sans_DOF[side][species] = evaluer_exp(coef)
    return equation_sans_DOF

# print(calc_DOF_a_un(calcul_coef_stoechio(dictionary_input_smiles)))
# Fonction qui detecte les DOF dans la resolution des coef. stoic. et, s'il en existe, les égalise a 1

def calculate_molecular_weight(
    dictionary_input_smiles: dict[str, list[str]],
) -> dict[str, dict[str, float]]:
    mol_weights_dict: dict[str, dict[str, float]] = {"reactants": {}, "products": {}}
    for key in ["reactants", "products"]:
        for individual_smiles in dictionary_input_smiles.get(key, []):
            molecule = Chem.MolFromSmiles(individual_smiles)
            mol_weight = Descriptors.MolWt(molecule)
            mol_weights_dict[key][individual_smiles] = mol_weight
    return mol_weights_dict

# print(calculate_molecular_weight(dictionary_input_smiles))
# Fonction qui calcule le poids moléculaire de chaque réactif et produit d'une réaction chimique apartir d'un dictionnaire contenant les listes de SMILES des réactifs et des produits, retournant un dictionnaire contenant les poids moléculaires correspondants.


def calcul_eco_atom_reactants_M(dict_coef_stoic_reactants: dict[str, int],mol_weight_reactants : dict[str,float]) -> float:
    # print(dict_coef_stoic_reactants)
    # print(mol_weight_reactants)
    eco_atom_reactants_M: float = 0.0
    for key_reac, coef in dict_coef_stoic_reactants.items():
        eco_atom_reactants_M += coef * mol_weight_reactants.get(key_reac, ())
    return eco_atom_reactants_M

# print(calcul_eco_atom_reactants_M(dictionary_input_smiles))
# Fonction qui calcule la partie des réactif pour l'economie d'atome en fonciton de la masse molaire

def calculate_eco_atm_M(
    dictionary_input_smiles: dict[str, list[str]],
) -> float:
    dict_coef_stoic : dict[str,dict[str,int | str]] = calcul_coef_stoechio(dictionary_input_smiles)
    dict_coef_stoic_sans_DOF = calc_DOF_a_un(dict_coef_stoic)
    dict_coef_stoic_sans_DOF_reac: dict[str,int] = dict_coef_stoic_sans_DOF["reactants"] 
    mol_weight_dict : dict[str,dict[str,float]] = calculate_molecular_weight(dictionary_input_smiles)
    mol_weight_reactants : dict[str,float] = mol_weight_dict["reactants"] 
    eco_atom_reactants_M : float = calcul_eco_atom_reactants_M(dict_coef_stoic_sans_DOF_reac,mol_weight_reactants)

    prod_list: list[str] = dictionary_input_smiles["products"]
    main_prod_smiles : str = prod_list[0]
    main_prod_molecule: str = Chem.MolFromSmiles(main_prod_smiles)
    main_prod_M : float = Descriptors.MolWt(main_prod_molecule)
    dict_coef_stoic_sans_DOF_prod: dict[str,int] = dict_coef_stoic_sans_DOF["products"]
    main_prod_coef_sans_DOF : int = dict_coef_stoic_sans_DOF_prod[main_prod_smiles]
    eco_atom_main_prod_M : float = main_prod_M * main_prod_coef_sans_DOF
    print(main_prod_coef_sans_DOF , main_prod_M)
    eco_atom_M : float = (eco_atom_main_prod_M*100)/eco_atom_reactants_M
    return eco_atom_M

print(calculate_eco_atm_M(dictionary_input_smiles))

# test:list[str]= ["ea", "eb", "ec", "ed", "ee"]
# print(test.index("ec"))
