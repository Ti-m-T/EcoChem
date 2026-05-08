import streamlit as st
<<<<<<< HEAD

=======
>>>>>>> 209e6781c4d592477e60f1ca0a6df3ca6def18d9
st.title("Reaction builder")

num_r = st.number_input("Number of reactants", min_value=1, max_value=10, value=2)
num_p = st.number_input("Number of products", min_value=1, max_value=10, value=1)

reactants = []
products = []

# Number of reactants and products

st.subheader("Reactants")
for i in range(num_r):
    val = st.text_input(f"Reactant {i + 1}", key=f"r{i}")
    if val:
        reactants.append(val)

# Reactant input

st.subheader("Products")
for i in range(num_p):
    val = st.text_input(f"Product {i + 1}", key=f"p{i}")
    if val:
        products.append(val)

# Product input

reaction_dict = {"reactants": reactants, "products": products}

st.write("Reaction input:", reaction_dict)

# Création du dictionnaire de réaction à partir des inputs

if st.button("Compute"):
    reac, prod = calcul_coef_stoechio(reaction_dict)
    result: str = "Les coefficients stœchiométriques des réactifs et des produits sont respectivement donnés par :"
    for key_reac, value_reac in reac.items():
        result += f" {key_reac} : {value_reac} ,"
    for key_prod, value_prod in prod.items():
        result += f" {key_prod} : {value_prod} ,"
    result = result[:-1]

    st.write(result)


# On importe la fonction de calcul des coefficients stoechiométriques et on l'applique au dictionnaire de réaction créé à partir des inputs, affichant les résultats dans l'application Streamlit.

