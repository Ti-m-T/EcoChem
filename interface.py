import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw


st.set_page_config(page_title="Green Chemistry Calculator", page_icon=":leaves:")
st.title("🌿 Green Chemistry Calculator")
st.header( "Welcome to Green Chemistry Calculator a powerful tool that can calculate the greenness of your reaction!")
st.divider()
#Settings of the Home page 
if "page_active" not in st.session_state: st.session_state.page_active="Home"
if "number_reagents" not in st.session_state: st.session_state.number_reagents = 1
if "number_solvent_and_catalyst" not in st.session_state: st.session_state.number_solvent_and_catalyst = 1
if "number_products" not in st.session_state: st.session_state.number_products = 1

if "reag_list" not in st.session_state: st.session_state.reag_list = []
if "solv_list" not in st.session_state: st.session_state.solv_list = []
if "prod_list" not in st.session_state: st.session_state.prod_list = []

def go_to(page_name):                            #Definition of a function to create buttons to redirect users to other pages
    st.session_state.page_active= page_name

st.sidebar.title("Menu")

if st.sidebar.button("🏠 Home"):
    go_to("Home")

if st.sidebar.button("🧪 Stoichiometry"):
    go_to("Stoichiometry")
    st.title("Stoichiometry")
    st.write(
        "Write the chemical formula of your reactants, products and solvent to get the correct stoichiometry of the reaction!"
    )

if st.sidebar.button("⚛️ Molecular visulization"):
    go_to("Molecular visualization")


if st.session_state.page_active == "Molecular visualization":
    st.title("Draw Molecule")
    st.write("Draw your molecule to get started.")
    tab1, tab2, tab3 = ["Reactants", "Solvent/Catalyst", "Products"]


    if "number_reagents" not in st.session_state:
        st.session_state.number_reagents = 1

    def add_reagents():
        if st.session_state.number_reagents < 5:
            st.session_state.number_reagents += 1

    def remove_reagents():
        if st.session_state.number_reagents > 1:
            st.session_state.number_reagents -= 1


    col1, col2 = st.columns(2)
    with col1:
        st.button("➕ Add Reagent", on_click=add_reagents)
    with col2:
        st.button("➖ Remove Reagent", on_click=remove_reagents)

    tabs_reag = st.tabs([f"Reagent {i + 1}" for i in range(st.session_state.number_reagents)])
    reag_list = []

    from streamlit_ketcher import st_ketcher

    for i, tab in enumerate(tabs_reag):
        with tab:
            smiles = st_ketcher(key=f"editor_reagente_{i}")

            if smiles and smiles.strip() != "":
                reag_list.append(smiles)
                st.success(f"**SMILES Reagent {i + 1}:** `{smiles}`")
            else:
                st.info(f"Please draw Reagent {i + 1} above.")
        st.session_state.reag_list = reag_list

    if "number_solvent_and_catalyst" not in st.session_state:
        st.session_state.number_solvent_and_catalyst = 1

    def add_solvent_and_catalyst():
        if st.session_state.number_solvent_and_catalyst < 5:
            st.session_state.number_solvent_and_catalyst += 1

    def remove_solvent_and_catalyst():
        if st.session_state.number_solvent_and_catalyst > 1:
            st.session_state.number_solvent_and_catalyst -= 1

    col3, col4 = st.columns(2)
    with col3:
        st.button("➕ Add Solvent and/or Catalyst", on_click=add_solvent_and_catalyst)
    with col4:
        st.button("➖ Remove Solvent and/or Catalyst", on_click=remove_solvent_and_catalyst)

    tabs_solv = st.tabs([f"Solvent {i + 1}"for i in range(st.session_state.number_solvent_and_catalyst)])
    solv_list = []

    from streamlit_ketcher import st_ketcher

    for i, tab in enumerate(tabs_solv):
        with tab:
            smiles = st_ketcher(key=f"editor_solvent_and_catalyst_{i}")

            if smiles and smiles.strip() != "":
                solv_list.append(smiles)
                st.success(f"**SMILES Solvent/Catalyst {i + 1}:** `{smiles}`")
            else:
                st.info(f"Please draw Solvent/Catalyst {i + 1} above.")
        st.session_state.solv_list = solv_list

    if "number_products" not in st.session_state:
        st.session_state.number_products = 1

    def add_products():
        if st.session_state.number_products < 5:
            st.session_state.number_products += 1

    def remove_products():
        if st.session_state.number_products > 1:
            st.session_state.number_products -= 1
    
    col5, col6 = st.columns(2)
    with col5: 
        st.button("➕ Add Product", on_click=add_products)
    with col6: 
        st.button("➖ Remove Product", on_click=remove_products)
    
    prod_labels = []
    for i in range(st.session_state.number_products):
        if i == 0:
            prod_labels.append("Main Product")
        else:
            prod_labels.append(f"Product {i + 1}")

    tabs_prod = st.tabs(prod_labels)
    prod_list = []
    for i, tab in enumerate(tabs_prod):
        with tab:
            smiles = st_ketcher(key=f"editor_product_{i}")
            if smiles: 
                prod_list.append(smiles)
            label = "Main Product" if i == 0 else f"Product {i + 1}"
            if smiles and smiles.strip() != "":
                st.success(f"**SMILES Product {i + 1}:** `{smiles}`")
            else:
                st.info(f"Please draw Product {i + 1} above.")
    st.session_state.prod_list = prod_list
    
    

    reag_str = ".".join([s for s in st.session_state.reag_list if s])
    solv_str = ".".join([s for s in st.session_state.solv_list if s])
    prod_str = ".".join([s for s in st.session_state.prod_list if s])
    full_reaction = f"{reag_str}>{solv_str}>{prod_str}"

    st.divider()
    if st.button("🚀 Generate Reaction SMILES"):
        if reag_str and prod_str:
            st.success("### Reaction SMILES Created!")
                
            try:
                reagents = [Chem.MolFromSmiles(s) for s in st.session_state.reag_list if s]
                solvents = [Chem.MolFromSmiles(s) for s in st.session_state.solv_list if s]
                products = [Chem.MolFromSmiles(s) for s in st.session_state.prod_list if s]
                rxn = AllChem.ChemicalReaction()
                    
            
                for m in reagents:
                    if m:
                        AllChem.Compute2DCoords(m)
                        rxn.AddReactantTemplate(m)
                for m in solvents:
                    if m:
                        AllChem.Compute2DCoords(m)
                        rxn.AddAgentTemplate(m)
                for m in products:
                    if m:
                        AllChem.Compute2DCoords(m)
                        rxn.AddProductTemplate(m)


                img = Draw.ReactionToImage(rxn, subImgSize=(400, 400),useSVG=False)

                st.image(img, use_container_width=True)
                    
                st.subheader("Complete Reaction Drawing!")
                    
            except Exception as e:
                st.error(f"Could not render reaction image{e}")

            st.info(f"**Equation:** {reag_str} ⎯⎯⎯({solv_str if solv_str else 'no catalyst'})⎯⎯→ {prod_str}")
            st.session_state.final_smiles = full_reaction
            all_smiles = st.session_state.reag_list + st.session_state.solv_list + st.session_state.prod_list
            string_list = ", ".join(all_smiles)
        
            st.write(f"**SMILES List:** {string_list}")

        else:
            st.warning("⚠️ Please draw at least one reactant and one product.")


