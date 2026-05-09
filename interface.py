import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import pubchempy as pcp

st.set_page_config(page_title="Green Chemistry Calculator", page_icon=":leaves:")
st.title("🌿 Green Chemistry Calculator")
st.header("Welcome to Green Chemistry Calculator a powerful tool that can calculate the greenness of your reaction!")
st.divider()
#Settings of the Home page 
if "page_active" not in st.session_state: st.session_state.page_active="Home"



def go_to(page_name):                            #Definition of a function to create buttons to redirect users to other pages
    st.session_state.page_active= page_name

st.sidebar.title("Menu")

if st.sidebar.button("🏠 Home"):
    go_to('Home')

if st.sidebar.button("🧪 Stoichiometry"):
    go_to('Stoichiometry')
    st.title("Stoichiometry")
    st.write("Write the chemical formula of your reactants, products and solvent to get the correct stoichiometry of the reaction!")
    
if st.sidebar.button("⚛️ Molecular visulization"):
     go_to("Molecular visualization")

if st.session_state.page_active == 'Molecular visualization':
        st.title("Draw Molecule")
        st.write("Draw the reactants, solvent/catalysts and products of your reaction to get the SMILES!")

    # Initializaion session state
        if 'reag_list' not in st.session_state: st.session_state.reag_list = []
        if 'solv_list' not in st.session_state: st.session_state.solv_list = []
        if 'prod_list' not in st.session_state: st.session_state.prod_list = []
        if "number_reagents" not in st.session_state: st.session_state.number_reagents = 1
        if "number_solvent_and_catalyst" not in st.session_state: st.session_state.number_solvent_and_catalyst = 1
        if "number_products" not in st.session_state: st.session_state.number_products = 1

    # Add/Remove functions
        def add_reagents(): st.session_state.number_reagents = min(st.session_state.number_reagents + 1, 5)
        def remove_reagents(): st.session_state.number_reagents = max(st.session_state.number_reagents - 1, 1)
        def add_solv(): st.session_state.number_solvent_and_catalyst = min(st.session_state.number_solvent_and_catalyst + 1, 5)
        def remove_solv(): st.session_state.number_solvent_and_catalyst = max(st.session_state.number_solvent_and_catalyst - 1, 1)
        def add_prod(): st.session_state.number_products = min(st.session_state.number_products + 1, 5)
        def remove_prod(): st.session_state.number_products = max(st.session_state.number_products - 1, 1)

        from streamlit_ketcher import st_ketcher

    # Reagents
        st.subheader("Reactants")
        c1, c2 = st.columns(2)
        c1.button("➕ Add Reagent", on_click=add_reagents, key="ar")
        c2.button("➖ Remove Reagent", on_click=remove_reagents, key="rr")
    
        tabs_r = st.tabs([f"Reagent {i+1}" for i in range(st.session_state.number_reagents)])
        current_reag = []
        for i, tab in enumerate(tabs_r):
            with tab:
                sm = st_ketcher(key=f"reag_tab_{i}")
                if sm and sm.strip():
                    st.success(f"SMILES: `{sm}`")
                    current_reag.append(sm)
        st.session_state.reag_list = current_reag

    # Solvents/Catalysts
        st.subheader("Solvents & Catalysts")
        c3, c4 = st.columns(2)
        c3.button("➕ Add Solvent", on_click=add_solv, key="as")
        c4.button("➖ Remove Solvent", on_click=remove_solv, key="rs")
    
        tabs_s = st.tabs([f"Solvent {i+1}" for i in range(st.session_state.number_solvent_and_catalyst)])
        current_solv = []
        for i, tab in enumerate(tabs_s):
            with tab:
                sm = st_ketcher(key=f"solv_tab_{i}")
                if sm and sm.strip():
                    st.success(f"SMILES: `{sm}`")
                    current_solv.append(sm)
        st.session_state.solv_list = current_solv

    # Products
        st.subheader("Products")
        c5, c6 = st.columns(2)
        c5.button("➕ Add Product", on_click=add_prod, key="ap")
        c6.button("➖ Remove Product", on_click=remove_prod, key="rp")

        tabs_names=["Main Product"] + [f"Product {i+1}" for i in range(1,st.session_state.number_products)]
        tabs_p = st.tabs(tabs_names)
        current_prod = []
        for i, tab in enumerate(tabs_p):
            with tab:
                sm = st_ketcher(key=f"prod_tab_{i}")
                if sm and sm.strip():
                    st.success(f"SMILES: `{sm}`")
                    current_prod.append(sm)
        st.session_state.prod_list = current_prod

    # Reaction
        reag_str = ".".join(st.session_state.reag_list)
        solv_str = ".".join(st.session_state.solv_list)
        prod_str = ".".join(st.session_state.prod_list)
        full_reaction = f"{reag_str}>{solv_str}>{prod_str}"

        st.divider()
        if st.button("🚀 Generate Reaction SMILES"):
            if reag_str and prod_str:
                st.success("### Reaction SMILES Created!")
                
                try:
                    reagents = [Chem.MolFromSmiles(s) for s in st.session_state.reag_list if s]
                    solvents = [Chem.MolFromSmiles(s) for s in st.session_state.solv_list if s]
                    products = [Chem.MolFromSmiles(s) for s in st.session_state.prod_list if s]
                    
                    if reagents and products:
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

                        img = Draw.ReactionToImage(
                    rxn, 
                    subImgSize=(400, 400),
                    useSVG=False)

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