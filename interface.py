import streamlit as st

st.set_page_config(page_title="Green Chemistry Calculator", page_icon=":leaves:")
st.title("🌿 Green Chemistry Calculator")
st.header(
    "Welcome to Green Chemistry Calculator a powerful tool that can calculate the greenness of your reaction!"
)
st.divider()
if "page_active" not in st.session_state:
    st.session_state.page_active = "Home"


def go_to(page_name):
    st.session_state.page_active = page_name


st.sidebar.title("Menu")

if st.sidebar.button("🏠 Home"):
    go_to("Home")

if st.sidebar.button("🧪 Stoichiometry"):
    go_to("Stoichiometry")
    st.title("Stoichiometry")
<<<<<<< HEAD
    st.write(
        "Write the chemical formula of your reactants, products and solvent to get the correct stoichiometry of the reaction!"
    )

=======
    st.write("Write the chemical formula of your reactants, products and solvent to get the correct stoichiometry of the reaction!")
    
>>>>>>> 209e6781c4d592477e60f1ca0a6df3ca6def18d9
if st.sidebar.button("⚛️ Molecular visulization"):
    go_to("Molecular visualization")

<<<<<<< HEAD
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

    tabs = st.tabs(
        [f"Reagent {i + 1}" for i in range(st.session_state.number_reagents)]
    )

    from streamlit_ketcher import st_ketcher

    for i, tab in enumerate(tabs):
        with tab:
            smiles = st_ketcher(key=f"editor_reagente_{i}")

            if smiles and smiles.strip() != "":
                st.success(f"**SMILES Reagent {i + 1}:** `{smiles}`")
            else:
                st.info(f"Please draw Reagent {i + 1} above.")

    if "number_solvent_and_catalyst" not in st.session_state:
        st.session_state.number_solvent_and_catalyst = 1

    def add_solvent_and_catalyst():
        if st.session_state.number_solvent_and_catlyst < 5:
            st.session_state.number_solvent_and_catalyst += 1

    def remove_solvent_and_catalyst():
        if st.session_state.number_solvent_and_catalyst > 1:
            st.session_state.number_solvent_and_catalyst -= 1

    col1, col2 = st.columns(2)
    with col1:
        st.button("➕ Add Solvent and/or Catalyst", on_click=add_solvent_and_catalyst)
    with col2:
        st.button(
            "➖ Remove Solvent and/or Catalyst", on_click=remove_solvent_and_catalyst
        )

    tabs = st.tabs(
        [
            f"Solvent {i + 1}"
            for i in range(st.session_state.number_solvent_and_catalyst)
        ]
    )

    from streamlit_ketcher import st_ketcher

    for i, tab in enumerate(tabs):
        with tab:
            smiles = st_ketcher(key=f"editor_solvent_and_catakyst_{i}")

            if smiles and smiles.strip() != "":
                st.success(f"**SMILES Reagent/Catalyst {i + 1}:** `{smiles}`")
            else:
                st.info(f"Please draw Reagent/Catalyst {i + 1} above.")

    if "number_products" not in st.session_state:
        st.session_state.number_products = 1

    def add_products():
        if st.session_state.number_products < 5:
            st.session_state.number_products += 1

    def remove_products():
        if st.session_state.number_products > 1:
            st.session_state.number_products -= 1
=======
if st.session_state.page_active == 'Molecular visualization':
        st.title("Draw Molecule")
        st.write("Draw your molecule to get started.")

    # 1. Initializaion session state
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
    
        tabs_p = st.tabs([f"Product {i+1}" for i in range(st.session_state.number_products)])
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
>>>>>>> 209e6781c4d592477e60f1ca0a6df3ca6def18d9

        st.divider()
        if st.button("🚀 Generate Reaction SMILES"):
            if reag_str and prod_str:
                st.success("### Reaction SMILES Created!")
                st.code(full_reaction)
                st.info(f"**Equation:** {reag_str} ⎯⎯⎯({solv_str if solv_str else 'no catalyst'})⎯⎯→ {prod_str}")
                st.session_state.final_smiles = full_reaction
                st.code(full_reaction, language="text")
                tutti_gli_smiles = st.session_state.reag_list + st.session_state.solv_list + st.session_state.prod_list
                lista_stringa = ", ".join(tutti_gli_smiles)
        
                st.write(f"**SMILES List:** {lista_stringa}")

<<<<<<< HEAD
    tabs = st.tabs(
        [f"Product {i + 1}" for i in range(st.session_state.number_reagents)]
    )

    from streamlit_ketcher import st_ketcher

    for i, tab in enumerate(tabs):
        with tab:
            smiles = st_ketcher(key=f"editor_product_{i}")

            if smiles and smiles.strip() != "":
                st.success(f"**SMILES Product {i + 1}:** `{smiles}`")
            else:
                st.info(f"Please draw Product {i + 1} above.")
=======
            else:
                st.warning("⚠️ Please draw at least one reactant and one product.")
>>>>>>> 209e6781c4d592477e60f1ca0a6df3ca6def18d9
