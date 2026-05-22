from altair import value
import streamlit as st
import sys
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import requests
import os 
import random
from PIL import Image
from rdkit.Chem import rdMolDescriptors
import re
import plotly.graph_objects as go
from dataclasses import dataclass, field
from thermo import chemical as density_finder
from experiment import Reaction
from experiment import Chemical
from experiment import ChemswithMass
#from experiment import LiquidChemical
from experiment import Solvent
from experiment import Extractant
from streamlit_ketcher import st_ketcher


st.set_page_config(page_title="EcoChem", page_icon=":leaves:")
st.image("EcoChem.jpg")
st.header("Welcome to EcoChem a powerful tool that can calculate the greenness of your reaction!")
st.divider()
st.write("""
With **EcoChem**, you can define how "green" your reaction is by calculating:
* **PMI** (Process Mass Intensity)
* **E-factor**
* **Atom Economy**, and more!

---

### 💡 Features:

* **Stoichiometry**: Don't know how to equilibrate your reaction? No problem! We will do it for you. Just head over to the **Stoichiometry** section.
* **Molecular Visualization**: If you want to visualize the reaction or you only have the molecular structures, go to the **Molecular Visualization** page.

Have fun using **EcoChem**! 🚀
""")


#Settings of the Home page 
if "page_active" not in st.session_state: st.session_state.page_active="Home"
if "number_reagents" not in st.session_state: st.session_state.number_reagents = 1
if "number_solvent_and_catalyst" not in st.session_state: st.session_state.number_solvent_and_catalyst = 1
if "number_products" not in st.session_state: st.session_state.number_products = 1

if "reag_list" not in st.session_state: st.session_state.reag_list = []
if "solv_list" not in st.session_state: st.session_state.solv_list = []
if "prod_list" not in st.session_state: st.session_state.prod_list = []

def convert_to_latex_subscripts(formula: str) -> str: #Function to write the molecular formulas in "latex" style
    return re.sub(r'(\d+)', r'_{\1}', formula)

def name_to_smiles(name: str) -> str: #Function that returns the smile given a chemical name
    try:
        chemical_name = name.strip()
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{chemical_name}/property/CanonicalSMILES/TXT"
        response = requests.get(url, timeout=5)

        if response.status_code == 200:
            return response.text.strip()
        else:
            return None
    except Exception:
        return None


def go_to(page_name):   #Definition of a function to create buttons to redirect users to other pages
    st.session_state.page_active = page_name

def get_formula(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
            latex_formula = re.sub(r'(\d+)', r'_{\1}', formula)
            return latex_formula
        return smiles
    except:
        return smiles


st.sidebar.title("Menu")

if st.sidebar.button("🏠 Home"):
    go_to("Home")


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

                st.image(img, width="stretch")
                    
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

all_smiles_dict = {"Reagent":[], "Solvent":[], "Product":[]}
if "reag_list" in st.session_state:
    all_smiles_dict["Reagent"] = [s for s in st.session_state.reag_list if s]

if "solv_list" in st.session_state:
    all_smiles_dict["Solvent"] = [s for s in st.session_state.solv_list if s]
        
if "prod_list" in st.session_state:
    all_smiles_dict["Product"] = [s for s in st.session_state.prod_list if s]
print(all_smiles_dict)

if "extr_list" not in st.session_state:
    st.session_state.extr_list = []


st.subheader("🧪 Add Extractants (for Work-up)")

extr_input = st.text_input("Extractant (Name or SMILES)", key="input_extr_name", placeholder="e.g., ethyl acetate, hexane, etc...")
extr_volume = st.number_input("Volume (mL)", min_value=0.0, step=10.0, key="input_extr_volume")
extr_density = st.number_input("Density (g/mL)", min_value=0.0, step=0.1, key="input_extr_density")
    
if st.button("➕ Add Extractant", key="btn_add_extractant_main"):
    if extr_input:
        calculated_smiles = name_to_smiles(extr_input)
        if not calculated_smiles:
            from rdkit import Chem
            if Chem.MolFromSmiles(extr_input):
                calculated_smiles = extr_input
        if calculated_smiles:
            st.session_state.extr_list.append({"smiles": calculated_smiles, "volume": extr_volume, "density": extr_density})
            st.success(f"Added successfully: `{calculated_smiles}`")
            st.rerun()
        else:
            st.error(f"❌ Could not resolve '{extr_input}'. Check the spelling.")
    else:
        st.error("Please enter a name or SMILES first.")

st.markdown("---")
if len(st.session_state.extr_list) > 0:
    st.write("### 📋 Current Extractants List:")
    with st.container():
        for i, e in enumerate(st.session_state.extr_list):

            col_text, col_btn = st.columns([0.8, 0.2])

            with col_text:

                st.markdown(f"**{i+1}.** `{e['smiles']}` — **{e['volume']} mL** — **{e['density']} g/mL**")

            with col_btn:
                unique_key = f"remove_btn_{i}"

                if st.button("❌ Remove", key=unique_key):
                    st.session_state.extr_list.pop(i)
                    st.rerun() 
else:
    st.info("No extractants added yet.")

if st.sidebar.button("🧪 Stoichiometry"):
    go_to("Stoichiometry")
    st.title("Stoichiometry")
    st.write("Write the chemical formula of your reactants, products and solvent to get the correct stoichiometry of the reaction!")

st.subheader("Yield of reaction")
yield_input = st.number_input("Reaction Yield [%]", value=100.0, step = 0.1)

if yield_input<=0 or yield_input>100:
    st.warning("⚠️ Please enter a valid yield value between 0 and 100 %.")
    st.stop()

else: 
    yield_fraction = yield_input / 100
    st.success(f"Reaction yield set to {yield_input:.2f}%.")



def build_reaction():
    reag_smiles_list = st.session_state.get("reag_list", [])
    prod_smiles_list = st.session_state.get("prod_list", [])
    
    reactants_objects = [Chemical(smiles=s) for s in reag_smiles_list]

    byproducts_objects = [Chemical(smiles=s) for s in prod_smiles_list[1:]]

    wanted_product_obj = ChemswithMass(
        smiles = prod_smiles_list[0],
        initial_mass = 1000 # For a mass of 1kg of the main product
    )

    return Reaction(
        reactants = reactants_objects,
        wanted_product = wanted_product_obj,
        byproducts = byproducts_objects,
        Catalysts = [ChemswithMass(smiles=s, initial_mass=m) for s, m in st.session_state.get("cat_list", [])],
        solvents = [Solvent(smiles=s, volume=v) for s, v in st.session_state.get("solv_list", [])],
        extractants =[Extractant(smiles=e["smiles"], volume=e["volume"], user_density=e["density"]) for e in st.session_state.get("extr_list", [])],
        Chosen_Yield = yield_fraction  
    )

if st.button("⚗️ Compute reaction stoichiometry, green metrics and risk assessment"):
    
    if len(st.session_state.get("reag_list", [])) == 0:
        st.warning("⚠️ Please draw at least one reactant to proceed.") # Checks for reactants
        st.stop()

    elif len(st.session_state.get("prod_list", [])) == 0:
        st.warning("⚠️ Please draw at least one product to proceed.") # Checks for products
        st.stop()

    for reactant in st.session_state.get("reag_list", []):
        if "." in reactant:
            index = reactant.index(".")
            st.warning(f"⚠️ The reactant {index} seems to contain multiple SMILES. Please separate them into different input fields.") # Checks for multiple SMILES in the same input for reactants.
            st.stop()

    for product in st.session_state.get("prod_list", []):
        if "." in product:
            index = product.index(".")
            st.warning(f"⚠️ The product {index} seems to contain multiple SMILES. Please separate them into different input fields.") # Checks for multiple SMILES in the same input for products.
            st.stop()

    experiment = build_reaction()

    reactant_set : set = {reactant.smiles for reactant in experiment.reactants}
    
    all_products_list = experiment.byproducts + [experiment.wanted_product]

    all_products_set : set = {products.smiles for products in all_products_list}

    if reactant_set & all_products_set:
        st.warning(f"⚠️ {reactant_set & all_products_set} appear both as reactants and products. Please check your input.") # Checks for the same molecule in reactants and products.
        st.stop()

    def atom_set_builder(st_list_input :list[Chemical]) -> set:

        atom_set = set()

        list_input = [Chemical(smiles=s) for s in st_list_input]

        for smiles_input in list_input:
            molecular_formula = smiles_input.mol_f
            atom = re.findall(r'[A-Z][a-z]?', molecular_formula) # Builds a set of the atoms present in the list of reactants or products.
            atom_set.update(atom)
        return atom_set
    
    atom_set_reactants = atom_set_builder(st.session_state.get("reag_list", []))

    atom_set_products = atom_set_builder(st.session_state.get("prod_list", []))

    if atom_set_reactants != atom_set_products:

        diff_in_atom = atom_set_reactants ^ atom_set_products
    
        diff_atom_list :list[str] = []

        for atom in diff_in_atom:

            diff_atom_list.append(atom)

            text_diff_atom = ", ".join(f'"{x}"' for x in diff_atom_list[:-1]) + f' and "{diff_atom_list[-1]}"'

        st.warning(f"⚠️ Atom {text_diff_atom} appear only in the reactants or only in the products. Please check your input.") # Checks for the same atoms in reactants and products.
        st.stop()

    if (
    len(experiment.reactants) > 0
    and all(r.mol for r in experiment.reactants)
    and experiment.wanted_product.mol
    ): 
        
        st.success("All reactants and main product are valid. Proceeding with calculations...")

        try:

            # 1. STOICHIOMETRY
            reactants_coeff, products_coeff = experiment.stoich_of_reaction()
            st.session_state.stoich_results = {
                "reactants": reactants_coeff,
                "products": products_coeff
            }

            # 2. GREEN METRICS
            ae_result = experiment.calcul_eco_atom()
            pmi_result = experiment.PMI()
            ef_result = experiment.e_factor()

            st.session_state.results = {
                "Atom Economy": ae_result,
                "PMI": pmi_result,
                "E-Factor": ef_result
            }

            # 3. DISPLAY STOICH
            st.subheader("⚖️ Balanced Reaction")

            def format_reaction_latex(experiment):
                def fmt(c, f):
                    return rf"{int(round(c))}\,{convert_to_latex_subscripts(f)}"

                reactants = " + ".join(fmt(r.coeff, r.mol_f) for r in experiment.reactants)

                products = " + ".join(fmt(p.coeff, p.mol_f) for p in experiment.byproducts + [experiment.wanted_product])

                return rf"{reactants} \;\longrightarrow\; {products}"
            
            st.latex(format_reaction_latex(experiment))

            st.success("Reaction Balanced!")

            st.subheader("🌿 Green Metrics")

            def show_skulls(num_skulls=20): #Definition of a function to pop skulls if the economy atom of the reaction is below 40
                skulls_html = f"""
                <div id="skulls-container-{random.randint(0, 1000)}" class="skulls-temporary">
                    {"".join([f'<span style="left:{random.uniform(5,95)}%; animation-delay:{random.uniform(0,2)}s;">💀</span>' for _ in range(num_skulls)])}
                </div>
                
                <style>
                .skulls-temporary {{
                    position: fixed;
                    top: 0; left: 0; width: 100%; height: 100%;
                    pointer-events: none; z-index: 9999;
                }}
                .skulls-temporary span {{
                    position: absolute;
                    bottom: -10%;
                    font-size: 2.5rem;
                    /* L'animazione viene eseguita 1 volta sola e si ferma alla fine (forwards) */
                    animation: floatUpOnce 4s ease-in forwards 1; 
                }}

                @keyframes floatUpOnce {{
                    0% {{ bottom: -10%; opacity: 0; transform: scale(0.5); }}
                    20% {{ opacity: 1; }}
                    80% {{ opacity: 1; }}
                    100% {{ bottom: 110%; opacity: 0; transform: scale(1.2); }}
                }}
                </style>
                """
                st.markdown(skulls_html, unsafe_allow_html=True)

            def display_linear_gauge(value, title="Atom Economy"): #Colored bar function
                pos = max(0, min(100, value)) #Def of min and max values
        
                if pos < 50: color = "#ff4b4b"   # Color of the text
                elif pos < 80: color = "#ffa500" 
                else: color = "#00cc96"          
        
                gauge_html = f"""
                <div style="font-family: sans-serif; margin: 20px 0; padding: 10px; background-color: #f9f9f9; border-radius: 10px;">
                    <div style="display: flex; justify-content: space-between; margin-bottom: 15px;">
                        <span style="font-weight: bold; font-size: 1.1rem; color: #333;">{title}</span>
                        <span style="font-weight: bold; font-size: 1.3rem; color: {color};">{value:.1f}%</span>
                    </div>
                    <div style="position: relative; height: 25px; width: 100%; border-radius: 12px; 
                                background: linear-gradient(to right, #ff4b4b 0%, #ffeb3b 50%, #00cc96 100%); 
                                box-shadow: inset 0 2px 4px rgba(0,0,0,0.1);">
                        <div style="position: absolute; left: calc({pos}% - 10px); top: -12px; 
                                    width: 0; height: 0; 
                                    border-left: 10px solid transparent;
                                    border-right: 10px solid transparent;
                                    border-top: 15px solid #333;">
                        </div>
                    </div>
                    <div style="display: flex; justify-content: space-between; font-size: 0.75rem; margin-top: 8px; color: #888;">
                        <span>Low Atom Economy</span>
                        <span>High Atom Economy</span>
                    </div>
                </div>
        """
                st.markdown(gauge_html, unsafe_allow_html=True)

            try :

                display_linear_gauge(ae_result, "Atom Economy")

                if ae_result > 90: 
                    st.success("Excellent! This reaction is very atom-efficient.")

                elif ae_result > 40: 
                    st.warning("⚠️ This reaction produces a significant amount of waste.")

                else : 
                    st.error(" 💀 This reaction is extremely inefficient in terms of atom economy.")
                    
            except Exception as e:
                st.error(f"Error during the calculation: {e}")

            def display_linear_gauge_pmi(value, title="PMI"): #Colored bar function
                pos = max(0, min(100, value)) #Def of min and max values
                
                if pos > 50: color = "#ff4b4b"   # Color of the text
                elif pos > 80: color = "#ffa500" 
                else: color = "#00cc96"          
                
                gauge_html = f"""
                <div style="font-family: sans-serif; margin: 20px 0; padding: 10px; background-color: #f9f9f9; border-radius: 10px;">
                    <div style="display: flex; justify-content: space-between; margin-bottom: 15px;">
                        <span style="font-weight: bold; font-size: 1.1rem; color: #333;">{title}</span>
                        <span style="font-weight: bold; font-size: 1.3rem; color: {color};">{value:.1f}%</span>
                    </div>
                    <div style="position: relative; height: 25px; width: 100%; border-radius: 12px; 
                                background: linear-gradient(to right, #00cc96 0%, #ffeb3b 50%, #ff4b4b 100%); 
                                box-shadow: inset 0 2px 4px rgba(0,0,0,0.1);">
                        <div style="position: absolute; left: calc({pos}% - 10px); top: -12px; 
                                    width: 0; height: 0; 
                                    border-left: 10px solid transparent;
                                    border-right: 10px solid transparent;
                                    border-top: 15px solid #333;">
                        </div>
                    </div>
                    <div style="display: flex; justify-content: space-between; font-size: 0.75rem; margin-top: 8px; color: #888;">
                        <span>Low PMI</span>
                        <span>High PMI</span>
                    </div>
                </div>
                """
                st.markdown(gauge_html, unsafe_allow_html=True)
            try :
                display_linear_gauge_pmi(pmi_result, "PMI")

                if pmi_result <= 10:    
                    st.success("Fantastic! Low material usage.")

                elif pmi_result <= 40:
                    st.warning("⚠️ This reaction uses a moderate amount of material relative to the target product.")

                elif pmi_result <= 70:
                    st.warning("⚠️ This reaction uses a significant amount of material relative to the target product.")

                else:
                    st.error("💀 This reaction uses a very high amount of material relative to the target product.")
            
            except Exception as e:

                st.error(f"Error during PMI calculation: {e}")

            def display_linear_gauge_efactor(value, title="E-Factor"): #Colored bar function
                pos = max(0, min(100, value)) #Def of min and max values
                
                if pos > 50: color = "#ff4b4b"   # Color of the text
                elif pos > 80: color = "#ffa500" 
                else: color = "#00cc96"          
                
                gauge_html = f"""
                <div style="font-family: sans-serif; margin: 20px 0; padding: 10px; background-color: #f9f9f9; border-radius: 10px;">
                    <div style="display: flex; justify-content: space-between; margin-bottom: 15px;">
                        <span style="font-weight: bold; font-size: 1.1rem; color: #333;">{title}</span>
                        <span style="font-weight: bold; font-size: 1.3rem; color: {color};">{value:.1f}%</span>
                    </div>
                    <div style="position: relative; height: 25px; width: 100%; border-radius: 12px; 
                                background: linear-gradient(to right, #00cc96 0%, #ffeb3b 50%, #ff4b4b 100%); 
                                box-shadow: inset 0 2px 4px rgba(0,0,0,0.1);">
                        <div style="position: absolute; left: calc({pos}% - 10px); top: -12px; 
                                    width: 0; height: 0; 
                                    border-left: 10px solid transparent;
                                    border-right: 10px solid transparent;
                                    border-top: 15px solid #333;">
                        </div>
                    </div>
                    <div style="display: flex; justify-content: space-between; font-size: 0.75rem; margin-top: 8px; color: #888;">
                        <span>Low E-Factor</span>
                        <span>High E-Factor</span>
                    </div>
                </div>
                """
                st.markdown(gauge_html, unsafe_allow_html=True)

            try:
                display_linear_gauge_efactor(ef_result, "E-Factor")

                if ef_result <= 10:
                    st.success("Excellent! This reaction generates very little waste.")

                elif ef_result < 50:
                    st.warning("⚠️ This reaction generates a significant amount of waste relative to the target product.")

                else :
                    st.error("💀 This reaction generates a very high amount of waste relative to the target product.")

            except Exception as e:
                st.error(f"Error during E-Factor calculation: {e}")

            #st.metric("PMI", f"{pmi_result:.1f}%")

            #st.metric("E-Factor", f"{ef_result:.1f}%")   # display of the E-Factor as a metric with a color code (green if low, red if high)

            # 5. EFFECTS

            if ae_result >= 90 and ef_result <= 10 and pmi_result <= 10:
                st.balloons()

            elif ae_result < 40 and ef_result >= 50 and pmi_result > 70:
                show_skulls()

            st.subheader("Risk Assessment")

            # To be implemented: risk assessment based on the GSK solvent selection guide and the hazard of the reactants and products.

        except Exception as e:

            if str(e) == "Failed to balance reaction":
                st.error("⚠️ Could not balance the reaction. Please check your input for any errors.")
            else:
                st.error(f"Error during computation: {e}")








