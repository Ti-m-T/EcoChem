import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import os 
import random
from PIL import Image
from Atom_economy_function import calcul_coef_stoechio, calculate_eco_atm_M
from rdkit.Chem import rdMolDescriptors
import re
import plotly.graph_objects as go


st.set_page_config(page_title="Green Chemistry Calculator", page_icon=":leaves:")
st.image(r"C:\Users\cvitt\Green_chem_project\EcoChem.jpg")
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


    from streamlit_ketcher import st_ketcher
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
if "solv_list" in st.session_state.solv_list:
    all_smiles_dict["Solvent"] = [s for s in st.session_state.solv_list if s]
        
if "prod_list" in st.session_state:
    all_smiles_dict["Product"] = [s for s in st.session_state.prod_list if s]
print(all_smiles_dict)



if st.sidebar.button("🧪 Stoichiometry"):
    go_to("Stoichiometry")
    st.title("Stoichiometry")
    st.write("Write the chemical formula of your reactants, products and solvent to get the correct stoichiometry of the reaction!")

input_data = {"reactants": st.session_state.get("reag_list", []), "products":st.session_state.get("prod_list", [])}

if st.button("⚖️ Calculate Coefficients"):
    if input_data["reactants"] and input_data["products"]:
        try:
            results = calcul_coef_stoechio(input_data)
            st.success("Reaction Balanced!")
            col1, col2 = st.columns(2)
            with col1:
                st.subheader("Reactants")
                for smiles, coef in results["reactants"].items():
                    formula = get_formula(smiles)
                    st.latex(rf"{coef} \text{{ }} {formula}")
            with col2:
                st.subheader("Products")
                for smiles, coef in results["products"].items():
                    formula = get_formula(smiles)
                    st.latex(rf"{coef} \text{{ }} {formula}")
            st.session_state.stoich_results = results
        except Exception as e:
            st.error(f"Error during balancing: {e}")
            st.info("Check if all atoms in reactants are also present in products.")
        else:
            if len(input_data["reactants"]) == 0 or len(input_data["products"]) == 0:
                st.warning("⚠️ Please add both Reagents and Products first!")

st.divider()
st.subheader("🌿 Green Metrics: Atom Economy")
    
input_data = {"reactants": st.session_state.get("reag_list", []),"products": st.session_state.get("prod_list", [])}

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
            <span>Low Efficiency</span>
            <span>Ideal Synthesis</span>
        </div>
    </div>
    """
    st.markdown(gauge_html, unsafe_allow_html=True)

def st_skulls(num_skulls=20): #Definition of a function to pop skulls if the economy atom of the reaction is below 40
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


if st.button("🌿 Green Metrics: Atom Economy"):
    if len(input_data["reactants"]) > 0 and len(input_data["products"]) > 0:
        try:
            
            ae_result = calculate_eco_atm_M(input_data)
            
            
            display_linear_gauge(ae_result, "Atom Economy")
            
            
            if ae_result > 90: #Congratulation message
                st.balloons()
                st.success("Excellent! This reaction is extremely atom-efficient.")
            elif ae_result < 40: #Warning message
                st_skulls()
                st.warning("Warning: this reaction produces too much waste.")
                
        except Exception as e:
            st.error(f"Error during the calculation: {e}")
    else:
        st.error("Please, add a reactant and a product.")






