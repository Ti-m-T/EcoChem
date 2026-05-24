import sys
import os
import random
import re
import requests
import streamlit as st
import plotly.graph_objects as go
from altair import value
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdMolDescriptors
from thermo import chemical as density_finder
from dataclasses import dataclass, field
from streamlit_ketcher import st_ketcher
import uuid
from DataClasses import Reaction
from DataClasses import Chemical
from DataClasses import ChemswithMass
from DataClasses import Solvent
from DataClasses import Extractant
from DataClasses import LiquidChemical


def make_id():
    return uuid.uuid4().hex[:8]


st.set_page_config(page_title="EcoChem", page_icon=":leaves:")


# -----------------------------
# Session State Initialization
# -----------------------------
if "page_active" not in st.session_state:
    st.session_state.page_active = "Home"

if "reag_list" not in st.session_state:
    st.session_state.reag_list = []  # [{"id":..., "smiles":...}]

if "prod_list" not in st.session_state:
    st.session_state.prod_list = []  # [{"id":..., "smiles":...}]

if "solvcat_list" not in st.session_state:
    st.session_state.solvcat_list = (
        []
    )  # [{"id":..., "role":..., "smiles":..., "volume":..., "density":..., "mass":...}]

if "extr_list" not in st.session_state:
    st.session_state.extr_list = []

if "wanted_product_mass" not in st.session_state:
    st.session_state.wanted_product_mass = 1.0


# -----------------------------
# Utility Functions
# -----------------------------
def convert_to_latex_subscripts(formula: str) -> str:
    return re.sub(r"(\d+)", r"_{\1}", formula)


def name_to_smiles(name: str) -> str:
    try:
        chemical_name = name.strip()
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{chemical_name}/property/CanonicalSMILES/TXT"
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            return response.text.strip()
        return None
    except Exception:
        return None


def go_to(page_name):
    st.session_state.page_active = page_name


def get_formula(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            formula = rdMolDescriptors.CalcMolFormula(mol)
            return re.sub(r"(\d+)", r"_{\1}", formula)
        return smiles
    except Exception:
        return smiles


def atom_set_builder(smiles_list: list[str]) -> set:
    atom_set = set()
    for s in smiles_list:
        mol = Chemical(smiles=s)
        molecular_formula = mol.mol_f
        atoms = re.findall(r"[A-Z][a-z]?", molecular_formula)
        atom_set.update(atoms)
    return atom_set


# -----------------------------
# Build Reaction Object
# -----------------------------
def build_reaction():
    reag_smiles = [
        item["smiles"] for item in st.session_state.reag_list if item["smiles"]
    ]
    prod_smiles = [
        item["smiles"] for item in st.session_state.prod_list if item["smiles"]
    ]

    reactants_objects = [Chemical(smiles=s) for s in reag_smiles]

    wanted_product_obj = ChemswithMass(
        smiles=prod_smiles[0] if prod_smiles else "",
        initial_mass=st.session_state.wanted_product_mass,
    )

    byproducts_objects = [Chemical(smiles=s) for s in prod_smiles[1:]]

    solvents = []
    catalysts = []

    for item in st.session_state.solvcat_list:
        if item["role"] == "solvent":
            solvents.append(
                Solvent(
                    smiles=item["smiles"],
                    volume=item.get("volume", 0.0),
                    user_density=item.get("density", 0.0),
                )
            )
        else:
            catalysts.append(
                ChemswithMass(
                    smiles=item["smiles"],
                    initial_mass=item.get("mass", 0.0),
                )
            )

    return Reaction(
        reactants=reactants_objects,
        wanted_product=wanted_product_obj,
        byproducts=byproducts_objects,
        Catalysts=catalysts,
        solvents=solvents,
        extractants=[
            Extractant(
                smiles=e["smiles"],
                volume=e["volume"],
                user_density=e["density"],
            )
            for e in st.session_state.extr_list
        ],
        Chosen_Yield=yield_fraction,
    )


# -----------------------------
# Reaction Latex Formatter
# -----------------------------
def format_reaction_latex(experiment):
    def fmt(c, f):
        return rf"{int(round(c))}\,{convert_to_latex_subscripts(f)}"

    reactants = " + ".join(fmt(r.coeff, r.mol_f) for r in experiment.reactants)
    products = " + ".join(
        fmt(p.coeff, p.mol_f)
        for p in experiment.byproducts + [experiment.wanted_product]
    )
    return rf"{reactants} \;\longrightarrow\; {products}"


# -----------------------------
# Gauges and Skulls
# -----------------------------
def show_skulls(num_skulls=20):
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


def display_linear_gauge(value, title="Atom Economy"):
    pos = max(0, min(100, value))
    if pos < 50:
        color = "#ff4b4b"
    elif pos < 80:
        color = "#ffa500"
    else:
        color = "#00cc96"

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


def display_linear_gauge_pmi(value, title="PMI"):
    pos = max(0, min(100, value))
    if pos > 80:
        color = "#ff4b4b"
    elif pos > 40:
        color = "#ffa500"
    else:
        color = "#00cc96"

    gauge_html = f"""
    <div style="font-family: sans-serif; margin: 20px 0; padding: 10px; background-color: #f9f9f9; border-radius: 10px;">
        <div style="display: flex; justify-content: space-between; margin-bottom: 15px;">
            <span style="font-weight: bold; font-size: 1.1rem; color: #333;">{title}</span>
            <span style="font-weight: bold; font-size: 1.3rem; color: {color};">{value:.1f}</span>
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


def display_linear_gauge_efactor(value, title="E-Factor"):
    pos = max(0, min(100, value))
    if pos > 80:
        color = "#ff4b4b"
    elif pos > 40:
        color = "#ffa500"
    else:
        color = "#00cc96"

    gauge_html = f"""
    <div style="font-family: sans-serif; margin: 20px 0; padding: 10px; background-color: #f9f9f9; border-radius: 10px;">
        <div style="display: flex; justify-content: space-between; margin-bottom: 15px;">
            <span style="font-weight: bold; font-size: 1.1rem; color: #333;">{title}</span>
            <span style="font-weight: bold; font-size: 1.3rem; color: {color};">{value:.1f}</span>
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
        <div style="display: flex; justify-content: space-between; font-size: 0.75rem; margin-top: 8px; color: #888%;">
            <span>Low E-Factor</span>
            <span>High E-Factor</span>
        </div>
    </div>
    """
    st.markdown(gauge_html, unsafe_allow_html=True)


# -----------------------------
# Toxcicity functions
# -----------------------------
def render_chemicals(title, chem_list):
    st.markdown(f"#### {title}")

    for chem in chem_list:
        chem.get_CID()

        if chem.CID is None:
            st.warning(f"No data was found for {chem.smiles}.")
            continue

        chem.get_GHS()
        chem.get_pictograms()

        list_picto = list(chem.pictograms)
        list_GHS = list(chem.GHS.keys())

        if not list_GHS:
            text_GHS = "No hazard statements available"
        elif len(list_GHS) == 1:
            text_GHS = f'"{list_GHS[0]}"'
        else:
            text_GHS = (
                ", ".join(f'"{x}"' for x in list_GHS[:-1]) + f' and "{list_GHS[-1]}"'
            )

        st.warning(f"⚠️ [{chem.smiles}] hazards: {text_GHS}")

        for i in range(0, len(list_picto), 3):
            row = list_picto[i : i + 3]
            cols = st.columns(3)

            for col, picto in zip(cols, row):
                img_path = GHS_pictograms.get(picto)
                if img_path:
                    col.image(img_path, width=100)


# -----------------------------
# UI: Sidebar Navigation
# -----------------------------
st.sidebar.title("Menu")

if st.sidebar.button("🏠 Home"):
    go_to("Home")

if st.sidebar.button("⚛️ Reaction Builder"):
    go_to("Reaction Builder")

if st.sidebar.button("🧪 Compute"):
    go_to("Compute")


# -----------------------------
# HOME PAGE
# -----------------------------
if st.session_state.page_active == "Home":
    st.title("🍃 EcoChem: Green Chemistry Calculator")
    st.markdown("""
    Welcome to **EcoChem**, your interactive assistant for designing reactions and calculating green chemistry metrics. 
    
    ### 🚀 Getting Started
    
    1. **⚛️ Reaction Builder**: Head over to the builder page to input or sketch your reactants, products, solvents, and catalysts using the interactive sketcher.
    2. **🧪 Compute**: Set your target parameters (like desired mass and expected yield), add extractions for your work-up, and compute your green metrics.
    
    ### 📊 Metrics Evaluated
    * **Atom Economy (AE):** Measures how many atoms from your starting materials end up in your desired product.
    * **Process Mass Intensity (PMI):** Evaluates the total mass of materials used per gram of product.
    * **E-Factor:** Quantifies the exact waste-to-product ratio of your process.
    """)

    st.info(
        "👈 Use the sidebar menu to navigate to the **Reaction Builder** and get started!"
    )


# -----------------------------
# REACTION BUILDER PAGE
# -----------------------------
elif st.session_state.page_active == "Reaction Builder":

    st.title("⚛️ Reaction Builder")

    # -----------------------------
    # Reagents Section
    # -----------------------------
    st.header("Reactants")

    if st.button("➕ Add Reagent"):
        st.session_state.reag_list.append({"id": make_id(), "smiles": ""})

    for i, item in enumerate(st.session_state.reag_list):
        rid = item["id"]
        with st.expander(f"Reagent {i+1}", expanded=True):

            typed = st.text_input(
                f"SMILES for Reagent {i+1}",
                key=f"typed_reag_{rid}",
                value=item["smiles"],
            )
            drawn = st_ketcher(key=f"drawn_reag_{rid}")

            if typed.strip():
                item["smiles"] = typed.strip()
            elif drawn and drawn.strip():
                item["smiles"] = drawn.strip()

            if item["smiles"]:
                st.success(f"Reagent {i+1} added: `{item['smiles']}`")

            if st.button("❌ Remove", key=f"remove_reag_{rid}"):
                st.session_state.reag_list.pop(i)
                st.session_state.pop(f"typed_reag_{rid}", None)
                st.session_state.pop(f"drawn_reag_{rid}", None)
                st.rerun()

    st.divider()

    # -----------------------------
    # Solvents & Catalysts Section
    # -----------------------------
    st.header("Solvents & Catalysts")

    if st.button("➕ Add Solvent/Catalyst"):
        st.session_state.solvcat_list.append(
            {
                "id": make_id(),
                "role": "solvent",  # default
                "smiles": "",
                "volume": 0.0,
                "density": 0.0,
                "mass": 0.0,
            }
        )

    for i, item in enumerate(st.session_state.solvcat_list):
        sid = item["id"]
        with st.expander(f"Item {i+1}", expanded=True):

            role = st.radio(
                "Role",
                ["solvent", "catalyst"],
                key=f"role_{sid}",
                index=0 if item["role"] == "solvent" else 1,
            )
            item["role"] = role

            typed = st.text_input(
                f"SMILES for {role.capitalize()} {i+1}",
                key=f"typed_solvcat_{sid}",
                value=item["smiles"],
            )
            drawn = st_ketcher(key=f"drawn_solvcat_{sid}")

            if typed.strip():
                item["smiles"] = typed.strip()
            elif drawn and drawn.strip():
                item["smiles"] = drawn.strip()

            if role == "solvent":
                item["volume"] = st.number_input(
                    "Volume (mL)",
                    min_value=0.0,
                    step=1.0,
                    key=f"vol_{sid}",
                    value=float(item["volume"]),
                )
                item["density"] = st.number_input(
                    "Density (g/mL)",
                    min_value=0.0,
                    step=0.1,
                    key=f"dens_{sid}",
                    value=float(item["density"]),
                )
            else:
                item["mass"] = st.number_input(
                    "Mass (g)",
                    min_value=0.0,
                    step=0.1,
                    key=f"mass_{sid}",
                    value=float(item["mass"]),
                )
            if item["smiles"]:
                st.success(f"{role.capitalize()} {i+1} added: `{item['smiles']}`")
            if st.button("❌ Remove", key=f"remove_solvcat_{sid}"):
                st.session_state.solvcat_list.pop(i)
                st.rerun()

    st.divider()

    # -----------------------------
    # Products Section
    # -----------------------------
    st.header("Products")

    if st.button("➕ Add Product"):
        st.session_state.prod_list.append({"id": make_id(), "smiles": ""})

    for i, item in enumerate(st.session_state.prod_list):
        pid = item["id"]
        label = "Main Product" if i == 0 else f"Byproduct {i}"
        with st.expander(label, expanded=True):

            typed = st.text_input(
                f"SMILES for {label}",
                key=f"typed_prod_{pid}",
                value=item["smiles"],
            )
            drawn = st_ketcher(key=f"drawn_prod_{pid}")

            if typed.strip():
                item["smiles"] = typed.strip()
            elif drawn and drawn.strip():
                item["smiles"] = drawn.strip()

            if item["smiles"]:
                st.success(f"{label} added: `{item['smiles']}`")

            if st.button("❌ Remove", key=f"remove_prod_{pid}"):
                st.session_state.prod_list.pop(i)
                st.rerun()

    st.divider()

    # -----------------------------
    # Generate Reaction SMILES
    # -----------------------------
    if st.button("🚀 Generate Reaction SMILES"):

        reag_str = ".".join(
            [item["smiles"] for item in st.session_state.reag_list if item["smiles"]]
        )
        solv_str = ".".join(
            [item["smiles"] for item in st.session_state.solvcat_list if item["smiles"]]
        )
        prod_str = ".".join(
            [item["smiles"] for item in st.session_state.prod_list if item["smiles"]]
        )

        full_reaction = f"{reag_str}>{solv_str}>{prod_str}"

        if not reag_str:
            st.warning("⚠️ Please add at least one reactant.")
            st.stop()

        if not prod_str:
            st.warning("⚠️ Please add at least one product.")
            st.stop()

        for item in st.session_state.reag_list:
            if "." in item["smiles"]:
                st.warning("⚠️ A reactant contains multiple SMILES.")
                st.stop()

        for item in st.session_state.prod_list:
            if "." in item["smiles"]:
                st.warning("⚠️ A product contains multiple SMILES.")
                st.stop()

        global yield_fraction
        yield_fraction = 1
        experiment = build_reaction()

        reactant_set = {r.smiles for r in experiment.reactants}
        product_set = {
            p.smiles for p in experiment.byproducts + [experiment.wanted_product]
        }

        if reactant_set & product_set:
            st.warning("⚠️ A molecule appears as both reactant and product.")
            st.stop()

        atom_reag = atom_set_builder(
            [item["smiles"] for item in st.session_state.reag_list]
        )
        atom_prod = atom_set_builder(
            [item["smiles"] for item in st.session_state.prod_list]
        )

        if atom_reag != atom_prod:
            diff = atom_reag ^ atom_prod
            st.warning(f"⚠️ Atom mismatch: {diff}")
            st.stop()

        try:
            reag_mols = [
                Chem.MolFromSmiles(item["smiles"])
                for item in st.session_state.reag_list
            ]
            solv_mols = [
                Chem.MolFromSmiles(item["smiles"])
                for item in st.session_state.solvcat_list
            ]
            prod_mols = [
                Chem.MolFromSmiles(item["smiles"])
                for item in st.session_state.prod_list
            ]

            rxn = AllChem.ChemicalReaction()

            for m in reag_mols:
                if m:
                    rxn.AddReactantTemplate(m)

            for m in solv_mols:
                if m:
                    rxn.AddAgentTemplate(m)

            for m in prod_mols:
                if m:
                    rxn.AddProductTemplate(m)

            img = Draw.ReactionToImage(rxn, subImgSize=(400, 400), useSVG=False)
            st.image(img, use_container_width=True)

        except Exception as e:
            st.error(f"Could not render reaction image: {e}")

        st.success("Reaction SMILES generated!")
        st.write(f"**Reaction SMILES:** `{full_reaction}`")
        st.session_state.final_smiles = full_reaction


# -----------------------------
# COMPUTE PAGE
# -----------------------------
elif st.session_state.page_active == "Compute":

    st.title("🧪 Compute Reaction Metrics")

    # -----------------------------
    # Extractants Section
    # -----------------------------
    st.subheader("Add Extractants (Work-up)")

    extr_input = st.text_input(
        "Extractant (Name or SMILES)",
        key="input_extr_name",
        placeholder="e.g., ethyl acetate, hexane, etc...",
    )
    extr_volume = st.number_input(
        "Volume (mL)", min_value=0.0, step=10.0, key="input_extr_volume"
    )
    extr_density = st.number_input(
        "Density (g/mL)", min_value=0.0, step=0.1, key="input_extr_density"
    )

    if st.button("➕ Add Extractant"):
        if extr_input:
            calculated_smiles = name_to_smiles(extr_input)
            if not calculated_smiles:
                if Chem.MolFromSmiles(extr_input):
                    calculated_smiles = extr_input

            if calculated_smiles:
                st.session_state.extr_list.append(
                    {
                        "smiles": calculated_smiles,
                        "volume": extr_volume,
                        "density": extr_density,
                    }
                )
                st.success(f"Added: `{calculated_smiles}`")
                st.rerun()
            else:
                st.error(f"Could not resolve '{extr_input}'.")
        else:
            st.error("Please enter a name or SMILES first.")

    st.markdown("---")

    if st.session_state.extr_list:
        st.write("### Current Extractants:")
        for i, e in enumerate(st.session_state.extr_list):
            col_text, col_btn = st.columns([0.8, 0.2])
            with col_text:
                st.markdown(
                    f"**{i+1}.** `{e['smiles']}` — {e['volume']} mL — {e['density']} g/mL"
                )
            with col_btn:
                if st.button("❌ Remove", key=f"remove_extr_{i}"):
                    st.session_state.extr_list.pop(i)
                    st.rerun()
    else:
        st.info("No extractants added yet.")

    st.divider()

    # -----------------------------
    # Reaction Parameters
    # -----------------------------
    st.subheader("Reaction Parameters")

    wanted_product_mass = st.number_input(
        "Desired Mass of Main Product (g)",
        min_value=0.0,
        value=st.session_state.wanted_product_mass,
        step=0.1,
    )
    st.session_state.wanted_product_mass = wanted_product_mass

    yield_input = st.number_input("Reaction Yield (%)", value=100.0)

    if yield_input <= 0 or yield_input > 100:
        st.warning("Yield must be between 0 and 100.")
        st.stop()
    else:
        yield_fraction = yield_input / 100
        st.success(f"Yield set to {yield_input:.1f}%")

    st.divider()

    # -----------------------------
    # Compute Button
    # -----------------------------
    if st.button("⚗️ Compute Reaction Stoichiometry and Green Metrics"):

        if not st.session_state.reag_list:
            st.warning("Please add at least one reactant.")
            st.stop()

        if not st.session_state.prod_list:
            st.warning("Please add at least one product.")
            st.stop()

        for item in st.session_state.reag_list:
            if "." in item["smiles"]:
                st.warning("A reactant contains multiple SMILES.")
                st.stop()

        for item in st.session_state.prod_list:
            if "." in item["smiles"]:
                st.warning("A product contains multiple SMILES.")
                st.stop()

        experiment = build_reaction()

        reactant_set = {r.smiles for r in experiment.reactants}
        product_set = {
            p.smiles for p in experiment.byproducts + [experiment.wanted_product]
        }

        if reactant_set & product_set:
            st.warning("A molecule appears as both reactant and product.")
            st.stop()

        atom_reag = atom_set_builder(
            [item["smiles"] for item in st.session_state.reag_list]
        )
        atom_prod = atom_set_builder(
            [item["smiles"] for item in st.session_state.prod_list]
        )

        if atom_reag != atom_prod:
            diff = atom_reag ^ atom_prod
            st.warning(f"Atom mismatch: {diff}")
            st.stop()

        try:
            reactants_coeff, products_coeff = experiment.stoich_of_reaction()
            st.session_state.stoich_results = {
                "reactants": reactants_coeff,
                "products": products_coeff,
            }

            ae_result = experiment.calcul_eco_atom()
            pmi_result = experiment.PMI()
            ef_result = experiment.e_factor()

            st.session_state.results = {
                "Atom Economy": ae_result,
                "PMI": pmi_result,
                "E-Factor": ef_result,
            }

            st.subheader("⚖️ Balanced Reaction")
            st.latex(format_reaction_latex(experiment))
            st.success("Reaction Balanced!")

            st.subheader("🌿 Green Metrics")

            # -----------------------------
            # Atom Economy Comments
            # -----------------------------
            display_linear_gauge(ae_result, "Atom Economy")

            if ae_result >= 90:
                st.success("Excellent! This reaction is very atom‑efficient.")
            elif ae_result >= 60:
                st.info(
                    "Good atom economy — there is some waste, but overall acceptable."
                )
            elif ae_result >= 40:
                st.warning("Moderate atom economy — significant waste is produced.")
            else:
                st.error(
                    "💀 Very poor atom economy — most atoms do not end up in the product."
                )

            # -----------------------------
            # PMI Comments
            # -----------------------------
            display_linear_gauge_pmi(pmi_result, "PMI")

            if pmi_result <= 10:
                st.success("Fantastic! Very low material usage.")
            elif pmi_result <= 25:
                st.info("Good PMI — material usage is reasonable.")
            elif pmi_result <= 50:
                st.warning("Moderate PMI — the process uses a lot of material.")
            else:
                st.error(
                    "💀 Very high PMI — the process is extremely material‑intensive."
                )

            # -----------------------------
            # E‑Factor Comments
            # -----------------------------
            display_linear_gauge_efactor(ef_result, "E-Factor")

            if ef_result <= 5:
                st.success("Excellent! Very little waste is generated.")
            elif ef_result <= 20:
                st.info("Good E‑Factor — waste generation is acceptable.")
            elif ef_result <= 50:
                st.warning("High E‑Factor — significant waste is generated.")
            else:
                st.error(
                    "💀 Extremely high E‑Factor — this process produces a lot of waste."
                )

            if ae_result >= 90 and ef_result <= 10 and pmi_result <= 10:
                st.balloons()
            elif ae_result < 40 and ef_result >= 50 and pmi_result > 70:
                show_skulls()

            st.subheader("Risk Assessment")
            set_CID: set = set()
            for chem in (
                experiment.reactants
                + experiment.byproducts
                + [experiment.wanted_product]
                + experiment.Catalysts
                + experiment.solvents
                + experiment.extractants
            ):
                chem.get_CID()
                set_CID.add(chem.CID)

            if set_CID == {None}:
                st.warning("No hazard data found.")

            from pathlib import Path

            BASE_DIR = Path(__file__).resolve().parents[2]
            ASSETS_DIR = BASE_DIR / "assets"
            st.write(ASSETS_DIR)
            set_pictograms = set()
            for chem in (
                experiment.reactants
                + experiment.byproducts
                + [experiment.wanted_product]
                + experiment.Catalysts
                + experiment.solvents
                + experiment.extractants
            ):
                chem.get_pictograms()
                for picto in chem.pictograms:
                    set_pictograms.add(picto)

            if not set_pictograms:
                st.success("No hazard pictograms found.")
            else:
                st.warning(
                    "⚠️ This reaction involves hazardous substances. Hazard information is shown below."
                )  # Warning of the reaciton GHS

                GHS_pictograms = {
                    "Exploding bomb": "../assets/GHS_Exploding_bomb.png",
                    "Flame": "../assets/GHS_Flame.png",
                    "Oxidizer (flame over circle)": "../assets/GHS_Oxidizer.png",
                    "Gas cylinder": "../assets/GHS_Gas_cylinder.png",
                    "Corrosion": "../assets/GHS_Corrosion.png",
                    "Skull and crossbones": "../assets/GHS_Skull.png",
                    "Exclamation mark": "../assets/GHS_Exclamation_mark.png",
                    "Health hazard": "../assets/GHS_Health_hazard.png",
                    "Environment": "../assets/GHS_Environment.png",
                }

                pictos = list(set_pictograms)

                for i in range(0, len(pictos), 3):
                    cols = st.columns(3)

                    for col, picto in zip(cols, pictos[i : i + 3]):
                        img_path = GHS_pictograms.get(picto)

                        if img_path:
                            col.image(img_path, width=100)

                        else:
                            col.write(picto)

                st.markdown("---")
                st.subheader("💀 Toxicity of specific species")

                render_chemicals("Reactants toxicity", experiment.reactants)

                render_chemicals(
                    "Products toxicity",
                    [experiment.wanted_product] + experiment.byproducts,
                )

                if experiment.solvents:
                    render_chemicals("Solvent toxicity", experiment.solvents)

                if experiment.Catalysts:
                    render_chemicals("Catalysts toxicity", experiment.Catalysts)

                if experiment.extractants:
                    render_chemicals("Extractants toxicity", experiment.extractants)

                st.info("The analysis of your reaction is done!")
                st.balloons()

        except Exception as e:
            st.error(f"Error during computation details: {e}")
