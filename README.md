
<p align="center">
  <img src="https://github.com/Ti-m-T/EcoChem/blob/main/EcoChem_loc.jpg" alt="EcoChem Logo" width="500">
</p>
<p align="center">
  <img src="https://img.shields.io/badge/coverage-76%25-green" alt="Coverage">&nbsp;
  <img src="https://img.shields.io/badge/Python-3.10-blue?logo=python&logoColor=white" alt="Python">&nbsp;
  <img src="https://img.shields.io/badge/Contributors-4-orange" alt="Contributors">&nbsp;
  <img src="https://img.shields.io/badge/License-MIT-red" alt="License">
</p>
EcoChem is a tool useful to calculate how green a reaction is by computing the atom economy, PMI and E-factor of the given reaction. It is a python based package that uses RDKIT, Streamlit, numpy and pandas.

---
## 📄*About the package*
With EcoChem it is possible to study how environmentally friendly a reaction is. It all starts by defining a reaction in terms of reagents, products, solvent/catalysts and extractants. Firstly, the SMILES of the participants are obtained: the reagents, solvent/catalyst and products are drawn while the SMILES of the extractants are directly searched on PubChem by typing the name of the compound. 

Once the SMILES are obtained the reaction is equilibrated and can be visualized. 

Finally, the atom economy, PMI and E-factor are calculated and the reaction is classified according to the three factors.

Because a reaction should be eco-friendly but also health-friendly, the toxicity of each participant is displayed to have a complete analysis of the safety of the reaction.

EcoChem is a package perfect for a first-hand chemist as well as an experienced chemist, as it permits to be more self-aware of what we do in the lab.

---
## 💻*Installation and set-up*
The package can be installed in different ways.
Install the package by running on the terminal
```bash
pip install ecochem
```
Install the package by using the URL of the GitHub repository of the project by running on the terminal
```bash
pip install git+https://github.com/Ti-m-T/EcoChem
```
You can also copy the repository and install the project in editable mode by following the following steps
```bash
git clone https://github.com/Ti-m-T/EcoChem
cd EcoChem
pip install -e .
```
Create a separate environement to use the package to prevent bugs
```bash
conda create -n environement_name python=3.10
conda activate
```
## 📦*Required packages*
EcoChem uses the following packages, make sure to have them
```bash
rdkit
pandas
numpy
matplotlib
requests
os
re
plotly
PIL
```
If you are not sure which package you have or not write on your terminal
```bash
conda list 
```
If a package is missing you can install it by running
```bash
pip install nam_of_the_package
```
## ⚡*Structure*
* Draw the reactans, solvent/catalyst and the products and click apply to visualize the smiles 
* Add extractants by writing the IUPAC name or directly the smiles and select the volume used
* Equilibrate the reaction and visualize it
* Calculate the atom economy, PMI and E-factor and visualize how ecofriendly the reaction is
* Obtain information about the toxicity of the compounds used

## *🖊️Autors*
This package was made by four chemistry students at EPFL:
* Sami Meghezzi [![GitHub](https://img.shields.io/badge/GitHub-Sami--Elevated-181717?style=flat&logo=github&logoColor=white)](https://github.com/Sami-Elevated)
* Tim Tellier [![GitHub](https://img.shields.io/badge/GitHub-Ti--m--T-181717?style=flat&logo=github&logoColor=white)](https://github.com/Ti-m-T)
* Julien Tchaplyguine [![GitHub](https://img.shields.io/badge/GitHub-J--Tchp-181717?style=flat&logo=github&logoColor=white)](https://github.com/J-Tchp)
* Claudia Vittorangeli [![GitHub](https://img.shields.io/badge/GitHub-ClauVitt-181717?style=flat&logo=github&logoColor=white)](https://github.com/ClauVitt)





