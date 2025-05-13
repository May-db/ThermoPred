# ThermoPred: Reaction thermodynamics Predictor
ThermoPred is an interactive Web App that allows chemists to predict chemical reactions and their thermodynamic favorability. By simply drawing two reactant molecules, the application:
1. Predicts possible reaction products
2. Calculates reaction energies using quantum mechanical methods
3. Determines if the reaction is thermodynamically favorable at 0K
4. Provides 3D visualization of all molecular structures

## Installation
## Prerequisites
* Python 3.8+
* RDKit
* xTB(Semienpirical Extended Tight-Binding Program)
* Streamlit, streamlit_ketcher, stmol, py3Dmol, mols2grid, scipy
* rxn-utils
* numpy
* pandas
* Morfeus
* Aizynthtrain

## Step 1: Environment Setup
```bash
#create and activate a virtual environment
python -m venv thermopred_env
# Activate it (Windows)
thermopred_env\Scripts\activate
# Activate  it (macOS/Linus)
source thermopred_env/bin/activate
```
## Install Dependencies
```bash
#Install RDKit
conda install -c conda-forge rdkit
# Install other Python dependencies
pip install streamlit py3Dmol morfeus rxnutils stmol streamlit-ketcher
#and the rest
```
## Step 3: Install xTB
Follow the [official installation instructions](https://xtb-python.readthedocs.io/en/latest/installation.html#conda-forge) for xTB.
Make sure the xtb executable is in your system PATH.
## Step 4: Clone and Install ThermoPred
```bash
#Clone the repository
git clone https://github.com/May-db/ThermoPred.git
cd ThermoPred
# Install the package
pip install -e .
```
## Usage
## Starting the Web App
```bash
cd path/to/ThermoPred
streamlit run app.py
```
Tha application will be available at http://localhost:8501 by default

## Troubleshooting
## Common Issues
1. __xTB not Found__: Ensure xTB is properly installed and added to your PATH
2. __RDKit Import Error__: Make sure RDKit is installed in your active environment
3. __3D Structure Generation Fails__: Some complex molecules may need manual adjustment
## Error Messages

* xTB(Semienpirical Extended Tight-Binding Program)
* Streamlit, streamlit_ketcher, stmol, py3Dmol, mols2grid, scipy
* rxn-utils
* numpy
* pandas
* Morfeus
* Aizynthtrain

## Step 1: Environment Setup
```bash
#create and activate a virtual environment
python -m venv thermopred_env
# Activate it (Windows)
thermopred_env\Scripts\activate
# Activate  it (macOS/Linus)
source thermopred_env/bin/activate
```
## Install Dependencies
```bash
#Install RDKit
conda install -c conda-forge rdkit
# Install other Python dependencies
pip install streamlit py3Dmol morfeus rxnutils stmol streamlit-ketcher
#and the rest
```
## Step 3: Install xTB
Follow the [official installation instructions](https://xtb-python.readthedocs.io/en/latest/installation.html#conda-forge) for xTB.
Make sure the xtb executable is in your system PATH.
## Step 4: Clone and Install ThermoPred
```bash
#Clone the repository
git clone https://github.com/May-db/ThermoPred.git
cd ThermoPred
# Install the package
pip install -e .
```
## Usage
## Starting the Web App
```bash
cd path/to/ThermoPred
streamlit run app.py
```
Tha application will be available at http://localhost:8501 by default

## Troubleshooting
## Common Issues
1. __xTB not Found__: Ensure xTB is properly installed and added to your PATH
2. __RDKit Import Error__: Make sure RDKit is installed in your active environment
3. __3D Structure Generation Fails__: Some complex molecules may need manual adjustment
## Error Messages
