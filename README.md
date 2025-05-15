# ThermoPred: Reaction thermodynamics Predictor
ThermoPred is an interactive web application that allows chemists to predict chemical reactions and their thermodynamic favorability. By simply drawing two reactant molecules, the application:
1. Predicts possible reaction products using template-based reaction prediction
2. Calculates reaction energies using molecular force fields
3. Determines if the reaction is thermodynamically favorable at 0K
4. Provides 3D visualization of all molecular structures

## Features
* __User-friendly drawing interface__ powered by Ketcher for creating molecular structures
* __Template-based reaction prediction__ supports common organic reactions including: 
  * Nucleophilic substitution (SN2)
  * Amide formation 
  * Ester formation
  * Carbon-carbon coupling reactions
  * Halogen exchange reactions
* __Real-time 3D visualization__ of reactants and products
* __Energy calculations__ using RDKit force fileds for thermodynamic analysis
* __Simple interpretation__ of reaction favorability with clear visual indicators

## Installation
### Prerequisites
* Python 3.8+
* RDKit
* xTB(Semienpirical Extended Tight-Binding Program)
* Streamlit, streamlit_ketcher, stmol, py3Dmol
* numpy
* pandas

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
pip install streamlit py3Dmol stmol streamlit-ketcher py3Dmol 
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
### Starting the Web App
```bash
cd path/to/ThermoPred
streamlit run app.py
```
1. __Draw your first reactant__ in the Molecule 1 sketcher
2. __Draw your second reactant__ in the Molecule 2 sketcher
3. __Select a reaction template__ or use "Auto-detect" to let the system find a matching template
4. __Click "Predict Product"__ to see the predicted reaction outcome
5. __View the 3D structures__ by expanding the 3D visualization panels
6. __Check thermodynamic analysis__ to see if the reaction is favorable at 0K

Tha application will be available at http://localhost:8501 by default

## Technical Details
### Reaction Prediction
ThermoPred uses a template-based approach for reaction prediction, similar to how any retrosynthesis systems work. The templates are defined using SMARTS patters with atom mapping to identify how atoms move from reactants to products.

Key reaction types supported: 
* SN2 reactions (various leaving groups and nucleophiles)
* Amide formation from carboxylic acids and amines
* Ester formation from carboxylic acids and alcohols
* Carbon-carbon coupling reactions (alkyl halides)
* Speail case handling for elemental halogen reactions

## Energy Calculations
the application uses RDKit's force fields(MMFF94) to : 
1. Generate optimized 3D conformations of molecules
2. Calculate energy values for reactants and products
3. Determine the energy difference (ΔE) of the reaction
4. Convert energy values to both Hartree and kcal/mol units


## Troubleshooting
### Common Issues
1. __Invalid SMILES Error__: Ensure you molecular structures are chemically valid
2. __No Templates Matched__: The current reaction may not be supported by existing templates
3. __3D Structure Generation Fails__: Some complex molecules may need simpler representations
4. __Energy Calculation Error__: Very large or complex might cause force field calculation failures
### Extending Reaction Templates
If you need to add support for new reaction types, you can extend the `REACTION_TEMPLATES` dictionary in `reaction_utils.py` with additional SMARTS patterns.

## License ?

