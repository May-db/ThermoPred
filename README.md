# ThermoPred: Reaction thermodynamics Predictor
ThermoPred is an interactive web application that allows chemists to predict chemical reactions and their thermodynamic favorability. By simply drawing two reactant molecules, users can explore reaction outcomes and energetics.

## Features
* __Reaction Prediction__ : Predicts products for common organic reactions using template-based matching
* __Thermodynamic Analysis__ : Calculates reaction energies tp determine if reactions are favorable at 0K
* __3D Visualization__ : Interactive 3D models of reactants and products
* __Chemical Drawing Interface__ : User-friendly molecule sketcher powered by Ketcher
Supported reaction types include:
  - Nucleophilic substitution (SN2)
  - Amide formation 
  - Ester formation
  - Carbon-carbon coupling 
  - Halogen exchange 
> the reaction that works best is the nucleophilic substition with halogens

## Installation
### Prerequisites
* Python 3.8+
* RDKit
* Streamlit, streamlit_ketcher, stmol, py3Dmol
* numpy
* pandas
* os
* pyarrow
* pathlib

## Step 1: Environment Setup
```bash
#create and activate a virtual environment
python -m venv thermopred_env
# Activate it (Windows)
thermopred_env\Scripts\activate
# Activate  it (macOS/Linus)
source thermopred_env/bin/activate
```
## Step 2: Install Dependencies
```bash
#Install RDKit
conda install -c conda-forge rdkit
# Install other Python dependencies
pip install streamlit py3Dmol stmol streamlit-ketcher py3Dmol 
#and the rest
```
## Step 3: Clone and Install ThermoPred
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
1. Draw you reactants molecules
2. Choose a reaction template or use "Auto-detect"
3. Click "Predict Product" to see the result and energy analysis

Tha application will be available at http://localhost:8501 by default

## Troubleshooting
### Common Issues
1. __Invalid SMILES Error__: Ensure you molecular structures are chemically valid
2. __No Templates Matched__: The current reaction may not be supported by existing templates
3. __3D Structure Generation Fails__: Some complex molecules may need simpler representations
4. __Energy Calculation Error__: Very large or complex might cause force field calculation failures
### Extending Reaction Templates
If you need to add support for new reaction types, you can extend the `REACTION_TEMPLATES` dictionary in `reaction_utils.py` with additional SMARTS patterns.

## License 
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

