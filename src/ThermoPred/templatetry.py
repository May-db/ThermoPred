



from rdkit import Chem
import pandas as pd
from pathlib import Path
from rxnutils.chem.reaction import ChemicalReaction

csv_file = Path("C:/Users/Maria/Desktop/projectppchem/ThermoPred/src/ThermoPred/rxn_data_reduced.csv")
def get_product(reactant1, reactant2):
    
    try:
        df = pd.read_csv(csv_file)
        print(f"Loaded {len(df)} reactions")
    except Exception as e:
        print(f"Error loading CSV: {e}")
        return None
    
    input_reactants = sorted([reactant1, reactant2])
    
    for _, row in df.iterrows():
        reaction_smiles = row['REACTION']
        if '>>' not in reaction_smiles:
            continue
            
        reactants_str, product_str = reaction_smiles.split('>>')
        
        db_reactants = sorted(reactants_str.split('.'))
        
        if db_reactants == input_reactants:
            print(f"Found matching reaction")
            
            #if 'CLASS' in row:
             #   print(f"Reaction class: {row['CLASS']}")
            
            mapped_reaction = row['MAPPED_REACTION']
            
            try:
                rxn = ChemicalReaction(mapped_reaction)
                print(f"Reactants from mapped reaction: {rxn.reactants_list}")
                print(f"Products from mapped reaction: {rxn.products_list}")
                rxn.generate_reaction_template(radius=1)
                print(f"Product SMILES: {product_str}")
                return product_str
            
            except Exception as e:
                print(f"Error processing with rxnutils: {e}")
                return product_str
    
    print("No matching reaction found")
    return None


#product = get_product(r1, r2)
    
   # if product:
    #    print(f"\nFinal product: {product}")
    #else:
      #  print("\nNo product found")