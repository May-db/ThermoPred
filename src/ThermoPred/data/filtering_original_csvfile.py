import pandas as pd
from pathlib import Path



csv_file= Path("C:/Users/Maria/Desktop/projectppchem/ThermoPred/src/ThermoPred/data/rxn_data_reduced.csv")
output=Path("C:/Users/Maria/Desktop/projectppchem/ThermoPred/src/ThermoPred/data/rxn_data.csv")

df = pd.read_csv(csv_file)


# Optionally filter by class
common_classes = ["SN2", "E1", "E2", "Addition", "Aldol"]
df_filtered = df[df["CLASS"].isin(common_classes)]

# Drop rows with badly formatted REACTION strings
df_filtered = df_filtered[df_filtered["REACTION"].str.contains(">>")]

# Save to new CSV
df_filtered.to_csv(output, index=False)
