import re
import pandas as pd

# Ask the user for input log file and output CSV file
log_file = input("Enter the path to the ChimeraX log file: ")
output_csv = input("Enter the path for the output CSV file: ")

# Read the log file
with open(log_file, "r") as f:
    log_content = f.readlines()

# Regex patterns to extract relevant data
match_pattern = re.compile(r"Matchmaker (\S+), chain \S+ .* with (\S+), chain \S+.* sequence alignment score = (\d+\.\d+)")
rmsd_pattern = re.compile(r"RMSD between (\d+) pruned atom pairs is (\d+\.\d+) angstroms; \(across all (\d+) pairs: (\d+\.\d+)\)")

# Extracted data
matches = []
current_match = None

for line in log_content:
    match = match_pattern.search(line)
    if match:
        current_match = {
            "Structure 1": match.group(1),
            "Structure 2": match.group(2),
            "Alignment Score": float(match.group(3)),
            "RMSD (pruned pairs)": None,
            "Pruned Atom Pairs": None,
            "RMSD (all pairs)": None
        }
        matches.append(current_match)
    
    rmsd_match = rmsd_pattern.search(line)
    if rmsd_match and current_match:
        current_match["Pruned Atom Pairs"] = int(rmsd_match.group(1))
        current_match["RMSD (pruned pairs)"] = float(rmsd_match.group(2))
        current_match["RMSD (all pairs)"] = float(rmsd_match.group(4))

# Convert to DataFrame
df = pd.DataFrame(matches)

# Save the table to a CSV file
df.to_csv(output_csv, index=False)

# Print the table
print(df)
