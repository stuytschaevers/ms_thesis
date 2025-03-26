# This script compares all models in a folder and saves the log for further analysis

#from chimerax.core.commands import run
#run(session, "command on commandline")

from chimerax.core.commands import run
import os

# Define the working directory
working_directory = "/Users/stuytschaevers/Desktop/Thesis/msa/homologs/all_models"
log_directory = "/Users/stuytschaevers/Desktop/Thesis/msa/homologs/matchmaker_output"

# Ensure the log directory exists
os.makedirs(log_directory, exist_ok=True)

# Get a list of all PDB files in the directory
pdb_files = [file for file in os.listdir(working_directory) if file.endswith(".pdb")]

print(pdb_files)

# Load each PDB file into ChimeraX
for pdb_file in pdb_files:
    full_path = os.path.join(working_directory, pdb_file)  # Ensure full file path
    run(session, f"open {full_path}")

# Compare each model to all other models
for i in range(1, len(pdb_files) + 1):
    for j in range(i + 1, len(pdb_files) + 1):
        run(session, f"match #{i} to #{j}")
        
# Save the log file
log_file_path = os.path.join(log_directory, "file.html")
run(session, f"log save {log_file_path}")

print(f"Comparison log saved to: {log_file_path}")
