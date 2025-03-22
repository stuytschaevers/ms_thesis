# This script prepares the homologs of TrwD
# It takes in the heximeric protein and saves the first subunit, A
# It also removes any ions and solvents
import os
from chimerax.core.commands import run

# Directory where the script is located
working_directory = os.getcwd()

# New directory for cleaned files
output_directory = os.path.join(working_directory, "cleaned_files")
os.makedirs(output_directory, exist_ok=True)  # Ensure directory exists

pdb_ids = ["1nlz", "2gza"]

# Process each PDB file
for pdb_id in pdb_ids:
    # Open PDB file
    run(session, f"open {pdb_id} fromDatabase pdb format mmcif")

    # Select solvent and delete
    run(session, "select solvent")
    run(session, "delete sel")

    # Select hydrogen atoms and delete
    run(session, "select H")
    run(session, "delete sel")

    # Select ions and delete
    run(session, "select ions")
    run(session, "delete sel")

    # Select everything except for the A chain and delete it
    run(session, "select /A")
    run(session, "select ~sel")
    run(session, "delete sel")

    # Save cleaned structure to a new PDB file in the output directory
    clean_pdb_file = os.path.join(output_directory, f"{pdb_id}_monomer.pdb")
    run(session, f"save {clean_pdb_file}")

    # Close the current model
    run(session, "close all")
