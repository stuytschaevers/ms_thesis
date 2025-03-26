#this script separates the clusters from SwissDock into individual files
#this script is run on the terminal

import os
import re

# Function to parse the PDB file and extract atom coordinates and cluster information
def parse_pdb(pdb_file):
    clusters = {}
    current_cluster = None
    
    with open(pdb_file, 'r') as file:
        for line in file:
            # Match REMARK lines
            if line.startswith('REMARK'):
                remark_match = re.match(r'REMARK\s+Cluster:\s+(\d+)', line)
                if remark_match:
                    current_cluster = int(remark_match.group(1))
                    clusters[current_cluster] = []
            # Match ATOM lines
            elif line.startswith('ATOM'):
                if current_cluster is not None:
                    clusters[current_cluster].append(line)
                    
    return clusters

# Function to create a new subfolder if it doesn't exist
def create_subfolder(folder_name):
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

# Function to write each cluster to a separate PDB file in a subfolder
def write_clusters_to_pdb(clusters, base_folder, base_filename, ligand_name):
    for cluster_id, cluster in clusters.items():
        cluster_filename = os.path.join(base_folder, f'{ligand_name}_{base_filename}_{cluster_id}.pdb')
        with open(cluster_filename, 'w') as file:
            for atom_line in cluster:
                file.write(atom_line)

def main():
    # Ask user for input
    pdb_file = input("Enter the path to the PDB file containing clusters: ") #  /Users/stuytschaevers/Desktop/Thesis/SwissDock/old/atp_new/TrwD_new_ATP/clusters.dock4.pdb
    output_folder_name = input("Enter the name for the output folder: ") #  /Users/stuytschaevers/Desktop/Thesis/SwissDock/old/atp_new/clusters_atp
    ligand_name = input("Enter the ligand name: ") #    atp

    # Validate input folder name
    if not output_folder_name:
        print("Error: Output folder name cannot be empty.")
        return

    # Create output folder
    output_folder = os.path.join(os.getcwd(), output_folder_name)
    create_subfolder(output_folder)

    # Base filename for output PDB files
    output_base_filename = "cluster"

    # Parse PDB file and write clusters to output folder
    clusters = parse_pdb(pdb_file)
    write_clusters_to_pdb(clusters, output_folder, output_base_filename, ligand_name)

    print(f"Clusters separated and saved in '{output_folder}'.")

if __name__ == "__main__":
    main()
