#!/usr/bin/env python
# encoding: utf-8

#md traj its 1 based and meld is 0 based 

# Importing the required libraries
import numpy as np
from meld.remd import ladder, adaptor, leader
from meld import comm, vault
from meld import system
from meld import parse
import meld as meld
from meld.system.scalers import LinearRamp
from meld.system.amber import *
import glob
from openmm import app as app
from openmm import unit as u
import openmm as omm
import mdtraj as md
import re

# Defining some variables
# Number of replicas
#N_REPLICAS = 30
N_REPLICAS = 2  # This is for testing purposes
# Number of steps (units of exchange period), total step in simulation
#N_STEPS = 20000
N_STEPS = 8000  # This is for testing purposes
# Controls frequency of output
BLOCK_SIZE = 100

# Define a cutoff for RG
RG_CUTOFF = 1.8 # nm

# File in this directory
# TrwD_atp_min_moved_to_origin.pdb, AMBER minimized structure that was moved to the origin
templates = glob.glob('/home/sjt3532/meld_atp/TEMPLATES/*.pdb')       # read the template file, can be multiple
print(f'Found templates: {templates}')

# Function to get max contact distance from .txt file
def get_max_dist(filename):
    # Open the file and read its contents
    with open(filename, "r") as file:  # Replace 'file.txt' with your actual file name
        data = file.read().strip()  # Read and remove any extra whitespace or newline

    # Convert the string to a float
    number = float(data)

    # Print the number
    print("The number in the file is:", number)
    return number

# Get the max contact distance
MAX_CONTACT_DIST = get_max_dist('/home/sjt3532/meld_atp/input_files/max_distance.txt')

# Load the first template to get the protein chain
prot_lig_system = md.load(templates[0])

# Select the protein and ligand chain
PROTEIN_CHAIN = prot_lig_system.topology.select('protein') # 1-based
LIGAND_CHAIN = prot_lig_system.topology.select('resname ATP') #1-based
#Protein chain [   0    1    2 ... 2822 2823 2824]
#Ligand chain [2825 2826 2827 2828 2829 2830 2831 2832 2833 2834 2835 2836 2837 2838 2839 2840 2841 2842 2843 2844 2845 2846 2847 2848 2849 2850 2851 2852 2853 2854 2855 2856 2857 2858 2859 2860 2861 2862 2863 2864 2865 2866 2867]
#PROTEIN_CHAIN = np.arange(1, 2827)  # 2827 because arange stops at (n-1)
#LIGAND_CHAIN = np.arange(2827, 2871)  # 2868 because arange stops at (n-1)
#PROTEIN_CHAIN = prot_lig_system.topology.select('chainid == 0')
#LIGAND_CHAIN = prot_lig_system.topology.select('chainid == 1')

# Residue numbers 0 based
#PROTEIN_CHAIN = prot_lig_system.topology.select('resid 1 to 358')
#LIGAND_CHAIN = prot_lig_system.topology.select('resid 358')

print("PROTEIN_CHAIN:", PROTEIN_CHAIN)
print("LIGAND_CHAIN:", LIGAND_CHAIN)

# Remove numbers < 2827 for ligand chain
#LIGAND_CHAIN = [i for i in LIGAND_CHAIN if i > 2827]

# Add other numbers to protein chain, up to 2826
#PROTEIN_CHAIN = [i for i in range(2827)]

# Print the protein and ligand chain
#print("PROTEIN_CHAIN:", PROTEIN_CHAIN)
#print("LIGAND_CHAIN:", LIGAND_CHAIN)

# Get length of each chain
#PROTEIN_LENGTH = len(set(prot_lig_system.topology.atom(i).residue.index for i in PROTEIN_CHAIN))
#LIGAND_LENGTH = len(set(prot_lig_system.topology.atom(i).residue.index for i in LIGAND_CHAIN))
protein_length = len(PROTEIN_CHAIN)
ligand_length = len(LIGAND_CHAIN)
traj_length = protein_length + ligand_length

# Function to get radius of gyration
def get_rg(traj, atom_chain):
    #print("Atom chain indices:", atom_chain, type(atom_chain))
    rg = md.compute_rg(traj.atom_slice(atom_chain)) * 0.1  # convert Å to nm
    return rg

# Get the RG of the protein and ligand
prot_rg = get_rg(prot_lig_system, PROTEIN_CHAIN)
lig_rg = get_rg(prot_lig_system, LIGAND_CHAIN)

def compute_com(traj, atom_indices):
    masses = np.array([atom.element.mass for atom in traj.topology.atoms])[atom_indices]
    com = np.sum(traj.xyz[:, atom_indices, :] * masses[:, np.newaxis], axis=1) / masses.sum()
    return com

protein_com = compute_com(prot_lig_system, PROTEIN_CHAIN)
sphere_center = protein_com
system_rg = prot_rg + lig_rg
sphere_radius = 2 * system_rg + RG_CUTOFF
# Convert sphere_radius to int
sphere_radius = float(sphere_radius)

# Function to get heavy atoms of ligand
def get_ligand_heavy_atoms_from_file(file_path):
    ligand_heavy_atoms = []
    try:
        with open(file_path, 'r') as file:
            for line in file:
                atoms = line.split()  # Handles both space- and newline-separated values
                ligand_heavy_atoms.extend(atoms)  # Add atoms to the list
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
    return ligand_heavy_atoms

# Get the heavy atoms of the ligand
LIGAND_HEAVY_ATOMS = get_ligand_heavy_atoms_from_file("/home/sjt3532/meld_atp/input_files/unique_atom_names_amber.txt")

# Function to keep ligand within a sphere
# Confinementrestraint needs to be based from origin
# The force constant in a way decides how steep is the potential. the higher the value the steeper the potential is, and the more the ligand will be "confined" to the sphere.
def make_sphere_restraints(s, scaler, ligand_residues, ligand_atoms_dict, sphere_radius, force_const=350):
    """
    Create sphere restraints by matching ligand heavy atoms to their respective residues.

    Parameters:
    - s: MELD system
    - scaler: MELD scaler
    - ligand_residues: List of ligand residue indices
    - ligand_atoms_dict: Dictionary mapping residue indices to their heavy atom names
    - sphere_radius: Restraint radius (nanometers)
    - force_const: Force constant (default 350 kJ/mol/nm²)
    
    Returns:
    - List of MELD restraints
    """
    dist = []  # List to store the restraints
    #heavy_atoms = ligand_atoms_dict
    #print(f"Heavy atoms: {heavy_atoms}")

    """for resid in ligand_residues:
        print(f"Residue {resid}")
        if resid in ligand_atoms_dict:  # Ensure residue exists in the mapping
            print(f"Residue {resid} has atoms: {ligand_atoms_dict[resid]}")
            for atom_name in ligand_atoms_dict[resid]:  # Match atoms to residue
                print(f"Atom {atom_name}")
                try:
                    atom_index = s.index.atom(resid, atom_name)  # Find atom index
                    print(f"Restraint added: Residue {resid}, Atom {atom_name}, Index {atom_index}")
                    
                    rest = meld.system.restraints.ConfinementRestraint(
                        s, scaler, 
                        #LinearRamp(0.0, 100.0, 0.0, 1.0),
                        LinearRamp(0.0, 500.0, 0.0, 1.0),
                        atom_index=atom_index, 
                        radius=sphere_radius * u.nanometer,
                        force_const=force_const * u.kilojoules_per_mole / (u.nanometer * u.nanometer)
                    )
                    dist.append(rest)
                except KeyError:
                    print(f"Warning: Atom {atom_name} not found in residue {resid}")

    return dist"""
    print(f"Ligand residues: {ligand_residues}")
    print(f"Ligand atoms: {ligand_atoms_dict}")

    for num in ligand_residues:
        resid = num - 1
        #print(f"Residue {resid}")
        for atom_name in ligand_atoms_dict:  # Match atoms to residue
            #print(f"Atom {atom_name}")
            try:
                atom_index = s.index.atom(resid, atom_name)  # Find atom index
                print(f"Restraint added: Residue {resid}, Atom {atom_name}, Index {atom_index}")
                
                rest = meld.system.restraints.ConfinementRestraint(
                    s, scaler, 
                    #LinearRamp(0.0, 100.0, 0.0, 1.0),
                    LinearRamp(0.0, 500.0, 0.0, 1.0),
                    atom_index=atom_index, 
                    radius=sphere_radius * u.nanometer,
                    force_const=force_const * u.kilojoules_per_mole / (u.nanometer * u.nanometer)
                )
                dist.append(rest)
            except KeyError:
                #print(f"Warning: Atom {atom_name} not found in residue {resid}")
                continue

    return dist


# Function from shared folder, edited and commented 
# Function to create the cartesian restraints, this is used to keep the protein from unfolding by restraining the CA atoms
# The function takes the system, the scaler, the residues to restrain, the delta and the force constant as input
def make_cartesian_collections(s, scaler, residues, delta=0.2*u.nanometers, k=250.):
    cart = []   # List to store the restraints
    backbone = ['CA']  # The atoms to restrain
    
    for i in residues:
        for b in backbone:
            try:
                atom_index = s.index.atom(i - 1, b)  # Try to get the atom index
                x, y, z = (s.template_coordinates[atom_index] / 10) * u.nanometers  # Convert coordinates

                rest = s.restraints.create_restraint(
                    'cartesian', scaler, LinearRamp(0.0, 15.0, 0.0, 1.0), atom_index=atom_index,
                    x=x, y=y, z=z, delta=delta, force_const=k*u.kilojoules_per_mole/(u.nanometer*u.nanometer)
                )

                cart.append(rest)  # Append the restraint to the list
            except KeyError:
                #print(f"Skipping residue {i}: No 'CA' atom found.")  # Debugging message
                pass  # Silently skip residues without a 'CA' atom

    return cart

# Function to create the distance restraints
# The function takes the filename, the system, the scaler, the distance parameters and the trust factor as input
# The function reads the file containing the distance restraints, extracts the residue and atom information and creates the distance restraints in the simulation system
def get_dist_restraints(filename, s, scaler, max_dist, r1=0.0, r2=0.0):
    """
    Reads a file containing distance restraints and applies them to the system.
    The file should contain lines in the format:
    residue_number_1, residue_name_1, atom_name_1, residue_number_2, residue_name_2, atom_name_2, distance

    Parameters:
        filename (str): Path to the file containing distance restraints.
        s (System): MELD system object.
        scaler (Scaler): Scaler for restraint strength.
        max_dist (float): Maximum allowed distance for a restraint.
        r1, r2 (float, optional): Additional distance parameters in nanometers.

    Returns:
        List of restraints added to the system.
    """

    restraints = []
    
    # Convert r1, r2 to nanometers
    r1_value = r1 * u.nanometers
    r2_value = r2 * u.nanometers

    # Define r3 and r4 using the max contact distance and adjustments
    r3_value = (max_dist * u.nanometers) + (0.15 * u.nanometers)
    r4_value = r3_value + (0.4 * u.nanometers)  # Prevents energy from blowing up/getting to infinity

    with open(filename, 'r') as f:
        for line in f:
            # Use regular expressions to extract the data from the line
            match = re.match(r"(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+([0-9\.]+)", line.strip())
            if match:
                # Extract the values
                res_num_1 = int(match.group(1))
                res_name_1 = match.group(2)
                atom_1_name = match.group(3)
                res_num_2 = int(match.group(4))
                res_name_2 = match.group(5)
                atom_2_name = match.group(6)
                distance = float(match.group(7))

                # Debugging
                #print(f"Residue 1: {res_num_1} {res_name_1} {atom_1_name}")
                #print(f"Residue 2: {res_num_2} {res_name_2} {atom_2_name}")

                try:
                    # Convert residue numbers to residue indices if needed
                    res_idx_1 = s.index.residue(res_num_1)  # Modify if MELD uses different indexing
                    res_idx_2 = s.index.residue(res_num_2)
                    #print(f"Residue 1 index: {res_idx_1}")
                    #print(f"Residue 2 index: {res_idx_2}")

                    # Get atom index for atom 1
                    atom_1_index = s.index.atom(res_idx_1, atom_1_name, one_based=True)
                    #atom_1_index = s.index.atom(res_num_1-1, atom_1_name)
                    #print(f"Atom 1 index: {atom_1_index}")

                    # Get atom index for atom 2
                    atom_2_index = s.index.atom(res_idx_2, atom_2_name, one_based=True)
                    #atom_2_index = s.index.atom(res_num_2-1, atom_2_name)
                    #print(f"Atom 2 index: {atom_2_index}")

                except Exception as e:
                    print(f"Error: Failed to find atoms. Details: {e}")
                    continue  # Skip this restraint if an error occurs

                # Create DistanceRestraint object
                rest = s.restraints.create_restraint(
                    'distance',
                    scaler=scaler, 
                    ramp=LinearRamp(0.0, 500.0, 0.0, 1.0), 
                    r1=r1_value, 
                    r2=r2_value, 
                    r3=r3_value, 
                    r4=r4_value, 
                    k=350 * u.kilojoules_per_mole / (u.nanometer * u.nanometer),  
                    atom1=atom_1_index, 
                    atom2=atom_2_index
                )
                restraints.append(rest)

    return restraints

# Function to generate tleap commands
def tleap_commands(subsystem, mol_id = 'ATP'):
    tleap_commands = []
    tleap_commands = subsystem.generate_tleap_input(mol_id)
    return tleap_commands

# Function to Prepare any inputs needed for tleap
def prepare_tleap_input(subsystem, mol_id = 'ATP'):
    subsystem.prepare_for_tleap(mol_id)

# Function to create and setup the system
def setup_system(templates, PROTEIN_CHAIN, MAX_CONTACT_DIST):
    print("Setting up the system")
    #templates = glob.glob('TEMPLATES/*.pdb')       # read the template file, can be multiple
    print(f'Found templates: {templates}')

    # Class to handle building a System from SubSystems.
    build_options = meld.AmberOptions(forcefield='ff14sbside', implicit_solvent_model='gbNeck2', use_bigger_timestep=True, cutoff=1.8*u.nanometer)

    p = meld.AmberSubSystemFromPdbFile(templates[0])         #build the system

    list_of_tleap_commands = tleap_commands(p)                                       #generate tleap commands
    print(list_of_tleap_commands)
    print("Tleap commands generated")
    prepare_tleap_input(p)                                  #prepare tleap input
    print('Tleap prepared')
    print('Tleap commands generated and input prepared')

    frcmod_file = '/home/sjt3532/meld_atp/input_files/atp.frcmod'               #add force field file
    p.add_frcmod_file(frcmod_file)               #add force field file
    print(f'Added frcmod file {frcmod_file}')
    lib_file = '/home/sjt3532/meld_atp/input_files/atp.lib'                     #add library file
    p.add_lib_file(lib_file)                     #add library file
    print(f'Added lib file {lib_file}')
    prep_file = '/home/sjt3532/meld_atp/input_files/atp.prepi'                  #add prep file
    p.add_prep_file(prep_file)                  #add prep file
    print(f'Added prep file {prep_file}')
    print(f'Added frcmod file {frcmod_file}, lib file {lib_file}, and prep file {prep_file}')

    b = meld.AmberSystemBuilder(build_options)           #use amber builder
    s = b.build_system([p]).finalize()                 #build the system
    print(f'System finalized with {s.n_atoms} atoms')

    # Set up the temperature ladder and scalers
    s.temperature_scaler = meld.system.temperature.GeometricTemperatureScaler(0.0, 0.3, 300*u.kelvin, 450*u.kelvin)    #setup temperature range 300K to 500K for replicas. 0 is for the first replcia and 0.4 is for 30*0.4= 12th replica i.e. we assign temperature from 300 to 500K on first 12 replicas and then contast 500K for rest. This temperature range is distributed geometrically over 12 replcias.
    protein_scalers = s.restraints.create_scaler('constant')    # defining a constant distance scaler i.e. it will keep restraint strength equal through the replica ladder
    # Adaptive replica
    lig_protein_scalers = s.restraints.create_scaler('nonlinear', alpha_min=0.2, alpha_max=0.8, factor=4.0) # Defining a nonlinear distance scaler. 1st to 12th replica will have maximum restraint strength and then from 12 to 30th it will decreas making 0 at the 30th
    sphere_scaler = s.restraints.create_scaler('constant')    # defining a constant distance scaler i.e. it will keep restraint strength equal through the replica ladder
    print('Temperature ladder and scalers set up')

    # Set up caresian coordiantes
    system_index = list(range(1, 1+(s.n_atoms))) # This is 1-5724, this includes all the atoms in the system (ligand and protein)
    protein_positions = make_cartesian_collections(s, protein_scalers, PROTEIN_CHAIN)
    s.restraints.add_as_always_active_list(protein_positions)
    print(f'Protein positions enforced with {len(protein_positions)} restraints')

    # Set up Protein Ligand contacts
    trust_in_group = 33 # Cluster 1 has 66 contacts, so enforcing 50% of the contacts would give us 33, since cluster 1 has the least number of contacts
    prot_ligand_restraints = get_dist_restraints('/home/sjt3532/meld_atp/input_files/final_contacts_updated.txt', s, lig_protein_scalers, MAX_CONTACT_DIST, r1 = 0.0, r2 = 0.0) # Update this
    s.restraints.add_selectively_active_collection(prot_ligand_restraints, trust_in_group)  # Trusting 50% of contacts in the shortest group in the restraint file providing flexibility to the receptor.
    print(f'Protein ligand contacts enforced with {len(prot_ligand_restraints)} restraints')

    # Print inpput for sphere restraints
    print("Sphere radius:", sphere_radius)
    print("Ligand heavy atoms:", LIGAND_HEAVY_ATOMS)
    print("Ligand chain:", LIGAND_CHAIN)
    print("Protein chain:", PROTEIN_CHAIN)

    # Set up the sphere restraints
    sphere_restraints = make_sphere_restraints(s, sphere_scaler, [359], LIGAND_HEAVY_ATOMS, sphere_radius)
    s.restraints.add_as_always_active_list(sphere_restraints)
    print(f'Sphere restraints enforced with {len(sphere_restraints)} restraints')

    # Options
    #We save 1 frame in each 11111 frames, i.e. every 50 ps
    options = meld.RunOptions(timesteps = 11111, minimize_steps = 20000)
    #options.implicit_solvent_model = 'gbNeck2'        #implicit solvent gbNeck2 model, might be automatically calling that
    remd = meld.setup_replica_exchange(s, n_replicas = N_REPLICAS, n_steps = N_STEPS)
    meld.setup_data_store(s, options, remd)
    print("Replica exchange setup complete")

    return s.n_atoms

# Set up system
n_atoms = setup_system(templates, PROTEIN_CHAIN, MAX_CONTACT_DIST)

print(f'Set up system with {n_atoms} atoms')