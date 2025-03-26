#!/usr/bin/env python
# coding: utf-8

#this script is new which includes amber14-all.xml over amber99sb.xml

import os
from sys import stdout
import numpy as np
from openmm import *
from openmm.app import *
from openmm.unit import *

# Get the directory of the script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Set the simulation parameters
simulation_steps = 250000000
output_prefix = os.path.join(script_dir, "TrwD_intensive")
reporter_interval = 10000  # Report every 10,000 steps for 250,000,000 steps simulation

# Create a force field for the simulation
forcefield = ForceField("amber14-all.xml", "implicit/gbn2.xml")

# Load the protein structure from the PDB file
pdb = PDBFile("TrwD_intensive.pdb")
print(pdb.topology)

# Create a Modeller object and add hydrogens
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)

# Create the simulation system with specified parameters
system = forcefield.createSystem(modeller.topology, nonbondedCutoff=3 * nanometer, constraints=HBonds)

# Initialize the Langevin integrator with specified parameters
integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond, 2 * femtoseconds)

simulation = Simulation(modeller.topology, system, integrator)

simulation.context.setPositions(modeller.positions)

# Minimize the energy of the system
simulation.minimizeEnergy()

# Get the state of the simulation (positions)
state = simulation.context.getState(getPositions=True)

# Save the minimized structure
output_topology = f"{output_prefix}_topology.pdb"
with open(output_topology, 'w') as f:
    PDBFile.writeFile(simulation.topology, state.getPositions(), f)
print(f"Topology PDB saved as {output_topology}")

# Configure and add reporters to the simulation
simulation.reporters = []

# DCDReporter writes the trajectory to a DCD file
output_trajectory = f"{output_prefix}_traj.dcd"
simulation.reporters.append(DCDReporter(output_trajectory, reporter_interval))

# StateDataReporter writes state data to a CSV file
output_csv = f"{output_prefix}_scalars.csv"
simulation.reporters.append(
    StateDataReporter(output_csv, reporter_interval, time=True, potentialEnergy=True, totalEnergy=True, temperature=True)
)

# Run the simulation for the specified number of steps
simulation.step(simulation_steps)

print(f"Simulation completed. Results saved as {output_topology}, {output_trajectory}, {output_csv}")
