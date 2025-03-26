# Creates .frcmod from .mol2
parmchk2 -i sustiva.mol2 -f mol2 -o sustiva.frcmod

# Runs tleap on the ligand
tleap -f tleap.in

# Runs tleap on the protein ligand system
tleap -f tleap2.in

# Run minimization on the system
sander -O -i /Users/stuytschaevers/Desktop/Thesis/AMBER/TrwD_atp_mg2/input_files/min.in -o complex_noH.out -p complex_noH.prmtop -c complex_noH.inpcrd  -r complex_noH_min.ncrst  &

# Converts .ncrst and .prmtop to PDB
ambpdb -p complex_noH.prmtop -c complex_noH_min.ncrst > complex_noH_min.pdb

