# Runs tleap on ligand
tleap -f tleap.in

# Runs tleap on the protein ligand system
tleap -f tleap2_pose50.in.txt

# Runs minimization of hte protein ligand system
sander -O -i min_pose50.in.txt -o min_pose50.out -p TrwD_pose_50.prmtop -c TrwD_pose_50.inpcrd  -r TrwD_pose_50_min.ncrst  &

# Converts .ncrst to PDB
ambpdb -p /Users/stuytschaevers/Desktop/Thesis/atp_mg2+/amber/TrwD_pose_50.prmtop -c /Users/stuytschaevers/Desktop/Thesis/atp_mg2+/amber/TrwD_pose_50_min.ncrst > /Users/stuytschaevers/Desktop/Thesis/atp_mg2+/amber/TrwD_pose_50_min.pdb
