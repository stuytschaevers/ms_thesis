#!/bin/bash

# Define the input directory containing the min.out files
min_out_dir="/Users/stuytschaevers/Desktop/Thesis/AMBER/TrwD_atp/min_out"

# Define the output directory where the energy and RMS values will be saved
output_dir="/Users/stuytschaevers/Desktop/Thesis/AMBER/TrwD_atp/min_energies"

# Ensure the output directory exists, or create it
mkdir -p "$output_dir"

# Iterate over each min.out file in the input directory
for min_file in "$min_out_dir"/*_min.out; do
    # Extract the base filename (e.g., TrwD_atp_cluster_24)
    base_filename=$(basename "$min_file" _min.out)
    
    # Define the output filename based on the base filename, but in the output directory
    output_file="${output_dir}/${base_filename}_energy_rms_values.txt"
    
    # Extract NSTEP, ENERGY, and RMS values and save them to the corresponding output file
    grep 'NSTEP' -A 1 "$min_file" | grep -v -e 'NSTEP' -e '--' | awk '{print $1, $2, $3}' > "$output_file"
    
    echo "NSTEP, ENERGY, and RMS values extracted from $min_file and saved to $output_file"
done

#grep Energy -A 1 min.out: searches for the word energy in your min.out file, -A 1 prints line 1 after the word energy
# |grep -v -e Energy -e ‘—‘: this I think filters out for extra characters like dashes

# awk{print}: prints the first two fields in a line. This gets the step number and the energy numerical value

# And writes it into a text file