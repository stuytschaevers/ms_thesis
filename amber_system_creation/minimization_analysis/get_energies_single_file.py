import re
import csv

def extract_data_from_file(file_path, output_csv_path):
    # Open the file and read it
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Store the extracted data in a list
    extracted_data = []

    # Loop through the lines and check for the pattern
    for line in lines:
        # Regular expression to match the desired line
        if re.match(r'^\d+,-?\d+\.\d+,\d+\.\d+', line):
            extracted_data.append(line.strip().split(','))  # Split line by commas

    # Write the extracted data to a CSV file
    with open(output_csv_path, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["Step", "Potential Energy (kJ/mole)", "Temperature (K)"])  # Write the header
        writer.writerows(extracted_data)  # Write the extracted data

# Example usage
file_path = '/Users/stuytschaevers/Desktop/Thesis/atp_mg2+/amber/test/from_cluster/output/pose_50_20089012.out'  # Input file path
output_csv_path = '/Users/stuytschaevers/Desktop/Thesis/atp_mg2+/amber/test/from_cluster/output/min_energies/pose_50_energies.csv'  # Output CSV file path
extract_data_from_file(file_path, output_csv_path)

print(f"Data has been extracted and saved to {output_csv_path}")


file_path = '/Users/stuytschaevers/Desktop/Thesis/atp_mg2+/amber/test/from_cluster/output/pose_241_20089013.out'  # Input file path
output_csv_path = '/Users/stuytschaevers/Desktop/Thesis/atp_mg2+/amber/test/from_cluster/output/min_energies/pose_241_energies.csv'  # Output CSV file path
extract_data_from_file(file_path, output_csv_path)

print(f"Data has been extracted and saved to {output_csv_path}")