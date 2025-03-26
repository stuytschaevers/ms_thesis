import re

def extract_data_from_file(file_path):
    # Open the file and read it
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Store the extracted data in a list
    extracted_data = []

    # Loop through the lines and check for the pattern
    for line in lines:
        # Regular expression to match the desired line
        if re.match(r'^\d+,-?\d+\.\d+,\d+\.\d+', line):
            extracted_data.append(line.strip())  # Add the line without leading/trailing spaces

    return extracted_data

# Example usage
file_path = 'output.txt'
data = extract_data_from_file(file_path)
for line in data:
    print(line)
