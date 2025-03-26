# Load the 'msa' library for multiple sequence alignment
library(msa)

# Load the 'Biostrings' library for sequence manipulation
library(Biostrings)

# Load the 'seqinr' library for converting MSA result to a seqinr alignment object
library(seqinr)

# Load the 'ape' library for phylogenetic tree construction
library(ape)

# Set the working directory to the desired path
setwd("/Users/stuytschaevers/Desktop/Thesis/msa")

# Read the amino acid sequences from the FASTA file (default assumes type = "AA")
fasta_file <- "/Users/stuytschaevers/Desktop/Thesis/msa/combined.fasta"
aa_sequences <- readAAStringSet(fasta_file)

# Default Alignment, clustalW
# Perform the multiple sequence alignment
msa_result <- msa(aa_sequences, type = "protein")
# Convert the MSA result to a seqinr alignment object
msa_seqinr  <- msaConvert(msa_result, type = "seqinr::alignment")
# Call msaPrettyPrint to generate the alignment in the PDF without limiting vertical positions
msaPrettyPrint(
  msa_result,
  output = "pdf",
  paperWidth=8,
  subset = 1:length(aa_sequences),
  showNames = "none",
  showLogo="top",
  logoColors="rasmol", 
  shadingMode="identical",
  consensusColor = "ColdHot",
  showLegend = TRUE,
  askForOverwrite = FALSE
)
# Compute the identity matrix
identity_matrix <- dist.alignment(msa_seqinr, matrix = "identity")
# Convert the identity matrix to a data frame
identity_df <- as.data.frame(as.matrix(identity_matrix))
rownames(identity_df) <- colnames(identity_matrix)
# Save the identity matrix data frame as a CSV file
csv_file <- "/Users/stuytschaevers/Desktop/Thesis/msa/identity_matrix.csv"
write.csv(identity_df, file = csv_file, row.names = TRUE)
# Print a message confirming the CSV file creation
cat("Identity matrix data frame saved as", csv_file, "\n")


#create phylogenetic try using neighbor-joining method
proteinTree <- nj(identity_matrix)
plot(proteinTree, main="Phylogenetic Tree of Protein Sequences Aligned with ClustalW")

# Check the version of the 'msa' library
print(packageVersion("msa"))
# Check the version of the 'Biostrings' library
print(packageVersion("Biostrings"))
# Check the version of the 'seqinr' library
print(packageVersion("seqinr"))
# Check the version of the 'ape' library
print(packageVersion("ape"))
