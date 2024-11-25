#!/bin/bash

# r base need to be install in linux system
# install.packages("bio3d", dependencies=TRUE) [Need to install in r/ add the command below before library(bio3d) at line no. 32]

# Automatically detect DCD and PDB files in the current directory
DCD_FILE=$(ls *.dcd 2>/dev/null | head -n 1)
PDB_FILE=$(ls *.pdb 2>/dev/null | head -n 1)

# Check if both files are found
if [ -z "$DCD_FILE" ] || [ -z "$PDB_FILE" ]; then
    echo "Error: Could not find both *.dcd and *.pdb files in the current directory."
    echo "Ensure that exactly one DCD and one PDB file exist in the directory."
    exit 1
fi

# Define output file paths
PC1_OUTPUT="pc1.pdb"
PC2_OUTPUT="pc2.pdb"
PCA_SCORES="PCA_scores.csv"
B_FACTORS_OUTPUT="b_factors.pdb"
B_FACTORS_CSV="b_factors.csv"

echo "Using DCD file: $DCD_FILE"
echo "Using PDB file: $PDB_FILE"

# Temporary R script file
R_SCRIPT=$(mktemp)

# Generate the R script dynamically
cat > "$R_SCRIPT" <<EOF

library(bio3d)

# Read input files
dcd <- read.dcd("$DCD_FILE")
pdb <- read.pdb("$PDB_FILE")

# Select CA atoms
ca.inds <- atom.select(pdb, elety="CA")

# Align trajectory to reference structure
xyz <- fit.xyz(fixed=pdb\$xyz, mobile=dcd,
               fixed.inds=ca.inds\$xyz,
               mobile.inds=ca.inds\$xyz)

# Perform PCA
pc <- pca.xyz(xyz[,ca.inds\$xyz])

# Plot PCA with a color gradient
plot(pc, col=bwr.colors(nrow(xyz)))

# Perform hierarchical clustering and color PCA plot by clusters
hc <- hclust(dist(pc\$z[,1:2]))
grps <- cutree(hc, k=3)
plot(pc, col=grps)

# Generate PCA trajectories
p1 <- mktrj.pca(pc, pc=1, b=pc\$au[,1], file="$PC1_OUTPUT")
p2 <- mktrj.pca(pc, pc=2, b=pc\$au[,2], file="$PC2_OUTPUT")

# Plot contribution of PCs to residue motion
plot.bio3d(pc\$au[,1], ylab="PC (A)", xlab="Residue Position", typ="l", resno=pdb)
points(pc\$au[,2], typ="l", col="blue")

# Calculate and plot DCCM
cij <- dccm(xyz[,ca.inds\$xyz])
plot(cij, resno=pdb)

# Save PCA scores
pc1 <- pc\$z[,1]
pc2 <- pc\$z[,2]
df <- data.frame(PC1=pc1, PC2=pc2)
write.csv(df, file="$PCA_SCORES")

# B-factor calculation
# Mean square fluctuations (MSFs) from selected PCs
msf <- rowSums(pc\$z[, 1:10]^2)  # Use first 10 PCs
b_factors <- 8 * pi^2 * msf      # Convert MSFs to B-factors

# Assign B-factors to the selected CA atoms only
pdb\$atom\$b[ca.inds\$atom] <- b_factors

# Update PDB file and save
write.pdb(pdb, file="$B_FACTORS_OUTPUT")

# Save B-factors as CSV for reference
write.csv(data.frame(Residue=ca.inds\$atom, B_Factor=b_factors), file="$B_FACTORS_CSV")

# Plot B-factors
plot(b_factors, type="h", xlab="Residue Index (CA only)", ylab="B-factor", main="B-factors from PCA")
EOF
EOF

# Run the R script
Rscript "$R_SCRIPT"

# Remove the temporary R script
rm -f "$R_SCRIPT"

# Check if output files are generated and notify the user
if [ -f "$PC1_OUTPUT" ] && [ -f "$PC2_OUTPUT" ] && [ -f "$PCA_SCORES" ] && [ -f "$B_FACTORS_OUTPUT" ] && [ -f "$B_FACTORS_CSV" ]; then
    echo "PCA analysis completed successfully!"
    echo "Output files:"
    echo "  PCA trajectory (PC1): $PC1_OUTPUT"
    echo "  PCA trajectory (PC2): $PC2_OUTPUT"
    echo "  PCA scores: $PCA_SCORES"
    echo "  B-factor PDB: $B_FACTORS_OUTPUT"
    echo "  B-factors CSV: $B_FACTORS_CSV"
else
    echo "Analysis failed. Check the R script output for details."
fi

