#!/bin/bash

# Command to run: "chmod +x combine_structures.sh" "./combine_structures.sh protein.gro ligand.gro"

# Check if the required input files are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 protein.gro ligand.gro"
    exit 1
fi

# Assign input file names
protein_gro="$1"
ligand_gro="$2"

# Extract the number of atoms from protein.gro and ligand.gro
protein_atoms=$(sed '2q;d' "$protein_gro" | awk '{print $1}')
ligand_atoms=$(sed '2q;d' "$ligand_gro" | awk '{print $1}')

# Calculate the total number of atoms
total_atoms=$((protein_atoms + ligand_atoms))

# Create complex.gro with the calculated total atoms
{
  head -n 1 "$protein_gro"
  echo "$total_atoms"
  sed '$d' "$protein_gro" | tail -n +3
  sed '$d' "$ligand_gro" | tail -n +3
} > complex.gro

echo "Complex structure created as complex.gro"
