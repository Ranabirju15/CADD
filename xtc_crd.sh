#!/bin/bash

# ********Convert Gromacst Topogoly and Trajectory files into Amber format*********

# Modify your input Gromacs file (structure file: .gro, trajectory: .xtc, topology file: .xtc) names.
# Also, copy the amber force field (ff) into the current directory and update the ff name in line no: 13. 
# Which can be found default: /usr/local/gromacs/share/gromacs/top/ 

# Define filenames (edit as necessary)
TOPOL="topol.top"
GRO="em.gro"
XTC="md_0_10_center.xtc"
FORCEFIELD="./amber99sb-ildn.ff"

# Progress function
show_progress() {
    local CURRENT=$1
    local TOTAL=$2
    local STEP_DESC=$3
    echo -ne "[$CURRENT/$TOTAL] $STEP_DESC...\r"
    sleep 1  # Optional: simulate time delay
}

TOTAL_STEPS=4  # Total number of steps

# Step 1: Update topology file to include the correct force field path
show_progress 1 $TOTAL_STEPS "Updating topology file"
sed -i "s|amber99sb-ildn.ff|$FORCEFIELD|g" $TOPOL

# Step 2: ParmEd conversion script
show_progress 2 $TOTAL_STEPS "Generating AMBER topology and coordinates"
cat > convert_gmx_to_amber.py <<EOF
import parmed as pmd

# Load GROMACS topology and GRO structure
gmx_top = pmd.load_file("$TOPOL", xyz="$GRO")

# Write AMBER prmtop and inpcrd files
gmx_top.save("system.prmtop", format="amber")
gmx_top.save("system.inpcrd", format="rst7")
EOF

# Run ParmEd conversion
python3 convert_gmx_to_amber.py || { echo "ParmEd conversion failed."; exit 1; }

# Step 3: Create cpptraj input to convert XTC to CRD
show_progress 3 $TOTAL_STEPS "Preparing cpptraj input for XTC to CRD conversion"
cat > convert.in <<EOF
parm system.prmtop
trajin $XTC
trajout system.crd crd
EOF

# Run cpptraj
cpptraj -i convert.in || { echo "cpptraj conversion failed."; exit 1; }

# Step 4: Finalize
show_progress 4 $TOTAL_STEPS "Finalizing"
echo -ne "\n"

# Completion message
echo "Conversion completed successfully!"
echo "Output file: ./system.crd"

