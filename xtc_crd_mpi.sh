#!/bin/bash

# ********Convert Gromacst Topogoly and Trajectory files into Amber format*********

# Modify your input Gromacs file (structure file: .gro, trajectory: .xtc, topology file: .xtc) names.
# Also, copy the amber force field field into the current directory. 
# Which can be found default: /usr/local/gromacs/share/gromacs/top/ 

# Define filenames
TOPOL="topol.top"
GRO="em.gro"
XTC="md_0_10_center.xtc"
FORCEFIELD="./amber99sb-ildn.ff"
PROCESSES=6  # Number of CPU threads for MPI

# Default interval (skip every n-th frame)
DEFAULT_INTERVAL=50

# Check for interval argument
if [ $# -gt 0 ]; then
    INTERVAL=$1
    echo "Using custom trajectory interval: $INTERVAL"
else
    INTERVAL=$DEFAULT_INTERVAL
    echo "Using default trajectory interval: $INTERVAL"
fi

# Progress function
show_progress() {
    local CURRENT=$1
    local TOTAL=$2
    local STEP_DESC=$3
    echo -ne "[$CURRENT/$TOTAL] $STEP_DESC...\r"
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
trajin $XTC 1 last $INTERVAL
trajout system.crd crd
EOF

# Run cpptraj.MPI with mpirun
show_progress 4 $TOTAL_STEPS "Running cpptraj.MPI with MPI"
mpirun -np $PROCESSES cpptraj.MPI -i convert.in || { echo "cpptraj.MPI conversion failed."; exit 1; }

# Finalize
echo -ne "\n"

# Completion message
echo "Conversion completed successfully!"
echo "Output file: ./system.crd"

