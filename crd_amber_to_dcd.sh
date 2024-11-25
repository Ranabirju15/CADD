#!/bin/bash

# ** Conversion of Amber trajectory to DCD trajectory and extracting structural PDB, "first_frame.pdb" **
# Amber topology (.prmtop), and trajectory (.crd) needed as input
# Modify the input files name and interval in below section

# Define input and output files
CRD_FILE="system.crd"
DCD_FILE="system.dcd"
PDB_FILE="first_frame.pdb"
TOPO_FILE="system.prmtop"
INTERVAL=1  # Default interval for skipping frames

# Check if the correct number of arguments are passed
if [ $# -gt 0 ]; then
    CRD_FILE=$1
    DCD_FILE=$2
    INTERVAL=$3
fi

# Progress function
show_progress() {
    local CURRENT=$1
    local TOTAL=$2
    local STEP_DESC=$3
    echo -ne "[$CURRENT/$TOTAL] $STEP_DESC...\r"
}

TOTAL_STEPS=4  # Total number of steps

# Step 1: Prepare cpptraj input for converting CRD to DCD and saving only the first frame as PDB
show_progress 1 $TOTAL_STEPS "Preparing cpptraj input for CRD to DCD conversion"
cat > convert_crd_to_dcd.in <<EOF
parm $TOPO_FILE
trajin $CRD_FILE 1 last $INTERVAL
# Remove water and ions
strip :WAT
strip :Na+
strip :Cl-
# Save the full trajectory in DCD format
trajout $DCD_FILE dcd
# Save only the first frame as a PDB file
#trajout $PDB_FILE pdb first
EOF


# Step 2: Run cpptraj to convert the trajectory from CRD to DCD, and generate PDB file for the first frame
show_progress 2 $TOTAL_STEPS "Running cpptraj for CRD to DCD conversion and generating PDB file"
cpptraj -i convert_crd_to_dcd.in || { echo "cpptraj conversion failed."; exit 1; }

# Step 3: Extract the first frame and save it to first_frame.pdb
show_progress 3 $TOTAL_STEPS "Extracting first frame PDB"
cat > first_pdb.in <<EOF
parm $TOPO_FILE
trajin $CRD_FILE 1 1 $INTERVAL
# Remove water and ions
strip :WAT
strip :Na+
strip :Cl-
# Save only the first frame as a PDB file
trajout $PDB_FILE
EOF
cpptraj -i first_pdb.in || { echo "cpptraj conversion failed."; exit 1; }
# Finalize
echo -ne "\n"

# Step 4: Completion message
echo "Conversion completed successfully!"
echo "Output files:"
echo "  DCD (full trajectory): $DCD_FILE"
echo "  PDB (first frame): $PDB_FILE"

