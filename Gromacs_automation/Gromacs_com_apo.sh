#!/bin/bash

# put 1)protein.pdb, 2)ions.mdp, 3)md.mdp, 4)min.mdp, 5)nvt.mdp, 6)npt.mdp  in same directory

# Command to run: "chmod +x Gromacs_com.sh" "./Gromacs_com.sh"

echo "#############Generate force field########################"
gmx pdb2gmx -f protein.pdb -o protein.gro


echo "Complex structure created as complex.gro"

echo "#############Generate Cubic box####################"
gmx editconf -f protein.gro -o newbox.gro -bt cubic -c -d 1.0

echo "##############Solvate the system########################"
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro 

echo "###########################Add ions##########################"
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral <<EOF
13
EOF

echo "########################Energy Minimization#######################"
gmx grompp -f min.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

echo "Equilibration NVT"
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt -v

echo "Equilibration NPT"
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -o npt.tpr
gmx mdrun -deffnm npt -v

echo "Generate md_run file md_0_10.tpr"
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_10.tpr
