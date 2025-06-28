#!/bin/bash

# put 1)protein.pdb, 2)ligand.gro, 3)ligand.itp, 4)combine_structures.sh, 5)ions.mdp, 6)md.mdp, 7)min.mdp, 8)nvt.mdp, 9)npt.mdp  in same directory

# Command to run: "chmod +x Gromacs_com.sh" "./Gromacs_com.sh"

echo "#############Generate force field########################"
gmx pdb2gmx -f protein.pdb -o protein.gro

./combine_structures.sh protein.gro ligand.gro

#topol.top update 

sed -i '/^; Include forcefield parameters$/ {
    n
    a #include "ligand.itp"
    n
}' topol.top

sed -i '/^; Compound        #mols$/ {
    n
    a lig                 1
}' topol.top

echo "Complex structure created as complex.gro"

echo "#############Generate Cubic box####################"
gmx editconf -f complex.gro -o newbox.gro -bt cubic -c -d 1.0

echo "##############Solvate the system########################"
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro 

echo "###########################Add ions##########################"
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral <<EOF
15
EOF

echo "########################Energy Minimization#######################"
gmx grompp -f min.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

echo "##########################Make Ligand-H in index##########################" 
gmx make_ndx -f ligand.gro -o index_ligand.ndx <<EOF
0 & ! a H*
q
EOF

echo "#########################Generate Position restrain for ligand###################"
gmx genrestr -f ligand.gro -n index_ligand.ndx -o posre_ligand.itp -fc 1000 1000 1000 <<EOF
3
EOF

sed -i '/^; Include forcefield parameters$/ {
    n
    n 
    a #ifdef POSRES
    a #include "posre_ligand.itp"
    a #endif
    n
}' topol.top

echo "Make Protein | Ligand index"  
gmx make_ndx -f em.gro -o index.ndx <<EOF
1 | 13
q
EOF

echo "Equilibration NVT"
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -n index.ndx
gmx mdrun -deffnm nvt -v

echo "Equilibration NPT"
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -o npt.tpr -n index.ndx
gmx mdrun -deffnm npt -v

echo "Generate md_run file md_0_10.tpr"
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_10.tpr -n index.ndx