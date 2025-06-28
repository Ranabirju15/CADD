import parmed as pmd
amber = pmd.load_file("complex.prmtop", "complex.inpcrd")
amber.save('topol.top')
amber.save('inpcrd.gro')
