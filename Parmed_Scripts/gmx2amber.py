import parmed as pmd

gmx_top = pmd.load_file('topol.top', xyz='inpcrd.gro')
gmx_top.save('pmaa.top', format='amber')
gmx_top.save('pmaa.crd', format='rst7')
