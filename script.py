'''
Run script
Dr. John T. Hwang
Sicheng He
May, 2016
'''

from __future__ import division
import numpy
import pylab
from utils import setup_problem, writeBDF, get_matrix, solve
import tecwrite


u1 = 0.005
u2 = 0.995

EA = 10.0e9

geom_file = 'CRM_AVL/wing_coarse.avl'
upp_file = 'airfoils/rae2822_upp.txt’
low_file = 'airfoils/rae2822_low.txt’
results_file = 'CRM_AVL/results_coarse.txt'

factor = 1

xyz, nodes, bars, constrained, forces, forcesArray = setup_problem(u1, u2, geom_file, upp_file, low_file, results_file,factor)

surfs = []
surfs.append(xyz[:, :,  0, :])
surfs.append(xyz[:, :, -1, :])
surfs.append(xyz[:,  0, :, :])
surfs.append(xyz[:, -1, :, :])
surfs.append(xyz[ 0, :, :, :])
surfs.append(xyz[-1, :, :, :])
tecwrite.write_surf_multi('surf.dat’, surfs)

writeBDF('mesh.bdf', nodes, bars+1)






mat, Kmat = get_matrix(EA, nodes, bars, constrained)

rhs = numpy.zeros(mat.shape[0])
rhs[:3*len(forces)] = forces.flatten()



sol = solve(mat, rhs)[:3*len(forces)]
sol_surf = sol.reshape(xyz.shape)
xyz += sol_surf

surfs = []
surfs.append(xyz[:, :,  0, :])
surfs.append(xyz[:, :, -1, :])
surfs.append(xyz[:,  0, :, :])
surfs.append(xyz[:, -1, :, :])
surfs.append(xyz[ 0, :, :, :])
surfs.append(xyz[-1, :, :, :])
tecwrite.write_surf_multi('surf2.dat’, surfs)

writeBDF('mesh2.bdf', nodes + sol.reshape(nodes.shape), bars+1)



numpy.savetxt('output/data_nodes.dat', nodes)
numpy.savetxt('output/data_elems.dat', bars)
numpy.savetxt('output/data_constraints.dat', constrained)
numpy.savetxt('output/data_forces.dat', forces)



import R
Rmat = numpy.matrix(R.getR(False)) # without the fixed pts
Rmat_fixpts = numpy.matrix(R.getR(True)) # with the fixed points
numpy.savetxt('output/data_Rmat_2x.dat', Rmat)
numpy.savetxt('output/data_Rmat_fixpts_2x.dat', Rmat_fixpts)



# following is a comparison between the direct assembly K and indirect (K=R EA/L R^T)

L = numpy.zeros(bars.shape[0])
for ielem in xrange(bars.shape[0]):
    xyz1 = nodes[bars[ielem, 0], :]
    xyz2 = nodes[bars[ielem, 1], :]
    L[ielem] = numpy.linalg.norm(xyz2 - xyz1)

EA_L_diag = numpy.diag(numpy.ones(Rmat.shape[1]) * EA / L)

Kmat2 = Rmat.dot(EA_L_diag).dot(Rmat.T)

numx, numy, numz, dim = xyz.shape
rhs2 = forcesArray[:,1:,:,:].flatten()
sol2 = numpy.linalg.solve(Kmat2, rhs2)
sol2 = numpy.reshape(sol2,(numx,numy-1,numz,3))
sol2 = sol2.tolist()

loc_mat = (numpy.zeros((numz,dim))).tolist()
for i in range(numx):
    sol2_x = sol2[i]
    sol2_x.insert(0,loc_mat)
    sol2[i] = sol2_x
sol2 = numpy.asarray(sol2)
writeBDF('mesh3.bdf', nodes + sol2.reshape(nodes.shape), bars+1)
Kmat1 = Kmat.todense()
