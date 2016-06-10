from __future__ import division
import numpy as numpy
import pylab
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
from scipy import *
from math import *
import sqlitedict


def plot(ax, sol, elems, nodes):
    for i in xrange(len(sol)):
        sol_loc = sol[i]

        ind_node_1 = elems[i][0]
        ind_node_2 = elems[i][1]

        x1 = nodes[ind_node_1][0]
        y1 = nodes[ind_node_1][1]
        z1 = nodes[ind_node_1][2]

        x2 = nodes[ind_node_2][0]
        y2 = nodes[ind_node_2][1]
        z2 = nodes[ind_node_2][2]

        x = [x1, x2]
        y = [y1, y2]
        z = [z1, z2]
        #print '++++++',z

        ax.plot(x,y,z,'-ro',
                #c=colorVal)
                markersize=0.1,
                linewidth=5e1*sol_loc)

        ax.plot(x,y,z,'-k',
                #c=colorVal)
                linewidth=5e-2)

        i += 1

    ax.view_init(elev=13., azim=-106)

    

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# s = s.replace(',', '')
sol = [line.strip() for line in open("Solution_3D.res", 'r')]
nodes = loadtxt('../output/data_nodes.dat')
elems = loadtxt('../output/data_elems.dat')

sol = sol[0].split(',')
sol = [1e-4*float(y) for y in sol]

db = sqlitedict.SqliteDict('data.db', 'openmdao')

index = 0
for case_name, case_data in db.iteritems():
    if "metadata" in case_name or "derivs" in case_name:
        continue # don't plot these cases

    sol2 = case_data['Unknowns']['areas']
    print case_data['Unknowns']['volume']

    if 1 and index % 10 == 0:
        plot(ax, sol2, elems, nodes)
        pylab.savefig('tmp2/frame%0i.png'%(index))
        print index
        ax.cla()

    print index
    index += 1


if 0:
    plot(ax, sol2, elems, nodes)
    pylab.savefig('comp0_stress.png')
    #plt.show()

