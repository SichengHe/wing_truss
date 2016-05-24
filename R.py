'''
Computes R matrix mapping internal forces to forces at the nodes
Sicheng He,
May, 2016
'''

from scipy import *
import numpy as np
import argparse
import sys
import math

def getR(totalFlag):

    # read in the data
    constraints = loadtxt('output/data_constraints.dat')
    elements = loadtxt('output/data_elems.dat')
    nodes = loadtxt('output/data_nodes.dat')

    # convert them to matrix form
    constraints =np. matrix(constraints)
    elements = np.matrix(elements)
    nodes = np.matrix(nodes)

    # convert the data type in constraints and elements to int
    constraints=constraints.astype(int)
    elements=elements.astype(int)

    # size of the problem
    NConstraints = matrix(constraints.shape)[0,1]
    NNode = matrix(nodes.shape)[0,0]
    NElement = matrix(elements.shape)[0,0]

    # create a connectivity list
    # connectivity[i][0]: # of nodes connected to ith node
    # connectivity[i][1][j]: index of connected node  with ith node
    # connectivity[i][2][j]: index of connected element with ith node

    connectivity = []

    for i in range(NNode):

        locCon = [[],[],[]]
        locCon[0].append(0)
        connectivity.append(locCon)


    for i in range(NElement):

        leftNode = elements[i,0]
        rightNode = elements[i,1]

        # add the connected node
        connectivity[leftNode][1].append(rightNode)
        connectivity[rightNode][1].append(leftNode)

        # add the connected element
        connectivity[leftNode][2].append(i)
        connectivity[rightNode][2].append(i)

        # update the # of neighbor
        connectivity[leftNode][0][0] += 1
        connectivity[rightNode][0][0] += 1

    # totalR matrix (including constrianed)

    R = zeros((3*NNode,NElement))

    for i in range(NNode):

        for j in range(connectivity[i][0][0]):

            locNode = connectivity[i][1][j]
            locElement = connectivity[i][2][j]


            distance = np.linalg.norm(nodes[locNode,:]-nodes[i,:])

            # x load

            R[3*i,locElement] = (nodes[locNode,0]-nodes[i,0])/distance

            # y load

            R[3*i+1,locElement] = (nodes[locNode,1]-nodes[i,1])/distance

            # z load

            R[3*i+2,locElement] = (nodes[locNode,2]-nodes[i,2])/distance


    # get rid of constrained ones
    # we assume the
    if (not totalFlag):

        for i in range(NConstraints):

            conNode = constraints[0,NConstraints-i-1]

            print(conNode)

            R=np.delete(R,3*conNode+2,0)
            R=np.delete(R,3*conNode+1,0)
            R=np.delete(R,3*conNode+0,0)

    print('eig00000000000',(np.linalg.eig(R.dot(np.transpose(R))))[0])


    return R
