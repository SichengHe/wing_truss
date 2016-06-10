'''
Computes finite element matrices for wing truss
Dr. John T. Hwang
November, 2015
'''

from __future__ import division
import numpy
import scipy.sparse
import scipy.sparse.linalg
from scipy import interpolate



def read_geom_file(file_name):
    ''' Reads AVL geometry file
    Gets root, Yehudi break, and tip values

    Outputs
    -------
    1. number of elements in the chord-wise direction
    2. number of elements in the span-wise direction
    3. section x values
    4. section y values
    5. section z values
    6. section c values
    7. section t values
    '''

    raw = open(file_name, 'r').readlines()

    sizes = raw[12].split()
    num_chord = int(sizes[0])
    num_span = int(sizes[2])

    data = raw[24].split()
    offset = numpy.array([float(data[ind]) for ind in xrange(3)])

    array = []

    data = raw[31].split()
    array.append([float(data[ind]) for ind in xrange(5)])

    data = raw[37].split()
    array.append([float(data[ind]) for ind in xrange(5)])

    data = raw[43].split()
    array.append([float(data[ind]) for ind in xrange(5)])

    array = numpy.array(array)

    for ind in xrange(3):
        array[ind, :3] += offset

    return num_chord, num_span, \
        array[:, 0], array[:, 1], array[:, 2], \
        array[:, 3], array[:, 4]

def read_airfoil_file(upp_name, low_name):
    ''' Reads and interpolates an airfoil's upper and lower surfaces '''

    upp = numpy.loadtxt(upp_name)
    low = numpy.loadtxt(low_name)

    upp = interpolate.interp1d(upp[:, 0], upp[:, 1])
    low = interpolate.interp1d(low[:, 0], low[:, 1])
    return upp, low

def read_results_file(file_name):
    ''' Reads pressures and strip areas from an AVL results file '''
    lines = open(file_name, 'r').readlines()

    pressure = []
    strip_area = []
    read = False
    for line in lines:
        terms = line.split()
        if len(terms) > 1 and terms[0] == 'I':
            read = True
        elif len(terms) > 1 and terms[0] == 'Yle':
            strip_area.append(float(terms[-1]))
        elif len(terms) != 7:
            read = False
        elif read:
            pressure.append(float(terms[-1]))

    return numpy.array(pressure), numpy.array(strip_area)



def interp_sec(t, tot_y, tot_q):
    ''' Utility function to get an intermediate value of a quantity
    at a point t along the span-wise direction '''
    y = t * tot_y[2]
    if y <= tot_y[1]:
        ind1 = 0
        ind2 = 1
    else:
        ind1 = 1
        ind2 = 2

    t2 = (y - tot_y[ind1]) / (tot_y[ind2] - tot_y[ind1])
    return tot_q[ind1] + (tot_q[ind2] - tot_q[ind1]) * t2



def read_files(u1, u2, geom_file, upp_file, low_file, results_file, factor):
    ''' '''
    num_chord, num_span, sec_x, sec_y, sec_z, sec_c, sec_t = read_geom_file(geom_file)
    upp, low = read_airfoil_file(upp_file, low_file)
    pressure, strip_area = read_results_file(results_file)
    pressure = pressure[:num_chord*num_span].reshape((num_chord, num_span), order='F')

    # Directions
    #  1 x u chord-wise
    #  2 y v span-wise
    #  3 z w vertical

    num1 = num_chord*factor + 1
    num2 = num_span*factor + 1
    num3 = 2

    uvw = numpy.zeros((num1, num2, num3, 3))
    for ind1 in xrange(num1):
        for ind2 in xrange(num2):
            for ind3 in xrange(num3):
                uvw[ind1, ind2, ind3, 0] = u1 + (u2-u1) * ind1 / (num1-1)
                uvw[ind1, ind2, ind3, 1] = ind2 / (num2-1)
                uvw[ind1, ind2, ind3, 2] = ind3

    axis_x = numpy.zeros(num2)
    axis_y = numpy.zeros(num2)
    axis_z = numpy.zeros(num2)
    axis_c = numpy.zeros(num2)
    axis_t = numpy.zeros(num2)
    for ind2 in xrange(num2):
        pos_y = ind2 / (num2-1)
        axis_x[ind2] = interp_sec(pos_y, sec_y, sec_x)
        axis_y[ind2] = interp_sec(pos_y, sec_y, sec_y)
        axis_z[ind2] = interp_sec(pos_y, sec_y, sec_z)
        axis_c[ind2] = interp_sec(pos_y, sec_y, sec_c)
        axis_t[ind2] = interp_sec(pos_y, sec_y, sec_t)

    xyz = numpy.zeros((num1, num2, num3, 3))

    # Start with the parametric coordinates
    xyz[:, :, :, :] = uvw[:, :, :, :]
    for ind2 in xrange(num2):
        xyz[:, ind2, 0, 2] = low(xyz[:, ind2, 0, 0])
        xyz[:, ind2, 1, 2] = upp(xyz[:, ind2, 1, 0])

    # Center about the quarter chord
    xyz[:, :, :, 0] -= 0.25

    # Rotate about the twist
    C = numpy.zeros((3, 3))
    C[1, 1] = 1.0
    for ind2 in xrange(num2):
        twist = axis_t[ind2] * numpy.pi/180.0
        C[0, 0] =  numpy.cos(twist)
        C[0, 2] =  numpy.sin(twist)
        C[2, 0] = -numpy.sin(twist)
        C[2, 2] =  numpy.cos(twist)
        for ind1 in xrange(num1):
            for ind3 in xrange(num3):
                xyz[ind1, ind2, ind3, :] = C.dot(xyz[ind1, ind2, ind3, :])

    # Make the LE the origin again
    xyz[:, :, :, 0] += 0.25

    # Scale by the chord
    for ind2 in xrange(num2):
        chord = axis_c[ind2]
        xyz[:, ind2, :, :] *= chord

    # Translate by the LE pos
    for ind2 in xrange(num2):
        xyz[:, ind2, :, 0] += axis_x[ind2]
        xyz[:, ind2, :, 1] += axis_y[ind2]
        xyz[:, ind2, :, 2] += axis_z[ind2]

    return xyz, pressure, strip_area
    #return xyz

def writeBDF(filename, nodes, bars):
    f = open(filename, 'w')

    def writeLine(line):
        write(line,r=80)
        write('\n')

    def write(line,l=0,r=0):
        if l is not 0:
            n = l - len(line)
            for i in range(n):
                line = ' ' + line
        if r is not 0:
            n = r - len(line)
            for i in range(n):
                line = line + ' '
        f.write(line)

    writeLine('$ Generated by ICEMCFD -  NASTRAN Interface Vers.  4.6.1')
    writeLine('$ Nastran input deck')
    writeLine('SOL 103')
    writeLine('CEND')
    writeLine('$')
    writeLine('BEGIN BULK')

    for k in xrange(nodes.shape[0]):
        ind = k + 1
        write('GRID*   ')
        write(str(ind),l=16)
        write('0',l=16)
        write('%.8E' % nodes[k,0],l=16)
        write('%.8E' % nodes[k,1],l=16)
        write('*')
        write(str(ind),l=7)
        write('\n')
        write('*')
        write(str(ind),l=7)
        write('%.8E' % nodes[k,2],l=16)
        write('0',l=16)
        write(' ',l=16)
        write('0',l=16)
        write(' ',l=8)
        write('\n')

    for k in xrange(bars.shape[0]):
        ind = k + 1
        write('CROD    ')
        write(str(ind),l=8)
        write(str(21),l=8)
        write('%8i' % bars[k,0],l=8)
        write('%8i' % bars[k,1],l=8)
        write('\n')


    writeLine('END BULK')

    f.close()




def setup_problem(u1, u2, geom_file, upp_file, low_file, results_file,factor):
    xyz, pressure, strip_area = read_files(u1, u2, geom_file, upp_file, low_file, results_file, factor)

    num_nodes = numpy.prod(xyz.shape[:3])
    nodes = xyz.reshape((num_nodes, 3))

    num1, num2, num3 = xyz.shape[:3]
    indices = numpy.arange(num_nodes).reshape(xyz.shape[:3])

    num = 0
    num += (num1-1) * (num2-1) * num3
    num += num1 * (num2-1) * num3
    num += num1 * (num2-1) * (num3-1)
    num += (num1-1) * (num2-1) * num3 * 2
    num += (num1-1) * (num2-1) * (num3-1) * 2
    num += num1 * (num2-1) * (num3-1) * 2
    bars = numpy.zeros((num, 2), int)
    counter = 0

    # x-edges
    for ind1 in xrange(num1-1):
        for ind2 in xrange(num2-1):
            for ind3 in xrange(num3):
                bars[counter, 0] = indices[ind1+0, ind2+1, ind3]
                bars[counter, 1] = indices[ind1+1, ind2+1, ind3]
                counter += 1

    # y-edges
    for ind1 in xrange(num1):
        for ind2 in xrange(num2-1):
            for ind3 in xrange(num3):
                bars[counter, 0] = indices[ind1, ind2+0, ind3]
                bars[counter, 1] = indices[ind1, ind2+1, ind3]
                counter += 1

    # z-edges
    for ind1 in xrange(num1):
        for ind2 in xrange(num2-1):
            for ind3 in xrange(num3-1):
                bars[counter, 0] = indices[ind1, ind2+1, ind3+0]
                bars[counter, 1] = indices[ind1, ind2+1, ind3+1]
                counter += 1

    # x-crosses
    for ind1 in xrange(num1):
        for ind2 in xrange(num2-1):
            for ind3 in xrange(num3-1):
                bars[counter, 0] = indices[ind1, ind2+0, ind3+0]
                bars[counter, 1] = indices[ind1, ind2+1, ind3+1]
                counter += 1
                bars[counter, 0] = indices[ind1, ind2+1, ind3+0]
                bars[counter, 1] = indices[ind1, ind2+0, ind3+1]
                counter += 1

    # y-crosses
    for ind1 in xrange(num1-1):
        for ind2 in xrange(num2-1):
            for ind3 in xrange(num3-1):
                bars[counter, 0] = indices[ind1+0, ind2+1, ind3+0]
                bars[counter, 1] = indices[ind1+1, ind2+1, ind3+1]
                counter += 1
                bars[counter, 0] = indices[ind1+0, ind2+1, ind3+1]
                bars[counter, 1] = indices[ind1+1, ind2+1, ind3+0]
                counter += 1

    # z-crosses
    for ind1 in xrange(num1-1):
        for ind2 in xrange(num2-1):
            for ind3 in xrange(num3):
                bars[counter, 0] = indices[ind1+0, ind2+0, ind3]
                bars[counter, 1] = indices[ind1+1, ind2+1, ind3]
                counter += 1
                bars[counter, 0] = indices[ind1+1, ind2+0, ind3]
                bars[counter, 1] = indices[ind1+0, ind2+1, ind3]
                counter += 1

    # constraints
    constrained = indices[:, 0, :].flatten()

    # loads
    Q = 0.5 * 0.4135 * 280**2

    num1_org = int((num1-1)/factor)
    num2_org = int((num2-1)/factor)

    for ind2 in xrange(num2_org-1):
        pressure[:, ind2] *= strip_area[ind2] / (num1_org-1) * Q

    forces = numpy.zeros((num1, num2, num3, 3))
    for i in range(num1_org):
        for j in range(num2_org):
            for k in range(num3):

                force_persurf_loc = pressure[i,j]/2.0
                force_distri_loc = force_persurf_loc/(factor**2)
                force_aux_mat = numpy.zeros((factor,factor,3))
                force_aux_mat[:,:,2] += force_distri_loc

                force_pernode_loc = numpy.zeros((factor+1,factor+1,3))

                force_pernode_loc[:-1, :-1, 2] += force_aux_mat[:,:,2]
                force_pernode_loc[ 1:, :-1, 2] += force_aux_mat[:,:,2]
                force_pernode_loc[:-1,  1:, 2] += force_aux_mat[:,:,2]
                force_pernode_loc[ 1:,  1:, 2] += force_aux_mat[:,:,2]

                for s in range(factor+1):
                    for t in range(factor+1):

                        indx_glo = i*(factor)+s
                        indy_glo = j*(factor)+t

                        forces[indx_glo,indy_glo,k,2] += force_pernode_loc[s,t,2]


    return xyz, nodes, bars, constrained, forces.reshape((num_nodes, 3)),forces 




def get_stiffness_matrix(EA, nodes, bars):
    num_nodes = nodes.shape[0]
    num_elems = bars.shape[0]

    Kdata = numpy.zeros(36 * num_elems)
    Krows = numpy.zeros(36 * num_elems, int)
    Kcols = numpy.zeros(36 * num_elems, int)

    Telem = numpy.zeros((2, 6))
    Kelem = numpy.array([[1, -1], [-1, 1]])

    data_6x6 = numpy.zeros((6, 6), dtype=numpy.int)
    rows_6x6 = numpy.zeros((6, 6), dtype=numpy.int)
    cols_6x6 = numpy.zeros((6, 6), dtype=numpy.int)

    for ind1 in xrange(6):
        for ind2 in xrange(6):
            rows_6x6[ind1, ind2] = ind1%3
            cols_6x6[ind1, ind2] = ind2%3

    data_36 = data_6x6.reshape(36)
    rows_36 = rows_6x6.reshape(36)
    cols_36 = cols_6x6.reshape(36)

    ones11 = numpy.zeros((6, 6), dtype=numpy.int)
    ones12 = numpy.zeros((6, 6), dtype=numpy.int)
    ones21 = numpy.zeros((6, 6), dtype=numpy.int)
    ones22 = numpy.zeros((6, 6), dtype=numpy.int)

    ones11[0:3, 0:3] = 1
    ones12[0:3, 3:6] = 1
    ones21[3:6, 0:3] = 1
    ones22[3:6, 3:6] = 1

    ones11 = ones11.reshape(36)
    ones12 = ones12.reshape(36)
    ones21 = ones21.reshape(36)
    ones22 = ones22.reshape(36)

    index = 0
    for ielem in xrange(bars.shape[0]):
        xyz1 = nodes[bars[ielem, 0], :]
        xyz2 = nodes[bars[ielem, 1], :]
        L = numpy.linalg.norm(xyz2 - xyz1)

        cos_xyz = (xyz2 - xyz1) / L
        Telem[0, 0:3] = cos_xyz
        Telem[1, 3:6] = cos_xyz

        data_6x6[:, :] = Telem.T.dot(Kelem).dot(Telem) * EA / L

        ind1 = 3 * bars[ielem, 0]
        ind2 = 3 * bars[ielem, 1]


        Kdata[index: index+36] += data_36
        Krows[index: index+36] += rows_36
        Kcols[index: index+36] += cols_36

        Krows[index: index+36] += ones11 * ind1
        Kcols[index: index+36] += ones11 * ind1

        Krows[index: index+36] += ones12 * ind1
        Kcols[index: index+36] += ones12 * ind2

        Krows[index: index+36] += ones21 * ind2
        Kcols[index: index+36] += ones21 * ind1

        Krows[index: index+36] += ones22 * ind2
        Kcols[index: index+36] += ones22 * ind2

        index += 36

    return scipy.sparse.csc_matrix((Kdata, (Krows, Kcols)),
                                   shape=(3 * num_nodes, 3 * num_nodes))


def get_constraint_matrix(nodes, constrained):
    num = len(constrained)
    num_nodes = len(nodes)

    data = numpy.ones(3*num)
    lins = numpy.arange(3*num)
    cons = numpy.zeros(3*num)

    for ind in xrange(3):
        cons[ind::3] += 3 * constrained + ind

    C1 = scipy.sparse.csc_matrix((data, (lins, cons)),
                                 shape=(3*num, 3*num_nodes))
    C2 = scipy.sparse.csc_matrix((data, (cons, lins)),
                                 shape=(3*num_nodes, 3*num))
    zmat = scipy.sparse.csc_matrix((3*num, 3*num))

    return C1, C2, zmat


def get_matrix(EA_L, nodes, bars, constrained):
    Kmat = get_stiffness_matrix(EA_L, nodes, bars)
    C1, C2, Zmat = get_constraint_matrix(nodes, constrained)

    mat = scipy.sparse.bmat(
        [
            [Kmat, C2],
            [C1, Zmat],
        ], format='csc')
    return mat, Kmat


def solve(mat, rhs):
    lu = scipy.sparse.linalg.splu(mat)
    sol = lu.solve(rhs)
    return sol
