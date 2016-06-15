# OpenMDAO components for the continuous 3-D wing-truss optimization problem
# Dr. John T. Hwang
# Sicheng He
# June, 2016

from __future__ import division
import numpy
import scipy.sparse
import time
from scipy.sparse.linalg import gmres, LinearOperator, aslinearoperator

from openmdao.api import Component

import lib


class SysDispAug(Component):

    def __init__(self, nodes, elements, loads, cons, E):
        '''
        nodes : ndarray[num_nodes, 3], float
        elements : ndarray[num_elems, 2], int
        loads : ndarray[num_nodes, 3], float
        cons : ndarray[num_cons], int
        E : float
        '''
        super(SysDispAug, self).__init__()

        num_nodes = nodes.shape[0]
        num_elem = elements.shape[0]
        num_cons = cons.shape[0]

        self.add_param('areas', val=numpy.zeros(num_elem))
        self.add_state('disp_aug', val=numpy.zeros(3*num_nodes + 3*num_cons))

        self.nodes = nodes
        self.elements = elements
        self.loads = loads
        self.cons = cons
        self.rhs = numpy.zeros(3*num_nodes + 3*num_cons)
        self.rhs[:3*num_nodes] = loads
        self.E = E


        
        nodes = self.nodes
        elems = self.elements
        cons = self.cons
        E = self.E

        num_nodes = nodes.shape[0]
        num_elems = elems.shape[0]
        num_cons = cons.shape[0]

        nnz = 36 * num_elems + 2 * 3 * num_cons

        data, self.rows, self.cols = lib.getmtx(num_nodes, num_elems, num_cons, nnz,
                                      E, nodes, elems+1, numpy.ones(num_elems), cons+1)



        nodes = self.nodes
        elems = self.elements
        cons = self.cons
        E = self.E

        num_nodes = nodes.shape[0]
        num_elems = elems.shape[0]
        num_cons = cons.shape[0]
        nnz = 36 * num_elems

        out = lib.getresder2(num_nodes, num_elems, nnz, E, nodes, elems+1)
        self.data2, self.rows2, self.cols2, self.ind_aug2 = out
        
#        self.deriv_options['type'] = 'cs'
#        self.deriv_options['form'] = 'central'

    def get_stiffness_matrix(self, areas):
        nodes = self.nodes
        bars = self.elements
        E = self.E

        num_nodes = nodes.shape[0]
        num_elems = bars.shape[0]

        nnz = 36 * num_elems
        
        Kdata = numpy.zeros(36 * num_elems)
        Krows = numpy.zeros(36 * num_elems, int)
        Kcols = numpy.zeros(36 * num_elems, int)

        Telem = numpy.zeros((2, 6))
        Kelem = numpy.array([[1, -1], [-1, 1]])

        data_6x6 = numpy.zeros((6, 6))
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

            data_6x6[:, :] = Telem.T.dot(Kelem).dot(Telem) * E * areas[ielem] / L

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
            
        Kmat = scipy.sparse.csc_matrix((Kdata, (Krows, Kcols)),
                                       shape=(3 * num_nodes, 3 * num_nodes))

        return Kmat

    def get_constraint_matrix(self):
        nodes = self.nodes
        constrained = self.cons
        
        num = len(constrained)
        num_nodes = len(nodes)

        data = numpy.ones(3*num)
        lins = numpy.arange(3*num)
        cons = numpy.zeros(3*num, int)

        for ind in xrange(3):
            cons[ind::3] += 3 * constrained + ind

        C1 = scipy.sparse.csc_matrix((data, (lins, cons)),
                                     shape=(3*num, 3*num_nodes))
        C2 = scipy.sparse.csc_matrix((data, (cons, lins)),
                                     shape=(3*num_nodes, 3*num))
        zmat = scipy.sparse.csc_matrix((3*num, 3*num))

        return C1, C2, zmat

    def get_matrix(self, areas):
        t0 = time.time()
        if 0:
            Kmat = self.get_stiffness_matrix(areas)
            C1, C2, Zmat = self.get_constraint_matrix()

            mat = scipy.sparse.bmat(
                [
                    [Kmat, C2],
                    [C1, Zmat],
                ], format='csc')
        else:
            nodes = self.nodes
            elems = self.elements
            cons = self.cons
            E = self.E

            num_nodes = nodes.shape[0]
            num_elems = elems.shape[0]
            num_cons = cons.shape[0]

            nnz = 36 * num_elems + 2 * 3 * num_cons

            data = lib.getmtx2(num_nodes, num_elems, nnz,
                               E, nodes, elems+1, areas)

            size = 3 * num_nodes + 3 * num_cons
            mat = scipy.sparse.csc_matrix((data, (self.rows, self.cols)),
                                          shape=(size, size))

        t1 = time.time()

        return mat

    def apply_nonlinear(self, params, unknowns, resids):
        self.mat = self.get_matrix(params['areas'])

        resids['disp_aug'] = self.mat.dot(unknowns['disp_aug']) - self.rhs

    def callback(self, res):
        print self.counter, numpy.linalg.norm(res)
        self.counter += 1

    def gmres(self, mat, rhs, lu):
        if 1:
            size = len(rhs)
            A = aslinearoperator(mat)
            M = LinearOperator((size, size), dtype=float, matvec=lu.solve)
            self.counter = 0
            sol, info = gmres(A, rhs, M=M, maxiter=10,
                              #callback=self.callback,
                              tol=1e-12)
            return sol
        else:
            return lu.solve(rhs)

    def solve_nonlinear(self, params, unknowns, resids):
        self.mat = self.get_matrix(params['areas'])
        self.lu = scipy.sparse.linalg.splu(self.mat)

        #unknowns['disp_aug'] = self.lu.solve(self.rhs)
        unknowns['disp_aug'] = self.gmres(self.mat, self.rhs, self.lu)

    def linearize(self, params, unknowns, resids):
        self.lu_T = scipy.sparse.linalg.splu(self.mat.T)

        if 1:
            nodes = self.nodes
            elems = self.elements
            cons = self.cons
            disp_aug = unknowns['disp_aug']
            E = self.E
            
            num_nodes = nodes.shape[0]
            num_elems = elems.shape[0]
            num_cons = cons.shape[0]
            num_aug = disp_aug.shape[0]
            nnz = 36 * num_elems
            
            data, rows, cols = lib.getresder(num_nodes, num_elems, num_aug, nnz,
                                             E, nodes, elems+1, disp_aug)
            
            size = 3 * num_nodes + 3 * num_cons
            mat = scipy.sparse.csc_matrix((data, (rows, cols)),
                                          shape=(size, num_elems))
            
        else:
            nodes = self.nodes
            elems = self.elements
            cons = self.cons
            disp_aug = unknowns['disp_aug']
        
            num_nodes = nodes.shape[0]
            num_elems = elems.shape[0]
            num_cons = cons.shape[0]
        
            data = self.data2 * disp_aug[self.ind_aug2]

            size = 3 * num_nodes + 3 * num_cons
            mat = scipy.sparse.csc_matrix((data, (self.rows2, self.cols2)),
                                          shape=(size, num_elems))
        
        jac = {} # self.alloc_jacobian()
        jac['disp_aug', 'disp_aug'] = self.mat
        jac['disp_aug', 'areas'] = mat

        return jac

    def solve_linear(self, dumat, drmat, vois, mode=None):
        if mode == 'fwd':
            sol_vec, rhs_vec = self.dumat, self.drmat
            t = 0
        else:
            sol_vec, rhs_vec = self.drmat, self.dumat
            t = 1

        for voi in vois:
            if mode == 'fwd':
                sol_vec[voi].vec[:] = self.lu.solve(rhs_vec[voi].vec)
            else:
                #sol_vec[voi].vec[:] = self.lu_T.solve(
                sol_vec[voi].vec[:] = self.gmres(self.mat.T, rhs_vec[voi].vec,
                                                 self.lu_T)



class SysDisplacements(Component):
    """ Selects displacements from augmented vector """

    def __init__(self, nodes, cons):
        super(SysDisplacements, self).__init__()

        n = nodes.shape[0]
        size = 3 * n + 3 * cons.shape[0]
        self.n = n

        self.add_param('disp_aug', val=numpy.zeros((size)))
        self.add_output('disp', val=numpy.zeros((n, 3)))

        #self.deriv_options['type'] = 'cs'
        #self.deriv_options['form'] = 'central'
        #self.deriv_options['extra_check_partials_form'] = "central"

        data = numpy.ones(3 * n)
        lins = numpy.arange(3 * n)
        self.mat = scipy.sparse.csc_matrix((data, (lins, lins)),
                                           shape=(3*n, size))

    def solve_nonlinear(self, params, unknowns, resids):
        n = self.n
        unknowns['disp'] = numpy.array(params['disp_aug'][:3*n].reshape((n, 3)))

    def linearize(self, params, unknowns, resids):
        jac = {}
        jac['disp', 'disp_aug'] = self.mat
        return jac



class SysCompliance(Component):

    def __init__(self, nodes, loads):
        super(SysCompliance, self).__init__()

        self.n = nodes.shape[0]
        self.loads = loads

        self.add_param('disp', val=numpy.zeros((self.n, 3)))
        self.add_output('compliance', val=0.)

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['compliance'] = numpy.sum(params['disp'] * self.loads)

    def linearize(self, params, unknowns, resids):
        jac = {}
        jac['compliance', 'disp'] = self.loads.reshape((1, 3*self.n))
        return jac



class SysVolume(Component):

    def __init__(self, elems, lengths):
        super(SysVolume, self).__init__()

        self.n = elems.shape[0]
        self.lengths = lengths

        self.add_param('areas', val=numpy.zeros(self.n))
        self.add_output('volume', val=0.)

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['volume'] = numpy.sum(params['areas'] * self.lengths)

    def linearize(self, params, unknowns, resids):
        jac = {}
        jac['volume', 'areas'] = self.lengths.reshape((1, self.n))
        return jac




class SysStress(Component):

    def __init__(self, nodes, elems, E, s0):
        super(SysStress, self).__init__()
        
        self.nodes = nodes
        self.elems = elems
        self.E = E
        self.s0 = s0

        self.add_param('disp', val=numpy.zeros((nodes.shape[0], 3)))
        self.add_output('stress', val=numpy.zeros(elems.shape[0]))
        
        nodes = self.nodes
        elems = self.elems
        E = self.E
        s0 = self.s0

        num_nodes = nodes.shape[0]
        num_elems = elems.shape[0]
        nnz = 2 * 3 * num_elems
        
        data, rows, cols = lib.getstressder(num_nodes, num_elems, nnz,
                                            E, nodes, elems+1)
        data /= s0
        self.mat = scipy.sparse.csc_matrix((data, (rows, cols)),
                                           shape=(num_elems, 3 * num_nodes))

    def solve_nonlinear(self, params, unknowns, resids):
        nodes = self.nodes
        elems = self.elems
        E = self.E
        s0 = self.s0

        num_nodes = nodes.shape[0]
        num_elems = elems.shape[0]
        
        #unknowns['stress'] = lib.getstresses(num_nodes, num_elems,
        #                                     E, nodes, elems+1, params['disp']) / s0
        unknowns['stress'] = self.mat.dot(params['disp'].flatten())
        
    def linearize(self, params, unknowns, resids):
        '''
        nodes = self.nodes
        elems = self.elems
        E = self.E
        s0 = self.s0

        num_nodes = nodes.shape[0]
        num_elems = elems.shape[0]
        nnz = 2 * 3 * num_elems
        
        data, rows, cols = lib.getstressder(num_nodes, num_elems, nnz,
                                            E, nodes, elems+1)
        data /= s0
        mat = scipy.sparse.csc_matrix((data, (rows, cols)),
                                      shape=(num_elems, 3 * num_nodes))
        '''
        
        jacs = {}
        jacs['stress', 'disp'] = self.mat
        
        return jacs



class SysKS(Component):
    """ Aggregates failure constraints from the structure """

    def __init__(self, elems, s0, rho=10):
        super(SysKS, self).__init__()

        self.n = elems.shape[0]
        self.s0 = s0
        self.rho = rho

        self.add_param('stress', val=numpy.zeros(self.n))
        self.add_output('minstress', val=0.)
        self.add_output('maxstress', val=0.)

    def solve_nonlinear(self, params, unknowns, resids):
        s0 = self.s0
        rho = self.rho
        stress = params['stress']

        fmin = -s0 - stress
        fmax =  stress - s0

        maxfmin = numpy.max(fmin)
        maxfmax = numpy.max(fmax)

        unknowns['minstress'] = maxfmin + 1 / rho * \
                                numpy.log(numpy.sum(numpy.exp(rho*(fmin - maxfmin))))
        unknowns['maxstress'] = maxfmax + 1 / rho * \
                                numpy.log(numpy.sum(numpy.exp(rho*(fmax - maxfmax))))

    def linearize(self, params, unknowns, resids):
        n = self.n
        s0 = self.s0
        rho = self.rho
        stress = params['stress']

        fmin = -s0 - stress
        fmax =  stress - s0
        sgnmin = -1.0
        sgnmax =  1.0

        indmin = numpy.argmax(fmin)
        indmax = numpy.argmax(fmax)
        maxfmin = fmin[indmin]
        maxfmax = fmax[indmax]

        dmaxfmin_ds = numpy.zeros(self.n)
        dmaxfmax_ds = numpy.zeros(self.n)
        dmaxfmin_ds[indmin] = sgnmin * 1.0
        dmaxfmax_ds[indmax] = sgnmax * 1.0

        datamin = dmaxfmin_ds + 1/rho * 1.0 /\
                  numpy.sum(numpy.exp(rho * (fmin - maxfmin))) * \
                  numpy.exp(rho * (fmin - maxfmin)) * (sgnmin * rho)
        datamax = dmaxfmax_ds + 1/rho * 1.0 /\
                  numpy.sum(numpy.exp(rho * (fmax - maxfmax))) * \
                  numpy.exp(rho * (fmax - maxfmax)) * (sgnmax * rho) 
        datamin[indmin] -= 1/rho * sgnmin * rho
        datamax[indmax] -= 1/rho * sgnmax * rho

        jacs = {}
        jacs['minstress', 'stress'] = datamin.reshape((1, self.n))
        jacs['maxstress', 'stress'] = datamax.reshape((1, self.n))
        return jacs
        
        
