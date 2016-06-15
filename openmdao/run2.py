# Run script for the continuous 3-D wing-truss optimization problem
# Dr. John T. Hwang
# Sicheng He
# June, 2016

from __future__ import division
import numpy
import time

from openmdao.api import IndepVarComp, Problem, Group, ScipyOptimizer, SqliteRecorder, pyOptSparseDriver
from stiffness import SysDispAug, SysDisplacements, SysCompliance, SysVolume, SysStress, SysKS
from utils import setup_problem, writeBDF
from openmdao.devtools.partition_tree_n2 import view_tree



E = 100.0e9
s0 = 0.215e9

u1 = 0.005
u2 = 0.995

geom_file = '../CRM_AVL/wing_coarse.avl'
upp_file = '../airfoils/rae2822_upp.txt'
low_file = '../airfoils/rae2822_low.txt'
results_file = '../CRM_AVL/results_coarse.txt'

factor = 1

xyz, nodes, elements, cons, forces, forcesArray = setup_problem(u1, u2, geom_file, upp_file, low_file, results_file,factor)
forces /= 1.e1

root = Group()
root.add('comp_areas',
         IndepVarComp([('areas', 4.e-3 * numpy.ones(elements.shape[0]))]),
         promotes=['*'])
root.add('sys_disp_aug',
         SysDispAug(nodes, elements, forces.flatten(), cons, E),
         promotes=['*'])
root.add('sys_displacements',
         SysDisplacements(nodes, cons),
         promotes=['*'])
root.add('sys_compliance',
         SysCompliance(nodes, forces),
         promotes=['*'])
root.add('sys_volume',
         SysVolume(elements, numpy.ones(len(elements))),
         promotes=['*'])
root.add('sys_stress',
         SysStress(nodes, elements, E, s0),
         promotes=['*'])
root.add('sys_ks',
         SysKS(elements, 1.0),
         promotes=['*'])

prob = Problem()
prob.root = root
#prob.root.deriv_options['type'] = 'fd'
prob.setup()

t0 = time.time()
prob.run()
t1 = time.time()

#print t1-t0

nodes0 = nodes
nodes1 = nodes + prob['disp']

writeBDF('jig.bdf', nodes0, elements+1)
writeBDF('deflected.bdf', nodes1, elements+1)

if 0:
    prob.check_partial_derivatives(compact_print=True)
    exit()

prob.driver = pyOptSparseDriver()
prob.driver.options['optimizer'] = "SNOPT"
prob.driver.opt_settings = {'Major optimality tolerance': 1.0e-7,
                            'Major feasibility tolerance': 1.0e-7,
                            'Iterations limit': int(1e6),
}

prob.driver.add_desvar('areas',lower=0.0001, upper=0.1, scaler=1e0) # test
prob.driver.add_objective('volume', scaler=1e0)
prob.driver.add_constraint('minstress', upper=0.)
prob.driver.add_constraint('maxstress', upper=0.)
# setup data recording
prob.driver.add_recorder(SqliteRecorder('data.db'))
prob.setup()
#view_tree(prob, outfile="aerostruct.html", show_browser=True)
prob.run()

writeBDF('optimized.bdf', nodes+prob['disp'], elements+1)
