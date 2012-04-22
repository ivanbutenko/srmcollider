"""
This file tests the functionality of the c_integrated module. 
"""

import unittest
from nose.plugins.attrib import attr
from nose.tools import nottest

import sys
sys.path.append( '..')
sys.path.append( '../external')
import collider

from test_shared import *
import test_shared 
import time
from Residues import Residues


try:
    import c_integrated
except ImportError:
    print "=" * 75, """
Module c_integrated is not available. Please compile it if you want to use it.
If you compiled it, please check the linking and compile errors.
""", "=" * 75

class Test_cintegrated(unittest.TestCase): 

    def setUp(self):
        self.transitions = transitions_def1
        self.collisions  = collisions_def1
        class Minimal: pass
        self.par = Minimal()
        self.par.q3_window = 4.0
        self.par.ppm = False
        self.MAX_UIS = 5
        self.q3_high = 1500
        self.q3_low = 300
        self.par.bions      =  True
        self.par.yions      =  True
        self.par.aions      =  False
        self.par.aMinusNH3  =  False
        self.par.bMinusH2O  =  False
        self.par.bMinusNH3  =  False
        self.par.bPlusH2O   =  False
        self.par.yMinusH2O  =  False
        self.par.yMinusNH3  =  False
        self.par.cions      =  False
        self.par.xions      =  False
        self.par.zions      =  False
        self.par.MMinusH2O  = False
        self.par.MMinusNH3  = False

        def returnrange(): return self.q3_high, self.q3_low
        self.par.get_q3range_collisions = returnrange

        import sys
        self.R = Residues('mono')

        self.acollider = collider.SRMcollider()
        self.aparamset = collider.testcase()

    @attr('slow') 
    def test_getMinNeededTransitions_1(self):
        pep = test_shared.runpep1
        transitions = test_shared.runtransitions1
        precursors = test_shared.runprecursors1
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        par = self.par
        q3_high = self.q3_high
        q3_low = self.q3_low
        R = self.R
        par.max_uis = 15

        collisions = list(collider.SRMcollider._get_all_collisions_calculate_sub(
                collider.SRMcollider(), precursors, par, R, q3_low, q3_high))
        m = self.acollider._getMinNeededTransitions(par, transitions, collisions)
        self.assertEqual(m, 8)

        m = c_integrated.getMinNeededTransitions(transitions, tuple(precursors), 
            par.max_uis, par.q3_window, par.ppm, par)
        self.assertEqual(m, 8)

        #now also test with lower q3 window
        par.q3_window = 1.0
        collisions = list(collider.SRMcollider._get_all_collisions_calculate_sub(
                collider.SRMcollider(), precursors, par, R, q3_low, q3_high))
        m = self.acollider._getMinNeededTransitions(par, transitions, collisions)
        self.assertEqual(m, 4)

        m = c_integrated.getMinNeededTransitions(transitions, tuple(precursors), 
            par.max_uis, par.q3_window, par.ppm, par)
        self.assertEqual(m, 4)

    @attr('slow') 
    def test_getMinNeededTransitions_2(self):
        pep = test_shared.runpep2
        transitions = test_shared.runtransitions2
        precursors = test_shared.runprecursors2
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        par = self.par
        q3_high = self.q3_high
        q3_low = self.q3_low
        R = self.R
        par.max_uis = len(transitions) +  1

        collisions = list(collider.SRMcollider._get_all_collisions_calculate_sub(
                collider.SRMcollider(), precursors, par, R, q3_low, q3_high))
        m = self.acollider._getMinNeededTransitions(par, transitions, collisions)
        #negative result, the transitions are not sufficient
        self.assertEqual(m, -1)

        m = c_integrated.getMinNeededTransitions(transitions, tuple(precursors), 
            par.max_uis, par.q3_window, par.ppm, par)
        # TODO fix this
        # self.assertEqual(m, 16)

        #now also test with lower q3 window
        par.q3_window = 1.0
        collisions = list(collider.SRMcollider._get_all_collisions_calculate_sub(
                collider.SRMcollider(), precursors, par, R, q3_low, q3_high))
        m = self.acollider._getMinNeededTransitions(par, transitions, collisions)
        self.assertEqual(m, 6)

        m = c_integrated.getMinNeededTransitions(transitions, tuple(precursors), 
            par.max_uis, par.q3_window, par.ppm, par)
        self.assertEqual(m, 6)

if __name__ == '__main__':
    unittest.main()

