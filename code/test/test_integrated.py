"""
This file tests the functionality of the c_integrated module. 
"""

import unittest
from nose.plugins.attrib import attr
from nose.tools import nottest

import sys
sys.path.append( '.')
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

    def test_integrated(self):
        pep = test_shared.runpep1
        transitions = test_shared.runtransitions1
        precursors = test_shared.runprecursors1
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        par = self.par
        q3_high = self.q3_high
        q3_low = self.q3_low
        R = self.R
        par.max_uis = 15

        alltuples = [ (p[1], p[2], p[2], p[3], p[0], 25, -1,-1, 0) for p in precursors]

        import c_rangetree
        r = c_rangetree.ExtendedRangetree_Q1_RT.create()
        r.new_rangetree()
        r.create_tree(tuple(alltuples))

        c_integrated.create_tree(tuple(alltuples))

        q1 = 450
        par.q1_window = 5
        par.isotopes_up_to = 2
        isotope_correction = 1
        ssrcalc_low = 0
        ssrcalc_high = 100
        MAX_UIS=5
        result = c_integrated.wrap_all_magic(transitions, q1 - par.q1_window, 
            ssrcalc_low, q1 + par.q1_window,  ssrcalc_high, -1,
            MAX_UIS, par.q3_window, par.ppm, par.isotopes_up_to, isotope_correction, par)
        self.assertEqual(result, [12, 66, 220, 495, 790])

        q1 = 450
        par.q1_window = 0.1
        par.q3_window = 0.1
        par.isotopes_up_to = 2
        isotope_correction = 1
        ssrcalc_low = 0
        ssrcalc_high = 100
        MAX_UIS=5
        result = c_integrated.wrap_all_magic(transitions, q1 - par.q1_window, 
            ssrcalc_low, q1 + par.q1_window,  ssrcalc_high, -1,
            MAX_UIS, par.q3_window, par.ppm, par.isotopes_up_to, isotope_correction, par)
        self.assertEqual(result, [12, 35, 20, 3, 0] )

        q1 = 450
        par.q1_window = 0.08
        par.q3_window = 0.01
        par.isotopes_up_to = 2
        isotope_correction = 1
        ssrcalc_low = 0
        ssrcalc_high = 100
        MAX_UIS=5
        result = c_integrated.wrap_all_magic(transitions, q1 - par.q1_window, 
            ssrcalc_low, q1 + par.q1_window,  ssrcalc_high, -1,
            MAX_UIS, par.q3_window, par.ppm, par.isotopes_up_to, isotope_correction, par)
        self.assertEqual(result, [10, 10, 5, 1, 0] )


if __name__ == '__main__':
    unittest.main()

