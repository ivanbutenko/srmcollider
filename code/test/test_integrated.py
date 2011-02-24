import unittest

import sys
sys.path.append( '..')
import collider

from test_shared import *
import test_shared 
import time


try:
    import c_integrated
except ImportError:
    print "=" * 75, """
Module c_integrated is not available. Please compile it if you want to use it.
""", "=" * 75

class Test_cintegrated(unittest.TestCase): 

    def setUp(self):
        #print "setup integrated"
        self.transitions = transitions_def1
        self.collisions  = collisions_def1
        class Minimal: pass
        self.par = Minimal()
        self.par.q3_window = 4.0
        self.par.ppm = False
        self.MAX_UIS = 5
        self.q3_high = 1500
        self.q3_low = 300

        def returnrange(): return self.q3_high, self.q3_low
        self.par.get_q3range_collisions = returnrange

        import sys
        sys.path.append( '/home/hroest/projects/' )
        sys.path.append( '/home/hroest/lib/' )
        import silver
        self.R = silver.Residues.Residues('mono')

        self.acollider = collider.SRMcollider()
        self.aparamset = collider.testcase()

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
                collider.SRMcollider(), precursors, R, q3_low, q3_high))
        m = self.acollider._getMinNeededTransitions(par, transitions, collisions)
        self.assertEqual(m, 8)

        m = c_integrated.getMinNeededTransitions(transitions, tuple(precursors), 
            par.max_uis, par.q3_window, par.ppm)
        self.assertEqual(m, 8)

        #now also test with lower q3 window
        par.q3_window = 1.0
        collisions = list(collider.SRMcollider._get_all_collisions_calculate_sub(
                collider.SRMcollider(), precursors, R, q3_low, q3_high))
        m = self.acollider._getMinNeededTransitions(par, transitions, collisions)
        self.assertEqual(m, 4)

        m = c_integrated.getMinNeededTransitions(transitions, tuple(precursors), 
            par.max_uis, par.q3_window, par.ppm)
        self.assertEqual(m, 4)


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
                collider.SRMcollider(), precursors, R, q3_low, q3_high))
        m = self.acollider._getMinNeededTransitions(par, transitions, collisions)
        self.assertEqual(m, 16)

        m = c_integrated.getMinNeededTransitions(transitions, tuple(precursors), 
            par.max_uis, par.q3_window, par.ppm)
        self.assertEqual(m, 16)

        #now also test with lower q3 window
        par.q3_window = 1.0
        collisions = list(collider.SRMcollider._get_all_collisions_calculate_sub(
                collider.SRMcollider(), precursors, R, q3_low, q3_high))
        m = self.acollider._getMinNeededTransitions(par, transitions, collisions)
        self.assertEqual(m, 6)

        m = c_integrated.getMinNeededTransitions(transitions, tuple(precursors), 
            par.max_uis, par.q3_window, par.ppm)
        self.assertEqual(m, 6)



if __name__ == '__main__':
    unittest.main()
