import unittest

import sys
sys.path.append( '..')
import collider

from test_shared import *
import test_shared 

class Test_collider_function(unittest.TestCase): 

    def setUp(self):
        self.transitions = transitions_def1
        self.collisions  = collisions_def1
        class Minimal: pass
        self.par = Minimal()
        self.par.q3_window = 4.0
        self.par.ppm = False
        self.MAX_UIS = 5

    def test_get_non_UIS_from_transitions1(self): 
        oldnon_uis = collider.get_non_UIS_from_transitions_old(self.transitions, self.collisions, 
                                         self.par, self.MAX_UIS)
        self.assertEqual([len(l) for l in oldnon_uis[1:]], test_shared.lennonuis1)
        self.assertEqual(oldnon_uis, test_shared.refnonuis1)

    def test_get_non_UIS_from_transitions2(self): 
        self.transitions = test_shared.transitions_def2
        self.collisions  = test_shared.collisions_def2
        oldnon_uis = collider.get_non_UIS_from_transitions_old(self.transitions,
            self.collisions, self.par, self.MAX_UIS)
        self.assertEqual([len(l) for l in oldnon_uis[1:]], test_shared.lennonuis2)
        self.assertEqual(oldnon_uis, test_shared.refnonuis2_sorted)

    def test_get_non_UIS_from_transitions2_unsorted(self): 
        #here we have the transitions in the wrong order
        #it should still work
        self.transitions = transitions_def2_unsorted
        self.collisions  = collisions_def2
        oldnon_uis = collider.get_non_UIS_from_transitions_old(self.transitions,
            self.collisions, self.par, self.MAX_UIS, unsorted=True)
        self.assertEqual([len(l) for l in oldnon_uis[1:]], test_shared.lennonuis2)
        self.assertEqual(oldnon_uis, test_shared.refnonuis2_unsorted)

    def test_get_non_UIS_from_transitions3(self): 
        self.transitions = transitions_def3
        self.collisions  = collisions_def3
        oldnon_uis = collider.get_non_UIS_from_transitions_old(self.transitions,
            self.collisions, self.par, self.MAX_UIS)
        self.assertEqual([len(l) for l in oldnon_uis[1:]], test_shared.lennonuis3)
        self.assertEqual(oldnon_uis, test_shared.refnonuis3)

    def test_get_non_UIS_from_transitions4(self): 
        self.transitions = transitions_def4
        self.collisions  = collisions_def4
        oldnon_uis = collider.get_non_UIS_from_transitions_old(self.transitions,
            self.collisions, self.par, self.MAX_UIS)
        self.assertEqual([len(l) for l in oldnon_uis[1:]], test_shared.lennonuis4)
        self.assertEqual(oldnon_uis, test_shared.refnonuis4)

    def test_get_uis(self):
        non_uis_list = collider.get_non_UIS_from_transitions_old(self.transitions,
                 self.collisions, self.par, self.MAX_UIS)
        srm_ids = [t[1] for t in self.transitions]
        rr = collider.get_uis(srm_ids, non_uis_list[2], 2)
        self.assertEqual(len(rr), 0)
        #
        self.transitions = transitions_def2
        self.collisions  = collisions_def2
        non_uis_list = collider.get_non_UIS_from_transitions_old(self.transitions,
                 self.collisions, self.par, self.MAX_UIS)
        srm_ids = [t[1] for t in self.transitions]
        rr = collider.get_uis(srm_ids, non_uis_list[2], 2)
        self.assertEqual(len(rr), 1)

    def test_get_UIS_from_transitions(self):
        res = collider.get_UIS_from_transitions(self.transitions, 
            self.collisions, self.par, self.MAX_UIS)
        self.assertEqual( [len(r) for r in res] , [0,0,0,0,0,0])

    def test_get_UIS_from_transitions2(self):
        self.transitions = transitions_def2
        self.collisions  = collisions_def2
        res = collider.get_UIS_from_transitions(self.transitions,
            self.collisions, self.par, self.MAX_UIS)
        self.assertEqual( [len(r) for r in res] , [0,0,1,2,1,0])

    def test_get_UIS_from_transitions3(self):
        self.transitions = transitions_def3
        self.collisions  = collisions_def3
        res = collider.get_UIS_from_transitions(self.transitions,
            self.collisions, self.par, self.MAX_UIS)
        self.assertEqual( [len(r) for r in res] , [0,0,0,0,0,0])

    def test_get_UIS_from_transitions4(self):
        self.transitions = transitions_def4
        self.collisions  = collisions_def4
        res = collider.get_UIS_from_transitions(self.transitions,
            self.collisions, self.par, self.MAX_UIS)
        self.assertEqual( [len(r) for r in res] , [0, 0, 7, 17, 15, 6])

if __name__ == '__main__':
    unittest.main()

