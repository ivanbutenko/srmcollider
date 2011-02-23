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
        self.q3_high = 1500
        self.q3_low = 300

        import sys
        sys.path.append( '/home/hroest/projects/' )
        sys.path.append( '/home/hroest/lib/' )
        import silver
        self.R = silver.Residues.Residues('mono')

        self.acollider = collider.SRMcollider()
        self.aparamset = collider.testcase()

    def test_getMinNeededTransitions(self):
        par = self.par
        par.max_uis = 5
        transitions = self.transitions
        collisions = self.collisions
        m = self.acollider._getMinNeededTransitions(par, transitions, collisions)
        self.assertEqual(m, 5)

        par.max_uis = 10 
        transitions = transitions_def3
        collisions = collisions_def3
        m = self.acollider._getMinNeededTransitions(par, transitions, collisions)
        self.assertEqual(m, 10)

        transitions = transitions_def4
        collisions = collisions_def4
        m = self.acollider._getMinNeededTransitions(par, transitions, collisions)
        self.assertEqual(m, 4)

        transitions = transitions_def5
        collisions = collisions_def5
        m = self.acollider._getMinNeededTransitions(par, transitions, collisions)
        self.assertEqual(m, 4)

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

    def test_calculate_collisions_1(self):
        pep = test_shared.runpep1
        transitions = test_shared.runtransitions1
        precursors = test_shared.runprecursors1
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        par = self.par
        q3_high = self.q3_high
        q3_low = self.q3_low
        R = self.R

        collisions = list(collider.SRMcollider._get_all_collisions_calculate_sub(
                collider.SRMcollider(), precursors, R, q3_low, q3_high))
        collisions_per_peptide = {}
        q3_window_used = par.q3_window
        for t in transitions:
            if par.ppm: q3_window_used = par.q3_window * 10**(-6) * t[0]
            for c in collisions:
                if abs( t[0] - c[0] ) <= q3_window_used:
                    #gets all collisions
                    if collisions_per_peptide.has_key(c[3]):
                        if not t[1] in collisions_per_peptide[c[3]]:
                            collisions_per_peptide[c[3]].append( t[1] )
                    else: 
                        collisions_per_peptide[c[3]] = [ t[1] ] ; 

        self.assertEqual(collisions_per_peptide, test_shared.collpepresult1)

    def test_calculate_collisions_2(self):
        pep = test_shared.runpep2
        transitions = test_shared.runtransitions2
        precursors = test_shared.runprecursors2
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        par = self.par
        q3_high = self.q3_high
        q3_low = self.q3_low
        R = self.R

        collisions = list(collider.SRMcollider._get_all_collisions_calculate_sub(
                collider.SRMcollider(), precursors, R, q3_low, q3_high))
        collisions_per_peptide = {}
        q3_window_used = par.q3_window
        for t in transitions:
            if par.ppm: q3_window_used = par.q3_window * 10**(-6) * t[0]
            for c in collisions:
                if abs( t[0] - c[0] ) <= q3_window_used:
                    #gets all collisions
                    if collisions_per_peptide.has_key(c[3]):
                        if not t[1] in collisions_per_peptide[c[3]]:
                            collisions_per_peptide[c[3]].append( t[1] )
                    else: 
                        collisions_per_peptide[c[3]] = [ t[1] ] ; 

        self.assertEqual(collisions_per_peptide, test_shared.collpepresult2)

    def test_get_uis_extra1(self):

        order =  1
        nonuis = set([(1627247L,), (1627240L,), (1627241L,), (1627242L,), (1627248L,), (1627243L,), (1627238L,), (1627249L,), (1627244L,), (1627239L,), (1627250L,), (1627251L,), (1627252L,)])
        srm_ids = [1627238L, 1627239L, 1627240L, 1627241L, 1627242L, 1627243L, 1627244L, 1627247L, 1627248L, 1627249L, 1627250L, 1627251L, 1627252L]
        uis_list = collider.get_uis(srm_ids, nonuis, order)
        self.assertEqual( len(uis_list), 0)







if __name__ == '__main__':
    unittest.main()

