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
        newnon_uis = collider.get_non_UIS_from_transitions(self.transitions, self.collisions, 
                                              self.par, self.MAX_UIS)
        newnon_uis = [set( newn.keys() ) for newn in newnon_uis]
        oldnon_uis = collider.get_non_UIS_from_transitions_old(self.transitions, self.collisions, 
                                         self.par, self.MAX_UIS)
        for non_uis in [newnon_uis, oldnon_uis]:
            self.assertEqual(non_uis[1], set([(2,), (3,), (1,)]) )
            self.assertEqual(non_uis[2], set([(1, 2), (1, 3), (2, 3)]) )
            self.assertEqual(non_uis[3], set([(1, 2, 3),]) )
            self.assertEqual(non_uis[4], set([]) )
            self.assertEqual(non_uis[5], set([]) )

        self.assertEqual([len(l) for l in oldnon_uis[1:]], test_shared.lennonuis1)
        self.assertEqual(oldnon_uis, test_shared.refnonuis1)

        self.assertEqual([len(l) for l in newnon_uis[1:]], test_shared.lennonuis1)
        self.assertEqual(newnon_uis, test_shared.refnonuis1)
       
    def test_get_non_UIS_from_transitions2(self): 
        self.transitions = transitions_def2
        self.collisions  = collisions_def2
        newnon_uis = collider.get_non_UIS_from_transitions(self.transitions, self.collisions, 
                                              self.par, self.MAX_UIS)
        newnon_uis = [set( newn.keys() ) for newn in newnon_uis]
        oldnon_uis = collider.get_non_UIS_from_transitions_old(self.transitions, self.collisions, 
                                         self.par, self.MAX_UIS)
        for non_uis in [newnon_uis, oldnon_uis]:
            self.assertEqual(non_uis[1], set([(1,), (2,), (3,), (4,)]) )
            self.assertEqual(non_uis[2], set([(1, 2), (2, 3),  (1, 3),\
                                              (2, 4), (3, 4)]) )
            self.assertEqual(non_uis[3], set([(1, 2, 3), (2, 3, 4), ]) )
            self.assertEqual(non_uis[4], set([]) )
            self.assertEqual(non_uis[5], set([]) )

            self.assertEqual(len(non_uis[1]),4)
            self.assertEqual(len(non_uis[2]),5)
            self.assertEqual(len(non_uis[3]),2)
            self.assertEqual(len(non_uis[4]),0)
            self.assertEqual(len(non_uis[5]),0)

        self.assertEqual([len(l) for l in oldnon_uis[1:]], test_shared.lennonuis2)
        self.assertEqual(oldnon_uis, test_shared.refnonuis2_sorted)

        self.assertEqual([len(l) for l in newnon_uis[1:]], test_shared.lennonuis2)
        self.assertEqual(newnon_uis, test_shared.refnonuis2_sorted)

    def test_get_non_UIS_from_transitions2_unsorted(self): 
        #here we have the transitions in the wrong order
        #it should still work
        self.transitions = transitions_def2_unsorted
        self.collisions  = collisions_def2
        newnon_uis = collider.get_non_UIS_from_transitions(self.transitions, self.collisions, 
                                              self.par, self.MAX_UIS)
        newnon_uis = [set( newn.keys() ) for newn in newnon_uis]
        oldnon_uis = collider.get_non_UIS_from_transitions_old(self.transitions, self.collisions, 
                                         self.par, self.MAX_UIS, unsorted=True)
        for non_uis in [newnon_uis, oldnon_uis]:
            self.assertEqual(non_uis[1], set([(1,), (2,), (3,), (4,)]) )
            self.assertEqual(non_uis[2], set([(1, 2), (2, 3),  (1, 3),\
                                              (2, 4), (4, 3)]) )
            self.assertEqual(non_uis[3], set([(1, 2, 3), (2, 4, 3), ]) )
            self.assertEqual(non_uis[4], set([]) )
            self.assertEqual(non_uis[5], set([]) )

            self.assertEqual(len(non_uis[1]),4)
            self.assertEqual(len(non_uis[2]),5)
            self.assertEqual(len(non_uis[3]),2)
            self.assertEqual(len(non_uis[4]),0)
            self.assertEqual(len(non_uis[5]),0)

        self.assertEqual([len(l) for l in oldnon_uis[1:]], test_shared.lennonuis2)
        self.assertEqual(oldnon_uis, test_shared.refnonuis2_unsorted)

        self.assertEqual([len(l) for l in newnon_uis[1:]], test_shared.lennonuis2)
        self.assertEqual(newnon_uis, test_shared.refnonuis2_unsorted)

    def test_get_non_UIS_from_transitions3(self): 
        #here we have the transitions in the wrong order
        #it should still work
        self.transitions = transitions_def3
        self.collisions  = collisions_def3
        newnon_uis = collider.get_non_UIS_from_transitions(self.transitions, self.collisions, 
                                              self.par, self.MAX_UIS)
        newnon_uis = [set( newn.keys() ) for newn in newnon_uis]
        oldnon_uis = collider.get_non_UIS_from_transitions_old(self.transitions, self.collisions, 
                                         self.par, self.MAX_UIS)
        for non_uis in [newnon_uis, oldnon_uis]:
            self.assertEqual(non_uis[1], set([(1,), (2,), (3,), (4,), (5,), (6,)]) )
            self.assertEqual(non_uis[2], set([
                (5,6),
                (4,6),
                (3,6),
                (2,6),
                (1,6),
                (4,5),
                (3,5),
                (2,5),
                (1,5),
                (3,4),
                (2,4),
                (1,4),
                (2,3),
                (1,3),
                (1,2),
                ]) )
            self.assertEqual(non_uis[3], set([(3, 4, 6), (2, 3, 5), (1, 2, 6),
                (2, 5, 6), (4, 5, 6), (2, 3, 6), (1, 3, 6), (2, 4, 6), (1, 4,
                5), (1, 2, 5), (1, 2, 3), (1, 3, 5), (3, 5, 6), (2, 4, 5), (1,
                3, 4), (3, 4, 5), (1, 4, 6), (1, 5, 6), (1, 2, 4), (2, 3, 4)]) )
            self.assertEqual(non_uis[4], set([
                # present     #absent
                (1, 2, 3, 4), # 5,6
                (1, 2, 3, 5), # 4,6
                (1, 2, 4, 5), # 3,6
                (1, 3, 4, 5), # 2,6
                (2, 3, 4, 5), # 1,6
                (1, 2, 3, 6), # 4,5 
                (1, 2, 4, 6), # 3,5
                (1, 3, 4, 6), # 2,5
                (2, 3, 4, 6), # 1,5
                (1, 2, 5, 6), # 3,4
                (1, 3, 5, 6), # 2,4
                (2, 3, 5, 6), # 1,4
                (1, 4, 5, 6), # 2,3
                (2, 4, 5, 6), # 1,3
                (3, 4, 5, 6), # 1,2
                ]) )

            self.assertEqual(non_uis[5], set([(1, 2, 3, 4, 5),
                (1, 2, 3, 4, 6),
                (1, 2, 3, 5, 6),
                (1, 2, 4, 5, 6),
                (1, 3, 4, 5, 6),
                (2, 3, 4, 5, 6),
                ]) )

            self.assertEqual(len(non_uis[1]),6)
            self.assertEqual(len(non_uis[2]),15)
            self.assertEqual(len(non_uis[3]),20)
            self.assertEqual(len(non_uis[4]),15)
            self.assertEqual(len(non_uis[5]),6)

        self.assertEqual([len(l) for l in oldnon_uis[1:]], test_shared.lennonuis3)
        self.assertEqual(oldnon_uis, test_shared.refnonuis3)

        self.assertEqual([len(l) for l in newnon_uis[1:]], test_shared.lennonuis3)
        self.assertEqual(newnon_uis, test_shared.refnonuis3)


    def test_get_non_UIS_from_transitions4(self): 
        #here we have the transitions in the wrong order
        #it should still work
        self.transitions = transitions_def4
        self.collisions  = collisions_def4
        newnon_uis = collider.get_non_UIS_from_transitions(self.transitions, self.collisions, 
                                              self.par, self.MAX_UIS)
        newnon_uis = [set( newn.keys() ) for newn in newnon_uis]
        oldnon_uis = collider.get_non_UIS_from_transitions_old(self.transitions, self.collisions, 
                                         self.par, self.MAX_UIS)
        for non_uis in [newnon_uis, oldnon_uis]:
            self.assertEqual(non_uis[1], set([(1,), (2,), (3,), (4,), (5,), (6,)]) )
            self.assertEqual(non_uis[2], set([
                #203
                (5,6),
                (4,6),
                (4,5),
                #202
                (3,4),
                (2,4),
                #201
                (2,3),
                (1,3),
                (1,2),
                ]) )
            self.assertEqual(non_uis[3], set([(1, 2, 3), (2, 3, 4), (4, 5, 6) ]) )
            self.assertEqual(non_uis[4], set([]))
            self.assertEqual(non_uis[5], set([]))

            self.assertEqual(len(non_uis[1]),6)
            self.assertEqual(len(non_uis[2]),8)
            self.assertEqual(len(non_uis[3]),3)
            self.assertEqual(len(non_uis[4]),0)
            self.assertEqual(len(non_uis[5]),0)

        self.assertEqual([len(l) for l in oldnon_uis[1:]], test_shared.lennonuis4)
        self.assertEqual(oldnon_uis, test_shared.refnonuis4)

        self.assertEqual([len(l) for l in newnon_uis[1:]], test_shared.lennonuis4)
        self.assertEqual(newnon_uis, test_shared.refnonuis4)


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

