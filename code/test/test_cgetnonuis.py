import unittest
import test_shared 
"""
This file tests the functionality of the c_getnonuis module. 
"""

import sys
sys.path.append( '..')
sys.path.append( '.')
import collider

from test_shared import ignoreImportError_rangetree, ignoreImportError_cget

try:
    import c_getnonuis
except ImportError:
    print "=" * 75, """
Module c_getnonuis is not available. Please compile it if you want to use it.
""", "=" * 75



def get_non_UIS_from_transitions(transitions, collisions, par, MAX_UIS, 
                                forceset=False):
    """ Get all combinations that are not UIS 
    
    Note that the new version returns a dictionary. To convert it to a set, one 
    needs to force the function to return a set.
    """
    import c_getnonuis
    non_uis_list = [{} for i in range(MAX_UIS+1)]
    collisions_per_peptide = c_getnonuis.getnonuis(
        transitions, collisions, par.q3_window, par.ppm)
    for order in range(1,MAX_UIS+1):
        non_uis_list[order] = c_getnonuis.get_non_uis(
            collisions_per_peptide, order)

    return non_uis_list



class Test_cgetnonuis(unittest.TestCase):
    def setUp(self):
            self.transitions = test_shared.transitions_def1
            self.collisions = test_shared.collisions_def1
            self.pep1 = test_shared.peptide1
            self.pep2 = test_shared.peptide2

            self.transitions_12_between300_1500 = test_shared.transitions_12_between300_1500
            self.pep1_yseries = test_shared.pep1_yseries
            self.pep1_bseries = test_shared.pep1_bseries

            tuples = []
            tuples.append(self.pep1)
            tuples.append(self.pep2)
            self.tuples = tuple( tuples)

            class Minimal: pass
            self.par = Minimal()
            self.par.q3_window = 4.0
            self.par.ppm = False
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


    def test_getnonuis(self):
            q3window = 1.0
            ppm = False
            #
            result = c_getnonuis.getnonuis( self.transitions, self.collisions, q3window, ppm)
            self.assertTrue( result[201] == [1,2] )
            self.assertTrue( result[202] == [1,3] )
            self.assertTrue( result[203] == [1,2,3] )
            self.assertEqual( result, test_shared.refcollperpep1)
            #Test 2
            transitions = test_shared.transitions_def2
            collisions = test_shared.collisions_def2
            result = c_getnonuis.getnonuis( transitions, collisions, q3window, ppm)
            self.assertTrue( result[201] == [1,2,3] )
            self.assertTrue( result[202] == [2,3,4] )

    def test_get_non_uis1(self):
        for order in range(1,6):
            res = c_getnonuis.get_non_uis(test_shared.refcollperpep1, order)
            res = set( res.keys() )
            self.assertEqual( res, test_shared.refnonuis1[order] )

    def test_get_non_uis2(self):
        for order in range(1,6):
            res = c_getnonuis.get_non_uis(test_shared.refcollperpep2, order)
            res = set( res.keys() )
            self.assertEqual( res, test_shared.refnonuis2_sorted[order] )

    def test_core_non_unique1(self):
            #collisions
            #q3, q1, srm_id, peptide_key
            #transitions
            #q3, srm_id
            #
            q3window = 1.0
            ppm = False
            result = c_getnonuis.core_non_unique( self.transitions, self.collisions, q3window, ppm)
            #
            self.assertTrue( abs(result[1] - 0.4) < 10**(-3) )
            self.assertTrue( abs(result[2] - 0.6) < 10**(-3) )
            self.assertTrue( abs(result[3] - 0.6) < 10**(-3) )

    def test_calculate_transitions_regular(self):
            trans = c_getnonuis.calculate_transitions( (self.pep1,), 300, 1500)
            self.assertEqual( len(trans), 10)
            for calc, ref in zip(trans, self.transitions_12_between300_1500):
                self.assertTrue(abs(calc[0] - ref) < 1e-3)

    def test_calculate_transitions_modifcation(self):
            trans = c_getnonuis.calculate_transitions( (self.pep2,), 300, 1500)
            #TODO check
            self.assertTrue(abs(trans[0][0]) - 909.333 < 1e-3)
            self.assertTrue(abs(trans[1][0]) - 780.290 < 1e-3)
            self.assertTrue(abs(trans[2][0]) - 683.238 < 1e-3)
            self.assertTrue(abs(trans[3][0]) - 523.207 < 1e-3)
            self.assertTrue(abs(trans[4][0]) - 410.123 < 1e-3)
            self.assertTrue(abs(trans[5][0]) - 330.112 < 1e-3)
            self.assertTrue(abs(trans[6][0]) - 490.143 < 1e-3)
            self.assertTrue(abs(trans[7][0]) - 603.227 < 1e-3)
            self.assertTrue(abs(trans[8][0]) - 718.254 < 1e-3)
            self.assertTrue(abs(trans[9][0]) - 865.289 < 1e-3)
            self.assertTrue(abs(trans[10][0]) - 455.170 < 1e-3) 
            self.assertTrue(abs(trans[11][0]) - 390.649 < 1e-3)
            self.assertTrue(abs(trans[12][0]) - 342.122 < 1e-3)
            self.assertTrue(abs(trans[13][0]) - 302.117 < 1e-3)
            self.assertTrue(abs(trans[14][0]) - 359.630 < 1e-3)
            self.assertTrue(abs(trans[15][0]) - 433.148 < 1e-3)

    def test_calculate_transitions_ch_regular(self):
            trans = c_getnonuis.calculate_transitions_ch( (self.pep1,), [1,2], 300, 1500)
            self.assertEqual( len(trans), 10)
            for calc, ref in zip(trans, self.transitions_12_between300_1500):
                self.assertTrue(abs(calc[0] - ref) < 1e-3)

    def test_calculate_transitions_inner(self):
            mypep = self.pep1

            transitions = c_getnonuis.calculate_transitions_inner(mypep, 2)
            self.assertEqual( len(transitions), 12)
            self.assertEqual( len(transitions), 2*len(mypep[1]) - 2 )

            y_series = transitions[:len(mypep[1])-1]
            b_series = transitions[len(mypep[1])-1:]
            self.assertEqual( len(y_series), 6)
            self.assertEqual( len(b_series), 6)

            for calc, ref in zip(y_series, self.pep1_yseries):
                self.assertTrue(abs(calc - ref) < 1e-3)
            for calc, ref in zip(b_series, self.pep1_bseries):
                self.assertTrue(abs(calc - ref) < 1e-3)

    def test_calculate_calculate_collisions_per_peptide_debug1(self):
            """
            Debug test, if there is something wrong rather not use the big ones.

            Is contained in the big test
            """
            pep = test_shared.runpep1
            transitions = test_shared.runtransitions1
            #precursors = test_shared.runprecursors1
            par = self.par
            q3_high = self.q3_high
            q3_low = self.q3_low

            precursors = ( 
                (449.72058221399999, 'SYVAWDR', 11498839L, 2),
            )
            transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
            collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide( 
                transitions, tuple(precursors), q3_low, q3_high, par.q3_window, par.ppm)
            #

            self.assertEqual(collisions_per_peptide,
                             { 11498839: [3]} )

    def test_calculate_calculate_collisions_per_peptide_debug2(self):
            """
            Debug test, if there is something wrong rather not use the big ones.

            Is contained in the big test
            """
            pep = test_shared.runpep1
            transitions = test_shared.runtransitions1
            #precursors = test_shared.runprecursors1
            par = self.par
            q3_high = self.q3_high
            q3_low = self.q3_low

            precursors = ( 
                (450.57777992699999, 'GPGPALAGEPAGSLR', 10682370L, 3),
            )
            transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
            collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide( 
                transitions, tuple(precursors), q3_low, q3_high, par.q3_window, par.ppm)
            #

            self.assertEqual(collisions_per_peptide,
                             {10682370: [0, 1, 2, 7, 8, 9, 11],
                             })

    def test_calculate_calculate_collisions_per_peptide_debug1and2(self):
            """
            Debug test, if there is something wrong rather not use the big ones.

            Is contained in the big test
            """
            pep = test_shared.runpep1
            transitions = test_shared.runtransitions1
            #precursors = test_shared.runprecursors1
            par = self.par
            q3_high = self.q3_high
            q3_low = self.q3_low

            precursors = ( 
                (450.57777992699999, 'GPGPALAGEPAGSLR', 10682370L, 3),
                (449.72058221399999, 'SYVAWDR', 11498839L, 2),
            )
            transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
            collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide( 
                transitions, tuple(precursors), q3_low, q3_high, par.q3_window, par.ppm)
            #

            self.assertEqual(collisions_per_peptide,
                             {10682370: [0, 1, 2, 7, 8, 9, 11],
                              11498839: [3]} )

    def test_calculate_calculate_collisions_per_peptide_1(self):
            pep = test_shared.runpep1
            transitions = test_shared.runtransitions1
            precursors = test_shared.runprecursors1
            par = self.par
            q3_high = self.q3_high
            q3_low = self.q3_low

            transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
            collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide( 
                transitions, tuple(precursors), q3_low, q3_high, par.q3_window, par.ppm)
            self.assertEqual(collisions_per_peptide, test_shared.collpepresult1)

    def test_calculate_calculate_collisions_per_peptide_2(self):
            pep = test_shared.runpep2
            transitions = test_shared.runtransitions2
            precursors = test_shared.runprecursors2
            par = self.par
            q3_high = self.q3_high
            q3_low = self.q3_low


            transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
            collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide( 
                transitions, tuple(precursors), q3_low, q3_high, par.q3_window, par.ppm)
            self.assertEqual(collisions_per_peptide, test_shared.collpepresult2)

    def test_calculate_calculate_collisions_per_peptide_other_ionseries_debug1(self):
            """
            Debug test, if there is something wrong rather not use the big ones.

            Is contained in the big test
            """
            pep = test_shared.runpep1
            transitions = test_shared.runtransitions1
            #precursors = test_shared.runprecursors1
            par = self.par
            q3_high = self.q3_high
            q3_low = self.q3_low

            precursors = ( 
                (449.72058221399999, 'SYVAWDR', 11498839L, 2),
            )
            transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
            collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide_other_ion_series( 
                transitions, tuple(precursors), q3_low, q3_high, par.q3_window, par.ppm, par)
            #

            self.assertEqual(collisions_per_peptide,
                             { 11498839: [3]} )

    def test_calculate_calculate_collisions_per_peptide_1_other(self):
            pep = test_shared.runpep1
            transitions = test_shared.runtransitions1
            precursors = test_shared.runprecursors1
            par = self.par
            q3_high = self.q3_high
            q3_low = self.q3_low

            transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
            collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide_other_ion_series( 
                transitions, tuple(precursors), q3_low, q3_high, par.q3_window, par.ppm, par)
            for key in collisions_per_peptide:
                self.assertEqual(collisions_per_peptide[key], test_shared.collpepresult1[key])
            self.assertEqual(len(collisions_per_peptide), len(test_shared.collpepresult1))
            self.assertEqual(collisions_per_peptide, test_shared.collpepresult1)

    def test_calculate_calculate_collisions_per_peptide_2_other(self):
            pep = test_shared.runpep2
            transitions = test_shared.runtransitions2
            precursors = test_shared.runprecursors2
            par = self.par
            q3_high = self.q3_high
            q3_low = self.q3_low

            transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
            collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide_other_ion_series( 
                transitions, tuple(precursors), q3_low, q3_high, par.q3_window, par.ppm, par)
            for key in collisions_per_peptide:
                self.assertEqual(collisions_per_peptide[key], test_shared.collpepresult2[key])
            self.assertEqual(len(collisions_per_peptide), len(test_shared.collpepresult2))
            self.assertEqual(collisions_per_peptide, test_shared.collpepresult2)

    def test_calculate_calculate_collisions_per_peptide_other_ionseries_debug_part1(self):
            """
            Debug test, if there is something wrong rather not use the big ones.

            Is contained in the big test
            """
            pep = test_shared.runpep1
            transitions = test_shared.runtransitions1
            #precursors = test_shared.runprecursors1
            par = self.par

            par.bions      =  True
            par.yions      =  True
            par.aions      =  True
            par.aMinusNH3  =  True
            par.bMinusH2O  =  False
            par.bMinusNH3  =  False
            par.bPlusH2O   =  False
            par.yMinusH2O  =  False
            par.yMinusNH3  =  False
            par.cions      =  False
            par.xions      =  True
            par.zions      =  True

            q3_high = self.q3_high
            q3_low = self.q3_low

            # there is a b3 350.17164  interfering with a y3 347.22949 close to transition 3
            # there is an a4 393.21384                      close to transition 8
            # there is an a4-NH3 393.21384 - 17 = 376.21384 close to transition 11
            # there is an x2 316.12575                      close to transition 9  
            # there is an z3 459.19925                      close to transition 2
            precursors = ( 
                (449.72058221399999, 'SYVAWDR', 11498839L, 2),
            )
            transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
            collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide_other_ion_series( 
                transitions, tuple(precursors), q3_low, q3_high, par.q3_window, par.ppm, par)
            #

            self.assertEqual(collisions_per_peptide,
                             { 11498839: [2, 3, 8, 9, 11]} )

    def test_calculate_calculate_collisions_per_peptide_other_ionseries_debug_part2(self):
            """
            Debug test, if there is something wrong rather not use the big ones.

            Is contained in the big test
            """
            pep = test_shared.runpep1
            transitions = test_shared.runtransitions1
            #precursors = test_shared.runprecursors1
            par = self.par

            par.bions      =  False #
            par.yions      =  False #
            par.aions      =  False #
            par.aMinusNH3  =  False #
            par.bMinusH2O  =  False
            par.bMinusNH3  =  False
            par.bPlusH2O   =  True
            par.yMinusH2O  =  True
            par.yMinusNH3  =  False
            par.cions      =  False
            par.xions      =  False #
            par.zions      =  False #

            q3_high = self.q3_high
            q3_low = self.q3_low

            precursors = ( 
                (449.72058221399999, 'SYVAWDR', 11498839L, 2),
            )
            transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
            collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide_other_ion_series( 
                transitions, tuple(precursors), q3_low, q3_high, par.q3_window, par.ppm, par)
            #

            self.assertEqual(collisions_per_peptide,
                             { 11498839: [1,2,4, 8, 9 ]} )

    def get_charged_mass(self):
            res = c_getnonuis.calculate_charged_mass( (0, 'VASHIPNLK', 1), 1)
            self.assertTrue( abs(res - 978.57368) < 1e-5)
            res = c_getnonuis.calculate_charged_mass( (0, 'VASHIPNLK', 1), 2)
            self.assertTrue( abs(res - 489.79078) < 1e-5)
            res = c_getnonuis.calculate_charged_mass( (0, 'VASHIPNLK', 1), 3)
            self.assertTrue( abs(res - 326.86314) < 1e-5)
            res = c_getnonuis.calculate_charged_mass( (0, 'VASHIPNLK', 1), 4)
            self.assertTrue( abs(res - 245.39932) < 1e-5)

class Test_cgetnonuis_get_non_UIS_from_transitions(unittest.TestCase):
    """ Tests the c_getnonuis module over the collider.

    By calling collider.get_non_UIS_from_transitions, we test the two functions

        * c_getnonuis.getnonuis 
        * c_getnonuis.get_non_uis

    in tandem.
    """

    def setUp(self):
        class Minimal: pass
        self.par = Minimal()
        self.par.q3_window = 4.0
        self.par.ppm = False
        self.MAX_UIS = 5

    def test_get_non_UIS_from_transitions1(self): 
            self.transitions = test_shared.transitions_def1
            self.collisions  = test_shared.collisions_def1
            newnon_uis = get_non_UIS_from_transitions(self.transitions, 
                self.collisions, self.par, self.MAX_UIS)
            newnon_uis = [set( newn.keys() ) for newn in newnon_uis]
            self.assertEqual([len(l) for l in newnon_uis[1:]], test_shared.lennonuis1)
            self.assertEqual(newnon_uis, test_shared.refnonuis1)

    def test_get_non_UIS_from_transitions2(self): 
            self.transitions = test_shared.transitions_def2
            self.collisions  = test_shared.collisions_def2
            newnon_uis = get_non_UIS_from_transitions(self.transitions, 
                self.collisions, self.par, self.MAX_UIS)
            newnon_uis = [set( newn.keys() ) for newn in newnon_uis]
            self.assertEqual([len(l) for l in newnon_uis[1:]], test_shared.lennonuis2)
            self.assertEqual(newnon_uis, test_shared.refnonuis2_sorted)

    def test_get_non_UIS_from_transitions2_unsorted(self): 
            #here we have the transitions in the wrong order
            #it should still work
            self.transitions = test_shared.transitions_def2_unsorted
            self.collisions  = test_shared.collisions_def2
            newnon_uis = get_non_UIS_from_transitions(self.transitions, 
                self.collisions, self.par, self.MAX_UIS)
            newnon_uis = [set( newn.keys() ) for newn in newnon_uis]
            self.assertEqual([len(l) for l in newnon_uis[1:]], test_shared.lennonuis2)
            self.assertEqual(newnon_uis, test_shared.refnonuis2_unsorted)

    def test_get_non_UIS_from_transitions3(self): 
            self.transitions = test_shared.transitions_def3
            self.collisions  = test_shared.collisions_def3
            newnon_uis = get_non_UIS_from_transitions(self.transitions, 
                self.collisions, self.par, self.MAX_UIS)
            newnon_uis = [set( newn.keys() ) for newn in newnon_uis]
            self.assertEqual([len(l) for l in newnon_uis[1:]], test_shared.lennonuis3)
            self.assertEqual(newnon_uis, test_shared.refnonuis3)

    def test_get_non_UIS_from_transitions4(self): 
            self.transitions = test_shared.transitions_def4
            self.collisions  = test_shared.collisions_def4
            newnon_uis = get_non_UIS_from_transitions(self.transitions, 
                self.collisions, self.par, self.MAX_UIS)
            newnon_uis = [set( newn.keys() ) for newn in newnon_uis]
            self.assertEqual([len(l) for l in newnon_uis[1:]], test_shared.lennonuis4)
            self.assertEqual(newnon_uis, test_shared.refnonuis4)



import inspect, types
for name, fn in inspect.getmembers(Test_cgetnonuis):
    if isinstance(fn, types.UnboundMethodType):
        setattr(Test_cgetnonuis, name, ignoreImportError_cget(fn))
for name, fn in inspect.getmembers(Test_cgetnonuis_get_non_UIS_from_transitions):
    if isinstance(fn, types.UnboundMethodType):
        setattr(Test_cgetnonuis_get_non_UIS_from_transitions, name, ignoreImportError_cget(fn))


if __name__ == '__main__':
    unittest.main()

