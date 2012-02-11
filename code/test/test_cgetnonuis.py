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
from test_shared import get_non_UIS_from_transitions, getnonuis

try:
    import c_getnonuis
except ImportError:
    print "=" * 75, """
Module c_getnonuis is not available. Please compile it if you want to use it.
""", "=" * 75

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
            self.par.MMinusH2O      =  False
            self.par.MMinusNH3      =  False


    def test_getnonuis(self):
            q3window = 1.0
            ppm = False
            #
            result = getnonuis(self.transitions, self.collisions, q3window, ppm)
            self.assertTrue( result[201] == [1,2] )
            self.assertTrue( result[202] == [1,3] )
            self.assertTrue( result[203] == [1,2,3] )
            self.assertEqual( result, test_shared.refcollperpep1)
            #Test 2
            transitions = test_shared.transitions_def2
            collisions = test_shared.collisions_def2
            result = getnonuis(transitions, collisions, q3window, ppm)
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

    def test_calculate_transitions_regular(self):
            trans = c_getnonuis.calculate_transitions_ch( (self.pep1,), [1,2], 300, 1500)
            self.assertEqual( len(trans), 10)
            for calc, ref in zip(trans, self.transitions_12_between300_1500):
                self.assertTrue(abs(calc[0] - ref) < 1e-3)

    def test_calculate_transitions_modifcation(self):
            trans = c_getnonuis.calculate_transitions_ch( (self.pep2,), [1,2], 300, 1500)
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

            transitions = c_getnonuis.calculate_transitions_ch( (mypep,), [2] , 0, 5000)
            transitions = [t[0] for t in transitions]
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
            #precursors = [Precursor(modified_sequence='SYVAWDR', q1=449.720582, q1_charge=2, transition_group = 11498839L, isotopically_modified=0) ]

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

            precursors = [Precursor(modified_sequence='SYVAWDR',
                transition_group=11498839L, isotopically_modified=0)] 
            transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
            collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide_other_ion_series( 
                transitions, precursors, par, q3_low, q3_high, par.q3_window, par.ppm, False)
            #

            self.assertEqual(collisions_per_peptide,
                             { 11498839: [3]} )

    def test_calculate_calculate_collisions_per_peptide_1_other(self):
            pep = test_shared.runpep1
            transitions = test_shared.runtransitions1
            precursors = test_shared.runprecursors_obj1
            par = self.par
            q3_high = self.q3_high
            q3_low = self.q3_low

            transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
            collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide_other_ion_series( 
                transitions, precursors, par, q3_low, q3_high, par.q3_window, par.ppm, False)
            for key in collisions_per_peptide:
                self.assertEqual(collisions_per_peptide[key], test_shared.collpepresult1[key])
            self.assertEqual(len(collisions_per_peptide), len(test_shared.collpepresult1))
            self.assertEqual(collisions_per_peptide, test_shared.collpepresult1)

    def test_calculate_calculate_collisions_per_peptide_2_other(self):
            pep = test_shared.runpep2
            transitions = test_shared.runtransitions2
            precursors = test_shared.runprecursors_obj2
            par = self.par
            q3_high = self.q3_high
            q3_low = self.q3_low

            transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
            collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide_other_ion_series( 
                transitions, precursors, par, q3_low, q3_high, par.q3_window, par.ppm, False)
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
            precursors = [Precursor(modified_sequence='SYVAWDR',
                transition_group=11498839L, isotopically_modified=0)] 
            transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
            collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide_other_ion_series( 
                transitions, precursors, par, q3_low, q3_high, par.q3_window, par.ppm, False)
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
            par.MMinusH2O      =  False
            par.MMinusNH3      =  False

            q3_high = self.q3_high
            q3_low = self.q3_low

            precursors = [Precursor(modified_sequence='SYVAWDR',
                transition_group=11498839L, isotopically_modified=0)] 
            transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
            collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide_other_ion_series( 
                transitions, precursors, par, q3_low, q3_high, par.q3_window, par.ppm, False)
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

    By calling collider.get_non_UIS_from_transitions, we test the function

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

from precursor import Precursor
import c_getnonuis
class Test_three_peptide_example(unittest.TestCase): 

    """
    The target is YYLLDYR with these transitions and numbers

      (842.4412197, 0), y6+
      (679.3778897, 1), y5+
      (566.2938297, 2), y4+
      (453.2097697, 3), y3+
      (440.2185450, 4), b3+
      (553.3026050, 5), b4+
      (668.3295450, 6), b5+
      (831.3928750, 7)  b6+ 

    The peptides GGLIVELGDK b5+ ion interferes with the targets b3+ ion which leads to 665: [4]

    The peptides NGTDGGLQVAIDAMR b9+ ion (842.4008) interferes with the targets y6+ ion
    and also the y11++ ion (565.8035) interferes with the targets y4+ ion which leads to 618: [0, 2].

    """
    def setUp(self):
      import sys

      self.acollider = collider.SRMcollider()
      self.aparamset = collider.testcase()
      self.EPSILON = 10**-5

      par = collider.SRM_parameters()
      par.q1_window = 25 / 2.0
      par.q3_window = 1 / 2.0
      par.ppm = False
      par.q3_low = 400
      par.q3_high = 1400

      par.q3_range = [par.q3_low, par.q3_high]
      par.set_default_vars()
      par.eval()
      self.real_parameters = par

      self.precursor = test_shared.ThreePeptideExample.precursor
      self.interfering_precursors = test_shared.ThreePeptideExample.interfering_precursors
      self.oldstyle_precursors = tuple([(p.q1, p.modified_sequence, p.transition_group, p.q1_charge, p.isotopically_modified) for p in self.interfering_precursors])

    def test_find_clashes_forall_other_series(self):
      """ Test how to calculate the transitions of the target
        nonunique = c_getnonuis._find_clashes_forall_other_series( 
            tuple(transitions), tuple(precursors), q3_low, q3_high, 
            par.q3_window, par.ppm, par, q1 - par.q1_window)

      """
      par = self.real_parameters
      q3_low, q3_high = self.real_parameters.get_q3range_transitions()
      precursor = self.precursor
      transitions = precursor.calculate_transitions(q3_low, q3_high)

    
      nonunique = c_getnonuis._find_clashes_forall_other_series( 
        tuple(transitions), self.interfering_precursors, par, q3_low, q3_high,
            par.q3_window, par.ppm, precursor.q1 - par.q1_window, False)

      self.assertEqual( len( nonunique ), 3)
      self.assertEqual( nonunique.keys(), [0,2,4] )

      self.assertTrue( abs(nonunique[0][0][0] - 842.4008) < self.EPSILON )
      self.assertEqual( nonunique[0][0][4], 'b')
      self.assertEqual( nonunique[0][0][5], 9)
      self.assertEqual( nonunique[0][0][-1], 1)

      self.assertTrue( abs(nonunique[2][0][0] - 565.8035) < self.EPSILON )
      self.assertEqual( nonunique[2][0][4], 'y')
      self.assertEqual( nonunique[2][0][5], 4)
      self.assertEqual( nonunique[2][0][-1], 2)

      self.assertTrue( abs(nonunique[4][0][0] - 440.287275) < self.EPSILON )
      self.assertEqual( nonunique[4][0][4], 'b')
      self.assertEqual( nonunique[4][0][5], 5)
      self.assertEqual( nonunique[4][0][-1], 1)

if __name__ == '__main__':
    unittest.main()

