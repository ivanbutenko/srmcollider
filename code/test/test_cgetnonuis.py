import unittest
import test_shared 

import sys
sys.path.append( '..')

try:
    import c_getnonuis
except ImportError:
    print "=" * 75, """
Module c_getnonuis is not available. Please compile it if you want to use it.
""", "=" * 75

class Test_cgetnonuis(unittest.TestCase):

    def setUp(self):
        try:
            import c_getnonuis
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

            #self.pep1 = (500, 'PEPTIDE', 1, 1)
            #self.pep2 = (400, 'CEPC[160]IDM[147]E',2,2)
            #for peptide 1
            #self.transitions_12_between300_1500 = [703.31500971599996,
            #                                       574.27241971600006,
            #                                       477.21965971600002,
            #                                       376.17197971600001,
            #                                       324.155935032,
            #                                       425.20361503200002,
            #                                       538.28767503200004,
            #                                       653.31461503200001,
            #                                       352.161417374,
            #                                       327.16122003200002]
            ##this is the double charged b and y series (peptide 1)
            #self.pep1_yseries = [352.161417374, 287.64012237400004,
            #                     239.11374237400003, 188.58990237400002,
            #                     132.04787237400001, 74.53440237400001]
            #self.pep1_bseries = [49.534205032000003, 114.055500032,
            #                     162.58188003200002, 213.10572003200002,
            #                     269.64775003200003, 327.16122003200002]

        except ImportError:
            pass

    def test_getnonuis(self):

        try:
            import c_getnonuis
            q3window = 1.0
            ppm = False
            #
            result = c_getnonuis.getnonuis( self.transitions, self.collisions, q3window, ppm)
            self.assertTrue( result[201] == [1,2] )
            self.assertTrue( result[202] == [1,3] )
            self.assertTrue( result[203] == [1,2,3] )
            #Test 2
            transitions = test_shared.transitions_def2
            collisions = test_shared.collisions_def2
            result = c_getnonuis.getnonuis( transitions, collisions, q3window, ppm)
            self.assertTrue( result[201] == [1,2,3] )
            self.assertTrue( result[202] == [2,3,4] )
        except ImportError: pass

    def test_get_non_uis(self):
        #non_uis_list[order] = c_getnonuis.get_non_uis(
        #                            collisions_per_peptide, order)

        pass

    def test_core_non_unique1(self):
    
        try:
            import c_getnonuis
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
        except ImportError: pass

    def test_calculate_transitions_regular(self):
        try:
            import c_getnonuis
            trans = c_getnonuis.calculate_transitions( (self.pep1,), 300, 1500)
            self.assertEqual( len(trans), 10)
            for calc, ref in zip(trans, self.transitions_12_between300_1500):
                self.assertTrue(abs(calc[0] - ref) < 1e-3)
        except ImportError: pass

    def test_calculate_transitions_modifcation(self):
        try:
            import c_getnonuis
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
        except ImportError: pass

    def test_calculate_transitions_ch_regular(self):
        try:
            import c_getnonuis
            trans = c_getnonuis.calculate_transitions_ch( (self.pep1,), [1,2], 300, 1500)
            self.assertEqual( len(trans), 10)
            for calc, ref in zip(trans, self.transitions_12_between300_1500):
                self.assertTrue(abs(calc[0] - ref) < 1e-3)
        except ImportError: pass

    def test_calculate_transitions_inner(self):
        try:
            import c_getnonuis
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
        except ImportError: pass


if __name__ == '__main__':
    unittest.main()

