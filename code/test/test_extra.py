import random
import unittest

import sys
sys.path.append( '..')
sys.path.append( '../external')
import collider
from Residues import Residues

import test_shared 

class Test_fragmentation(unittest.TestCase):

    def setUp(self):

        import sys, os, time

        import DDB 
        self.DDB = DDB
        self.R = Residues('mono')

        self.pep1 = (500, 'PEPTIDE', 1, 1)
        self.pep2 = (400, 'CEPC[160]IDM[147]E',2,2)
        self.charge = 2

        self.q3_low = 300
        self.q3_high = 1500

        self.transitions_12_between300_1500 = test_shared.transitions_12_between300_1500
        self.pep1_yseries = test_shared.pep1_yseries
        self.pep1_bseries = test_shared.pep1_bseries

    def test_fragmentation_ddbpep(self):
        mypep = self.pep1

        R = self.R
        peptide = self.DDB.Peptide()
        peptide.set_sequence(mypep[1])
        peptide.charge = self.charge
        peptide.create_fragmentation_pattern(R)

        ch = peptide.charge
        y_series = [( pred + (ch -1)*R.mass_H)/ch for pred in peptide.y_series]
        b_series = [( pred + (ch -1)*R.mass_H)/ch for pred in peptide.b_series]

        self.assertEqual( len(y_series), 6)
        self.assertEqual( len(b_series), 6)
        for calc, ref in zip(y_series, self.pep1_yseries):
            self.assertTrue(abs(calc - ref) < 1e-3)
        for calc, ref in zip(b_series, self.pep1_bseries):
            self.assertTrue(abs(calc - ref) < 1e-3)

    def test_fragmentation_ddbpep2(self):
        R = self.R
        peptide = self.DDB.Peptide()
        peptide.set_sequence(self.pep1[1])
        peptide.charge = self.charge
        peptide.create_fragmentation_pattern(R)

        b_series = peptide.b_series
        y_series = peptide.y_series
        res = []
        ch = self.charge
        q1 = self.pep1[0]
        peptide_key = 445
        q3_low = self.q3_low
        q3_high = self.q3_high

        for ch in [1,2]:
            for pred in y_series:
                q3 = ( pred + (ch -1)*R.mass_H)/ch
                if q3 < q3_low or q3 > q3_high: continue
                res.append( (q3, 8, 0, 8) )
            for pred in b_series:
                q3 = ( pred + (ch -1)*R.mass_H)/ch
                if q3 < q3_low or q3 > q3_high: continue
                res.append( (q3, 8, 0, 8) )



        self.assertEqual( len(res), 10)
        for calc, ref in zip(res, self.transitions_12_between300_1500):
            self.assertTrue(abs(calc[0] - ref) < 1e-3)

if __name__ == '__main__':
    unittest.main()
