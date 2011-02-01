import random
import unittest

import sys
sys.path.append( '..')
import collider

from test_shared import *
import test_shared 

class Test_collider_mysql(unittest.TestCase):

    def setUp(self):

        try:
            import MySQLdb
        except ImportError:
            print """Module MySQLdb not available.
            
            Please install it if you want to use it.
            Use the following command (on Ubuntu systems):
                sudo apt-get install python-mysqldb
            """

        self.db = MySQLdb.connect(read_default_file="~/.my.cnf")

    def test_1(self):

        par  = collider.testcase()
        mycollider = collider.SRMcollider()
        mycollider.find_clashes_small(self.db, par) 
        self.assertTrue( abs(sum( mycollider.allpeps.values() ) - 975.6326447245566) < 10**(-3) )
        self.assertEqual( len(mycollider.allpeps ), 978 )
        self.assertEqual( len([v for v in mycollider.allpeps.values() if v == 1.0 ] ), 960 )

        assert abs(mycollider.allpeps[111] - 0.888888888889) < 1e-5
        assert abs(mycollider.allpeps[191] - 0.7) < 1e-5
        assert abs(mycollider.allpeps[375] - 0.941176470588) < 1e-5
        assert abs(mycollider.allpeps[585] - 0.875) < 1e-5
        assert abs(mycollider.allpeps[609] - 0.888888888889) < 1e-5
        assert abs(mycollider.allpeps[929] - 0.857142857143) < 1e-5
        assert abs(mycollider.allpeps[1089] - 0.7) < 1e-5
        assert abs(mycollider.allpeps[1101] - 0.916666666667) < 1e-5
        assert abs(mycollider.allpeps[1177] - 0.928571428571) < 1e-5
        assert abs(mycollider.allpeps[1493] - 0.928571428571) < 1e-5
        assert abs(mycollider.allpeps[1537] - 0.857142857143) < 1e-5
        assert abs(mycollider.allpeps[1585] - 0.933333333333) < 1e-5
        assert abs(mycollider.allpeps[1663] - 0.928571428571) < 1e-5
        assert abs(mycollider.allpeps[1763] - 0.8125) < 1e-5
        assert abs(mycollider.allpeps[1925] - 0.928571428571) < 1e-5
        assert abs(mycollider.allpeps[1957] - 0.833333333333) < 1e-5
        assert abs(mycollider.allpeps[1967] - 0.857142857143) < 1e-5
        assert abs(mycollider.allpeps[1991] - 0.857142857143) < 1e-5

        self.assertEqual( mycollider.non_unique_count,  26 )
        self.assertEqual( mycollider.total_count, 12502 )
        self.assertTrue( abs(mycollider.allpeps[1585] - 0.93333) < 10**(-3) )

    def test_2(self):
        #verify that with toptrans=False we get the same results
        par  = collider.testcase()
        mycollider = collider.SRMcollider()
        mycollider.find_clashes(self.db, par, toptrans=False) 
        self.assertTrue (abs(sum( mycollider.allpeps.values() ) - 975.6326447245566) < 10**(-3) )
        self.assertEqual( mycollider.non_unique_count,  26 )
        self.assertEqual( mycollider.total_count, 12502 )
        self.assertTrue (abs(mycollider.allpeps[1585] - 0.93333) < 10**(-3) )
        #self.assertEqual( len( mycollider.q3min_distr ),  26)
        #self.assertEqual( len( mycollider.q1min_distr ),  26)
        self.assertEqual( len( mycollider.found3good ),  978)

    def test_3(self):
        #now test with isotopes enabled
        par  = collider.testcase()
        par.considerIsotopes = True
        par.eval()
        mycollider = collider.SRMcollider()
        #mycollider.find_clashes_small(db, par) 
        mycollider.find_clashes(self.db, par, toptrans=False) 
        self.assertTrue( abs(sum( mycollider.allpeps.values() ) - 971.2154365242601) < 10**(-3) )
        self.assertEqual( mycollider.non_unique_count,  71 )
        self.assertEqual( mycollider.total_count,  12502 )
        self.assertTrue( abs(mycollider.allpeps[1585] - 0.93333) < 10**(-3) )
        #self.assertEqual( len( mycollider.q3min_distr ) == 71 )
        #self.assertEqual( len( mycollider.q1min_distr ) == 71 )
        self.assertEqual( len( mycollider.found3good ), 978 )

class Test_collider_sqlite(unittest.TestCase):

    def setUp(self):

        try:
            import sqlite
            from sqlite import DatabaseError
        except ImportError:
            print """Module sqlite not available.
            
            Please install it if you want to use it.
            Use the following command (on Ubuntu systems):
                sudo apt-get install python-sqlite
            """

        try:
            #the database file must exist and the databases must be created
            self.db = sqlite.connect('/tmp/testdb')
            cursor = self.db.cursor()
            cursor.execute("select count(*) from srmPeptides_test")
            #print cursor.fetchall()
        except DatabaseError:
            print "=" * 75
            print """The sqlite database is not available.
            
            Please run the sqlite_setupdb.py first."""


    def test_1(self):

        par = collider.testcase()
        par.transition_table = 'srmTransitions_test'
        par.peptide_table = 'srmPeptides_test'
        cursor = self.db.cursor()
        #for historic reasons, we only select a subset of peptides
        exclude_pepids = [183, 203, 319, 345, 355, 365, 367, 385, 393, 425, 1227, 1233, 1297, 1299, 1303, 1305, 1307, 1309, 1311, 1509, 1681, 1683]
        mycollider = collider.SRMcollider()
        pepids = mycollider._get_unique_pepids(par, cursor, ignore_genomeoccurence=True)
        pepids = [p for p in pepids if p['parent_id'] not in exclude_pepids]

        mycollider.find_clashes_small(self.db, par, pepids=pepids) 

        self.assertTrue( abs(sum( mycollider.allpeps.values() ) - 975.6326447245566) < 10**(-3) )
        self.assertEqual( len([v for v in mycollider.allpeps.values() if v == 1.0 ] ), 960 )
        [(k,v) for k,v in mycollider.allpeps.iteritems() if v != 1.0 ] 


        assert abs(mycollider.allpeps[111] - 0.888888888889) < 1e-5
        assert abs(mycollider.allpeps[191] - 0.7) < 1e-5
        assert abs(mycollider.allpeps[375] - 0.941176470588) < 1e-5
        assert abs(mycollider.allpeps[585] - 0.875) < 1e-5
        assert abs(mycollider.allpeps[609] - 0.888888888889) < 1e-5
        assert abs(mycollider.allpeps[929] - 0.857142857143) < 1e-5
        assert abs(mycollider.allpeps[1089] - 0.7) < 1e-5
        assert abs(mycollider.allpeps[1101] - 0.916666666667) < 1e-5
        assert abs(mycollider.allpeps[1177] - 0.928571428571) < 1e-5
        assert abs(mycollider.allpeps[1493] - 0.928571428571) < 1e-5
        assert abs(mycollider.allpeps[1537] - 0.857142857143) < 1e-5
        assert abs(mycollider.allpeps[1585] - 0.933333333333) < 1e-5
        assert abs(mycollider.allpeps[1663] - 0.928571428571) < 1e-5
        assert abs(mycollider.allpeps[1763] - 0.8125) < 1e-5
        assert abs(mycollider.allpeps[1925] - 0.928571428571) < 1e-5
        assert abs(mycollider.allpeps[1957] - 0.833333333333) < 1e-5
        assert abs(mycollider.allpeps[1967] - 0.857142857143) < 1e-5
        assert abs(mycollider.allpeps[1991] - 0.857142857143) < 1e-5

        self.assertEqual( mycollider.non_unique_count,  26 )
        self.assertEqual( mycollider.total_count, 12502 )
        self.assertTrue( abs(mycollider.allpeps[1585] - 0.93333) < 10**(-3) )

if __name__ == '__main__':
    unittest.main()


