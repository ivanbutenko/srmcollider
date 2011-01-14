import random
import unittest

import collider
class Test_cgetnonuis(unittest.TestCase):

    def setUp(self):
        try:
            import c_getnonuis
            self.transitions = ( (500.0,1), 
                            (600,2), 
                            (700,3), 
                          )
            self.collisions = (  (500.4,400.0,101,201),
                            (600.6,400.0,102,201), 
                            (500.6,401.0,103,202),
                            (700.6,401.0,104,202), 
                            (500.7,400.0,105,203),
                            (600.7,400.0,106,203),
                            (700.7,400.0,107,203), 
                         )
        except ImportError:
            print """Module c_getnonuis is not available.
            
            Please compile it if you want to use it."""

    def test_getnonuis(self):

        try:
            import c_getnonuis
            #collisions
            #q3, q1, srm_id, peptide_key
            #transitions
            #q3, srm_id
            #
            q3window = 1.0
            ppm = False
            #
            result = c_getnonuis.getnonuis( self.transitions, self.collisions, q3window, ppm)
            self.assertTrue( result[201] == [1,2] )
            self.assertTrue( result[202] == [1,3] )
            self.assertTrue( result[203] == [1,2,3] )
            #Test 2
            transitions = ( (500.0,1), 
                            (600,2), 
                            (700,3), 
                            (800,4), 
                          )
            collisions = (  (500.4,400.0,101,201),
                            (600.6,400.0,102,201), 
                            (700.6,401.0,103,201),
                            (600.6,401.0,104,202), 
                            (700.7,400.0,105,202),
                            (800.7,400.0,106,202),
                         )
            result = c_getnonuis.getnonuis( transitions, collisions, q3window, ppm)
            self.assertTrue( result[201] == [1,2,3] )
            self.assertTrue( result[202] == [2,3,4] )
        except ImportError: pass

    def test_get_non_uis(self):
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

    def test_calculate_transitions(self):
        pass


class Test_crangetree(unittest.TestCase):

    def setUp(self):
        try:
            import c_rangetree
            self.parent_id = 101
            self.q1 = 501.0
            self.ssrcalc = 24
            self.mytuple1 = (
                ('PEPTIDE', 1, self.parent_id, 2, self.q1, self.ssrcalc),
            )
        except ImportError:
            print """Module c_rangetree is not available.

            Please compile it if you want to use it."""

    def test_rangetree(self):
        try:
            import c_rangetree
            c_rangetree.create_tree( self.mytuple1 )

            #we get our peptide out again with a large window
            res = c_rangetree.query_tree( self.q1 - 1, self.ssrcalc -1, 
                                         self.q1 + 1,  self.ssrcalc + 1 ) 
            self.assertEqual( len(res), 1)
            self.assertEqual( res[0][0], 101)

            #same result when lower boundary equals the value
            res = c_rangetree.query_tree( self.q1 , self.ssrcalc ,
                                         self.q1 + 1,  self.ssrcalc + 1 ) 
            self.assertEqual( len(res), 1)
            self.assertEqual( res[0][0], 101)

            #no result when upper boundary equals the value
            res = c_rangetree.query_tree( self.q1 - 1, self.ssrcalc -1, 
                                         self.q1,  self.ssrcalc ) 
            self.assertEqual( len(res), 0)

        except ImportError: pass

class Test_collider(unittest.TestCase):

    def setUp(self):

        import MySQLdb
        self.db = MySQLdb.connect(read_default_file="~/.my.cnf")

    def test_1(self):
        par  = collider.testcase()
        mycollider = collider.SRMcollider()
        mycollider.find_clashes_small(self.db, par) 
        print "I ran the testcase in %ss" % mycollider.total_time
        self.assertTrue( abs(sum( mycollider.allpeps.values() ) - 975.6326447245566) < 10**(-3) )
        self.assertEqual( mycollider.non_unique_count,  26 )
        self.assertEqual( mycollider.total_count, 12502 )
        self.assertTrue( abs(mycollider.allpeps[1585] - 0.93333) < 10**(-3) )

    def test_2(self):
        #verify that with toptrans=False we get the same results
        par  = collider.testcase()
        mycollider = collider.SRMcollider()
        mycollider.find_clashes(self.db, par, toptrans=False) 
        print "I ran the testcase in %ss" % mycollider.total_time
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
        print "I ran the testcase in %ss" % mycollider.total_time
        self.assertTrue( abs(sum( mycollider.allpeps.values() ) - 971.2154365242601) < 10**(-3) )
        self.assertEqual( mycollider.non_unique_count,  71 )
        self.assertEqual( mycollider.total_count,  12502 )
        self.assertTrue( abs(mycollider.allpeps[1585] - 0.93333) < 10**(-3) )
        #self.assertEqual( len( mycollider.q3min_distr ) == 71 )
        #self.assertEqual( len( mycollider.q1min_distr ) == 71 )
        self.assertEqual( len( mycollider.found3good ), 978 )


if __name__ == '__main__':
    unittest.main()


