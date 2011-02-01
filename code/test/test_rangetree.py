import random
import unittest


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

if __name__ == '__main__':
    unittest.main()
