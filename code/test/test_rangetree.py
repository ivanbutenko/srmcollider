import random
import unittest
"""
This file tests the functionality of the c_rangetree module. 
"""

import sys
sys.path.extend(['.', '..', '../external/', 'external/'])

from test_shared import ignoreImportError_rangetree
from test_shared import check_crangetree_availability

try:
    import c_rangetree
except ImportError:
    print "=" * 75, """
Module c_rangetree is not available. Please compile it if you want to use it.
""", "=" * 75

@check_crangetree_availability
class Test_crangetree(unittest.TestCase):

    def setUp(self):
            self.parent_id = 101
            self.q1 = 501.0
            self.ssrcalc = 24
            self.mytuple1 = (
                ('PEPTIDE', 1, self.parent_id, 2, self.q1, self.ssrcalc),
            )
            self.mytuple2 = (
                ('PEPTIDE', 1, self.parent_id+1, 2, self.q1 + 10, self.ssrcalc),
            )

    def test_rangetree(self):
            c_rangetree.create_tree( self.mytuple1 )

            #we get our peptide out again with a large window
            res = c_rangetree.query_tree( self.q1 - 1, self.ssrcalc -1, 
                                         self.q1 + 1,  self.ssrcalc + 1, 1, 0) 
            self.assertEqual( len(res), 1)
            self.assertEqual( res[0][0], 101)

            #same result when lower boundary equals the value
            res = c_rangetree.query_tree( self.q1 , self.ssrcalc ,
                                         self.q1 + 1,  self.ssrcalc + 1, 1, 0) 
            self.assertEqual( len(res), 1)
            self.assertEqual( res[0][0], 101)

            #no result when upper boundary equals the value
            res = c_rangetree.query_tree( self.q1 - 1, self.ssrcalc -1, 
                                         self.q1,  self.ssrcalc, 1, 0) 
            self.assertEqual( len(res), 0)

    def test_rangetree_object_empty(self):
            mytree = c_rangetree.Rangetree_Q1_RT.create()
            mytree.new_rangetree()

            #we can create a new tree that is empty
            mytree.new_rangetree()
            res = mytree.query_tree( self.q1 - 1, self.ssrcalc -1, 
                                         self.q1 + 1,  self.ssrcalc + 1, 1, 0) 
            self.assertEqual( len(res), 0)

    def test_rangetree_object(self):
            mytree = c_rangetree.Rangetree_Q1_RT.create()
            mytree.new_rangetree()
            mytree.create_tree( self.mytuple1 )

            #we get our peptide out again with a large window
            res = mytree.query_tree( self.q1 - 1, self.ssrcalc -1, 
                                         self.q1 + 1,  self.ssrcalc + 1, 1, 0) 
            self.assertEqual( len(res), 1)
            self.assertEqual( res[0][0], 101)

            #same result when lower boundary equals the value
            res = mytree.query_tree( self.q1 , self.ssrcalc ,
                                         self.q1 + 1,  self.ssrcalc + 1, 1, 0) 
            self.assertEqual( len(res), 1)
            self.assertEqual( res[0][0], 101)

            #no result when upper boundary equals the value
            res = mytree.query_tree( self.q1 - 1, self.ssrcalc -1, 
                                         self.q1,  self.ssrcalc, 1, 0) 
            self.assertEqual( len(res), 0)

            #we can create a new tree that is empty
            mytree.new_rangetree()
            res = mytree.query_tree( self.q1 - 1, self.ssrcalc -1, 
                                         self.q1 + 1,  self.ssrcalc + 1, 1, 0) 
            self.assertEqual( len(res), 0)

            mytree2 = c_rangetree.Rangetree_Q1_RT.create()
            mytree2.new_rangetree()
            res = mytree2.query_tree( self.q1 - 1, self.ssrcalc -1, 
                                         self.q1 + 1,  self.ssrcalc + 1, 1, 0) 
            self.assertEqual( len(res), 0)

    def test_rangetree_object_two_trees(self):
            mytree = c_rangetree.Rangetree_Q1_RT.create()
            mytree.new_rangetree()
            mytree.create_tree( self.mytuple1 )

            mytree2 = c_rangetree.Rangetree_Q1_RT.create()
            mytree2.new_rangetree()
            mytree2.create_tree( self.mytuple2 )

            #we get our peptide out again from the first tree
            res = mytree.query_tree( self.q1 - 1, self.ssrcalc -1, 
                                     self.q1 + 1,  self.ssrcalc + 1, 1, 0) 
            self.assertEqual( len(res), 1)
            self.assertEqual( res[0][0], 101)

            #we get our peptide out again from the second tree
            res = mytree2.query_tree( self.q1 + 10 - 1, self.ssrcalc -1, 
                                      self.q1 + 10 + 1,  self.ssrcalc + 1, 1, 0) 
            self.assertEqual( len(res), 1)
            self.assertEqual( res[0][0], 102)

            # but not the other way round
            res = mytree2.query_tree( self.q1 - 1, self.ssrcalc -1, 
                                      self.q1 + 1,  self.ssrcalc + 1, 1, 0) 
            self.assertEqual( len(res), 0)

if __name__ == '__main__':
    unittest.main()
