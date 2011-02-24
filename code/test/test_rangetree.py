import random
import unittest

import sys
sys.path.append( '..')

from test_shared import ignoreImportError_rangetree

try:
    import c_rangetree
except ImportError:
    print "=" * 75, """
Module c_rangetree is not available. Please compile it if you want to use it.
""", "=" * 75

class Test_crangetree(unittest.TestCase):

    def setUp(self):
            import c_rangetree
            self.parent_id = 101
            self.q1 = 501.0
            self.ssrcalc = 24
            self.mytuple1 = (
                ('PEPTIDE', 1, self.parent_id, 2, self.q1, self.ssrcalc),
            )

    def test_rangetree(self):
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

import inspect, types
for name, fn in inspect.getmembers(Test_crangetree):
    if isinstance(fn, types.UnboundMethodType):
        setattr(Test_crangetree, name, ignoreImportError_rangetree(fn))

if __name__ == '__main__':
    unittest.main()
