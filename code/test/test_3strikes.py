import unittest
"""
This file tests the functionality of the collider.py module. 
"""

import sys
sys.path.extend(['.', '..', '../external/', 'external/'])

import uis_functions
from collider import thisthirdstrike
from test_shared import check_cgetnonuis_availability

try: import c_getnonuis
except ImportError: pass

class Test_3strikes(unittest.TestCase): 

    def setUp(self):
        self.ssrcalcvalues_two_example  = [  
            [ 1,   1.1, 1.8, 1.9], 
            [2.5, 2.6, 2.7, 2.9], 
        ]
        self.ssrcalcvalues_four_example  = [  
            [ 1,   5,   8,   8.5], 
            [ 1.5,       8.2], 
            [ 2.1, 5.1,  8.1], 
            [1.0] 
        ]

    def test_example_tworows(self):
        # So with 0.5 we find that 2.5 - 1.9 = 0.6 is too far away
        strike3_ssrcalcwindow = 0.5
        ssrcalcvalues = self.ssrcalcvalues_two_example
        N = [len(v) for v in ssrcalcvalues]
        res = thisthirdstrike(N, ssrcalcvalues, strike3_ssrcalcwindow)
        newdic = dict([(i,list(v)) for i,v in enumerate(res)])
        expanded = uis_functions.get_nonuis_list(newdic, 2)

        self.assertEqual(len(expanded[2]), 0 )

    def test_example_tworows_larger_window(self):
        # So with 1.0 we find that 2.5 - 1.9 = 0.6 is close enough
        strike3_ssrcalcwindow = 1.0
        ssrcalcvalues = self.ssrcalcvalues_two_example
        N = [len(v) for v in ssrcalcvalues]
        res = thisthirdstrike(N, ssrcalcvalues, strike3_ssrcalcwindow)
        newdic = dict([(i,list(v)) for i,v in enumerate(res)])
        expanded = uis_functions.get_nonuis_list(newdic, 2)

        self.assertEqual(len(expanded[2]), 1 )

    def test_example_four_transitions(self):
        strike3_ssrcalcwindow = 1.0
        ssrcalcvalues  = self.ssrcalcvalues_four_example  
        
        N = [len(v) for v in ssrcalcvalues]
        res = thisthirdstrike(N, ssrcalcvalues, strike3_ssrcalcwindow)
        newdic = dict([(i,list(v)) for i,v in enumerate(res)])
        expanded = uis_functions.get_nonuis_list(newdic, 5)

        # There are five combinations of length 2 that are not eUIS: 
        #   (0,1), (0,2), (0,3), (1,2) and (1,3)
        # There are two combinations of length 3 that are not eUIS: 
        #   (0,1,3), (0,1,2)
        self.assertEqual(len(expanded[2]), 5 )
        self.assertEqual(len(expanded[3]), 2 )
        self.assertEqual(len(expanded[4]), 0 )

        self.assertTrue((0,1) in expanded[2])
        self.assertTrue((0,2) in expanded[2])
        self.assertTrue((0,3) in expanded[2])
        self.assertTrue((1,2) in expanded[2])
        self.assertTrue((1,3) in expanded[2])

        self.assertTrue((0,1,3) in expanded[3])
        self.assertTrue((0,1,2) in expanded[3])

    def test_example_four_transitions_small_window(self):
        strike3_ssrcalcwindow = 0.3
        ssrcalcvalues  = self.ssrcalcvalues_four_example  
        
        N = [len(v) for v in ssrcalcvalues]
        res = thisthirdstrike(N, ssrcalcvalues, strike3_ssrcalcwindow)
        newdic = dict([(i,list(v)) for i,v in enumerate(res)])
        expanded = uis_functions.get_nonuis_list(newdic, 5)

        # Now we loose (1,3) because 1.5 and 1.0 is too big a distance
        self.assertEqual(len(expanded[2]), 4 )
        self.assertEqual(len(expanded[3]), 1 )
        self.assertEqual(len(expanded[4]), 0 )

        self.assertTrue((1,3) not in expanded[2])
        self.assertTrue((0,1,2) in expanded[3])

    @check_cgetnonuis_availability
    def test_cpp_implementation(self):
        strike3_ssrcalcwindow = 0.3
        ssrcalcvalues  = self.ssrcalcvalues_four_example  
        N = [len(v) for v in ssrcalcvalues]

        cpp_result = c_getnonuis.calculate_eUIS(N, ssrcalcvalues, strike3_ssrcalcwindow)
        python_result = thisthirdstrike(N, ssrcalcvalues, strike3_ssrcalcwindow)
        cpp_result_compare = {}
        for c in cpp_result:
            cpp_result_compare[ tuple(c)] = 0
        self.assertEqual(cpp_result_compare, python_result)

        #print "second cpp impl"
        strike3_ssrcalcwindow = 1.0
        cpp_result = c_getnonuis.calculate_eUIS(N, ssrcalcvalues, strike3_ssrcalcwindow)
        python_result = thisthirdstrike(N, ssrcalcvalues, strike3_ssrcalcwindow)
        cpp_result_compare = {}
        for c in cpp_result:
            cpp_result_compare[ tuple(c)] = 0
        self.assertEqual(cpp_result_compare, python_result)

if __name__ == '__main__':
    unittest.main()

