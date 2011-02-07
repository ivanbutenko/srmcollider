import test_all
import unittest

verbosity = 1

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=verbosity).run(test_all.independent_tests)
