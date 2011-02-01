import unittest
import test_rangetree, test_cgetnonuis, test_collider, test_db

verbosity = 1

cgetnonuis_suite = unittest.TestLoader().loadTestsFromTestCase(test_cgetnonuis.Test_cgetnonuis)
cgetnonuis_tr_suite = unittest.TestLoader().loadTestsFromTestCase(test_cgetnonuis.Test_cgetnonuis_get_non_UIS_from_transitions)
crangetree_suite = unittest.TestLoader().loadTestsFromTestCase(test_rangetree.Test_crangetree)
collider_sqlite_suite = unittest.TestLoader().loadTestsFromTestCase(test_db.Test_collider_sqlite)
collider_function_suite = unittest.TestLoader().loadTestsFromTestCase(test_collider.Test_collider_function)

unittest.TextTestRunner(verbosity=verbosity).run(cgetnonuis_suite)
unittest.TextTestRunner(verbosity=verbosity).run(cgetnonuis_tr_suite)
unittest.TextTestRunner(verbosity=verbosity).run(crangetree_suite)
unittest.TextTestRunner(verbosity=verbosity).run(collider_function_suite)
unittest.TextTestRunner(verbosity=verbosity).run(collider_sqlite_suite)




