import unittest
import test_rangetree, test_cgetnonuis, test_collider, test_db
import test_extra, test_integrated



class Test_dummy(unittest.TestCase): pass


verbosity = 1

sharedtestloader = unittest.TestLoader()
alltests = sharedtestloader.loadTestsFromTestCase(Test_dummy)
independent_tests = sharedtestloader.loadTestsFromTestCase(Test_dummy)
db_tests = sharedtestloader.loadTestsFromTestCase(Test_dummy)

independent_tests.addTests( sharedtestloader.loadTestsFromTestCase(test_cgetnonuis.Test_cgetnonuis) )
independent_tests.addTests( sharedtestloader.loadTestsFromTestCase(test_cgetnonuis.Test_cgetnonuis_get_non_UIS_from_transitions))
independent_tests.addTests( sharedtestloader.loadTestsFromTestCase(test_rangetree.Test_crangetree))
independent_tests.addTests( sharedtestloader.loadTestsFromTestCase(test_collider.Test_collider_function))

db_tests.addTests( sharedtestloader.loadTestsFromTestCase(test_db.Test_collider_mysql))
db_tests.addTests( sharedtestloader.loadTestsFromTestCase(test_db.Test_collider_sqlite))

alltests.addTests( sharedtestloader.loadTestsFromTestCase(test_extra.Test_fragmentation))
alltests.addTests( sharedtestloader.loadTestsFromTestCase(test_integrated.Test_cintegrated))
alltests.addTests( independent_tests )
alltests.addTests( db_tests )


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=verbosity).run(alltests)






