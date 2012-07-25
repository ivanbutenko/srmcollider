"""
This file tests the interoperability with the MySQL / SQLite databases.
"""

from nose.tools import nottest
from nose.plugins.attrib import attr
import random
import unittest
import sys
sys.path.append( '..')
sys.path.append( '../external')
import collider
from Residues import Residues

from test_shared import *
import test_shared 
from test_shared import check_cgetnonuis_availability

test_database = 'hroest'
mysql_conf_file = "~/.my.cnf.srmcollider"

def psort_old(x,y):
    if(x[1] != y[1]): return -cmp(x[1], y[1]) 
    else: return cmp( x[0], y[0] )

def psort(x,y):
    if(x.transition_group != y.transition_group): return -cmp(x.transition_group, y.transition_group) 
    else: return cmp( x.modified_sequence, y.modified_sequence )

def _get_unique_pepids(par, cursor, ignore_genomeoccurence=False):
    query = """
    select parent_id, q1, q1_charge, ssrcalc, peptide.id, modified_sequence, transition_group
     from %s
     inner join
     ddb.peptide on peptide.id = %s.peptide_key
     inner join ddb.peptideOrganism on peptide.id = peptideOrganism.peptide_key 
     where genome_occurence = 1
     %s
    """ % (par.peptide_tables[0], par.peptide_tables[0], par.query_add )
    if ignore_genomeoccurence:
        query = """
        select parent_id, q1, q1_charge, ssrcalc, peptide_key, modified_sequence, transition_group
         from %s
         where 4 = 4
         %s
        """ % (par.peptide_tables[0], par.query_add )
    if par.print_query: print query
    cursor.execute( query )
    res = cursor.fetchall()
    return [
        {
            'parent_id' :  r[0],
            'q1' :         r[1],
            'q1_charge' :  r[2],
            'ssrcalc' :    r[3],
            'peptide_key' :r[4],
            'mod_sequence':r[5],
            'transition_group':r[6],
        }
        for r in res
    ]

def find_clashes_small(self, mycollider, cursor, par, pepids):
        swath_mode = False
        mycollider.allpeps = {}
        mycollider.non_unique_count = 0
        mycollider.total_count = 0
        MAX_UIS = 5
        for kk, pep in enumerate(pepids):
            precursor = Precursor(parent_id=pep['parent_id'], q1=pep['q1'], 
                ssrcalc=pep['ssrcalc'], modified_sequence = pep['mod_sequence'], 
                transition_group = pep['transition_group'])

            q3_low, q3_high = par.get_q3range_transitions()
            transitions = precursor.calculate_transitions(q3_low, q3_high)
            nr_transitions = len(transitions)

            ### if use_db and not swath_mode:
            precursors = mycollider._get_all_precursors(par, precursor, cursor)
            collisions_per_peptide = collider.get_coll_per_peptide_from_precursors(mycollider,
                    transitions, precursors, par, precursor, forceNonCpp=False)
            collisions_per_peptide_python = collider.get_coll_per_peptide_from_precursors_obj_wrapper(mycollider,
                    transitions, precursors, par, precursor, forceNonCpp=True)
            self.assertEqual(collisions_per_peptide, collisions_per_peptide_python)
            non_uis_list = collider.get_nonuis_list(collisions_per_peptide, MAX_UIS)
            # 
            # here we count how many are locally clean, e.g. look at UIS of order 1
            mycollider.allpeps[precursor.parent_id] = 1.0 - len(non_uis_list[1]) * 1.0  / nr_transitions
            mycollider.non_unique_count += len(non_uis_list[1])
            mycollider.total_count += nr_transitions

def do_check_complete(self, mycollider):

        nonzero = [(k,v) for k,v in mycollider.allpeps.iteritems() if v < 1.0]
        self.assertTrue(abs(mycollider.allpeps[49] - 1.0) < 1e-5)
        self.assertTrue(abs(mycollider.allpeps[81] - 1.0) < 1e-5)
        self.assertEqual( len(nonzero), 18)

        self.assertTrue(abs(mycollider.allpeps[111] - 0.888888888889) < 1e-5  )
        self.assertTrue(abs(mycollider.allpeps[191] - 0.7) < 1e-5             )
        self.assertTrue(abs(mycollider.allpeps[375] - 0.941176470588) < 1e-5  )
        self.assertTrue(abs(mycollider.allpeps[585] - 0.875) < 1e-5           )
        self.assertTrue(abs(mycollider.allpeps[609] - 0.888888888889) < 1e-5  )
        self.assertTrue(abs(mycollider.allpeps[929] - 0.857142857143) < 1e-5  )
        self.assertTrue(abs(mycollider.allpeps[1089] - 0.7) < 1e-5            )
        self.assertTrue(abs(mycollider.allpeps[1101] - 0.916666666667) < 1e-5 )
        self.assertTrue(abs(mycollider.allpeps[1177] - 0.928571428571) < 1e-5 )
        self.assertTrue(abs(mycollider.allpeps[1493] - 0.928571428571) < 1e-5 )
        self.assertTrue(abs(mycollider.allpeps[1537] - 0.857142857143) < 1e-5 )
        self.assertTrue(abs(mycollider.allpeps[1585] - 0.933333333333) < 1e-5 )
        self.assertTrue(abs(mycollider.allpeps[1663] - 0.928571428571) < 1e-5 )
        self.assertTrue(abs(mycollider.allpeps[1763] - 0.8125) < 1e-5         )
        self.assertTrue(abs(mycollider.allpeps[1925] - 0.928571428571) < 1e-5 )
        self.assertTrue(abs(mycollider.allpeps[1957] - 0.833333333333) < 1e-5 )
        self.assertTrue(abs(mycollider.allpeps[1967] - 0.857142857143) < 1e-5 )
        self.assertTrue(abs(mycollider.allpeps[1991] - 0.857142857143) < 1e-5 )

        self.assertTrue( abs(sum( mycollider.allpeps.values() ) - 975.6326447245566) < 10**(-3) )
        self.assertEqual( len([v for v in mycollider.allpeps.values() if v == 1.0 ] ), 960 )
        self.assertEqual( mycollider.non_unique_count, 26 )
        self.assertEqual( mycollider.total_count, 12502 )

class Test_collider_mysql(unittest.TestCase):

    def setUp(self):

        try:
            import MySQLdb
            self.database_available = True
        except ImportError:
            print """Module MySQLdb not available.
            
            Please install it if you want to use it.
            Use the following command (on Ubuntu systems):
                sudo apt-get install python-mysqldb
            """
            self.database_available = False

        try: 
          self.db = MySQLdb.connect(read_default_file=mysql_conf_file)
          self.database_available = True
        except MySQLdb.OperationalError as e:
            print "Could not connect to database: Please check the configuration in test/test_db.py!\n", e
            self.database_available = False

        class Minimal: 
            def get_q3_window_transitions(self, q3):
                if self.ppm: return [q3 - self.q3_window * 10**(-6) *q3, q3 + self.q3_window * 10**(-6) * q3]
                else: return [q3 - self.q3_window, q3 + self.q3_window]

        self.par = Minimal()
        self.par.q3_window = 4.0
        self.par.ppm = False
        self.MAX_UIS = 5
        self.q3_high = 1500
        self.q3_low = 300

        self.par.q1_window = 1.2 / 2.0
        self.par.ssrcalc_window = 9999
        self.par.query2_add = ' and isotope_nr = 0 '
        self.par.peptide_tables = [test_database + '.srmPeptides_human']
        self.par.transition_table = test_database + '.srmTransitions_human'
        self.par.print_query = False
        self.par.select_floor = False
        self.par.isotopes_up_to = 3

        self.par.parent_charges      =  [2,3]

        self.par.bions      =  True
        self.par.yions      =  True
        self.par.aions      =  False
        self.par.aMinusNH3  =  False
        self.par.bMinusH2O  =  False
        self.par.bMinusNH3  =  False
        self.par.bPlusH2O   =  False
        self.par.yMinusH2O  =  False
        self.par.yMinusNH3  =  False
        self.par.cions      =  False
        self.par.xions      =  False
        self.par.zions      =  False
        self.par.MMinusH2O  =  False
        self.par.MMinusNH3  =  False

        def returntrue(): return True
        self.par.do_b_y_only = returntrue

        def returnrange(): return self.q3_low, self.q3_high
        self.par.get_q3range_collisions = returnrange
        #self.par.get_q3_window_transitions = get_q3_window_transitions
        #self.par.get_q3_window_transitions = collider.SRM_parameters.get_q3_window_transitions

        import sys
        self.R = Residues('mono')

        self.acollider = collider.SRMcollider()
        self.aparamset = collider.testcase(testdatabase=test_database)

    def _reducecollisionstoperpep(self, transitions, collisions, par):

        if not self.database_available: return

        collisions_per_peptide = {}
        q3_window_used = par.q3_window
        for t in transitions:
            if par.ppm: q3_window_used = par.q3_window * 10**(-6) * t[0]
            for c in collisions:
                if abs( t[0] - c[0] ) <= q3_window_used:
                    #gets all collisions
                    if collisions_per_peptide.has_key(c[3]):
                        if not t[1] in collisions_per_peptide[c[3]]:
                            collisions_per_peptide[c[3]].append( t[1] )
                    else: 
                        collisions_per_peptide[c[3]] = [ t[1] ] ; 
        return collisions_per_peptide

    @attr('slow') 
    def test_get_all_precursors1(self): 

        if not self.database_available: return

        pep = test_shared.runpep1
        pep = test_shared.runpep_obj1
        transitions = test_shared.runtransitions1
        precursors = test_shared.runprecursors1

        # filter out the isotopes
        pprec = []
        for p in precursors:
          if(p[1] not in [pp[1] for pp in pprec]): pprec.append(p)
        precursors = pprec
        runprecursors_obj2 = []
        for p in precursors:
            runprecursors_obj2.append(Precursor( q1 = p[0], q1_charge = p[3],
            modified_sequence = p[1], transition_group = p[2], isotopically_modified = p[4]))
        precursors = runprecursors_obj2

        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        par = self.par

        cursor = self.db.cursor()
        myprec = collider.SRMcollider()._get_all_precursors(par, pep, cursor)
        self.assertEqual( len(myprec), len(precursors))

        precursors.sort(psort)
        myprec.sort(psort)

        for i in range(len(myprec)):
            self.assertEqual( myprec[i].modified_sequence, precursors[i].modified_sequence)
            self.assertEqual( myprec[i].q1_charge, precursors[i].q1_charge)

    @attr('slow') 
    def test_get_all_precursors2(self): 

        if not self.database_available: return

        pep = test_shared.runpep_obj2
        transitions = test_shared.runtransitions2
        precursors = test_shared.runprecursors2

        # filter out the isotopes
        pprec = []
        for p in precursors:
          if(p[1] not in [pp[1] for pp in pprec]): pprec.append(p)
        precursors = pprec
        runprecursors_obj2 = []
        for p in precursors:
            runprecursors_obj2.append(Precursor( q1 = p[0], q1_charge = p[3],
            modified_sequence = p[1], transition_group = p[2], isotopically_modified = p[4]))
        precursors = runprecursors_obj2

        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        par = self.par

        cursor = self.db.cursor()
        myprec = collider.SRMcollider()._get_all_precursors(par, pep, cursor)
        self.assertEqual( len(myprec), len(precursors))

        precursors.sort(psort)
        myprec.sort(psort)

        for i in range(len(myprec)):
            self.assertEqual( myprec[i].modified_sequence, precursors[i].modified_sequence)
            self.assertEqual( myprec[i].q1_charge, precursors[i].q1_charge)

    @attr('slow') 
    def test_get_coll_per_peptide(self): 

        if not self.database_available: return

        pep = test_shared.runpep_obj1
        transitions = test_shared.runtransitions1
        precursors = test_shared.runprecursors1
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        par = self.par
        q3_high = self.q3_high
        q3_low = self.q3_low
        R = self.R
        cursor = self.db.cursor()

        collisions_per_peptide = collider.get_coll_per_peptide(self.acollider, 
           transitions, par, pep, cursor)
        self.assertEqual(len(collisions_per_peptide), len(test_shared.collpepresult1))
        for k in collisions_per_peptide.keys():
          self.assertEqual(collisions_per_peptide[k], test_shared.collpepresult1[k])
        self.assertEqual(collisions_per_peptide, test_shared.collpepresult1)

        collisions_per_peptide = collider.get_coll_per_peptide(self.acollider, 
           transitions, par, pep, cursor, forceNonCpp=True)
        self.assertEqual(collisions_per_peptide, test_shared.collpepresult1)

    @attr('slow') 
    def test_get_coll_per_peptide2(self): 

        if not self.database_available: return

        pep = test_shared.runpep_obj2
        transitions = test_shared.runtransitions2
        precursors = test_shared.runprecursors2
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])

        par = self.par
        q3_high = self.q3_high
        q3_low = self.q3_low
        R = self.R
        cursor = self.db.cursor()

        collisions_per_peptide = collider.get_coll_per_peptide(self.acollider, 
           transitions, par, pep, cursor)
        self.assertEqual(len(collisions_per_peptide), len(test_shared.collpepresult2))
        for k in collisions_per_peptide.keys():
          self.assertEqual(collisions_per_peptide[k], test_shared.collpepresult2[k])
        self.assertEqual(collisions_per_peptide, test_shared.collpepresult2)

        collisions_per_peptide = collider.get_coll_per_peptide(self.acollider, 
           transitions, par, pep, cursor, forceNonCpp=True)
        self.assertEqual(len(collisions_per_peptide), len(test_shared.collpepresult2))
        self.assertEqual(collisions_per_peptide, test_shared.collpepresult2)

    def test_get_coll_per_peptide_debug(self): 
        """This is a small example with very few interferences, 
        good to debug"""

        if not self.database_available: return

        pep = test_shared.runpep_obj2
        transitions = test_shared.runtransitions2
        precursors = test_shared.runprecursors2
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])

        self.par.q3_window = 0.1 / 2.0
        self.par.q1_window = 0.1 / 2.0
        self.par.ssrcalc_window = 2.0 / 2.0 
        self.par.query2_add = ''
        self.par.peptide_tables = [test_database + '.srmPeptides_human']
        self.par.query2_add = ' and isotope_nr = 0'  # necessary because of old style tables
        self.par.print_query = False

        # Our q1 from peptide 'ELNQLEDK' is 494.751462374, the q1 of the
        # interfering NISVAEVEK is 494.769652374
        #
        # Our 5 transitions from peptide 'ELNQLEDK' are
        # [3, 5, 6, 9, 14]
        # [(504.26693971600002, 3), (357.17740503200002, 5), (485.23598503200003,
        #   6), (842.38957503200004, 9), (421.69870003200003, 14)]
        # The 5 transitions of the interfering peptide are
        #
        # 21406312 |   359414 |  1 | 504.266939716 | y4 | 
        # 21406339 |   359414 |  2 | 357.195595032 | b7 | 
        # 21406321 |   359414 |  1 | 485.272365032 | b5 | 
        # 21406324 |   359414 |  1 | 842.425955032 | b8 | 
        # 21406340 |   359414 |  2 | 421.716890032 | b8 | 
        #

        par = self.par
        q3_high = self.q3_high
        q3_low = self.q3_low
        R = self.R
        cursor = self.db.cursor()
        collisions_per_peptide = collider.get_coll_per_peptide(self.acollider, 
           transitions, par, pep, cursor)

        self.assertEqual(len(collisions_per_peptide[359414]), 5)
        self.assertEqual(collisions_per_peptide[359414], [3, 5, 6, 9, 14])
        
        self.assertEqual(59, len(collisions_per_peptide))
        self.assertEqual(19, len([k for k,v in collisions_per_peptide.iteritems() if len(v) == 1]))
        self.assertEqual(25, len([k for k,v in collisions_per_peptide.iteritems() if len(v) == 2]))
        self.assertEqual(6, len([k for k,v in collisions_per_peptide.iteritems() if len(v) == 3]))
        self.assertEqual(8, len([k for k,v in collisions_per_peptide.iteritems() if len(v) == 4]))
        self.assertEqual(1, len([k for k,v in collisions_per_peptide.iteritems() if len(v) == 5]))
        self.assertEqual(0, len([k for k,v in collisions_per_peptide.iteritems() if len(v) == 6]))

    @check_cgetnonuis_availability
    @attr('slow') 
    def test_complete_without_isotopes(self):

        if not self.database_available: return

        par = collider.testcase(testdatabase=test_database)
        par.query2_add = ' and isotope_nr = 0 ' # still necessary, old style tables
        par.quiet = True
        mycollider = collider.SRMcollider()
        cursor = self.db.cursor()

        #for historic reasons, we only select a subset of peptides
        exclude_pepids = [183, 203, 319, 345, 355, 365, 367, 385, 393, 425, 1227, 1233, 1297, 1299, 1303, 1305, 1307, 1309, 1311, 1509, 1681, 1683]
        mycollider = collider.SRMcollider()
        pepids = _get_unique_pepids(par, cursor, ignore_genomeoccurence=True)
        pepids = [p for p in pepids if p['parent_id'] not in exclude_pepids]

        find_clashes_small(self, mycollider, cursor, par, pepids)
        do_check_complete(self, mycollider)

    @check_cgetnonuis_availability
    @attr('slow') 
    def test_complete_with_isotopes(self):

        if not self.database_available: return

        #now test with isotopes enabled
        par  = collider.testcase(testdatabase=test_database)
        par.quiet = True
        par.isotopes_up_to = 3
        par.eval()
        par.query2_add = ' and isotope_nr = 0 ' # need to do this because we also have other isotopes in there!
        mycollider = collider.SRMcollider()
        cursor = self.db.cursor()

        #for historic reasons, we only select a subset of peptides
        exclude_pepids = [183, 203, 319, 345, 355, 365, 367, 385, 393, 425, 1227, 1233, 1297, 1299, 1303, 1305, 1307, 1309, 1311, 1509, 1681, 1683]
        mycollider = collider.SRMcollider()
        pepids = _get_unique_pepids(par, cursor, ignore_genomeoccurence=True)
        pepids = [p for p in pepids if p['parent_id'] not in exclude_pepids]

        find_clashes_small(self, mycollider, cursor, par, pepids) 

        assert abs(mycollider.allpeps[111] - 0.888888888889) < 1e-5
        assert abs(mycollider.allpeps[191] - 0.7) < 1e-5
        assert abs(mycollider.allpeps[375] - 0.941176470588) < 1e-5
        assert abs(mycollider.allpeps[585] - 0.875) < 1e-5
        assert abs(mycollider.allpeps[609] - 0.6666666) < 1e-5
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

        self.assertTrue( abs(sum( mycollider.allpeps.values() ) - 971.215436) < 10**(-3) )
        self.assertEqual( len([v for v in mycollider.allpeps.values() if v == 1.0 ] ), 922 )

        self.assertEqual( mycollider.non_unique_count,  71)
        self.assertEqual( mycollider.total_count, 12502 )

class Test_collider_sqlite(unittest.TestCase):

    def setUp(self):

        try:
            import sqlite
            from sqlite import DatabaseError
            self.database_available = True
        except ImportError:
            print """Module sqlite not available.
            
            Please install it if you want to use it.
            Use the following command (on Ubuntu systems):
                sudo apt-get install python-sqlite
            """
            self.database_available = False

        try:
            #the database file must exist and the databases must be created
            self.db = sqlite.connect(test_shared.SQLITE_DATABASE_LOCATION)
            cursor = self.db.cursor()
            cursor.execute("select count(*) from srmPeptides_test")
            #print cursor.fetchall()
        except DatabaseError:
            print 
            print "=" * 75
            print """The sqlite database is not available.
            
            Please run the sqlite_setupdb.py script first."""
            print "=" * 75
            self.database_available = False

    def tearDown(self):
      self.db.close()

    @check_cgetnonuis_availability
    def test_1(self):

        if not self.database_available: return
        verbose = False
        par = collider.testcase(testdatabase=test_database)
        #par.query2_add = ' and isotope_nr = 0 ' # still necessary, old style tables
        par.quiet = True
        par.transition_table = 'srmTransitions_test'
        par.peptide_tables = ['srmPeptides_test']
        cursor = self.db.cursor()

        #for historic reasons, we only select a subset of peptides
        exclude_pepids = [183, 203, 319, 345, 355, 365, 367, 385, 393, 425, 1227, 1233, 1297, 1299, 1303, 1305, 1307, 1309, 1311, 1509, 1681, 1683]
        mycollider = collider.SRMcollider()
        pepids = _get_unique_pepids(par, cursor, ignore_genomeoccurence=True)
        pepids = [p for p in pepids if p['parent_id'] not in exclude_pepids]

        find_clashes_small(self, mycollider, cursor, par, pepids)
        do_check_complete(self, mycollider)

if __name__ == '__main__':
    unittest.main()


