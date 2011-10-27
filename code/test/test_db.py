import random
import unittest
"""
This file tests the interoperability with the MySQL / SQLite databases.
"""

import sys
sys.path.append( '..')
sys.path.append( '../external')
import collider
from Residues import Residues

from test_shared import *
import test_shared 

test_database = 'hroest'

def psort(x,y):
    if(x[1] != y[1]): return -cmp(x[1], y[1]) 
    else: return cmp( x[0], y[0] )


def _get_unique_pepids(par, cursor, ignore_genomeoccurence=False):
    query = """
    select parent_id, q1, q1_charge, ssrcalc, peptide.id, modified_sequence
     from %s
     inner join
     ddb.peptide on peptide.id = %s.peptide_key
     inner join ddb.peptideOrganism on peptide.id = peptideOrganism.peptide_key 
     where genome_occurence = 1
     %s
    """ % (par.peptide_table, par.peptide_table, par.query_add )
    if ignore_genomeoccurence:
        query = """
        select parent_id, q1, q1_charge, ssrcalc, peptide_key, modified_sequence
         from %s
         where 4 = 4
         %s
        """ % (par.peptide_table, par.query_add )
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
        }
        for r in res
    ]

def _get_all_transitions(self, par, pep, cursor, values = "q3, srm_id"):
    q3_low, q3_high = par.get_q3range_transitions()
    query1 = """
    select %(values)s
    from %(pep)s
    inner join %(trans)s
      on %(pep)s.transition_group = %(trans)s.group_id
    where parent_id = %(parent_id)s
    and q3 > %(q3_low)s and q3 < %(q3_high)s         
    %(query_add)s
    """ % { 'parent_id' : pep['parent_id'], 'q3_low' : q3_low,
           'q3_high' : q3_high, 'query_add' : par.query1_add,
           'pep' : par.peptide_table, 'trans' : par.transition_table, 
           'values' : values}
    if par.print_query: print query1
    cursor.execute( query1 )
    return cursor.fetchall()

def _get_all_collisions(self, par, pep, cursor, 
            values="q3, q1, srm_id, peptide_key", transitions=None,
            bysequence=False):
    q3_low, q3_high = par.get_q3range_collisions()
    vdict = { 'q1' : pep['q1'], 'ssrcalc' : pep['ssrcalc'],
            'peptide_key' : pep['peptide_key'], 'q1_window' : par.q1_window,
            'query_add' : par.query2_add, 'ssr_window' : par.ssrcalc_window,
            'pep' : par.peptide_table, 'values' : values, 
            'pepseq' : pep['mod_sequence'], 
            'q3_low':q3_low,'q3_high':q3_high, 
            'trans' : par.transition_table,
            }
    #we compare the parent ion against 4 different parent ions
    #thus we need to take the PEPTIDE key here
    if bysequence: selectby = "and %(pep)s.modified_sequence != '%(pepseq)s'" % vdict
    else: selectby = "and %(pep)s.peptide_key != %(peptide_key)d" % vdict
    vdict['selectby'] = selectby
    query2 = """
    select %(values)s
    from %(pep)s
    inner join %(trans)s
      on %(pep)s.transition_group = %(trans)s.group_id
    where ssrcalc > %(ssrcalc)s - %(ssr_window)s 
        and ssrcalc < %(ssrcalc)s + %(ssr_window)s
    and q1 > %(q1)s - %(q1_window)s and q1 < %(q1)s + %(q1_window)s
    %(selectby)s
    and q3 > %(q3_low)s and q3 < %(q3_high)s
    %(query_add)s
    """ % vdict
    if transitions is None: 
        #if par.select_floor: query2 += "\n GROUP BY FLOOR(q3)"
        if par.print_query: print query2
        cursor.execute( query2 )
    #Here we add N conditions to the SQL query where N = 2*len(transitions)
    #we only select those transitions that also could possible interact
    #We gain about a factor two to three [the smaller the window, the better]
    else:
        txt = 'and ('
        q3_window_used = par.q3_window 
        for q3, id in transitions:
            if par.ppm: q3_window_used = par.q3_window * 10**(-6) * q3
            txt += '(q3 > %s and q3 < %s) or\n' % (q3 - q3_window_used, q3 + q3_window_used)
        txt = txt[:-3] + ')'
        #if par.select_floor: txt += "\n GROUP BY FLOOR(q3)"
        if par.print_query: print query2 + txt
        cursor.execute( query2 + txt )
    return cursor.fetchall()

def find_clashes_small(self, db, par, use_per_transition=False,
                       UIS_only=False, exp_key=None, calc_collisions=False, 
                       pepids=None):
    #make sure we only get unique peptides
    cursor = db.cursor()
    if pepids is not None: self.pepids = pepids
    else: self.pepids = _get_unique_pepids(par, cursor)
    MAX_UIS = par.max_uis
    common_filename = par.get_common_filename()
    self.allpeps = {}
    self.non_unique_count = 0
    self.total_count = 0
    for ii, pep in enumerate(self.pepids):
        p_id = pep['parent_id']
        q1 = pep['q1']
        non_unique = {}
        non_unique_q1 = {}
        non_unique_ppm = {}
        non_unique_clash = {}
        transitions = _get_all_transitions(self,par, pep, cursor)
        nr_transitions = len( transitions )
        if nr_transitions == 0: continue #no transitions in this window
        if use_per_transition: collisions = self._get_all_collisions_per_transition(par, pep, transitions, cursor)
        elif calc_collisions: 
            collisions = list(self._get_all_collisions_calculate(par, pep, cursor))
        else: collisions = _get_all_collisions(self,par, pep, cursor, transitions = transitions)
        #
        if UIS_only:
            non_uis_list = get_non_UIS_from_transitions(transitions, 
                                        collisions, par, MAX_UIS)
            for i in range(1,min(MAX_UIS+1, nr_transitions+1)): 
                restable = 'hroest.result_srmuis'
                cursor.execute( "insert into %s (non_useable_UIS, total_UIS, \
                               parent_key, uisorder, exp_key) values " % restable +\
                               "(%s,%s,%s,%s,%s)" % (len(non_uis_list[i]) , 
                                   choose(nr_transitions, i), p_id , i, exp_key) )
            continue
        #here we loop through all possible combinations of transitions and
        #potential collisions and check whether we really have a collision
        q3_window_used = par.q3_window
        for t in transitions:
            if par.ppm: q3_window_used = par.q3_window * 10**(-6) * t[0]
            this_min = q3_window_used
            for c in collisions:
                if abs( t[0] - c[0] ) <= this_min:
                    non_unique[ t[1] ] = t[0] - c[0]
        self.allpeps[ p_id ] = 1.0 - len( non_unique ) * 1.0  / nr_transitions
        self.non_unique_count += len( non_unique )
        self.total_count += len( transitions)



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

        self.db = MySQLdb.connect(read_default_file="~/.my.cnf.srmcollider")

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
        self.par.peptide_table = test_database + '.srmPeptides_human'
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

    def test_get_all_precursors1(self): 
        pep = test_shared.runpep1
        transitions = test_shared.runtransitions1
        precursors = test_shared.runprecursors1
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        par = self.par

        cursor = self.db.cursor()
        myprec = collider.SRMcollider()._get_all_precursors(par, pep, cursor)
        myprec = list(myprec)
        precursors.sort(psort)
        myprec.sort(psort)

        pprec = []
        for p in precursors:
          if(p[1] not in [pp[1] for pp in pprec]): pprec.append(p)
        precursors = pprec

        for i in range(len(myprec)):
            self.assertEqual( myprec[i][1], precursors[i][1])
            self.assertEqual( myprec[i][2], precursors[i][2])
            self.assertEqual( myprec[i][3], precursors[i][3])
        self.assertEqual( len(myprec), len(precursors))
        #self.assertEqual( myprec, precursors)

    def test_get_all_precursors2(self): 
        pep = test_shared.runpep2
        transitions = test_shared.runtransitions2
        precursors = test_shared.runprecursors2
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        par = self.par

        cursor = self.db.cursor()
        myprec = collider.SRMcollider()._get_all_precursors(par, pep, cursor)
        myprec = list(myprec)
        precursors.sort(psort)
        myprec.sort(psort)

        pprec = []
        for p in precursors:
          if(p[1] not in [pp[1] for pp in pprec]): pprec.append(p)

        precursors = pprec

        for i in range(len(myprec)):
            self.assertEqual( myprec[i][1], precursors[i][1])
            self.assertEqual( myprec[i][2], precursors[i][2])
            self.assertEqual( myprec[i][3], precursors[i][3])
        self.assertEqual( len(myprec), len(precursors))
        #self.assertEqual( myprec, precursors)

    def disabled_test_get_all_collisions1(self):
        pep = test_shared.runpep1
        transitions = test_shared.runtransitions1
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        par = self.par
        cursor = self.db.cursor()

        collisions = list(_get_all_collisions(
                collider.SRMcollider(), par, pep, cursor))
        collisions_per_peptide = self._reducecollisionstoperpep(transitions, collisions, par)
        self.assertEqual(len(collisions_per_peptide), len(test_shared.collpepresult1))
        self.assertEqual(collisions_per_peptide, test_shared.collpepresult1)
        #now test the with transitions fxn
        collisions = list(_get_all_collisions(
                collider.SRMcollider(), par, pep, cursor, transitions=transitions))
        collisions_per_peptide = self._reducecollisionstoperpep(transitions, collisions, par)
        self.assertEqual(collisions_per_peptide, test_shared.collpepresult1)

    def disabled_test_get_all_collisions2(self):
        pep = test_shared.runpep2
        transitions = test_shared.runtransitions2
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        par = self.par
        cursor = self.db.cursor()

        collisions = list(_get_all_collisions(
                collider.SRMcollider(), par, pep, cursor))
        collisions_per_peptide = self._reducecollisionstoperpep(transitions, collisions, par)
        self.assertEqual(len(collisions_per_peptide), len(test_shared.collpepresult2) )
        self.assertEqual(collisions_per_peptide, test_shared.collpepresult2)
        #now test the with transitions fxn
        collisions = list(_get_all_collisions(
                collider.SRMcollider(), par, pep, cursor, transitions=transitions))
        collisions_per_peptide = self._reducecollisionstoperpep(transitions, collisions, par)
        self.assertEqual(collisions_per_peptide, test_shared.collpepresult2)

    def test_get_coll_per_peptide(self): 
        pep = test_shared.runpep1
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
        self.assertEqual(collisions_per_peptide, test_shared.collpepresult1)
        #collisions_per_peptide = collider.get_coll_per_peptide(self.acollider, 
        #   transitions, par, pep, cursor, do_not_calculate=True)
        self.assertEqual(collisions_per_peptide, test_shared.collpepresult1)
        collisions_per_peptide = collider.get_coll_per_peptide(self.acollider, 
           transitions, par, pep, cursor, forceNonCpp=True)
        self.assertEqual(collisions_per_peptide, test_shared.collpepresult1)

    def test_get_coll_per_peptide2(self): 
        pep = test_shared.runpep2
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
        self.assertEqual(collisions_per_peptide, test_shared.collpepresult2)
        #collisions_per_peptide = collider.get_coll_per_peptide(self.acollider, 
        #transitions, par, pep, cursor, do_not_calculate=True)
        self.assertEqual(collisions_per_peptide, test_shared.collpepresult2)
        collisions_per_peptide = collider.get_coll_per_peptide(self.acollider, 
           transitions, par, pep, cursor, forceNonCpp=True)
        self.assertEqual(collisions_per_peptide, test_shared.collpepresult2)

    def disabled_test_get_coll_per_peptide_debug(self): 
        """This is a small example with very few interferences, 
        good to debug"""
        pep = test_shared.runpep2
        transitions = test_shared.runtransitions2
        precursors = test_shared.runprecursors2
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])

        self.par.q3_window = 0.1 / 2.0
        self.par.q1_window = 0.1 / 2.0
        self.par.ssrcalc_window = 2.0 / 2.0 
        self.par.query2_add = ''
        self.par.peptide_table = test_database + '.srmPeptides_yeast'
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
        self.assertEqual(59, len(collisions_per_peptide))
        self.assertEqual(19, len([k for k,v in collisions_per_peptide.iteritems() if len(v) == 1]))
        self.assertEqual(25, len([k for k,v in collisions_per_peptide.iteritems() if len(v) == 2]))
        self.assertEqual(6, len([k for k,v in collisions_per_peptide.iteritems() if len(v) == 3]))
        self.assertEqual(8, len([k for k,v in collisions_per_peptide.iteritems() if len(v) == 4]))
        self.assertEqual(1, len([k for k,v in collisions_per_peptide.iteritems() if len(v) == 5]))
        self.assertEqual(0, len([k for k,v in collisions_per_peptide.iteritems() if len(v) == 6]))

    #TODO also test uis option of clashes_small
    def test_1(self):

        par  = collider.testcase(testdatabase=test_database)
        par.quiet = True
        mycollider = collider.SRMcollider()
        find_clashes_small(mycollider, self.db, par) 
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
        par  = collider.testcase(testdatabase=test_database)
        par.quiet = True
        mycollider = collider.SRMcollider()
        find_clashes_small(mycollider, self.db, par) 
        self.assertTrue (abs(sum( mycollider.allpeps.values() ) - 975.6326447245566) < 10**(-3) )
        self.assertEqual( mycollider.non_unique_count,  26 )
        self.assertEqual( mycollider.total_count, 12502 )
        self.assertTrue (abs(mycollider.allpeps[1585] - 0.93333) < 10**(-3) )
        #self.assertEqual( len( mycollider.q3min_distr ),  26)
        #self.assertEqual( len( mycollider.q1min_distr ),  26)
        #self.assertEqual( len( mycollider.found3good ),  978)

    def test_3(self):
        #now test with isotopes enabled
        par  = collider.testcase(testdatabase=test_database)
        par.quiet = True
        par.isotopes_up_to = 3
        par.eval()
        mycollider = collider.SRMcollider()
        find_clashes_small(mycollider, self.db, par) 
        self.assertTrue( abs(sum( mycollider.allpeps.values() ) - 975.6326447245566) < 10**(-3) )
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
            self.db = sqlite.connect('/tmp/testdb')
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


    def test_1(self):

        if not self.database_available: return
        par = collider.testcase(testdatabase=test_database)
        par.quiet = True
        par.transition_table = 'srmTransitions_test'
        par.peptide_table = 'srmPeptides_test'
        cursor = self.db.cursor()
        #for historic reasons, we only select a subset of peptides
        exclude_pepids = [183, 203, 319, 345, 355, 365, 367, 385, 393, 425, 1227, 1233, 1297, 1299, 1303, 1305, 1307, 1309, 1311, 1509, 1681, 1683]
        mycollider = collider.SRMcollider()
        pepids = _get_unique_pepids(par, cursor, ignore_genomeoccurence=True)
        pepids = [p for p in pepids if p['parent_id'] not in exclude_pepids]

        find_clashes_small(mycollider, self.db, par, pepids=pepids) 

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


