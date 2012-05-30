"""
This file tests the speed of different implementations, mostly C++ vs Python.
"""

from nose.tools import nottest
from nose.plugins.attrib import attr
import random
import unittest
import time

import sys
sys.path.append( '.')
sys.path.append( '..')
sys.path.append( '../test/')
sys.path.append( './test/')
sys.path.append( '../external/')
import collider
import DDB 
import test_shared
from Residues import Residues
from precursor import Precursors

try:
    import c_getnonuis
except ImportError:
    print "=" * 75, """
Module c_getnonuis is not available. Please compile it if you want to use it.
""", "=" * 75

try:
    import c_rangetree
except ImportError:
    print "=" * 75, """
Module c_rangetree is not available. Please compile it if you want to use it.
""", "=" * 75

try:
    import c_integrated
except ImportError:
    print "=" * 75, """
Module c_integrated is not available. Please compile it if you want to use it.
""", "=" * 75

def mysort(x,y):
    if cmp(x[3], y[3]) == 0:
        return cmp(x[0], y[0])
    else: return cmp(x[3], y[3])

from test_shared import _get_unique_pepids
from test_shared import _get_all_collisions
from test_shared import get_non_UIS_from_transitions_old
from test_shared import get_non_UIS_from_transitions

from precursor import Precursor
def runcpp(self, pep, transitions, precursors):
    #first we run the C++ version
    st = time.time()
    pre_obj = [Precursor(modified_sequence=p[1], q1_charge=2, q1=p[0], transition_group=p[2], isotopically_modified=0) for p in precursors]
    for kk in range(10):
        tmp = c_getnonuis.calculate_collisions_per_peptide_other_ion_series( 
            transitions, pre_obj, self.par, self.q3_low, self.q3_high, self.par.q3_window, self.par.ppm, False)
    ctime = (time.time() - st)/10.0
    return tmp, ctime

def runpy(self, pep, transitions, precursors):
    #now we run the pure Python version
    st = time.time()
    collisions = list(collider.SRMcollider._get_all_collisions_calculate_sub(
            collider.SRMcollider(), precursors, self.par, self.R, self.q3_low, self.q3_high))
    collisions_per_peptide = {}
    q3_window_used = self.par.q3_window
    for t in transitions:
        if self.par.ppm: q3_window_used = self.par.q3_window * 10**(-6) * t[0]
        for c in collisions:
            if abs( t[0] - c[0] ) <= q3_window_used:
                #gets all collisions
                if collisions_per_peptide.has_key(c[3]):
                    if not t[1] in collisions_per_peptide[c[3]]:
                        collisions_per_peptide[c[3]].append( t[1] )
                else: 
                    collisions_per_peptide[c[3]] = [ t[1] ] ; 
    oldtime = time.time() - st
    return collisions_per_peptide, oldtime

class Test_speed(unittest.TestCase):

    def setUp(self):
        try:
            import c_getnonuis
            self.transitions = test_shared.transitions_def1
            self.collisions = test_shared.collisions_def1
            self.pep1 = test_shared.peptide1
            self.pep2 = test_shared.peptide2

            self.transitions_12_between300_1500 = test_shared.transitions_12_between300_1500
            self.pep1_yseries = test_shared.pep1_yseries
            self.pep1_bseries = test_shared.pep1_bseries

            tuples = []
            tuples.append(self.pep1)
            tuples.append(self.pep2)
            self.tuples = tuple( tuples)

            self.R = Residues('mono')

            class Minimal: pass
            self.par = Minimal()
            self.par.q3_window = 4.0
            self.par.ppm = False
            self.q3_high = 1500
            self.q3_low = 300
            self.MAX_UIS = 5

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

        except ImportError:
            pass

    def test_calculatetrans1(self):
        print '\nTesting calculate_transitions 1'
        pep = test_shared.runpep1
        transitions = test_shared.runtransitions1
        precursors = test_shared.runprecursors1[240:300]
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])

        st = time.time()
        collisions = list(collider.SRMcollider._get_all_collisions_calculate_sub(
                collider.SRMcollider(), precursors, self.par, self.R, self.q3_low, self.q3_high))

        oldtime = time.time() - st
        st = time.time()
        tr_new = c_getnonuis.calculate_transitions_ch(tuple(precursors), [1,2], self.q3_low, self.q3_high)
        ctime = time.time() - st

        tr_new.sort(mysort)
        collisions.sort(mysort)

        self.assertEqual(len(tr_new), len(collisions))
        for o,n in zip(collisions, tr_new):
            for oo,nn in zip(o,n):
                self.assertTrue(oo - nn < 1e-5)

        print ctime, oldtime
        print "Speedup:", oldtime / ctime

    def test_calculatetrans2(self):
        print '\nTesting calculate_transitions 2'
        pep = test_shared.runpep2
        transitions = test_shared.runtransitions2
        precursors = test_shared.runprecursors2
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])

        st = time.time()
        collisions = list(collider.SRMcollider._get_all_collisions_calculate_sub(
                collider.SRMcollider(), precursors, self.par, self.R, self.q3_low, self.q3_high))

        oldtime = time.time() - st
        st = time.time()
        tr_new = c_getnonuis.calculate_transitions_ch(tuple(precursors), [1,2], self.q3_low, self.q3_high)
        ctime = time.time() - st

        tr_new.sort(mysort)
        collisions.sort(mysort)

        self.assertEqual(len(tr_new), len(collisions))
        for o,n in zip(collisions, tr_new):
            for oo,nn in zip(o,n):
                self.assertTrue(oo - nn < 1e-5)

        print ctime, oldtime
        print "Speedup:", oldtime / ctime

    def test_collperpeptide1(self):
        print '\nTesting collisions_per_peptide'

        pep = test_shared.runpep1
        transitions = test_shared.runtransitions1
        precursors = test_shared.runprecursors1
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])

        tmp, ctime = runcpp(self, pep, transitions, precursors)
        collisions_per_peptide, oldtime = runpy(self, pep, transitions, precursors)
        self.assertEqual(collisions_per_peptide, tmp )

        print ctime, oldtime
        print "Speedup:", oldtime / ctime

    def test_collperpeptide2(self):
        print '\nTesting collisions_per_peptide'

        pep = test_shared.runpep2
        transitions = test_shared.runtransitions2
        precursors = test_shared.runprecursors2
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])

        tmp, ctime = runcpp(self, pep, transitions, precursors)
        collisions_per_peptide, oldtime = runpy(self, pep, transitions, precursors)
        self.assertEqual(collisions_per_peptide, tmp )

        print ctime, oldtime
        print "Speedup:", oldtime / ctime
        
    def test_get_non_uis1(self):
        print '\nTesting non_uis_list'

        pep = test_shared.runpep1
        transitions = test_shared.runtransitions1
        precursors = test_shared.runprecursors1
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        collisions_per_peptide, ctime = runcpp(self, pep, transitions, precursors)

        MAX_UIS = self.MAX_UIS

        non_uis_list = [set() for i in range(MAX_UIS+1)]
        #here we calculate the UIS for this peptide with the given RT-range
        st = time.time()
        for pepc in collisions_per_peptide.values():
            for i in range(1,MAX_UIS+1):
                collider.get_non_uis(pepc, non_uis_list[i], i)
        oldtime = time.time() - st

        st = time.time()
        for kk in range(10):
            non_uis_list_new = [{} for i in range(MAX_UIS+1)]
            for order in range(1,MAX_UIS+1):
                non_uis_list_new[order] = c_getnonuis.get_non_uis(
                    collisions_per_peptide, order)

        non_uis_list_new = [set( res.keys() ) for res in non_uis_list_new]
        ctime = (time.time() - st)/10

        self.assertEqual(non_uis_list_new, non_uis_list)
        print ctime, oldtime
        print "Speedup:", oldtime / ctime

    def test_get_non_uis2(self):
        print '\nTesting non_uis_list'

        pep = test_shared.runpep2
        transitions = test_shared.runtransitions2
        precursors = test_shared.runprecursors2
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        collisions_per_peptide, ctime = runcpp(self, pep, transitions, precursors)

        MAX_UIS = self.MAX_UIS

        non_uis_list = [set() for i in range(MAX_UIS+1)]
        #here we calculate the UIS for this peptide with the given RT-range
        st = time.time()
        for pepc in collisions_per_peptide.values():
            for i in range(1,MAX_UIS+1):
                collider.get_non_uis(pepc, non_uis_list[i], i)
        oldtime = time.time() - st

        st = time.time()
        for kk in range(10):
            non_uis_list_new = [{} for i in range(MAX_UIS+1)]
            for order in range(1,MAX_UIS+1):
                non_uis_list_new[order] = c_getnonuis.get_non_uis(
                    collisions_per_peptide, order)

        non_uis_list_new = [set( res.keys() ) for res in non_uis_list_new]
        ctime = (time.time() - st)/10

        self.assertEqual(non_uis_list_new, non_uis_list)
        print ctime, oldtime
        print "Speedup:", oldtime / ctime

    def test_calculatetrans(self):
        print '\nTesting calculate_transitions'

        myprecursors = ((500, 'PEPTIDE', 1, 1, 0), (400, 'CEPC[160]IDM[147]E', 2, 2, 0))
        st = time.time()
        for i in range(100000):
            tr_new = c_getnonuis.calculate_transitions_ch(myprecursors, [1,2], self.q3_low, self.q3_high)
        ctime = time.time() - st
        st = time.time()
        coll = collider.SRMcollider()
        for i in range(100000):
            tr_old = list(collider.SRMcollider._get_all_collisions_calculate_sub(
                        coll, myprecursors, self.par, self.R, self.q3_low, self.q3_high))
        oldtime = time.time() - st

        tr_new.sort(mysort)
        tr_old.sort(mysort)

        self.assertEqual(len(tr_new), len(tr_old) )
        self.assertEqual(tr_new, tr_old)
        print ctime, oldtime
        print "Speedup:", oldtime / ctime

    def test_getmintrans1(self):
        print '\nTesting _getMinNeededTransitions'
        pep = test_shared.runpep1
        transitions = test_shared.runtransitions1
        precursors = test_shared.runprecursors1
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        par = self.par
        q3_high = self.q3_high
        q3_low = self.q3_low
        R = self.R
        par.max_uis = 15

        collisions = list(collider.SRMcollider._get_all_collisions_calculate_sub(
                collider.SRMcollider(), precursors, par, R, q3_low, q3_high))
        st = time.time()
        m = collider.SRMcollider()._getMinNeededTransitions(par, transitions, collisions)
        oldtime = time.time() - st

        import c_integrated
        st = time.time()
        m = c_integrated.getMinNeededTransitions(transitions, tuple(precursors), 
            par.max_uis, par.q3_window, par.ppm, par)
        ctime = time.time() - st

        print ctime, oldtime
        print "Speedup:", oldtime / ctime

import time,sys
from random import shuffle
class Test_speed_integrated(unittest.TestCase):

    def setUp(self):
        self.limit = 100
        self.limit_large = 100
        self.limit = 300
        self.limit_large = 600
        try:
            import MySQLdb
            #db_l = MySQLdb.connect(read_default_file="~/.my.cnf.local")
            db = MySQLdb.connect(read_default_file="~/.my.cnf.orl")
            cursor = db.cursor()
            self.cursor = cursor


            par = test_shared.get_default_setup_parameters()
            par.use_sqlite = True
            par.q1_window = 1.2 / 2.0 
            par.q3_window = 2.0 / 2.0
            par.ssrcalc_window = 4 / 2.0 
            par.ssrcalc_window = 9999 / 2.0 
            par.peptide_table = 'hroest.srmPeptides_yeast'
            par.transition_table = 'hroest.srmTransitions_yeast'
            par.isotopes_up_to = 3
            self.mycollider = collider.SRMcollider()

            self.par = par
            self.min_q1 = 440
            self.max_q1 = 450

            ## For debugging
            ##self.max_q1 = 440.18
            ##par.q1_window = 0.009 / 2.0 
            ##par.q3_window = 2.0 / 2.0
            ##par.isotopes_up_to = 1

            import Residues
            R = Residues.Residues('mono')
            isotope_correction = par.isotopes_up_to * R.mass_diffC13 / min(par.parent_charges)

            start = time.time()
            query =  """
            select modified_sequence, transition_group, parent_id, q1_charge, q1, ssrcalc, isotope_nr, 0, 0
            from %s where q1 between %s and %s
            #and isotope_nr = 0
                           """ % (par.peptide_table, self.min_q1 - par.q1_window, self.max_q1 + par.q1_window) 
            cursor.execute(query)
            self.alltuples =  list(cursor.fetchall() )
            #print "len alltuples", len(self.alltuples)

            query =  """
            select modified_sequence, transition_group, parent_id, q1_charge, q1, ssrcalc, isotope_nr, 0, 0
            from %s where q1 between %s - %s and %s
            and isotope_nr = 0
                           """ % (par.peptide_table, self.min_q1 - par.q1_window, isotope_correction, self.max_q1 + par.q1_window) 
            cursor.execute(query)
            self.alltuples_isotope_correction =  list(cursor.fetchall())
            #print "len alltuples zero", len(self.alltuples_isotope_correction)

            self.myprecursors = Precursors()

            # myprecursors.getFromDB(par, db.cursor(), self.min_q1 - par.q1_window, self.max_q1 + par.q1_window)
            ##### LEGACY getFromDB -- have isotope_nr = 0 in there!

            # Get all precursors from the DB within a window of Q1
            lower_q1 = self.min_q1 - par.q1_window
            upper_q1 =  self.max_q1 + par.q1_window
            self.myprecursors.precursors = []
            isotope_correction = par.isotopes_up_to * R.mass_diffC13 / min(par.parent_charges)
            q =  """
            select modified_sequence, transition_group, parent_id, q1_charge, q1, ssrcalc, modifications, missed_cleavages, isotopically_modified
            from %(peptide_table)s where q1 between %(lowq1)s - %(isotope_correction)s and %(highq1)s
            and isotope_nr = 0
            """ % {'peptide_table' : par.peptide_table, 
                          'lowq1'  : lower_q1,  # min_q1 - par.q1_window
                          'highq1' : upper_q1, # max_q1 + par.q1_window,
                          'isotope_correction' : isotope_correction
                  } 
            cursor.execute(q)
            for res in cursor.fetchall():
              p = Precursor()
              p.initialize(*res)
              # Only include those precursors that are actually have isotopes in the specified range
              if(p.included_in_isotopic_range(lower_q1, upper_q1, par, R) ): 
                self.myprecursors.precursors.append(p)

            ##### END LEGACY getFromDB

            self.precursors_to_evaluate = self.myprecursors.getPrecursorsToEvaluate(self.min_q1, self.max_q1)

        except Exception as inst:
            print "something went wrong"
            print inst

    @attr('comparison') 
    def test_integrated_cpp(self):
            """ Compare the fully integrated vs the mixed C++ rangetree / Python solution.

            Here we are comparing the fully integrated solution of storing all
            precursors in a C++ range tree and never passing them to Python vs
            the solution where we store the precursors in the rangetree, pass
            them to Python and then evaluate them.
            """

            verbose = True
            verbose = False
            print '\n'*1
            print "Comparing fully integrated solution (c_integrated.wrap_all)"
            par = self.par
            cursor = self.cursor

            all_precursors = self.precursors_to_evaluate
            ## shuffle(all_precursors)
            all_precursors = all_precursors[:self.limit_large]

            self.myprecursors.build_parent_id_lookup()
            testrange = self.myprecursors.build_rangetree()
            import c_rangetree
            r = c_rangetree.ExtendedRangetree_Q1_RT.create()
            r.new_rangetree()
            r.create_tree(tuple(self.alltuples_isotope_correction))
            #c_integrated.create_tree(tuple(self.alltuples_isotope_correction))

            MAX_UIS = 5
            newtime = 0; oldtime = 0; ctime = 0
            oldsql = 0; newsql = 0
            alllocaltime = 0
            localprecursor = 0
            transitiontime = 0
            c_fromprecursortime = 0
            prepare  = []
            self._cursor = False
            print "Running experiment ", par.get_common_filename()
            print "calc_trans. = time to calculate the transitions of the precursor"
            print "old = use rangetree to get precursors, use C++ code to get collperpep"
            print "new = use rangetree to get precursors, use single C++ code to get collperpep"
            print "i\tcalc_trans.\tnew\t\told\t\t>>\tspeedup"
            for kk, precursor in enumerate(all_precursors):
                ii = kk + 1

                q1 =       precursor.q1
                ssrcalc =  precursor.ssrcalc
                sequence = precursor.modified_sequence
                peptide_key = precursor.transition_group
                p_id = precursor.parent_id

                q3_low, q3_high = par.get_q3range_transitions()

                #new way to calculate the precursors
                mystart = time.time()
                transitions = c_getnonuis.calculate_transitions_ch(
                    ((q1, sequence, p_id),), [1,2], q3_low, q3_high)
                nr_transitions = len( transitions )
                #fake some srm_id for the transitions
                transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])

                q1_low = q1 - par.q1_window 
                q1_high = q1 + par.q1_window
                innerstart = time.time()
                #correct rounding errors, s.t. we get the same results as before!
                ssrcalc_low = ssrcalc - par.ssrcalc_window + 0.001
                ssrcalc_high = ssrcalc + par.ssrcalc_window - 0.001
                transitiontime += time.time() - mystart

                isotope_correction = par.isotopes_up_to * Residues.mass_diffC13 / min(par.parent_charges)

                ###################################
                # New way to calculate non_uislist 
                #  - start out from transitions, get non_uislist

                mystart = time.time()
                newresult = c_integrated.wrap_all_magic(transitions, q1_low, ssrcalc_low,
                    q1_high,  ssrcalc_high, peptide_key,
                    min(MAX_UIS,nr_transitions) , par.q3_window, #q3_low, q3_high,
                    par.ppm, par.isotopes_up_to, isotope_correction, par, r)
                newtime += time.time() - mystart

                ###################################
                # Old way to calculate non_uislist:
                #  - get collisions per peptide
                #  - get non_uislist

                mystart = time.time()
                collisions_per_peptide_obj = self.myprecursors.get_collisions_per_peptide_from_rangetree(
                    precursor, precursor.q1 - par.q1_window, precursor.q1 + par.q1_window, 
                    transitions, par)

                ## # if False:
                ## #     precursor_ids = tuple(c_rangetree.query_tree( q1_low, ssrcalc_low, 
                ## #         q1_high,  ssrcalc_high, par.isotopes_up_to, isotope_correction)) 
                ## #     precursors = tuple([parentid_lookup[myid[0]] for myid in precursor_ids
                ## #                         #dont select myself 
                ## #                        if parentid_lookup[myid[0]][2]  != peptide_key])
                ## #     collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide( 
                ## #         transitions, precursors, q3_low, q3_high, par.q3_window, par.ppm)

                non_uis_list = [{} for i in range(MAX_UIS+1)]
                for order in range(1,MAX_UIS+1):
                    non_uis_list[order] = c_getnonuis.get_non_uis(
                        collisions_per_peptide_obj, order)

                oldtime += time.time() - mystart

                non_uis_list_old_way = [set(kk.keys()) for kk in non_uis_list]
                non_uis_list_len = [len(kk) for kk in non_uis_list_old_way[1:]]

                ###################################
                # Assert equality, print out result
                self.assertEqual( newresult , non_uis_list_len) 

                mys =  "%s\t%0.4fms\t%0.2fms\t\t%0.2fms\t\t>>\t%0.2f" % \
                (ii, transitiontime *1000/ ii, 
                        newtime*1000/ii, oldtime*1000/ii,
                oldtime *1.0 / newtime)
                #start a new line for each output?
                #mys += '\t%s\t%s' % ( len(precursors), len(precursor_new) )
                if False: mys += '\n'

                self.ESC = chr(27)
                sys.stdout.write(self.ESC + '[2K')
                if self._cursor:
                    sys.stdout.write(self.ESC + '[u')
                self._cursor = True
                sys.stdout.write(self.ESC + '[s')
                sys.stdout.write(mys)
                sys.stdout.flush()


            # except Exception as inst:
            #     print "something went wrong"
            #     print inst

    #@attr('comparison') 
    def test_two_table_mysql(self):
            """Compare the two table vs the one table MySQL approach
            
            Here we are comparing the old (querying the transitions database as
            well as the precursor database) and the new way (only query the
            precursor database and calculate the transitions on the fly) way of
            calculating the collisions.
            """

            print '\n'*1
            print "Comparing one vs two table MySQL solution"
            par = self.par
            cursor = self.cursor

            mycollider = collider.SRMcollider()
            mypepids = _get_unique_pepids(par, cursor)
            self.mycollider.pepids = mypepids
            self.mycollider.calcinner = 0
            shuffle( self.mycollider.pepids )
            self.mycollider.pepids = self.mycollider.pepids[:self.limit]

            MAX_UIS = 5
            c_newuistime = 0; oldtime = 0; c_fromprecursortime = 0
            oldsql = 0; newsql = 0
            oldcalctime = 0; localsql = 0
            self._cursor = False
            print "oldtime = get UIS from collisions and transitions (getting all collisions from the transitions db)"
            print "cuis + oldsql = as oldtime but calculate UIS in C++"
            print "py+newsql = only get the precursors from the db and calculate collisions in Python"
            print "ctime + newsql = only get the precursors from the db and calculate collisions in C++"
            print "new = use fast SQL and C++ code"
            print "old = use slow SQL and Python code"
            print "i\toldtime\tcuis+oldsql\tpy+newsql\tctime+newsql\t>>>\toldsql\tnewsql\t...\t...\tspeedup"
            for kk, pep in enumerate(self.mycollider.pepids):
                ii = kk + 1
                p_id = pep['parent_id']
                q1 = pep['q1']
                q3_low, q3_high = par.get_q3range_transitions()
                precursor = Precursor(q1=pep['q1'], transition_group=pep['transition_group'], parent_id = pep['parent_id'], modified_sequence=pep['mod_sequence'], ssrcalc=pep['ssrcalc'])
                transitions = collider.calculate_transitions_ch(
                    ((q1, pep['mod_sequence'], p_id),), [1], q3_low, q3_high)
                #fake some srm_id for the transitions
                transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
                ##### transitions = self.mycollider._get_all_transitions(par, pep, cursor)
                nr_transitions = len( transitions )
                if nr_transitions == 0: continue #no transitions in this window
                #
                mystart = time.time()
                collisions = list(self.mycollider._get_all_collisions_calculate_new(par, precursor, cursor))
                oldcolllen = len(collisions)
                oldcalctime += time.time() - mystart
                #
                mystart = time.time()
                collisions = _get_all_collisions(mycollider, par, pep, cursor, transitions = transitions)
                oldcsqllen = len(collisions)
                oldsql += time.time() - mystart
                #
                par.query2_add = ' and isotope_nr = 0 ' # due to the new handling of isotopes
                mystart = time.time()
                self.mycollider.mysqlnewtime= 0
                precursors = self.mycollider._get_all_precursors(par, precursor, cursor)
                newsql += time.time() - mystart
                #
                mystart = time.time()
                #precursors = self.mycollider._get_all_precursors(par, pep, cursor_l)
                localsql += time.time() - mystart
                par.query2_add = '' # due to the new handling of isotopes
                #
                mystart = time.time()
                non_uis_list = get_non_UIS_from_transitions(transitions, 
                                            tuple(collisions), par, MAX_UIS)
                cnewuis = non_uis_list
                c_newuistime += time.time() - mystart
                #
                mystart = time.time()
                non_uis_list = get_non_UIS_from_transitions_old(transitions, 
                                            collisions, par, MAX_UIS)
                oldnonuislist = non_uis_list
                oldtime += time.time() - mystart
                #
                mystart = time.time()
                q3_low, q3_high = par.get_q3range_transitions()
                collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide_other_ion_series( 
                    transitions, precursors, par, q3_low, q3_high, par.q3_window, par.ppm, False)
                non_uis_list = [ {} for i in range(MAX_UIS+1)]
                for order in range(1,MAX_UIS+1):
                    non_uis_list[order] = c_getnonuis.get_non_uis(
                        collisions_per_peptide, order)
                c_fromprecursortime += time.time() - mystart

                newl = [len(n) for n in non_uis_list]
                self.assertEqual(newl, [len(o) for o in cnewuis])
                self.assertEqual(newl, [len(o) for o in oldnonuislist])

                non_uis_list = [set(k.keys()) for k in non_uis_list]
                cnewuis = [set(k.keys()) for k in cnewuis]

                self.assertEqual(non_uis_list, cnewuis)
                self.assertEqual(non_uis_list, oldnonuislist)

                mys =  "%s\t%0.fms\t%0.fms\t\t%0.fms\t\t%0.2fms\t\t>>>\t%0.fms\t%0.2fms\t...\t...\t%0.2f" % \
                 (ii,  #i
                 (oldtime + oldsql)*1000/ii,  #oldtime
                 (c_newuistime+oldsql)*1000/ii, #cuis + oldsql
                 (oldcalctime + newsql + oldtime)*1000/ii,  #newsql
                 (c_fromprecursortime + newsql)*1000/ii,  #ctime+newsql
                 #(c_fromprecursortime + localsql)*1000/ii,

                 oldsql*1000/ii, #newsql
                 newsql*1000/ii, #oldsql
                 #localsql*1000/ii,
                 #oldtime / c_newuistime
                 (oldtime + oldsql) / (c_fromprecursortime + newsql)
                )

                self.ESC = chr(27)
                sys.stdout.write(self.ESC + '[2K')
                if self._cursor:
                    sys.stdout.write(self.ESC + '[u')
                self._cursor = True
                sys.stdout.write(self.ESC + '[s')
                sys.stdout.write(mys)
                sys.stdout.flush()

    @attr('comparison') 
    def test_mysql_vs_integrated(self):
            """Compare the one table MySQL approach vs the fully integrated Cpp approach
            
            Here we are comparing the old (querying the transitions database as
            well as the precursor database) and the new way (only query the
            precursor database and calculate the transitions on the fly) way of
            calculating the collisions.
            """

            print '\n'*1
            print "Comparing one table MySQL solution vs integrated solution"
            par = self.par
            cursor = self.cursor

            mypepids = [
                        {
                            'mod_sequence'  :  r[0],
                            'peptide_key' :r[1],
                            'transition_group' :r[1],
                            'parent_id' :  r[2],
                            'q1_charge' :  r[3],
                            'q1' :         r[4],
                            'ssrcalc' :    r[5],
                        }
                        for r in self.alltuples
                if r[3] == 2 #charge is 2
                and r[6] == 0 #isotope is 0
                and r[4] >= self.min_q1
                and r[4] < self.max_q1
            ]

            mycollider = collider.SRMcollider()
            #mypepids = _get_unique_pepids(par, cursor)
            self.mycollider.pepids = mypepids
            self.mycollider.calcinner = 0
            shuffle( self.mycollider.pepids )
            self.mycollider.pepids = self.mycollider.pepids[:self.limit]

            import c_rangetree
            r = c_rangetree.ExtendedRangetree_Q1_RT.create()
            r.new_rangetree()
            r.create_tree(tuple(self.alltuples_isotope_correction))
            #c_integrated.create_tree(tuple(self.alltuples_isotope_correction))

            MAX_UIS = 5
            c_newuistime = 0; oldtime = 0; c_fromprecursortime = 0
            oldsql = 0; newsql = 0
            newtime = 0
            oldcalctime = 0; localsql = 0
            self._cursor = False
            print "i\toldtime\t\tnewtime\t>>\tspeedup"
            for kk, pep in enumerate(self.mycollider.pepids):
                ii = kk + 1
                p_id = pep['parent_id']
                q1 = pep['q1']
                q3_low, q3_high = par.get_q3range_transitions()
                q1_low = q1 - par.q1_window 
                q1_high = q1 + par.q1_window
                ssrcalc = pep['ssrcalc']
                peptide_key = pep['peptide_key']

                #correct rounding errors, s.t. we get the same results as before!
                ssrcalc_low = ssrcalc - par.ssrcalc_window + 0.001
                ssrcalc_high = ssrcalc + par.ssrcalc_window - 0.001
                isotope_correction = par.isotopes_up_to * Residues.mass_diffC13 / min(par.parent_charges)

                precursor = Precursor(q1=pep['q1'], transition_group=pep['transition_group'], parent_id = pep['parent_id'], modified_sequence=pep['mod_sequence'], ssrcalc=pep['ssrcalc'])

                transitions = collider.calculate_transitions_ch(
                    ((q1, pep['mod_sequence'], p_id),), [1], q3_low, q3_high)
                #fake some srm_id for the transitions
                transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
                ##### transitions = self.mycollider._get_all_transitions(par, pep, cursor)
                nr_transitions = len( transitions )
                if nr_transitions == 0: continue #no transitions in this window

                ###################################
                # Old way to calculate non_uislist 
                #  - get all precursors from the SQL database
                #  - calculate collisions per peptide in C++

                par.query2_add = ' and isotope_nr = 0 ' # due to the new handling of isotopes
                mystart = time.time()
                self.mycollider.mysqlnewtime= 0
                precursors = self.mycollider._get_all_precursors(par, precursor, cursor)
                newsql += time.time() - mystart

                mystart = time.time()
                q3_low, q3_high = par.get_q3range_transitions()
                collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide_other_ion_series( 
                    transitions, precursors, par, q3_low, q3_high, par.q3_window, par.ppm, False)
                non_uis_list = [ {} for i in range(MAX_UIS+1)]
                for order in range(1,MAX_UIS+1):
                    non_uis_list[order] = c_getnonuis.get_non_uis(
                        collisions_per_peptide, order)
                c_fromprecursortime += time.time() - mystart

                newl = [len(n) for n in non_uis_list]
                non_uis_list_old_way = [set(kk.keys()) for kk in non_uis_list]
                non_uis_list_len = [len(kk) for kk in non_uis_list_old_way[1:]]

                ###################################
                # New way to calculate non_uislist 
                #  - start out from transitions, get non_uislist
                mystart = time.time()
                newresult = c_integrated.wrap_all_magic(transitions, q1_low, ssrcalc_low,
                    q1_high,  ssrcalc_high, peptide_key,
                    min(MAX_UIS,nr_transitions) , par.q3_window, #q3_low, q3_high,
                    par.ppm, par.isotopes_up_to, isotope_correction, par, r)
                newtime += time.time() - mystart

                ###################################
                # Assert equality, print out result
                self.assertEqual( newresult , non_uis_list_len) 

                mys =  "%s\t%0.1fms\t\t%0.2fms\t>>>\t%0.1f" % \
                 (ii,  #i
                 (c_fromprecursortime + newsql)*1000/ii,  # oldtime
                 (newtime)*1000/ii, # newtime
                 (c_fromprecursortime + newsql) / (newtime), # speedup
                )

                self.ESC = chr(27)
                sys.stdout.write(self.ESC + '[2K')
                if self._cursor:
                    sys.stdout.write(self.ESC + '[u')
                self._cursor = True
                sys.stdout.write(self.ESC + '[s')
                sys.stdout.write(mys)
                sys.stdout.flush()

if __name__ == '__main__':
    unittest.main()

