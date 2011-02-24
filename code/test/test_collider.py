import unittest

import sys
sys.path.append( '..')
import collider

from test_shared import *
import test_shared 

#
# inc     means it is included in another test
# nc      means non critical (print fxn etc)
# db      means tested in test_db
# OK      means tested here
#
###################################
###################################
#         def eval(self):
#         def get_copy(self):
#         def get_q3range_transitions(self):
#         def get_q3range_collisions(self):
#         def get_q3_window_transitions(self, q3):
#         def get_common_filename(self):
#         def transition_db(self): return self.transition_table.split('.')[0]
#         def transition_tbl(self): return self.transition_table.split('.')[1]
#         def peptide_db(self): return self.peptide_table.split('.')[0]
#         def peptide_tbl(self): 
#         def peptide_tbl_identifier(self): return self.peptide_tbl[12:] #cut off 'srmPeptides'
#     
#     
# db      def find_clashes_small(self, db, par, use_per_transition=False,
# db      def _get_all_precursors(self, par, pep, cursor, 
# db      def _get_all_collisions_calculate(self, par, pep, cursor, 
# OK      def _get_all_collisions_calculate_sub(self, precursors, R, q3_low, q3_high):
# db      def find_clashes(self, db, par, toptrans=False, pepids=None, 
#         def find_clashes_toptrans_paola(self, db, par,
# OK      def _getMinNeededTransitions(self, par, transitions, collisions):
#         def find_clashes_toptrans_3strike(self, db, par, pepids=None, 
#         def _get_unique_pepids(self, par, cursor, ignore_genomeoccurence=False):
#         def _get_unique_pepids_toptransitions(self, par, cursor):
#         def _get_all_transitions_toptransitions(self, par, pep, cursor, values = 'q3, m.id'):
#         def _get_all_transitions(self, par, pep, cursor, values = "q3, srm_id"):
# db      def _get_all_collisions_per_transition(self, par, pep, transitions, cursor):
# db      def _get_all_collisions(self, par, pep, cursor, 
# incl     def _get_collisions_per_transition(self, par, pep, q3, cursor, 
# nc      def store_object(self, par):
# nc      def store_in_file(self, par):
# nc      def load_from_file(self, par, directory):
# nc      def print_unique_histogram(self, par):
# nc      def print_cumm_unique(self, par):
# nc      def print_cumm_unique_all(self, par, cursor):
# nc          #def print_cumm_unique(self, par, mydist, filename):
# nc      def print_q3min(self, par):
# nc      def print_q3min_ppm(self, par):
# nc      def print_q1all(self, par, bars = 50):
# nc      def print_q3all_ppm(self, par, bars = 50):
# nc      def print_stats(self):
#     def get_cum_dist(original):
# OK  def get_non_UIS_from_transitions(transitions, collisions, par, MAX_UIS, 
# db def get_coll_per_peptide(self, transitions, par, pep):
# OK  def get_non_UIS_from_transitions_old(transitions, collisions, par, MAX_UIS, unsorted=False):
# OK  def get_UIS_from_transitions(transitions, collisions, par, MAX_UIS):
# nc  def get_peptide_from_table(t, row):
# ??  def insert_peptide_in_db(self, db, peptide_table, transition_group):
# ??  def get_actual_mass(self):
# ??  def insert_in_db(self, db, fragment_charge, transition_table):
# ??  def fast_insert_in_db(self, db, fragment_charge, transition_table, transition_group):
# ??  def all_calculate_clashes_in_series_insert_db( S, S2, pairs_dic, 
# ??  def calculate_clashes_in_series_insert_db(S, S2, charge1, charge2, pairs_dic, 
# ??  def reset_pairs_unique(mass_bins):
# ??  def calculate_clashes_in_series(S, S2, charge1, charge2, pairs_dic, 
# inc def get_non_uis(pepc, non_uis, order):
# inc def get_non_uis_unsorted(pepc, non_uis, order):
# OK  def choose(i,r):
# OK  def get_uis(srm_ids, non_uis, order):
#     def permutations(iterable, r=None):
#     def _permutations(iterable, r=None):
# inc def _combinations(N, M):
# OK  def combinations(iterable, r):
#     
class Test_collider_function(unittest.TestCase): 

    def setUp(self):
        self.transitions = transitions_def1
        self.collisions  = collisions_def1
        class Minimal: pass
        self.par = Minimal()
        self.par.q3_window = 4.0
        self.par.ppm = False
        self.MAX_UIS = 5
        self.q3_high = 1500
        self.q3_low = 300

        def returnrange(): return self.q3_high, self.q3_low
        self.par.get_q3range_collisions = returnrange

        import sys
        sys.path.append( '/home/hroest/projects/' )
        sys.path.append( '/home/hroest/lib/' )
        import silver
        self.R = silver.Residues.Residues('mono')

        self.acollider = collider.SRMcollider()
        self.aparamset = collider.testcase()

    def test_getMinNeededTransitions(self):
        par = self.par
        par.max_uis = 5
        transitions = self.transitions
        collisions = self.collisions
        m = self.acollider._getMinNeededTransitions(par, transitions, collisions)
        self.assertEqual(m, 5)

        par.max_uis = 10 
        transitions = transitions_def3
        collisions = collisions_def3
        m = self.acollider._getMinNeededTransitions(par, transitions, collisions)
        self.assertEqual(m, 10)

        transitions = transitions_def4
        collisions = collisions_def4
        m = self.acollider._getMinNeededTransitions(par, transitions, collisions)
        self.assertEqual(m, 4)

        transitions = transitions_def5
        collisions = collisions_def5
        m = self.acollider._getMinNeededTransitions(par, transitions, collisions)
        self.assertEqual(m, 4)

    def test_get_non_uis(self):
        test = set()
        collider.get_non_uis( [1,2,3], test,2 )
        self.assertEqual(test, set([(1, 2), (1, 3), (2, 3)]) )
        test = set()
        collider.get_non_uis( [1,2,3,4], test,2 )
        self.assertEqual(test, set([(1, 2), (1, 3), (1, 4), (2, 3), (3, 4), (2, 4)]))

    def test_get_non_UIS_from_transitions1(self): 
        oldnon_uis = collider.get_non_UIS_from_transitions_old(self.transitions, self.collisions, 
                                         self.par, self.MAX_UIS)
        self.assertEqual([len(l) for l in oldnon_uis[1:]], test_shared.lennonuis1)
        self.assertEqual(oldnon_uis, test_shared.refnonuis1)

    def test_get_non_UIS_from_transitions2(self): 
        self.transitions = test_shared.transitions_def2
        self.collisions  = test_shared.collisions_def2
        oldnon_uis = collider.get_non_UIS_from_transitions_old(self.transitions,
            self.collisions, self.par, self.MAX_UIS)
        self.assertEqual([len(l) for l in oldnon_uis[1:]], test_shared.lennonuis2)
        self.assertEqual(oldnon_uis, test_shared.refnonuis2_sorted)

    def test_get_non_UIS_from_transitions2_unsorted(self): 
        #here we have the transitions in the wrong order
        #it should still work
        self.transitions = transitions_def2_unsorted
        self.collisions  = collisions_def2
        oldnon_uis = collider.get_non_UIS_from_transitions_old(self.transitions,
            self.collisions, self.par, self.MAX_UIS, unsorted=True)
        self.assertEqual([len(l) for l in oldnon_uis[1:]], test_shared.lennonuis2)
        self.assertEqual(oldnon_uis, test_shared.refnonuis2_unsorted)

    def test_get_non_UIS_from_transitions3(self): 
        self.transitions = transitions_def3
        self.collisions  = collisions_def3
        oldnon_uis = collider.get_non_UIS_from_transitions_old(self.transitions,
            self.collisions, self.par, self.MAX_UIS)
        self.assertEqual([len(l) for l in oldnon_uis[1:]], test_shared.lennonuis3)
        self.assertEqual(oldnon_uis, test_shared.refnonuis3)

    def test_get_non_UIS_from_transitions4(self): 
        self.transitions = transitions_def4
        self.collisions  = collisions_def4
        oldnon_uis = collider.get_non_UIS_from_transitions_old(self.transitions,
            self.collisions, self.par, self.MAX_UIS)
        self.assertEqual([len(l) for l in oldnon_uis[1:]], test_shared.lennonuis4)
        self.assertEqual(oldnon_uis, test_shared.refnonuis4)

    def test_get_uis(self):
        non_uis_list = collider.get_non_UIS_from_transitions_old(self.transitions,
                 self.collisions, self.par, self.MAX_UIS)
        srm_ids = [t[1] for t in self.transitions]
        rr = collider.get_uis(srm_ids, non_uis_list[2], 2)
        self.assertEqual(len(rr), 0)
        #
        self.transitions = transitions_def2
        self.collisions  = collisions_def2
        non_uis_list = collider.get_non_UIS_from_transitions_old(self.transitions,
                 self.collisions, self.par, self.MAX_UIS)
        srm_ids = [t[1] for t in self.transitions]
        rr = collider.get_uis(srm_ids, non_uis_list[2], 2)
        self.assertEqual(len(rr), 1)

    def test_get_UIS_from_transitions(self):
        res = collider.get_UIS_from_transitions(self.transitions, 
            self.collisions, self.par, self.MAX_UIS)
        self.assertEqual( [len(r) for r in res] , [0,0,0,0,0,0])

    def test_get_UIS_from_transitions2(self):
        self.transitions = transitions_def2
        self.collisions  = collisions_def2
        res = collider.get_UIS_from_transitions(self.transitions,
            self.collisions, self.par, self.MAX_UIS)
        self.assertEqual( [len(r) for r in res] , [0,0,1,2,1,0])

    def test_get_UIS_from_transitions3(self):
        self.transitions = transitions_def3
        self.collisions  = collisions_def3
        res = collider.get_UIS_from_transitions(self.transitions,
            self.collisions, self.par, self.MAX_UIS)
        self.assertEqual( [len(r) for r in res] , [0,0,0,0,0,0])

    def test_get_UIS_from_transitions4(self):
        self.transitions = transitions_def4
        self.collisions  = collisions_def4
        res = collider.get_UIS_from_transitions(self.transitions,
            self.collisions, self.par, self.MAX_UIS)
        self.assertEqual( [len(r) for r in res] , [0, 0, 7, 17, 15, 6])

    def test_calculate_collisions_1(self):
        pep = test_shared.runpep1
        transitions = test_shared.runtransitions1
        precursors = test_shared.runprecursors1
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        par = self.par
        q3_high = self.q3_high
        q3_low = self.q3_low
        R = self.R

        collisions = list(collider.SRMcollider._get_all_collisions_calculate_sub(
                collider.SRMcollider(), precursors, R, q3_low, q3_high))
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

        self.assertEqual(collisions_per_peptide, test_shared.collpepresult1)

    def test_calculate_collisions_2(self):
        pep = test_shared.runpep2
        transitions = test_shared.runtransitions2
        precursors = test_shared.runprecursors2
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        par = self.par
        q3_high = self.q3_high
        q3_low = self.q3_low
        R = self.R

        collisions = list(collider.SRMcollider._get_all_collisions_calculate_sub(
                collider.SRMcollider(), precursors, R, q3_low, q3_high))
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

        self.assertEqual(collisions_per_peptide, test_shared.collpepresult2)

    def test_get_uis_extra1(self):

        order =  1
        nonuis = set([(1627247L,), (1627240L,), (1627241L,), (1627242L,), (1627248L,), (1627243L,), (1627238L,), (1627249L,), (1627244L,), (1627239L,), (1627250L,), (1627251L,), (1627252L,)])
        srm_ids = [1627238L, 1627239L, 1627240L, 1627241L, 1627242L, 1627243L, 1627244L, 1627247L, 1627248L, 1627249L, 1627250L, 1627251L, 1627252L]
        uis_list = collider.get_uis(srm_ids, nonuis, order)
        self.assertEqual( len(uis_list), 0)

    def test_choose(self): 
        self.assertEqual(10,  collider.choose( 5,2) )
        self.assertEqual(10,  collider.choose( 5,3) )
        self.assertEqual(45,  collider.choose(10,2) )
        self.assertEqual(120, collider.choose(10,3) )
        self.assertEqual(210, collider.choose(10,4) )

    def test_combinations(self):
        comb52 = collider.combinations( range(5), 2 ) 
        self.assertEqual(comb52, [
            (0, 1),
            (0, 2),
            (0, 3),
            (0, 4),
            (1, 2),
            (1, 3),
            (1, 4),
            (2, 3),
            (2, 4),
            (3, 4)
        ])
        comb53 = collider.combinations( range(5), 3 ) 
        self.assertEqual( comb53, [
            (0, 1, 2),
            (0, 1, 3),
            (0, 1, 4),
            (0, 2, 3), 
            (0, 2, 4),
            (0, 3, 4),
            (1, 2, 3),
            (1, 2, 4),
            (1, 3, 4),
            (2, 3, 4)
        ])



if __name__ == '__main__':
    unittest.main()

