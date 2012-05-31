import unittest
"""
This file tests the functionality of the collider.py module. 
"""

import sys
sys.path.extend(['.', '..','../..' '../external/', 'external/'])
sys.path.extend(['test/'])
import collider

import test

from test_shared import *
import test_shared 
import time
from Residues import Residues

import precursor
from precursor import Precursors

class Test_integration_run_uis(unittest.TestCase): 

  def setUp(self):

    self.transitions = transitions_def1
    self.collisions  = collisions_def1

    self.EPSILON = 10**-5

    self.min_q1 = 400
    self.max_q1 = 1500

    par = collider.SRM_parameters()
    par.q1_window = 1 / 2.0
    par.q3_window = 1 / 2.0
    par.ssrcalc_window = 10 / 2.0
    par.ppm = False
    par.isotopes_up_to = 3
    par.q3_low = 400
    par.q3_high = 1400
    par.max_uis = 5
    par.peptide_table = 'srmPeptides_test'
    par.mysql_config = '~/.my.cnf'
    par.sqlite_database = test_shared.SQLITE_DATABASE_LOCATION
    par.use_sqlite = True
    par.quiet = False

    par.bions      =  True
    par.yions      =  True
    par.aions      =  False
    par.aMinusNH3  =  False
    par.bMinusH2O  =  False
    par.bMinusNH3  =  False
    par.bPlusH2O   =  False
    par.yMinusH2O  =  False
    par.yMinusNH3  =  False
    par.cions      =  False
    par.xions      =  False
    par.zions      =  False
    par.MMinusH2O  =  False
    par.MMinusNH3  =  False
    par.q3_range = [par.q3_low, par.q3_high]
    par.set_default_vars()
    par.eval()

    self.par = par
    self.R = Residues('mono')

    self.acollider = collider.SRMcollider()
    self.aparamset = collider.testcase()

    self.db = par.get_db()

    # Get the precursors
    ###########################################################################
    myprecursors = Precursors()
    cursor = self.db.cursor()
    myprecursors.getFromDB(par, cursor, self.min_q1 - par.q1_window, self.max_q1 + par.q1_window)
    testrange = myprecursors.build_rangetree()
    self.precursors_to_evaluate = myprecursors.getPrecursorsToEvaluate(self.min_q1, self.max_q1)
    myprecursors.build_parent_id_lookup()
    myprecursors.build_transition_group_lookup()
    self.myprecursors = myprecursors
    cursor.close()

  def tearDown(self):
    self.db.close()

  def get_final_report(self, par, prepare):
    final_report = {}
    for order in range(1,par.max_uis+1):
        sum_all = sum([p[0]*1.0/p[1] for p in prepare if p[3] == order]) 
        nr_peptides = len([p for p in prepare if p[3] == order])
        if not nr_peptides ==0: print "Order %s, Average non useable UIS %s" % (order, sum_all *1.0/ nr_peptides)
        final_report[order] = sum_all*1.0/nr_peptides
    return final_report

  def test_runuis_nonswath(self):
    self.assertEqual(len(self.precursors_to_evaluate), 905)
    par = self.par
    cursor = self.db.cursor()
    prepare = []

    for precursor in self.precursors_to_evaluate:

      q3_low, q3_high = par.get_q3range_transitions()
      transitions = precursor.calculate_transitions(q3_low, q3_high)
      nr_transitions = len(transitions)

      precursors_obj = self.acollider._get_all_precursors(par, precursor, cursor)
      collisions_per_peptide = collider.get_coll_per_peptide_from_precursors(self.acollider, 
                transitions, precursors_obj, par, precursor)
      non_uis_list = collider.get_nonuis_list(collisions_per_peptide, par.max_uis)

      for order in range(1,min(par.max_uis+1, nr_transitions+1)): 
        prepare.append( (len(non_uis_list[order]), collider.choose(nr_transitions, 
          order), precursor.parent_id , order, -1) )

    self.assertEqual(len(prepare), 905*par.max_uis)
    self.assertEqual(prepare[0], (0, 17.0, 1, 1, -1)) 

    final_report = self.get_final_report(par, prepare)
    self.check_final_report_nonswath(final_report)

  def test_runuis_nonswath_rangetree(self):
    self.assertEqual(len(self.precursors_to_evaluate), 905)
    par = self.par
    cursor = self.db.cursor()

    # If we dont use the DB, we use the rangetree to query and get our list of
    # precursors that are interfering. In SWATH we dont include a +/- q1_window
    # around our range or precursors because the precursor window is fixed to
    # (min_q1,max_q1) and no other precursors are considered.
    self.myprecursors.getFromDB(par, cursor, self.min_q1 - par.q1_window, self.max_q1 + par.q1_window)
    rtree = self.myprecursors.build_rangetree()

    prepare = []

    for precursor in self.precursors_to_evaluate:

      q3_low, q3_high = par.get_q3range_transitions()
      transitions = precursor.calculate_transitions(q3_low, q3_high)
      nr_transitions = len(transitions)

      # Use the rangetree, whether it is swath or not
      collisions_per_peptide = self.myprecursors.get_collisions_per_peptide_from_rangetree(precursor, precursor.q1 - par.q1_window, precursor.q1 + par.q1_window, transitions, par, rtree)
      non_uis_list = collider.get_nonuis_list(collisions_per_peptide, par.max_uis)

      for order in range(1,min(par.max_uis+1, nr_transitions+1)): 
        prepare.append( (len(non_uis_list[order]), collider.choose(nr_transitions, 
          order), precursor.parent_id , order, -1) )

    self.assertEqual(len(prepare), 905*par.max_uis)
    self.assertEqual(prepare[0], (0, 17.0, 1, 1, -1)) 
    final_report = self.get_final_report(par, prepare)
    self.check_final_report_nonswath(final_report)

  def test_runuis_swath(self):

    self.assertEqual(len(self.precursors_to_evaluate), 905)
    swath_mode = False
    par = self.par
    R = self.R
    cursor = self.db.cursor()
    prepare = []

    self.min_q1 = 500
    self.max_q1 = 525

    # Get the precursors (now for 500-525 instead of the full range)
    ###########################################################################
    myprecursors = Precursors()
    cursor = self.db.cursor()
    myprecursors.getFromDB(par, cursor, self.min_q1 - par.q1_window, self.max_q1 + par.q1_window)
    rtree = myprecursors.build_rangetree()
    self.precursors_to_evaluate = myprecursors.getPrecursorsToEvaluate(self.min_q1, self.max_q1)
    self.assertEqual(len(self.precursors_to_evaluate), 39)

    isotope_correction = par.isotopes_up_to * R.mass_diffC13 / min(par.parent_charges)
    temp_precursors = Precursors()
    temp_precursors.getFromDB(par, self.db.cursor(), self.min_q1 - isotope_correction, self.max_q1)
    all_swath_precursors = []
    for p in temp_precursors.precursors:
      if(p.included_in_isotopic_range(self.min_q1, self.max_q1, par, R) ): 
        all_swath_precursors.append(p)

    for precursor in self.precursors_to_evaluate:

      q3_low, q3_high = par.get_q3range_transitions()
      transitions = precursor.calculate_transitions(q3_low, q3_high)
      nr_transitions = len(transitions)


      if par.ssrcalc_window > 1000:
          precursors_obj = [p for p in all_swath_precursors if p.transition_group != precursor.transition_group]
      else:
          ssrcalc_low =  precursor.ssrcalc - par.ssrcalc_window 
          ssrcalc_high = precursor.ssrcalc + par.ssrcalc_window 
          precursors_obj = [p for p in all_swath_precursors if p.transition_group != precursor.transition_group
                       and p.ssrcalc > ssrcalc_low and p.ssrcalc < ssrcalc_high ]
      collisions_per_peptide = collider.get_coll_per_peptide_from_precursors(self.acollider, 
              transitions, precursors_obj, par, precursor)


      non_uis_list = collider.get_nonuis_list(collisions_per_peptide, par.max_uis)

      for order in range(1,min(par.max_uis+1, nr_transitions+1)): 
        prepare.append( (len(non_uis_list[order]), collider.choose(nr_transitions, 
          order), precursor.parent_id , order, -1) )

    self.assertEqual(len(prepare), 39*par.max_uis)
    self.assertEqual(prepare[0], (5, 8.0, 69, 1, -1) )

    final_report = self.get_final_report(par, prepare)
    self.check_final_report_swath(final_report)

  def test_runuis_swath_rangetree(self):

    self.assertEqual(len(self.precursors_to_evaluate), 905)
    swath_mode = False
    par = self.par
    R = self.R
    cursor = self.db.cursor()
    prepare = []

    self.min_q1 = 500
    self.max_q1 = 525

    # Get the precursors (now for 500-525 instead of the full range)
    ###########################################################################
    myprecursors = Precursors()
    cursor = self.db.cursor()
    myprecursors.getFromDB(par, cursor, self.min_q1 - par.q1_window, self.max_q1 + par.q1_window)
    rtree = myprecursors.build_rangetree()
    self.precursors_to_evaluate = myprecursors.getPrecursorsToEvaluate(self.min_q1, self.max_q1)
    self.assertEqual(len(self.precursors_to_evaluate), 39)

    # If we dont use the DB, we use the rangetree to query and get our list of
    # precursors that are interfering. In SWATH we dont include a +/- q1_window
    # around our range or precursors because the precursor window is fixed to
    # (min_q1,max_q1) and no other precursors are considered.
    myprecursors.getFromDB(par, cursor, self.min_q1, self.max_q1)
    rtree = myprecursors.build_rangetree()

    for precursor in self.precursors_to_evaluate:

      q3_low, q3_high = par.get_q3range_transitions()
      transitions = precursor.calculate_transitions(q3_low, q3_high)
      nr_transitions = len(transitions)

      # Use the rangetree, whether it is swath or not
      collisions_per_peptide = self.myprecursors.get_collisions_per_peptide_from_rangetree(precursor, self.min_q1, self.max_q1, transitions, par, rtree)
      non_uis_list = collider.get_nonuis_list(collisions_per_peptide, par.max_uis)

      for order in range(1,min(par.max_uis+1, nr_transitions+1)): 
        prepare.append( (len(non_uis_list[order]), collider.choose(nr_transitions, 
          order), precursor.parent_id , order, -1) )

    self.assertEqual(len(prepare), 39*par.max_uis)
    # self.assertEqual(prepare[0], (5, 8.0, 69, 1, -1) )

    final_report = self.get_final_report(par, prepare)
    self.check_final_report_swath(final_report)

  def check_final_report_swath(self, final_report):

    self.assertTrue(abs(final_report[1] - 0.478557767019   ) < self.EPSILON )
    self.assertTrue(abs(final_report[2] - 0.0377468685161  ) < self.EPSILON )
    self.assertTrue(abs(final_report[3] - 0.00329392829393 ) < self.EPSILON )
    self.assertTrue(abs(final_report[4] - 0.0002035002035  ) < self.EPSILON )
    self.assertTrue(abs(final_report[5] - 0.0              ) < self.EPSILON )

  def check_final_report_nonswath(self, final_report):
    self.assertTrue(abs(final_report[1] - 0.055404336098     ) < self.EPSILON )
    self.assertTrue(abs(final_report[2] - 0.00377300709314   ) < self.EPSILON )
    self.assertTrue(abs(final_report[3] - 0.000427182046582  ) < self.EPSILON )
    self.assertTrue(abs(final_report[4] - 4.72421355715e-05  ) < self.EPSILON )
    self.assertTrue(abs(final_report[5] - 3.86353977514e-06  ) < self.EPSILON )

if __name__ == '__main__':
    unittest.main()

