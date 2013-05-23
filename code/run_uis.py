#!/usr/bin/python
# -*- coding: utf-8  -*-
"""
 *
 * Program       : SRMCollider
 * Author        : Hannes Roest <roest@imsb.biol.ethz.ch>
 * Date          : 05.02.2011 
 *
 *
 * Copyright (C) 2011 - 2012 Hannes Roest
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307, USA
 *
"""

"""
This program will process all peptides of a given proteome and compare it to
the background in that proteome. It will output the number of UIS and the
number of total combinations for each order up to the specified limit for each
precursor.

Expected results on the test-database (see README and test-folder on how to set up the sqlite-testdatabase):

python run_uis.py 123456789 400 1500 --peptide_table=srmPeptides_test --max_uis 5 -i 3 --q1_window=1 --q3_window=1 --ssrcalc_window=10 --sqlite_database=/tmp/srmcollider_testdb 
[------------------------------------------------------------>] 100%  9519.7 peptides/sec (eta 0s)
It took 0s
Analyzed 905 peptides
Order 1, Average non useable UIS 0.055404336098
Order 2, Average non useable UIS 0.00377300709314
Order 3, Average non useable UIS 0.000427182046582
Order 4, Average non useable UIS 4.72421355715e-05
Order 5, Average non useable UIS 3.86353977514e-06

python run_uis.py 123456789 500 525 --peptide_table=srmPeptides_test --max_uis 5 -i 3 --q1_window=25 --q3_window=1 --ssrcalc_window=10 --sqlite_database=/tmp/srmcollider_testdb --swath_mode
[------------------------------------------------------------>] 100%  4886.1 peptides/sec (eta 0s)
It took 0s
Analyzed 39 peptides
Order 1, Average non useable UIS 0.478557767019
Order 2, Average non useable UIS 0.0377468685161
Order 3, Average non useable UIS 0.00329392829393
Order 4, Average non useable UIS 0.0002035002035
Order 5, Average non useable UIS 0.0

"""

import sys 
import collider, progress
from copy import copy

# count the number of interfering peptides
count_avg_transitions = False

from optparse import OptionParser, OptionGroup
usage = "usage: %prog experiment_key startQ1 endQ1 [options]"
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "Run uis Options")
group.add_option("--swath_mode", action='store_true', dest="swath_mode", default=False,
                  help="SWATH mode enabled (use fixed window)")
group.add_option("--use_db", action='store_true', dest="use_db", default=False,
                  help="Use db instead of rangetree")
group.add_option("--restable", dest="restable", default='srmcollider.result_srmuis', type="str",
                  help="MySQL result table" + " - Defaults to result_srmuis") 
group.add_option("--insert",
                  action="store_true", dest="insert_mysql", default=False,
                  help="Insert into mysql experiments table")
group.add_option("--query_peptide_table", type="str", help="Peptide table to get query peptides from")
parser.add_option_group(group)

# Run the collider
###########################################################################
# Parse options
mycollider = collider.SRMcollider()
par = collider.SRM_parameters()
par.parse_cmdl_args(parser)
options, args = parser.parse_args(sys.argv[1:])
par.parse_options(options)
par.eval()

if len(sys.argv) < 4: 
    raise Exception("Wrong number of arguments (provide at least experiment_key startQ1 endQ1")
    sys.exit()

#local arguments
exp_key = sys.argv[1]
min_q1 = float(sys.argv[2])
max_q1 = float(sys.argv[3])
use_db = options.use_db
swath_mode = options.swath_mode
restable = options.restable

########
# Sanity check for SWATH : the provided window needs to be the SWATH window
if swath_mode:
    if not par.q1_window*2 >= max_q1 - min_q1:
        raise Exception("Your Q1 window needs to be at least as large as the min/max q1 you chose.")
        sys.exit()

###########################################################################
# Prepare the collider

db = par.get_db()
cursor = db.cursor()

if options.insert_mysql:
    assert False # you have to implement this yourself

# Get the precursors
###########################################################################
from precursor import Precursors
myprecursors = Precursors()
myprecursors.getFromDB(par, db.cursor(), min_q1 - par.q1_window, max_q1 + par.q1_window)
if not options.query_peptide_table is None and not options.query_peptide_table == "":
  print "Using a different table for the query peptides than for the background peptides!"
  print "Will use table %s " % options.query_peptide_table
  query_precursors = Precursors()
  query_par = copy(par)
  query_par.peptide_tables = [options.query_peptide_table]
  query_precursors.getFromDB(query_par, db.cursor(), min_q1 - par.q1_window, max_q1 + par.q1_window)
  precursors_to_evaluate = query_precursors.getPrecursorsToEvaluate(min_q1, max_q1)
else:
  precursors_to_evaluate = myprecursors.getPrecursorsToEvaluate(min_q1, max_q1)
myprecursors.build_parent_id_lookup()
myprecursors.build_transition_group_lookup()

# If we dont use the DB, we use the rangetree to query and get our list of
# precursors that are interfering. In SWATH we dont include a +/- q1_window
# around our range or precursors because the precursor window is fixed to
# (min_q1,max_q1) and no other precursors are considered.
if not use_db and swath_mode: 
  myprecursors.getFromDB(par, cursor, min_q1, max_q1)
  rtree = myprecursors.build_rangetree()
elif not use_db:
  #myprecursors.getFromDB(par, cursor, min_q1 - par.q1_window, max_q1 + par.q1_window)
  rtree = myprecursors.build_rangetree()

# In SWATH mode, select all precursors that are relevant for the background at
# once. Select all precursors between min_q1 - correction and max_q1 and then
# only append those to all_swath_precursors that are actually included in the
# range (min_q1,max_q1) with at least one isotope.
par.query2_add = ''
if swath_mode and use_db: 
    temp_precursors = Precursors()
    temp_precursors.getFromDB(par, db.cursor(), min_q1 - par.calculate_isotope_correction(), max_q1)
    all_swath_precursors = []
    for p in temp_precursors.precursors:
      if(p.included_in_isotopic_range(min_q1, max_q1, par) ): 
        all_swath_precursors.append(p)

allintertr = []
MAX_UIS = par.max_uis
progressm = progress.ProgressMeter(total=len(precursors_to_evaluate), unit='peptides')
prepare  = []
print "analyzing %s peptides" % len(precursors_to_evaluate)
for precursor in precursors_to_evaluate:

    transitions = precursor.calculate_transitions_from_param(par)
    nr_transitions = len(transitions)

    if use_db and not swath_mode:
        # Case 1: regular SRMCollider, get transitions from the DB
        precursors_obj = mycollider._get_all_precursors(par, precursor, cursor)
        collisions_per_peptide = collider.get_coll_per_peptide_from_precursors(mycollider, 
                transitions, precursors_obj, par, precursor)
    elif use_db and swath_mode:
        # Case 2: SRMCollider in SWATH mode, get transitions from the DB (e.g. the memory in this case)
        if par.ssrcalc_window > 1000:
            precursors_obj = [p for p in all_swath_precursors if p.transition_group != precursor.transition_group]
        else:
            ssrcalc_low =  precursor.ssrcalc - par.ssrcalc_window 
            ssrcalc_high = precursor.ssrcalc + par.ssrcalc_window 
            precursors_obj = [p for p in all_swath_precursors if p.transition_group != precursor.transition_group
                         and p.ssrcalc > ssrcalc_low and p.ssrcalc < ssrcalc_high ]
        collisions_per_peptide = collider.get_coll_per_peptide_from_precursors(mycollider, 
                transitions, precursors_obj, par, precursor)
    elif not use_db:
        if swath_mode:
        # Case 3: SRMCollider in SWATH mode, get transitions from the rangetree
          collisions_per_peptide = myprecursors.get_collisions_per_peptide_from_rangetree(
              precursor, min_q1, max_q1, transitions, par, rtree)
        else:
          # Case 4: regular SRMCollider, get transitions from the rangetree
          collisions_per_peptide = myprecursors.get_collisions_per_peptide_from_rangetree(
              precursor, precursor.q1 - par.q1_window, precursor.q1 + par.q1_window, 
              transitions, par, rtree)

    non_uis_list = collider.get_nonuis_list(collisions_per_peptide, MAX_UIS)
    ## 
    ## Lets count the number of peptides that interfere
    if count_avg_transitions:
        tr_arr = [0 for i in range(nr_transitions)]
        for v in collisions_per_peptide.values():
            for vv in v:
                tr_arr[vv] += 1
        allintertr.extend( tr_arr )
    ##
    #
    for order in range(1,min(MAX_UIS+1, nr_transitions+1)): 
        prepare.append( (len(non_uis_list[order]), collider.choose(nr_transitions, 
            order), precursor.parent_id , order, exp_key) )
    progressm.update(1)

if count_avg_transitions:
    print "\n"
    print "found %s transitions" % len(allintertr)#
    print "found max of %s interferences" % max(allintertr)
    print "found average of %s interferences" % ( sum(allintertr) * 1.0 / len(allintertr) )

# if any problems with the packet/buffer length occur, try this:
## set global max_allowed_packet=1000000000;
## set global net_buffer_length=1000000;
# cursor.executemany('insert into %s' % restable + ' (non_useable_UIS, total_UIS, \
#                   parent_key, uisorder, exp_key) values (%s,%s,%s,%s,%s)' , prepare)

print "Analyzed %s peptides" % len(precursors_to_evaluate)
for order in range(1,MAX_UIS+1):
    sum_all = sum([p[0]*1.0/p[1] for p in prepare if p[3] == order]) 
    nr_peptides = len([p for p in prepare if p[3] == order])
    if not par.quiet and not nr_peptides ==0:
      print "Order %s, Average non useable UIS %s" % (order, sum_all *1.0/ nr_peptides)
    # cursor.execute("insert into hroest.result_completegraph_aggr (sum_nonUIS, nr_peptides, uisorder, experiment) VALUES (%s,%s,%s,'%s')" % (sum_all, nr_peptides, order, exp_key))

