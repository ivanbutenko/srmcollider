#!/usr/bin/python
# -*- coding: utf-8  -*-
"""
 *
 * Program       : SRMCollider
 * Author        : Hannes Roest <roest@imsb.biol.ethz.ch>
 * Date          : 05.02.2011 
 *
 *
 * Copyright (C) 2011 Hannes Roest
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

$ python run_uis.py 123456789 400 1500 --peptide_table=srmPeptides_test --max_uis 5 -i 3 --q1_window=1 --q3_window=1 --ssrcalc_window=10 --sqlite_database=/tmp/testdb 
[------------------------------------------------------------>] 100%  9519.7 peptides/sec (eta 0s)
It took 0s
Analyzed 905 peptides
Order 1, Average non useable UIS 0.055404336098
Order 2, Average non useable UIS 0.00377300709314
Order 3, Average non useable UIS 0.000427182046582
Order 4, Average non useable UIS 4.72421355715e-05
Order 5, Average non useable UIS 3.86353977514e-06

"""

import sys 
import collider
import progress

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
#group.add_option("--dry_run", action='store_true', dest="dry_run", default=False,
#                  help="Only a dry run, do not start processing (but create experiment)")
group.add_option("--restable", dest="restable", default='srmcollider.result_srmuis', type="str",
                  help="MySQL result table" + " - Defaults to result_srmuis") 
group.add_option("--insert",
                  action="store_true", dest="insert_mysql", default=False,
                  help="Insert into mysql experiments table")
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
    print "wrong number of arguments"
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
    if not par.q1_window >= max_q1 - min_q1:
        raise Exception("Your Q1 window needs to be at least as large as the min/max q1 you chose.")
        sys.exit()

###########################################################################
# Prepare the collider

import Residues
R = Residues.Residues('mono')
db = par.get_db()
cursor = db.cursor()

if options.insert_mysql:
    common_filename = par.get_common_filename()
    superkey = 31
    if common_filename.split('_')[0] == 'human': superkey = 34
    query = """
    insert into srmcollider.experiment  (name, short_description,
    description, comment1, comment2,comment3,  super_experiment_key, ddb_experiment_key)
    VALUES (
        'ludovic_swath', '%s', '%s', '%s', '%s','%s', %s, 0
    )
    """ %( common_filename + '_' + par.peptide_table.split('.')[1], 
          par.experiment_type, par.peptide_table, par.transition_table,
          #comment3
          'q1: %s; q3 %s ; isotopes 0,1' % (par.q1_window, par.q3_window),
          superkey)
    cursor.execute(query)
    exp_key = db.insert_id()
    print "Inserted into mysql db with id ", exp_key


# Get the precursors
###########################################################################
from precursor import Precursors
myprecursors = Precursors()
myprecursors.getFromDB(par, db.cursor(), min_q1 - par.q1_window, max_q1 + par.q1_window)
testrange = myprecursors.build_rangetree()
precursors_to_evaluate = myprecursors.getPrecursorsToEvaluate(min_q1, max_q1)
myprecursors.build_parent_id_lookup()
myprecursors.build_transition_group_lookup()

print "analyzing %s peptides" % len(precursors_to_evaluate)

# If we dont use the DB, we use the rangetree to query and get our list of
# precursors that are interfering. In SWATH we dont include a +/- q1_window
# around our range or precursors because the precursor window is fixed to
# (min_q1,max_q1) and no other precursors are considered.
if not use_db and swath_mode: 
  myprecursors.getFromDB(par, cursor, min_q1, max_q1)
  testrange = myprecursors.build_rangetree()
elif not use_db:
  myprecursors.getFromDB(par, cursor, min_q1 - par.q1_window, max_q1 + par.q1_window)
  testrange = myprecursors.build_rangetree()

# In SWATH mode, select all precursors that are relevant for the background at
# once. Select all precursors between min_q1 - correction and max_q1 and then
# only append those to all_swath_precursors that are actually included in the
# range (min_q1,max_q1) with at least one isotope.
par.query2_add = ''
if use_db and swath_mode: 
    isotope_correction = par.isotopes_up_to * R.mass_diffC13 / min(par.parent_charges)
    temp_precursors = Precursors()
    temp_precursors.getFromDB(par, db.cursor(), min_q1 - isotope_correction, max_q1)
    all_swath_precursors = []
    for p in temp_precursors.precursors:
      if(p.included_in_isotopic_range(min_q1, max_q1, par, R) ): 
        all_swath_precursors.append(p)

# nr_interfering_prec = [ [] for i in range(7)]

allintertr = []
MAX_UIS = par.max_uis
progressm = progress.ProgressMeter(total=len(precursors_to_evaluate), unit='peptides')
prepare  = []
for precursor in precursors_to_evaluate:

    q3_low, q3_high = par.get_q3range_transitions()
    transitions = precursor.calculate_transitions(q3_low, q3_high)
    nr_transitions = len(transitions)

    if use_db and not swath_mode:
        precursors_obj = mycollider._get_all_precursors_obj(par, precursor, cursor)
        collisions_per_peptide = collider.get_coll_per_peptide_from_precursors_obj_wrapper(mycollider, 
                transitions, precursors_obj, par, precursor)
    elif use_db and swath_mode:
        if par.ssrcalc_window > 1000:
            precursors_obj = [p for p in all_swath_precursors if p.transition_group != precursor.transition_group]
        else:
            ssrcalc_low =  precursor.ssrcalc - par.ssrcalc_window 
            ssrcalc_high = precursor.ssrcalc + par.ssrcalc_window 
            precursors_obj = [p for p in all_swath_precursors if p.transition_group != precursor.transition_group
                         and p.ssrcalc > ssrcalc_low and p.ssrcalc < ssrcalc_high ]
        collisions_per_peptide = collider.get_coll_per_peptide_from_precursors_obj_wrapper(mycollider, 
                transitions, precursors_obj, par, precursor)
    elif not use_db:
        # Use the rangetree, whether it is swath or not
        collisions_per_peptide = myprecursors.get_collisions_per_peptide_from_rangetree(precursor, transitions, par)

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

for order in range(1,6):
    sum_all = sum([p[0]*1.0/p[1] for p in prepare if p[3] == order]) 
    nr_peptides = len([p for p in prepare if p[3] == order])
    if not par.quiet: print sum_all *1.0/ nr_peptides
    cursor.execute("insert into hroest.result_completegraph_aggr (sum_nonUIS, nr_peptides, uisorder, experiment) VALUES (%s,%s,%s,'%s')" % (sum_all, nr_peptides, order, exp_key))

