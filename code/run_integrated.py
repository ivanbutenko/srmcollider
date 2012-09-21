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

python run_integrated.py 123456789 400 1500 --peptide_table=srmPeptides_test --max_uis 5 -i 3 --q1_window=1 --q3_window=1 --ssrcalc_window=10 --sqlite_database=/tmp/srmcollider_testdb 
analyzing 905 peptides
[------------------------------------------------------------>] 100%  9828.8 peptides/sec (eta 0s)
It took 0s
Analyzed 905 peptides
Order 1, Average non useable UIS 0.055404336098
Order 2, Average non useable UIS 0.00377300709314
Order 3, Average non useable UIS 0.000427182046582
Order 4, Average non useable UIS 4.72421355715e-05
Order 5, Average non useable UIS 3.86353977514e-06

"""

import sys 
import c_integrated, collider, progress
from precursor import Precursors
from copy import copy

from optparse import OptionParser, OptionGroup
usage = "usage: %prog experiment_key startQ1 endQ1 [options]"
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "Run integrated Options")
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
print par.get_common_filename()

if len(sys.argv) < 4: 
    print "wrong number of arguments"
    sys.exit()

#local arguments
exp_key = sys.argv[1]
min_q1 = float(sys.argv[2])
max_q1 = float(sys.argv[3])
db = par.get_db()
cursor = db.cursor()

if par.max_uis ==0: 
    print "Please change --max_uis option, 0 does not make sense here"
    sys.exit()

# Get the precursors
###########################################################################
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
isotope_correction = par.calculate_isotope_correction()
r_tree = myprecursors.build_extended_rangetree ()

print "Will evaluate %s precursors" % len(precursors_to_evaluate)
progressm = progress.ProgressMeter(total=len(precursors_to_evaluate), unit='peptides')
prepare  = []
for precursor in precursors_to_evaluate:
    transitions = precursor.calculate_transitions_from_param(par)
    #correct rounding errors, s.t. we get the same results as before!
    ssrcalc_low = precursor.ssrcalc - par.ssrcalc_window + 0.001
    ssrcalc_high = precursor.ssrcalc + par.ssrcalc_window - 0.001
    try:
        result = c_integrated.wrap_all_bitwise(transitions,
            precursor.q1 - par.q1_window, ssrcalc_low, precursor.q1 + par.q1_window,  ssrcalc_high,
            precursor.transition_group, min(par.max_uis,len(transitions)), par.q3_window, par.ppm,
            par.isotopes_up_to, isotope_correction, par, r_tree)
    except ValueError: 
      print "Too many transitions for", precursor.modification
      continue

    for order in range(1,min(par.max_uis+1, len(transitions)+1)): 
        prepare.append( (result[order-1], collider.choose(len(transitions), 
            order), precursor.parent_id , order, exp_key)  )
    #//break;
    progressm.update(1)

for order in range(1,6):
    sum_all = sum([p[0]*1.0/p[1] for p in prepare if p[3] == order]) 
    nr_peptides = len([p for p in prepare if p[3] == order])
    if not par.quiet and not nr_peptides ==0: print "Order %s, Average non useable UIS %s" % (order, sum_all *1.0/ nr_peptides)
    #cursor.execute("insert into hroest.result_completegraph_aggr (sum_nonUIS, nr_peptides, uisorder, experiment) VALUES (%s,%s,%s,'%s')" % (sum_all, nr_peptides, order, exp_key))

"""

create table hroest.result_completegraph (
exp_key          int(11),
parent_key       int(11),
non_useable_UIS  int(11),
total_UIS        int(11),
uisorder         int(4) 
)


"""
