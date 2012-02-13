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
import c_integrated
import collider
import progress

from optparse import OptionParser, OptionGroup
usage = "usage: %prog experiment_key startQ1 endQ1 [options]"
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "Run integrated Options")
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
from precursor import Precursors
myprecursors = Precursors()
myprecursors.getFromDB(par, db.cursor(), min_q1 - par.q1_window, max_q1 + par.q1_window)
precursors_to_evaluate = myprecursors.getPrecursorsToEvaluate(min_q1, max_q1)

import Residues
R = Residues.Residues('mono')
isotope_correction = par.isotopes_up_to * R.mass_diffC13 / min(par.parent_charges)
alltuples = [ (p.modified_sequence, p.transition_group, p.parent_id, p.q1_charge, p.q1, p.ssrcalc,0,0,p.isotopically_modified) for p in myprecursors.precursors]
c_integrated.create_tree(tuple(alltuples))

MAX_UIS = par.max_uis
progressm = progress.ProgressMeter(total=len(precursors_to_evaluate), unit='peptides')
prepare  = []
for precursor in precursors_to_evaluate:

    ssrcalc = precursor.ssrcalc
    q1 = precursor.q1
    p_id = precursor.parent_id

    q3_low, q3_high = par.get_q3range_transitions()
    transitions = precursor.calculate_transitions(q3_low, q3_high)
    nr_transitions = len(transitions)

    #correct rounding errors, s.t. we get the same results as before!
    ssrcalc_low = ssrcalc - par.ssrcalc_window + 0.001
    ssrcalc_high = ssrcalc + par.ssrcalc_window - 0.001
    try:
        result = c_integrated.wrap_all_magic(transitions, q1 - par.q1_window, 
            ssrcalc_low, q1 + par.q1_window,  ssrcalc_high, precursor.transition_group,  
            min(MAX_UIS,nr_transitions) , par.q3_window, par.ppm, par.isotopes_up_to, isotope_correction, par)
    except ValueError: 
      print "Too many transitions for", precursor.modification
      continue

    for order in range(1,min(MAX_UIS+1, nr_transitions+1)): 
        prepare.append( (result[order-1], collider.choose(nr_transitions, 
            order), p_id , order, exp_key)  )
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
