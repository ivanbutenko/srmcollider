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
"""

import MySQLdb, sys 
import c_integrated, c_rangetree, c_getnonuis
sys.path.append( '/home/hroest/lib/hlib/' )
sys.path.append( '/home/hroest/projects' )
sys.path.append( '/home/hroest/projects/hlib' )
import collider
import progress

from optparse import OptionParser, OptionGroup
usage = "usage: %prog experiment_key startQ1 endQ1 [options]"
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "Run uis Options",
                    "None yet")
parser.add_option_group(group)

#Run the collider
###########################################################################
#Parse options
mycollider = collider.SRMcollider()
par = collider.SRM_parameters()
par.parse_cmdl_args(parser)
options, args = parser.parse_args(sys.argv[1:])
par.parse_options(options)
par.dontdo2p2f = False #also look at 2+ parent / 2+ fragment ions
par.eval()
print par.get_common_filename()

if len(sys.argv) < 4: 
    print "wrong number of arguments"
    sys.exit()

#local arguments
exp_key = sys.argv[1]
min_q1 = float(sys.argv[2])
max_q1 = float(sys.argv[3])

if False: #use_sqlite:
    import sqlite
    db = sqlite.connect(sqlite_database)
    cursor = db.cursor()
else:
    db = MySQLdb.connect(read_default_file=par.mysql_config)
    cursor = db.cursor()

if par.max_uis ==0: 
    print "Please change --max_uis option, 0 does not make sense here"
    sys.exit()

import Residues
R = Residues.Residues('mono')
isotope_correction = par.isotopes_up_to * R.mass_diffC13 / min(par.parent_charges)
q =  """
select modified_sequence, transition_group, parent_id, q1_charge, q1, ssrcalc, modifications, missed_cleavages, isotopically_modified
from %(peptide_table)s where q1 between %(lowq1)s - %(isotope_correction)s and %(highq1)s
""" % {'peptide_table' : par.peptide_table, 
              'lowq1'  : min_q1 - par.q1_window, 
              'highq1' : max_q1 + par.q1_window,
               'isotope_correction' : isotope_correction
      } 
#print q
cursor.execute(q)

alltuples =  list(cursor.fetchall() )
mypepids = [
            {
                'sequence'  :  r[0],
                'transition_group' :r[1],
                'parent_id' :  r[2],
                'q1_charge' :  r[3],
                'q1' :         r[4],
                'ssrcalc' :    r[5],
                'isotope_mod': r[8],
            }
            for r in alltuples
    if r[3] == 2  # charge is 2
    and r[6] == 0 # no modification 
    and r[7] == 0 # no missed cleavages 
    and r[4] >= min_q1
    and r[4] < max_q1
]

print "building tree with %s Nodes" % len(alltuples)
c_integrated.create_tree(tuple(alltuples))

self = mycollider
self.mysqlnewtime = 0
self.pepids = mypepids
MAX_UIS = par.max_uis
progressm = progress.ProgressMeter(total=len(self.pepids), unit='peptides')
prepare  = []
for pep in self.pepids:
    p_id = pep['parent_id']
    q1 = pep['q1']
    ssrcalc = pep['ssrcalc']
    q3_low, q3_high = par.get_q3range_transitions()
    #
    transitions = c_getnonuis.calculate_transitions_ch(
        ((q1, pep['sequence'], p_id),), [1], q3_low, q3_high)
    transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
    nr_transitions = len(transitions)
    #correct rounding errors, s.t. we get the same results as before!
    ssrcalc_low = ssrcalc - par.ssrcalc_window + 0.001
    ssrcalc_high = ssrcalc + par.ssrcalc_window - 0.001
    try:
        result = c_integrated.wrap_all_magic(transitions, q1 - par.q1_window, 
            ssrcalc_low, q1 + par.q1_window,  ssrcalc_high, pep['transition_group'],  
            min(MAX_UIS,nr_transitions) , par.q3_window, par.ppm, par.isotopes_up_to, isotope_correction, par)
    except ValueError: print "Too many transitions for", pep['sequence']; continue
    for order in range(1,min(MAX_UIS+1, nr_transitions+1)): 
        prepare.append( (result[order-1], collider.choose(nr_transitions, 
            order), p_id , order, exp_key)  )
    progressm.update(1)

for order in range(1,6):
    sum_all = sum([p[0]*1.0/p[1] for p in prepare if p[3] == order]) 
    nr_peptides = len([p for p in prepare if p[3] == order])
    if not par.quiet: print sum_all *1.0/ nr_peptides
    cursor.execute("insert into hroest.result_completegraph_aggr (sum_nonUIS, nr_peptides, uisorder, experiment) VALUES (%s,%s,%s,'%s')" % (sum_all, nr_peptides, order, exp_key))

"""
create table hroest.result_completegraph (
exp_key          int(11),
parent_key       int(11),
non_useable_UIS  int(11),
total_UIS        int(11),
uisorder         int(4) 
)

"""
