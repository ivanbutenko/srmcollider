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

import sys 
import collider
import progress

# count the number of interfering peptides
count_avg_transitions = False

from optparse import OptionParser, OptionGroup
usage = "usage: %prog experiment_key startQ1 endQ1 [options]"
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "Run uis Options",
                    "None yet")
group.add_option("--swath_mode", action='store_true', dest="swath_mode", default=False,
                  help="SWATH mode enabled (use fixed window)")
group.add_option("--use_db", action='store_true', dest="use_db", default=False,
                  help="Use db instead of rangetree")
group.add_option("--dry_run", action='store_true', dest="dry_run", default=False,
                  help="Only a dry run, do not start processing (but create experiment)")
group.add_option("--restable", dest="restable", default='srmcollider.result_srmuis', type="str",
                  help="MySQL result table" + 
                  "Defaults to result_srmuis") 
group.add_option("--sqlite_database", dest="sqlite_database", default='',
                  help="Use specified sqlite database instead of MySQL database" )
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
par.dontdo2p2f = False #also look at 2+ parent / 2+ fragment ions
par.eval()
if not par.quiet: print par.get_common_filename()
if not par.quiet: print "only b y ", par.do_b_y_only()
if not par.quiet: print par.aions

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
sqlite_database = options.sqlite_database
if sqlite_database != '': use_sqlite = True
else: use_sqlite = False

########
# Sanity check for SWATH : the provided window needs to be the SWATH window
if swath_mode:
    if not par.q1_window >= max_q1 - min_q1:
        raise Exception("Your Q1 window needs to be at least as large as the min/max q1 you chose.")
        sys.exit()

###########################################################################
# Prepare the collider

if use_sqlite:
    import sqlite
    db = sqlite.connect(sqlite_database)
    cursor = db.cursor()
else:
    import MySQLdb
    db = MySQLdb.connect(read_default_file=par.mysql_config)
    cursor = db.cursor()

if not par.quiet: print 'isotopes' , par.isotopes_up_to

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
myprecursors.getFromDB(par, db.cursor(), min_q1, max_q1)
testrange = myprecursors.build_rangetree()

precursors_to_evaluate = [p for p in myprecursors.precursors 
                         if p.q1_charge == 2 
                         and p.modifications == 0
                         and p.missed_cleavages == 0 
                         and p.q1 >= min_q1
                         and p.q1 <= max_q1
                         ]

myprecursors.build_parent_id_lookup()
myprecursors.build_transition_group_lookup()

print "analyzing %s peptides" % len(precursors_to_evaluate)

# use the rangetree?
if not use_db:
    if True:
        if swath_mode: this_min = min_q1; this_max = max_q1
        else: this_min = min_q1 - par.q1_window; this_max = max_q1 + par.q1_window
        import Residues
        R = Residues.Residues('mono')
        isotope_correction = par.isotopes_up_to * R.mass_diffC13 / min(par.parent_charges)
        q = """
        select modified_sequence, transition_group, parent_id, q1_charge, q1, ssrcalc, 0 as isotope_nr, isotopically_modified
        from %(peptide_table)s where q1 between %(lowq1)s - %(isotope_correction)s and %(highq1)s
        """ % {'peptide_table' : par.peptide_table, 
                      'lowq1'  : this_min,
                      'highq1' : this_max,
                      'isotope_correction' : isotope_correction
              }
        if not par.quiet: print q
        cursor.execute( q )
        if not swath_mode:
            alltuples = tuple(cursor.fetchall() )
        else:
            # filter out wrong isotopes
            result = cursor.fetchall() 
            new_result = []
            for r in result:
              append = False
              ch = r[3]
              for iso in range(par.isotopes_up_to+1):
                if (r[4] + (R.mass_diffC13 * iso)/ch > min_q1 and 
                    r[4] + (R.mass_diffC13 * iso)/ch < max_q1): append=True
              if(append): new_result.append(r)
            alltuples = tuple(new_result)
    import c_rangetree
    # parent_id : q1, sequence, trans_group, q1_charge
    parentid_lookup = [ [ r[2], (r[4], r[0], r[1], r[3], r[7]) ] 
            for r in alltuples ]
    parentid_lookup  = dict(parentid_lookup)
    print "building tree with %s Nodes" % len(alltuples)
    if use_sqlite: alltuples = [tuple(t) for t in alltuples]
    # print [p for p in alltuples if p[1] == 554689]
    c_rangetree.create_tree(tuple(alltuples))

# in SWATH mode, select all precursors that are relevant for the background at
# once
par.query2_add = ''
if use_db and swath_mode: 
    import Residues
    R = Residues.Residues('mono')
    isotope_correction = par.isotopes_up_to * R.mass_diffC13 / min(par.parent_charges)
    values="q1, modified_sequence, transition_group, q1_charge, isotopically_modified, ssrcalc"
    q1_low = min_q1; q1_high = max_q1
    query2 = """
    select %(values)s
    from %(pep)s
    where q1 >= %(q1_low)s - %(correction)s and q1 <= %(q1_high)s
    """ % { 'q1_high' : q1_high, 'q1_low'  : q1_low,
           'correction' : isotope_correction,
           'pep' : par.peptide_table,
           'values' : values}
    if not par.quiet: print 'swath ' , query2
    cursor.execute( query2 )
    result = cursor.fetchall() 
    # filter out wrong isotopes
    new_result = []
    for r in result:
      append = False
      ch = r[3]
      for iso in range(par.isotopes_up_to+1):
        if (r[0] + (R.mass_diffC13 * iso)/ch > min_q1 and 
            r[0] + (R.mass_diffC13 * iso)/ch < max_q1): append=True
      if(append): new_result.append(r)
    allprecursors = new_result

    temp_precursors = Precursors()
    temp_precursors.getFromDB(par, db.cursor(), min_q1 - isotope_correction, max_q1)
    all_swath_precursors = []
    for p in temp_precursors.precursors:
      if(p.included_in_isotopic_range(min_q1, max_q1, par, R) ): 
        all_swath_precursors.append(p)





nr_interfering_prec = [ [] for i in range(7)]

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

