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

import MySQLdb, time, sys 
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

old_tables_with_isotopes = False

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

if use_sqlite:
    import sqlite
    db = sqlite.connect(sqlite_database)
    cursor = db.cursor()
else:
    db = MySQLdb.connect(read_default_file=options.mysql_config)
    cursor = db.cursor()

print 'isotopes' , par.isotopes_up_to

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

start = time.time()
q = """
select modified_sequence, peptide_key, parent_id, q1_charge, q1, ssrcalc, 0 as isotope_nr
from %(peptide_table)s where q1 between %(lowq1)s and %(highq1)s
""" % {'peptide_table' : par.peptide_table, 
              'lowq1'  : min_q1 - par.q1_window, 
              'highq1' : max_q1 + par.q1_window,
      }
if old_tables_with_isotopes: 
    # only of the old tables are used
    q = """
    select modified_sequence, peptide_key, parent_id, q1_charge, q1, ssrcalc, isotope_nr
    from %(peptide_table)s where q1 between %(lowq1)s and %(highq1)s
    and isotope_nr <= %(nr_isotopes)s
    """ % {'peptide_table' : par.peptide_table, 
                  'lowq1'  : min_q1 - par.q1_window, 
                  'highq1' : max_q1 + par.q1_window,
                  'highq1' : max_q1 + par.q1_window,
                  'highq1' : max_q1 + par.q1_window,
                   'nr_isotopes' :  par.isotopes_up_to
          }
print q
cursor.execute( q )
alltuples =  list(cursor.fetchall() )

mypepids = [
            {
                'mod_sequence'  :  r[0],
                'peptide_key' :r[1],
                'parent_id' :  r[2],
                'q1_charge' :  r[3],
                'q1' :         r[4],
                'ssrcalc' :    r[5],
            }
            for r in alltuples
    if r[3] == 2 #charge is 2
    and r[6] == 0 #isotope is 0
    and r[4] >= min_q1
    and r[4] < max_q1
]

print "analyzing %s peptides" % len(mypepids)

# use the rangetree?
if not use_db:
    if not old_tables_with_isotopes: 
        if swath_mode: this_min = min_q1; this_max = max_q1
        else: this_min = min_q1 - par.q1_window; this_max = max_q1 + par.q1_window
        import Residues
        R = Residues.Residues('mono')
        isotope_correction = par.isotopes_up_to * R.mass_diffC13 / min(par.parent_charges)
        q = """
        select modified_sequence, peptide_key, parent_id, q1_charge, q1, ssrcalc, 0 as isotope_nr
        from %(peptide_table)s where q1 between %(lowq1)s - %(isotope_correction)s and %(highq1)s
        """ % {'peptide_table' : par.peptide_table, 
                      'lowq1'  : this_min,
                      'highq1' : this_max,
                      'isotope_correction' : isotope_correction
              }
        print q
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
    parentid_lookup = [ [ r[2], (r[4], r[0], r[1], r[3]) ] 
            for r in alltuples ]
    parentid_lookup  = dict(parentid_lookup)
    print "building tree with %s Nodes" % len(alltuples)
    if use_sqlite: alltuples = [tuple(t) for t in alltuples]
    c_rangetree.create_tree(tuple(alltuples))

# in SWATH mode, select all precursors that are relevant for the background at
# once
par.query2_add = ''
if use_db and swath_mode: 
    import Residues
    R = Residues.Residues('mono')
    isotope_correction = par.isotopes_up_to * R.mass_diffC13 / min(par.parent_charges)
    values="q1, modified_sequence, peptide_key, q1_charge, ssrcalc"
    q1_low = min_q1; q1_high = max_q1
    query2 = """
    select %(values)s
    from %(pep)s
    where q1 >= %(q1_low)s - %(correction)s and q1 <= %(q1_high)s
    """ % { 'q1_high' : q1_high, 'q1_low'  : q1_low,
           'correction' : isotope_correction,
           'pep' : par.peptide_table,
           'values' : values}
    print 'swath ' , query2
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

allintertr = []
self = mycollider
self.mysqlnewtime = 0
self.pepids = mypepids
MAX_UIS = par.max_uis
progressm = progress.ProgressMeter(total=len(self.pepids), unit='peptides')
prepare  = []
for kk, pep in enumerate(self.pepids):
    p_id = pep['parent_id']
    q1 = pep['q1']
    ssrcalc = pep['ssrcalc']
    q3_low, q3_high = par.get_q3range_transitions()
    #
    transitions = collider.calculate_transitions_ch(
        ((q1, pep['mod_sequence'], p_id),), [1], q3_low, q3_high)
    #fake some srm_id for the transitions
    transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
    nr_transitions = len(transitions)
    #
    # in SWATH mode, we have a fixed window independent of q1
    # in regular mode, we have to account for the isotope correction
    import Residues
    R = Residues.Residues('mono')
    isotope_correction = par.isotopes_up_to * R.mass_diffC13 / min(par.parent_charges)
    if swath_mode: q1_low = min_q1; q1_high = max_q1; q1_low -= isotope_correction
    else: q1_low = q1 - par.q1_window -isotope_correction ; q1_high = q1 + par.q1_window
    #
    if use_db and not swath_mode:
        precursors = self._get_all_precursors(par, pep, cursor)
    elif use_db and swath_mode:
        if par.ssrcalc_window > 1000:
            precursors = [p for p in allprecursors if p[2] != pep['peptide_key'] ]
        else:
            ssrcalc_low = ssrcalc - par.ssrcalc_window 
            ssrcalc_high = ssrcalc + par.ssrcalc_window 
            precursors = [p for p in allprecursors if p[2] != pep['peptide_key'] 
                         and p[4] > ssrcalc_low and p[4] < ssrcalc_high ]
    elif not use_db:
        # Use the rangetree, whether it is swath or not
        #correct rounding errors, s.t. we get the same results as before!
        ssrcalc_low = ssrcalc - par.ssrcalc_window + 0.001
        ssrcalc_high = ssrcalc + par.ssrcalc_window - 0.001
        precursor_ids = tuple(c_rangetree.query_tree( q1_low, ssrcalc_low, 
                                                     q1_high,  ssrcalc_high )  )
        if swath_mode:
            precursors = tuple([parentid_lookup[myid[0]] for myid in precursor_ids
                                #dont select myself 
                               if parentid_lookup[myid[0]][2]  != pep['peptide_key']])
        else:
            precursors = []
            # filter out wrong isotopes
            for myid in precursor_ids:
              append = False
              r = parentid_lookup[myid[0]]
              ch = r[3]
              for iso in range(par.isotopes_up_to+1):
                if (r[0] + (R.mass_diffC13 * iso)/ch > q1 - par.q1_window and 
                    r[0] + (R.mass_diffC13 * iso)/ch < q1 + par.q1_window): append=True
              if(append and r[2]  != pep['peptide_key']): precursors.append(r)

    import c_getnonuis
    precursors = tuple(precursors)
    collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide( 
            transitions, precursors, q3_low, q3_high, par.q3_window, par.ppm)
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
            order), p_id , order, exp_key) )
    progressm.update(1)

if count_avg_transitions:
    print "\n"
    print "found %s transitions" % len(allintertr)#
    print "found max of %s interferences" % max(allintertr)
    print "found average of %s interferences" % ( sum(allintertr) * 1.0 / len(allintertr) )

# if any problems with the packet/buffer length occur, try this:
## set global max_allowed_packet=1000000000;
## set global net_buffer_length=1000000;
cursor.executemany('insert into %s' % restable + ' (non_useable_UIS, total_UIS, \
                  parent_key, uisorder, exp_key) values (%s,%s,%s,%s,%s)' , prepare)


