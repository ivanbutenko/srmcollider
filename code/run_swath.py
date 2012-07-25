#!/usr/bin/python
# -*- coding: utf-8  -*-
# vim:set fdm=marker:
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
select count(distinct parent_key) from result_srmuis where exp_key = 104;
select count(*) from srmPeptides_mouse where isotope_nr = 0 and q1_charge = 2 and q1 between 400 and 1400;

python prepare_rangetree.py 104 60000 9999 hroest.srmPeptides_mouse
sh /tmp/tmp.sh


"""
import MySQLdb
import time
import sys 
import c_rangetree, c_getnonuis
sys.path.append( '/home/hroest/lib/hlib/' )
sys.path.append( '/home/hroest/projects' )
sys.path.append( '/home/hroest/projects/hlib' )
#db_l = MySQLdb.connect(read_default_file="~/.my.cnf.local")
db = MySQLdb.connect(read_default_file="~/.my.cnf")
cursor = db.cursor()
#cursor_l = db_l.cursor()
#cursor = cursor_l
import collider
import progress

exp_key = -1234
min_q1 = 1047
max_q1 = 1205
ssrcalcwin = 0.25

#11234556 1047 1205 0.25 hroest.srmPeptides_yeast

from optparse import OptionParser, OptionGroup
usage = "usage: %prog experiment_key startQ1 endQ1 [options]"
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "Run uis Options",
                    "None yet")
group.add_option("--swath_mode", action='store_true', dest="swath_mode", default=False,
                  help="SWATH mode enabled (use fixed window)")
group.add_option("--use_db", action='store_true', dest="use_db", default=False,
                  help="Use db instead of rangetree (slower, but necessary for large Q1 windows)")
group.add_option("--dry_run", action='store_true', dest="dry_run", default=False,
                  help="Only a dry run, do not start processing (but create experiment)")
group.add_option("--restable", dest="restable", default='hroest.result_srmuis', type="str",
                  help="MySQL result table" + 
                  "Defaults to result_srmuis") 
group.add_option("--insert",
                  action="store_true", dest="insert_mysql", default=False,
                  help="Insert into mysql experiments table")
parser.add_option_group(group)


#Run the collider
###########################################################################
#Parse options
db = MySQLdb.connect(read_default_file="~/.my.cnf")
cursor = db.cursor()
par = collider.SRM_parameters()
par.parse_cmdl_args(parser)
options, args = parser.parse_args(sys.argv[1:])

#local arguments
exp_key = sys.argv[1]
min_q1 = float(sys.argv[2])
max_q1 = float(sys.argv[3])
use_db = options.use_db
swath_mode = options.swath_mode
restable = options.restable
dry_run = options.dry_run
#ssrcalcwin = float(sys.argv[4])
par.__dict__.update( options.__dict__ )
par.q3_range = [options.q3_low, options.q3_high]
par.q1_window /= 2.0
par.q3_window /= 2.0
par.ssrcalc_window /= 2.0
if par.ppm == 'True': par.ppm = True
else: par.ppm = False
#par.q1_window = 1.2 / 2.0 #UIS paper = 1.2
#par.q3_window = 2.0 / 2.0 #UIS paper = 2.0
#par.ssrcalc_window = ssrcalcwin / 2.0
#par.ppm = False
#par.considerIsotopes = True
#par.max_uis = 5 #turn off uis if set to 0
par.dontdo2p2f = False #also look at 2+ parent / 2+ fragment ions
#par.q3_range = [400, 1400]
par.eval()
#print par.experiment_type
print par.get_common_filename()
mycollider = collider.SRMcollider()

#print 'par', par.ppm, type(par.ppm)
print 'isotopes' , par.considerIsotopes

qadd = ''
if not par.considerIsotopes: qadd = 'and isotope_nr = 0'

qadd += ' and isotope_nr in (0,1) '
par.query_add += ' and isotope_nr in (0,1) '
par.query2_add += ' and isotope_nr in (0,1) '


print "swath ", options.swath_mode
print "usedb ", use_db
print "add ", par.query_add
print "add 1", par.query1_add
print "add 2", par.query2_add


assert(len(par.peptide_tables) == 1)

if options.insert_mysql:
    common_filename = par.get_common_filename()
    superkey = 31
    if common_filename.split('_')[0] == 'human': superkey = 34
    query = """
    insert into hroest.experiment  (name, short_description,
    description, comment1, comment2,comment3,  super_experiment_key, ddb_experiment_key)
    VALUES (
        'ludovic_swath', '%s', '%s', '%s', '%s','%s', %s, 0
    )
    """ %( common_filename + '_' + par.peptide_tables[0].split('.')[1], 
          par.experiment_type, par.peptide_tables[0], par.transition_table,
          #comment3
          'q1: %s; q3 %s ; isotopes 0,1' % (par.q1_window, par.q3_window),
          superkey)
    cursor.execute(query)
    exp_key = db.insert_id()
    print "Inserted into mysql db with id ", exp_key

start = time.time()
q =  """
select modified_sequence, peptide_key, parent_id, q1_charge, q1, ssrcalc, isotope_nr
from %(peptide_table)s where q1 between %(lowq1)s and %(highq1)s
%(qadd)s
""" % {'peptide_table' : par.peptide_tables[0], 
              'lowq1'  : min_q1 ,
              'highq1' : max_q1,
              'qadd'   : qadd
} 
print q
cursor.execute(q)

#print "finished query execution: ", time.time() - start
alltuples =  list(cursor.fetchall() )

#print "finished query fetch: ", time.time() - start
mypepids = [
            {
                'sequence'  :  r[0],
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

# we use a rangetree to select the precursors
if not use_db:
    parentid_lookup = [ [ r[2], (r[4], r[0], r[1]) ] 
            for r in alltuples ]
    parentid_lookup  = dict(parentid_lookup)
    print "building tree with %s Nodes" % len(alltuples)
    c_rangetree.create_tree(tuple(alltuples))



# in SWATH mode, select all precursors that are relevant for the background at
# once
if swath_mode: 
    values="q1, modified_sequence, peptide_key, q1_charge, ssrcalc"
    q1_low = min_q1; q1_high = max_q1
    query2 = """
    select %(values)s
    from %(pep)s
    where q1 >= %(q1_low)s and q1 <= %(q1_high)s
    %(query_add)s
    """ % { 'q1_high' : q1_high, 'q1_low'  : q1_low,
           'query_add' : par.query2_add,
           'pep' : par.peptide_tables[0],
           'values' : values}
    print 'swath ' , query2
    cursor.execute( query2 )
    allprecursors =  cursor.fetchall()

if dry_run:
    sys.exit()

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
    #new way to calculate the precursors
    transitions = c_getnonuis.calculate_transitions_ch(
        ((q1, pep['sequence'], p_id),), [1], q3_low, q3_high)
    #fake some srm_id for the transitions
    transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
    nr_transitions = len(transitions)
    #
    #in SWATH mode, we have a fixed window independent of q1
    if swath_mode: q1_low = min_q1; q1_high = max_q1
    else: q1_low = q1 - par.q1_window ; q1_high = q1 + par.q1_window
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
        # Use the rangetree
        # correct rounding errors, s.t. we get the same results as before!
        # use the rangetree
        ssrcalc_low = ssrcalc - par.ssrcalc_window + 0.001
        ssrcalc_high = ssrcalc + par.ssrcalc_window - 0.001
        precursor_ids = tuple(c_rangetree.query_tree( q1_low, ssrcalc_low, 
                                                     q1_high,  ssrcalc_high )  )
        precursors = [parentid_lookup[myid[0]] for myid in precursor_ids
                            #dont select myself 
                           if parentid_lookup[myid[0]][2]  != pep['peptide_key']]
    #now calculate collisions in this area of the space
    colldensity = c_getnonuis.calculate_density( 
        transitions, precursors, q3_low, q3_high, par.q3_window, par.ppm)
    assert MAX_UIS == 1
    nonuseable = len( [c for c in colldensity if c >0 ] )
    prepare.append( (nonuseable, nr_transitions, p_id , 1, exp_key) )
    progressm.update(1)


cursor.executemany('insert into %s' % restable + ' (non_useable_UIS, total_UIS, \
                  parent_key, uisorder, exp_key) values (%s,%s,%s,%s,%s)' , prepare)


