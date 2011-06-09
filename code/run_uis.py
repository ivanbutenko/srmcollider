#!/usr/bin/python
# -*- coding: utf-8  -*-
# vim:set fdm=marker:
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
sys.path.append( '/home/hroest/srm_clashes/code' )
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
group.add_option("--restable", dest="restable", default='hroest.result_srmuis', type="str",
                  help="MySQL result table" + 
                  "Defaults to result_srmuis") 
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
#ssrcalcwin = float(sys.argv[4])
#peptide_table = sys.argv[5]
par.__dict__.update( options.__dict__ )
par.q3_range = [options.q3_low, options.q3_high]
par.q1_window /= 2.0
par.q3_window /= 2.0
par.ssrcalc_window /= 2.0
if par.ppm == 'True': par.ppm = True
elif par.ppm == 'False': par.ppm = False
else: 'wrong arg for ppm'; assert False
par.dontdo2p2f = False #do not look at 2+ parent / 2+ fragment ions
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
#par.peptide_table = peptide_table
print par.get_common_filename()
mycollider = collider.SRMcollider()

#print 'par', par.ppm, type(par.ppm)
print 'isotopes' , par.considerIsotopes

qadd = ''
if not par.considerIsotopes: qadd = 'and isotope_nr = 0'

start = time.time()
cursor.execute( """
select modified_sequence, peptide_key, parent_id, q1_charge, q1, ssrcalc, isotope_nr
from %(peptide_table)s where q1 between %(lowq1)s and %(highq1)s
%(qadd)s
""" % {'peptide_table' : par.peptide_table, 
              'lowq1'  : min_q1 - par.q1_window, 
              'highq1' : max_q1 + par.q1_window,
              'qadd'   : qadd
      } )

#print "finished query execution: ", time.time() - start
alltuples =  list(cursor.fetchall() )
#from random import shuffle
#shuffle(alltuples)

#print "finished query fetch: ", time.time() - start
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
if not use_db:
    parentid_lookup = [ [ r[2], (r[4], r[0], r[1]) ] 
            for r in alltuples ]
    parentid_lookup  = dict(parentid_lookup)
    print "building tree with %s Nodes" % len(alltuples)
    c_rangetree.create_tree(tuple(alltuples))


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
    #new way to calculate the precursors
    transitions = c_getnonuis.calculate_transitions_ch(
        ((q1, pep['mod_sequence'], p_id),), [1], q3_low, q3_high)
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
        values="q1, modified_sequence, peptide_key, q1_charge"
        #we compare the parent ion against 4 different parent ions
        #thus we need to take the PEPTIDE key here
        query2 = """
        select %(values)s
        from %(pep)s
        where ssrcalc > %(ssrcalc)s - %(ssr_window)s 
            and ssrcalc < %(ssrcalc)s + %(ssr_window)s
        and q1 > %(q1_low)s and q1 < %(q1_high)s
        and %(pep)s.peptide_key != %(peptide_key)d
        %(query_add)s
        """ % { 'q1_high' : q1_high, 'q1_low'  : q1_low,
               'ssrcalc' : pep['ssrcalc'], 
                'peptide_key' : pep['peptide_key'],
               'query_add' : par.query2_add, 'ssr_window' : par.ssrcalc_window,
               'pep' : par.peptide_table,
               'values' : values}
        cursor.execute( query2 )
        precursors=  cursor.fetchall()
    elif not use_db:
        #correct rounding errors, s.t. we get the same results as before!
        ssrcalc_low = ssrcalc - par.ssrcalc_window + 0.001
        ssrcalc_high = ssrcalc + par.ssrcalc_window - 0.001
        precursor_ids = tuple(c_rangetree.query_tree( q1_low, ssrcalc_low, 
                                                     q1_high,  ssrcalc_high )  )
        precursors = tuple([parentid_lookup[myid[0]] for myid in precursor_ids
                            #dont select myself 
                           if parentid_lookup[myid[0]][2]  != pep['peptide_key']])
    #
    #now calculate coll per peptide the new way
    collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide( 
        transitions, precursors, q3_low, q3_high, par.q3_window, par.ppm)
    ## 
    ## Lets count the number of peptides that interfere
    if count_avg_transitions:
        tr_arr = [0 for i in range(nr_transitions)]
        for v in collisions_per_peptide.values():
            for vv in v:
                tr_arr[vv] += 1
        allintertr.extend( tr_arr )
    ##
    ## Get the nonuis list from the collperpep
    non_uis_list = [{} for i in range(MAX_UIS+1)]
    for order in range(1,MAX_UIS+1):
        non_uis_list[order] = c_getnonuis.get_non_uis(
            collisions_per_peptide, order)
    #
    for order in range(1,min(MAX_UIS+1, nr_transitions+1)): 
        prepare.append( (len(non_uis_list[order]), collider.choose(nr_transitions, 
            order), p_id , order, exp_key) )
    progressm.update(1)


print "\n"
print "found %s transitions" % len(allintertr)#
print "found max of %s interferences" % max(allintertr)
print "found average of %s interferences" % ( sum(allintertr) * 1.0 / len(allintertr) )

cursor.executemany('insert into %s' % restable + ' (non_useable_UIS, total_UIS, \
                  parent_key, uisorder, exp_key) values (%s,%s,%s,%s,%s)' , prepare)


