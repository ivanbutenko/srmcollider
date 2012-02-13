#!/usr/bin/python
# -*- coding: utf-8  -*-
# vim:set fdm=marker:
import MySQLdb
import time
import sys 
sys.path.append( '/home/hroest/srm_clashes/code' )
sys.path.append( '/home/hroest/lib/' )
db = MySQLdb.connect(read_default_file="~/.my.cnf")
cursor = db.cursor()
import collider
import hlib
import copy
import progress

from optparse import OptionParser, OptionGroup
usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)

group = OptionGroup(parser, "Spectral library Options", "")
group.add_option("--background", dest="background", default='',
                  help="A different background " )
parser.add_option_group(group)
par = collider.SRM_parameters()
par.parse_cmdl_args(parser)
options, args = parser.parse_args(sys.argv[1:])
#local arguments
background = options.background
par.__dict__.update( options.__dict__ )
par.parse_options(options)
if par.max_uis ==0: 
    print "Please change --max_uis option, 0 does not make sense here"
    sys.exit()


use_randomized_transitions = True
use_above_precursor = False



"""
background='yeast_pepatlas' 
par.ssrcalc_window=4/2.0 
par.max_uis=10 
par.peptide_table='hroest.srmPeptides_yeast' 
par.transition_table='hroest.srmTransitions_yeast'
"""

par.dontdo2p2f = False
par.do_1vs = False
par.eval()
print par.experiment_type
print par.get_common_filename()
par.query1_add += ' order by Intensity DESC' #only for the toptrans workflow

mycollider = collider.SRMcollider()
cmadd = par.peptide_table.split('.')[1]
#we need a different background for the peptide atlas
#because there we have differnt peptide_keys and cannot relate
#those with the MRMlink table
#only for peptide atlas since those peptide keys map to experiment 3452 and not
#to 3131 as we would need to then find the transition in the MRM table
#we only use the parbg for the calculations of the transitions. We still want 
if background != '':
    parbg = copy.deepcopy(par)
    parbg.transition_table = 'hroest.srmTransitions_' + background
    parbg.peptide_table = 'hroest.srmPeptides_' + background
    background = True
    cmadd = parbg.peptide_table.split('.')[1]
    bgpar = parbg

else: bgpar = par


self = mycollider
self.pepids = self._get_unique_pepids_toptransitions(par, cursor)
MAX_UIS = par.max_uis
common_filename = par.get_common_filename()
self.min_transitions = []
start = time.time()
progressm = progress.ProgressMeter(total=len(self.pepids), unit='peptides')
for i, pep in enumerate(self.pepids):
    p_id = pep['parent_id']
    q1 = pep['q1']
    if not use_randomized_transitions and not use_above_precursor: 
        transitions = self._get_all_transitions_toptransitions(par, pep, cursor)
    elif use_above_precursor: 
        import c_getnonuis
        transitions = c_getnonuis.calculate_transitions_ch(
            ((pep['q1'], pep['mod_sequence'], p_id),), [1], 0, 20000)
        #get all singly charged y ions above the precursor
        y_series = [(t[0], 0) for t in transitions[:len(transitions)/2] if t[0] > q1]
        y_series.reverse()
        y_series = y_series[:5]
        transitions = y_series
    elif use_randomized_transitions: 
        import c_getnonuis
        transitions = c_getnonuis.calculate_transitions_ch(
            ((pep['q1'], pep['mod_sequence'], p_id),), [1], 0, 20000)
        #get random fragment ions
        from random import shuffle
        shuffle(transitions)
        transitions = [ (t[0], 0) for t in transitions[:10] ]

    nr_transitions = len( transitions )
    if nr_transitions == 0: continue #no transitions in this window
    try:
        import c_integrated
        if background: precursors = self._get_all_precursors(bgpar, pep, cursor, bysequence=True)
        else: precursors = self._get_all_precursors(bgpar, pep, cursor)
        min_needed = c_integrated.getMinNeededTransitions(tuple(transitions), tuple(precursors), 
            par.max_uis, par.q3_window, par.ppm)
    except ImportError:
        if background: collisions = self._get_all_collisions(bgpar, pep, cursor, bysequence=True, transitions=transitions)
        else: collisions = self._get_all_collisions(bgpar, pep, cursor, transitions=transitions)
        min_needed = self._getMinNeededTransitions(par, transitions, collisions)

    self.min_transitions.append( [p_id, min_needed] )
    end = time.time()
    if not par.quiet: progressm.update(1)
self.total_time = end - start










## #calculate for precursor, peptide and protein the min nr transitions necessary
## #to measure them without ambiguity
## #create a new hroest.experiment
## common_filename = par.get_common_filename()
## query = """
## insert into hroest.experiment  (name, short_description,
## description, comment1, comment2, comment3, super_experiment_key, ddb_experiment_key)
## VALUES (
##     'paola 1Da/1Da', '%s', '%s', '%s', '%s', '%s', 3, 0
## )
## """ %( common_filename + '_' + cmadd, 
##       par.experiment_type, par.peptide_table, par.transition_table, 
##      'using new MRMAtlas_qtrap_final_no_pyroGlu')
## cursor.execute(query)
## myid = db.insert_id()

cursor.execute("insert into hroest.result_randomexp (dummy) values (0)")
myid = db.insert_id()


prepare = [ [p[0], p[1], myid]  for p in mycollider.min_transitions]
#save our result ("the min nr transitions per precursor") linked with 
#our new experiment key
cursor.executemany(
""" insert into hroest.result_randomtransitions (parent_key, min_transitions, exp_key) 
    VALUES (%s,%s,%s) """, prepare
)


