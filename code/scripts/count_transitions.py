#!/usr/bin/python
# -*- coding: utf-8  -*-
# vim:set fdm=marker:
""" This script will go through all peptides and count the number of
transitions and assays that are present for each peptide.
"""
import MySQLdb,time, sys, c_rangetree, c_getnonuis, collider, progress
sys.path.append( '/home/hroest/lib/hlib/' )
sys.path.append( '/home/hroest/projects' )
sys.path.append( '/home/hroest/projects/hlib' )
db = MySQLdb.connect(read_default_file="~/.my.cnf")
cursor = db.cursor()

min_q1 = 400
max_q1 = 800

peptide_table = 'hroest.srmPeptides_human'
qadd = ''
#qadd = 'and peptide_key = 10673050 '
qadd = 'and isotope_nr in (0) and q1_charge = 2'

#Run 
###########################################################################
#Parse options
db = MySQLdb.connect(read_default_file="~/.my.cnf")
cursor = db.cursor()

start = time.time()
q =  """
select modified_sequence, peptide_key, parent_id, q1_charge, q1, ssrcalc, isotope_nr
from %(peptide_table)s where q1 between %(lowq1)s and %(highq1)s
%(qadd)s
""" % {'peptide_table' : peptide_table, 
              'lowq1'  : min_q1 ,
              'highq1' : max_q1,
              'qadd'   : qadd
} 
print q
cursor.execute(q)
print "Obtained all peptides"

alltuples =  list(cursor.fetchall() )

progressm = progress.ProgressMeter(total=len(alltuples), unit='peptides')
prepare  = []
total_transitions = 0
total_assays = 0
for kk, pep in enumerate(alltuples):
    q3_low, q3_high = [400, 1200]
    q3_low, q3_high = [400, 1400]
    q3charges = [1,2]
    q3charges = [1]
    #q3_low, q3_high = [0, 12000]
    #
    #new way to calculate the precursors
    transitions = c_getnonuis.calculate_transitions_ch(
        ((-2, pep[0], -1),), q3charges, q3_low, q3_high)
    #print pep, len(pep[0]), len(transitions)
    total_transitions += len(transitions)
    total_assays += sum([ collider.choose(len(transitions), i) for i in range(1,6) 
                        if len(transitions) >= i ] )
    progressm.update(1)


print "Total number of precursors is ", len(alltuples)
print "Total number of transitions is ", total_transitions
print "Total number of assays is ", total_assays


