#!/usr/bin/python

"""
Store all parent ions and their transitions (srmPeptides and srmTransitions)
"""


import sys
from optparse import OptionParser, OptionGroup
usage = "usage: %prog spectrallibrary backgroundorganism [options]\n" 
parser = OptionParser(usage=usage)

group = OptionGroup(parser, "Spectral library Options", "") 
group.add_option("--exp_key", dest="exp_key", default='',
                  help="Experiment Key(s)"  , metavar='(3475, 3474)' ) 
group.add_option("--peptide_table", dest="peptide_table", default='hroest.srmPeptides_test',
                  help="MySQL table containing the peptides" )
group.add_option("--transition_table", dest="transition_table", default='hroest.srmTransitions_test',
                  help="MySQL table containing the transitions" )
group.add_option("--dotransitions", dest="dotransitions", default=False, action="store_true",
                  help="Also fill the transitions table (not always necessary)")
parser.add_option_group(group)

#parameters evaluation
#parameters = collider.SRM_parameters()
#parameters.parse_cmdl_args(parser)
options, args = parser.parse_args(sys.argv[1:])
peptide_table = options.peptide_table
transition_table = options.transition_table
exp_key = options.exp_key
dotransitions = options.dotransitions

##exp_key = 3061  #human (800 - 5000 Da)
##exp_key = 3130  #human (all)
###exp_key = 3120  #yeast
##exp_key = 3131  #yeast (all)
###exp_key = 3352  #yeast, 1200 peptides
###exp_key = 3445  #mouse (all)
##exp_key = '(3475, 3474)' #all TB

import MySQLdb
import sys 
sys.path.append( '/home/hroest/lib/' )
sys.path.append( '/home/hroest/srm_clashes/code' ) #Collider
sys.path.append( '/home/hroest/msa/code/tppGhost' ) #DDB
sys.path.append( '/home/hroest/lib/hlib' ) #utils
import time
from hlib import utils
db = MySQLdb.connect(read_default_file="~/.my.cnf")
c = db.cursor()
c2 = db.cursor()
import silver
import DDB 
import csv 
R = silver.Residues.Residues('mono')
import numpy
import progress
import collider

all_peptide_query  = """
select distinct peptide.sequence, molecular_weight, ssrcalc,
genome_occurence, 
peptide.id as peptide_key
from peptide 
inner join peptideOrganism on peptide.id = peptideOrganism.peptide_key
inner join compep.ssrcalc_prediction on ssrcalc_prediction.sequence =
peptide.sequence
where experiment_key in %s
and length( peptide.sequence ) > 1
""" % exp_key

#human go from parent_id = 1 to parent_id = 1177958
#R.recalculate_monisotopic_data_for_N15()

###################################
# A) store the transitions
###################################
#read in the data from DB and put into bins
# 75 peptides/sec
insert_db = True
modify_cysteins = False
read_from_db = False
c.execute( 'use ddb;' )
t = utils.db_table( c2 )
t.read( all_peptide_query )
print "executed the query"
rows = t.c.fetchall()

c.execute( " drop table if exists %(table)s;" % { 'table' : peptide_table} )
c.execute( """
create table %(table)s(
    parent_id INT PRIMARY KEY AUTO_INCREMENT,
    peptide_key INT,
    modified_sequence VARCHAR(255),
    q1_charge TINYINT,
    q1 DOUBLE,
    ssrcalc DOUBLE,
    isotope_nr TINYINT,
    transition_group INT
);
""" % { 'table' : peptide_table} )
c.execute("alter table %(table)s add index(peptide_key);" % { 'table' : peptide_table} )
c.execute("alter table %(table)s add index(q1);" % { 'table' : peptide_table} )
c.execute("alter table %(table)s add index(ssrcalc);" % { 'table' : peptide_table} )
c.execute("alter table %(table)s add index(transition_group);" % { 'table' : peptide_table} )
#c.execute('truncate table ' + peptide_table)
if dotransitions: c.execute('truncate table ' + transition_table)
mass_bins = [ []  for i in range(0, 10000) ]
rt_bins = [ []  for i in range(-100, 500) ]
tmp_c  = db.cursor()
progressm = progress.ProgressMeter(total=len(rows), unit='peptides')
transition_group = 0
for row in rows:
    progressm.update(1)
    transition_group += 1 
    #if transition_group >= 1001: break #####FOR TESTING ONLY
    for mycharge in [2,3]:  #precursor charge 2 and 3
        peptide = collider.get_peptide_from_table(t, row)
        peptide.charge = mycharge
        if modify_cysteins: peptide.modify_cysteins()
        peptide.create_fragmentation_pattern(R)
        if insert_db:
            #insert peptide into db
            collider.insert_peptide_in_db(peptide, db, peptide_table,
                                          transition_group=transition_group)
    #we want to insert the fragments only once per peptide
    if insert_db and dotransitions:
        #insert fragment charge 1 and 2 into database
        peptide.mass_H = R.mass_H
        collider.fast_insert_in_db(peptide, db, 1, transition_table,
                                   transition_group=transition_group)
        collider.fast_insert_in_db(peptide, db, 2, transition_table,
                                   transition_group=transition_group)


#A2 create the additional parent ions for the isotope patterns
PATTERNS_UP_TO = 3 #up to how many isotopes should be considered
cursor = c
vals = "peptide_key, q1_charge, q1, modified_sequence, ssrcalc, isotope_nr, transition_group"
query = "SELECT parent_id, %s FROM %s" % (vals, peptide_table)
cursor.execute( query )
allpeptides =  cursor.fetchall()
prepared = []
for p in allpeptides:
    q1_charge = p[2]
    for i in range(1,1+PATTERNS_UP_TO):
        prepared.append( [p[1], p[2], p[3] + (R.mass_diffC13 * i* 1.0) / q1_charge,
                              p[4], p[5], i, p[7]] )

q = "INSERT INTO %s (%s)" % (peptide_table, vals)  + \
     " VALUES (" + "%s," *6 + "%s)"
c.executemany( q, prepared)

