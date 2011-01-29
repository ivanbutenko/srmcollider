#!/usr/bin/python

"""
Store all parent ions and their transitions (srmPeptides and srmTransitions)
"""

#safety is on, remove to acutally run the script
from sys import exit
print "I didnt do anything, exiting now"
exit()

import MySQLdb
import sys 
sys.path.append( '/home/hroest/lib/' )
sys.path.append( '/home/hroest/srm_clashes/code' ) #Collider
sys.path.append( '/home/hroest/msa/code/tppGhost' ) #DDB
sys.path.append( '/home/hroest/lib/hlib' ) #utils
import time
from hlib import utils
#import pipeline 
#db = MySQLdb.connect(read_default_file="~/hroest/.my.cnf")
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

exp_key = 3061  #human (800 - 5000 Da)
exp_key = 3130  #human (all)
#exp_key = 3120  #yeast
#exp_key = 3131  #yeast (all)
#exp_key = 3352  #yeast, 1200 peptides
#exp_key = 3445  #mouse (all)
all_peptide_query  = """
select distinct peptide.sequence, molecular_weight, ssrcalc,
1 as genome_occurence, 
peptide.id as peptide_key
from peptide 
#inner join peptideOrganism on peptide.id = peptideOrganism.peptide_key
inner join compep.ssrcalc_prediction on ssrcalc_prediction.sequence =
peptide.sequence
where experiment_key = %s
and length( peptide.sequence ) > 1
order by molecular_weight
""" % exp_key

#human go from parent_id = 1 to parent_id = 1177958
#R.monoisotopic_data
#R.recalculate_monisotopic_data_for_N15()

for i in range(400, 1800, 100):
    db = MySQLdb.connect(read_default_file="~/.my.cnf")
    c = db.cursor()
    mytable = 'hroest.srmTransitions_perf_test4_' + str(i)
    c.execute("drop table if exists %(table)s " % {  'table' : mytable } )
    #except Exception: pass
    create_q = """
    create table %(table)s (
        srm_id INT PRIMARY KEY,
        group_id INT,
        q3_charge TINYINT ,
        q3 DOUBLE,
        type VARCHAR(8),
        fragment_number TINYINT
    ) ;
    alter table %(table)s add index(srm_id);
    alter table %(table)s add index(group_id);
    alter table %(table)s add index(q3);
    """ % {  'table' : mytable  }
    c.execute(create_q)

for i in range(400, 1800, 100):
    db = MySQLdb.connect(read_default_file="~/.my.cnf")
    c = db.cursor()
    mytable = 'hroest.srmPeptides_perf_test4_' + str(i)
    c.execute("drop table if exists %(table)s " % {  'table' : mytable } )
    #except Exception: pass
    create_q = """
    create table %(table)s(
        parent_id INT PRIMARY KEY,
        peptide_key INT,
        modified_sequence VARCHAR(255),
        q1_charge TINYINT,
        q1 DOUBLE,
        ssrcalc DOUBLE,
        isotope_nr TINYINT,
        transition_group INT
    ) ;
    alter table %(table)s add index(peptide_key);
    alter table %(table)s add index(q1);
    alter table %(table)s add index(ssrcalc);
    alter table %(table)s add index(transition_group);
        """ % {  'table' : mytable  }
    c.execute(create_q)


###################################
# A) store the transitions
###################################
#read in the data from DB and put into bins
# 75 peptides/sec
db = MySQLdb.connect(read_default_file="~/.my.cnf")
c = db.cursor()
c2 = db.cursor()
insert_db = True
modify_cysteins = True
read_from_db = False
c.execute( 'use ddb;' )
t = utils.db_table( c2 )
t.read( all_peptide_query )
print "executed the query"
rows = t.c.fetchall()

db = MySQLdb.connect(read_default_file="~/.my.cnf")
c = db.cursor()
#1000 entries / 9 s with fast (executemany), 10.5 with index
#1000 entries / 31 s with normal (execute), 34 with index
transition_table = 'hroest.srmTransitions_perf_test4_400'
base_transition_table = 'hroest.srmTransitions_perf_test4_'
peptide_table = 'hroest.srmPeptides_perf_test4_400'
base_peptide_table = 'hroest.srmPeptides_perf_test4_'
c.execute('truncate table ' + peptide_table)
c.execute('truncate table ' + transition_table)
mass_bins = [ []  for i in range(0, 10000) ]
rt_bins = [ []  for i in range(-100, 500) ]
tmp_c  = db.cursor()
progressm = progress.ProgressMeter(total=len(rows), unit='peptides')
transition_group = 0
parent_id = 0
srm_id  = 0

mapping = {}
overlap = 20 
for i in range(400, 1800, 100):
    i
    mapping[ (i-overlap, i+100+overlap)  ] = i

def determine_tables(weight):
    results = []
    for k,v in mapping.iteritems():
        if weight > k[0] and weight < k[1]: results.append( v)
    return results

def fast_insert_in_db(self, db, fragment_charge, transition_table, transition_group):
    c = db.cursor()
    srm_id = peptide.srm_id
    ch = fragment_charge
    vals = "srm_id, type, fragment_number, group_id, q3_charge, q3 "
    q = "insert into %s (%s)" % (transition_table, vals)  + " VALUES (%s, %s,%s,%s,%s,%s)" 
    tr = len(self.y_series)
    charged_y =  [ ( pred + (ch -1)*self.mass_H)/ch for pred in self.y_series ]
    charged_b =  [ ( pred + (ch -1)*self.mass_H)/ch for pred in self.b_series ]
    many = [ [srm_id + i, 'y', i+1, transition_group, ch, q3] for i, q3 in enumerate(reversed(charged_y))] 
    srm_id += len(charged_y)
    manyb = [ [srm_id + i, 'b', i+1, transition_group, ch, q3] for i, q3 in enumerate(charged_b) ]
    many.extend( manyb )
    c.executemany( q, many)
    srm_id += len(charged_b)
    peptide.srm_id = srm_id
    #
    #here we could recover the inserted ids but since they map to several 
    #peptides, this is not helpful
    return
    first_id = db.insert_id()
    self.fragment_ids[ ch ] = {
    'y' : [i for i in range(first_id, first_id + len(charged_y )) ] , 
    'b' : [i for i in range(first_id + len(charged_y ),  first_id + 2*len(charged_y )) ]
    }

def insert_peptide_in_db(self, db, peptide_table, transition_group, parent_id):
            c = db.cursor()
            #insert peptide into db
            vals = "parent_id, peptide_key, q1_charge, q1, ssrcalc, modified_sequence, isotope_nr, transition_group"
            q = "insert into %s (%s) VALUES (%s, %s,%s,%s,%s,'%s', %s, %s)" % (
                peptide_table,
                vals, 
                parent_id,
                self.id, self.charge, 
                self.charged_mass, self.ssr_calc, 
                self.get_modified_sequence(),
                0, #we only have the 0th isotope (0 C13 atoms)
                transition_group
            )
            c.execute(q)
            self.parent_id = parent_id

for row in rows:
    #print transition_table
    progressm.update(1)
    transition_group += 1 
    for mycharge in [2,3]:  #precursor charge 2 and 3
        peptide = collider.get_peptide_from_table(t, row)
        peptide.charge = mycharge
        if modify_cysteins: peptide.modify_cysteins()
        peptide.create_fragmentation_pattern(R)
        #
        weight = peptide.charged_mass
        tables_to_use = determine_tables(weight)
        #
        if insert_db:
            #insert peptide into db
            for table_to_use in tables_to_use:
                transition_table = base_transition_table + str(table_to_use)
                peptide_table = base_peptide_table + str(table_to_use)
                insert_peptide_in_db(peptide, db, peptide_table,
                              transition_group=transition_group, parent_id=parent_id)
                parent_id += 1
        #we want to insert the fragments only once per peptide
        if insert_db:
            #insert fragment charge 1 and 2 into database
            peptide.mass_H = R.mass_H
            peptide.srm_id = srm_id
            for table_to_use in tables_to_use:
                transition_table = base_transition_table + str(table_to_use)
                peptide_table = base_peptide_table + str(table_to_use)
                fast_insert_in_db(peptide, db, 1, transition_table,
                                           transition_group=transition_group)
                fast_insert_in_db(peptide, db, 2, transition_table,
                                           transition_group=transition_group)
            srm_id = peptide.srm_id


new_parent_id = parent_id
#with numbering of the primary key oursevles
for i in range(400, 1800, 100):
    db = MySQLdb.connect(read_default_file="~/.my.cnf")
    c = db.cursor()
    peptide_table = base_peptide_table + str(i)
    #except Exception: pass
    #A2 create the additional parent ions for the isotope patterns
    PATTERNS_UP_TO = 3 #up to how many isotopes should be considered
    cursor = c
    vals = "peptide_key, q1_charge, q1, modified_sequence, ssrcalc, isotope_nr, transition_group, parent_id"
    query = "SELECT parent_id, %s FROM %s" % (vals, peptide_table)
    cursor.execute( query )
    allpeptides =  cursor.fetchall()
    prepared = []
    for p in allpeptides:
        q1_charge = p[2]
        for i in range(1,1+PATTERNS_UP_TO):
            new_parent_id += 1
            prepared.append( [p[1], p[2], p[3] + (R.mass_diffC13 * i* 1.0) / q1_charge,
                                  p[4], p[5], i, p[7], new_parent_id] )
    q = "INSERT INTO %s (%s)" % (peptide_table, vals)  + \
         " VALUES (" + "%s," *7 + "%s)"
    c.executemany( q, prepared)

