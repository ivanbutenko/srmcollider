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
 *
 * Set up the database tables containing the peptide precursors.  If you start
 * out from a FASTA file, you should use the trypsinize.py script in the script
 * folder. Afterwards, use SSRCalc to calculate hydrophobicity and you can use
 * the input from the SSRCalc program directly with this script.
 *
"""

import sys
sys.path.append('external/')
import MySQLdb
import Residues
import DDB 
import progress
import collider

print "Script is deactivated, please edit if you want to use it."
# Since this drops tables, we done want to run it by accident.
sys.exit() #remove this if you want to use this script. 

from optparse import OptionParser, OptionGroup
usage = "usage: %prog [options]\n" 
parser = OptionParser(usage=usage)

group = OptionGroup(parser, "Create db tables Options", "") 
group.add_option("--exp_key", dest="exp_key", default='',
                  help="Experiment Key(s)"  , metavar='(3475, 3474)' ) 
group.add_option("--peptide_table", dest="peptide_table", default='hroest.srmPeptides_test',
                  help="MySQL table containing the peptides" )
group.add_option("--transition_table", dest="transition_table", default='hroest.srmTransitions_test',
                  help="MySQL table containing the transitions" )
group.add_option("--dotransitions", dest="dotransitions", default=False, action="store_true",
                  help="Also fill the transitions table (not always necessary)")
group.add_option("--tsv_file", dest="tsv_file", default='',
                  help="Take TSV file created by SSRCalc as input" )
group.add_option("--sqlite_database", dest="sqlite_database", default='',
                  help="Use specified sqlite database instead of MySQL database" )
group.add_option("--nr_isotopes", dest="nr_isotopes", default='3',
                  help="Number of isotopes of the precursor to consider (default 3)" )
group.add_option("--mysql_config", dest="mysql_config", default='~/.my.cnf',
                  help="Location of mysql config (.my.cnf) file" )
parser.add_option_group(group)

options, args = parser.parse_args(sys.argv[1:])
peptide_table = options.peptide_table
transition_table = options.transition_table
exp_key = options.exp_key
dotransitions = options.dotransitions
tsv_file = options.tsv_file
sqlite_database = options.sqlite_database
if tsv_file != '': use_tsv = True
else: use_tsv = False
if sqlite_database != '': use_sqlite = True
else: use_sqlite = False


#up to how many isotopes should be considered
PATTERNS_UP_TO = int(options.nr_isotopes) 

if use_sqlite:
    import sqlite
    conn = sqlite.connect(sqlite_database)
    db = conn
    c = conn.cursor()
    c2 = conn.cursor()

else:
    db = MySQLdb.connect(read_default_file=options.mysql_config)
    c = db.cursor()
    c2 = db.cursor()



## fxn
def get_peptide_from_table(t, row):
    peptide = DDB.Peptide()
    peptide.set_sequence( t.row(row, 'sequence')  )
    peptide.genome_occurence = t.row( row, 'genome_occurence' )
    peptide.ssr_calc         = t.row( row, 'ssrcalc' )
    #peptide.gene_id          = t.row( row, 'gene_id' )
    peptide.mw               = t.row( row, 'molecular_weight' )
    peptide.id               = t.row( row, 'peptide_key' )
    peptide.pairs            = []
    peptide.non_unique = {}
    return peptide

def insert_peptide_in_db(self, db, peptide_table, transition_group):
    c = db.cursor()
    #insert peptide into db
    vals = "peptide_key, q1_charge, q1, ssrcalc, modified_sequence, isotope_nr, transition_group"
    q = "insert into %s (%s) VALUES (%s,%s,%s,%s,'%s', %s, %s)" % (
        peptide_table,
        vals, 
        self.id, self.charge, 
        self.charged_mass, self.ssr_calc, 
        self.get_modified_sequence(),
        0, #we only have the 0th isotope (0 C13 atoms)
        transition_group
    )
    c.execute(q)
    self.parent_id = db.insert_id()


##exp_key = 3061  #human (800 - 5000 Da)
##exp_key = 3130  #human (all)
###exp_key = 3120  #yeast
##exp_key = 3131  #yeast (all)
###exp_key = 3352  #yeast, 1200 peptides
###exp_key = 3445  #mouse (all)
##exp_key = '(3475, 3474)' #all TB


residues = Residues.Residues('mono')

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

insert_db = True
modify_cysteins = False

# how to get the peptides 
import csv
rows = []
#f = open('ssrcalc.out')
f = open(tsv_file)
reader = csv.reader(f, delimiter='\t')
for id, line in enumerate(reader):
    if len(line[0]) < 2: continue
    peptide = DDB.Peptide()
    peptide.set_sequence( line[0] )
    peptide.ssr_calc = line[2] 
    peptide.id = id
    rows.append(peptide)



# where to store the peptides (in MySQL or sqlite)
if use_sqlite:
    try: c.execute( " drop table %(table)s;" % { 'table' : peptide_table} )
    except Exception: pass
    c.execute(
    """
    create table %(table)s(
        parent_id INT PRIMARY KEY,
        peptide_key INT,
        modified_sequence VARCHAR(255),
        q1_charge TINYINT,
        q1 DOUBLE,
        ssrcalc DOUBLE,
        isotope_nr TINYINT,
        transition_group INT
    );
    create index %(table)spepkey   on %(table)s (peptide_key);
    create index %(table)sq1       on %(table)s (q1);
    create index %(table)sssrcalc  on %(table)s (ssrcalc);
    create index %(table)strgroup  on %(table)s (transition_group);
        """ % {'table' : peptide_table}
    )

else:
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
    for mycharge in [2,3]:  #precursor charge 2 and 3
        if not use_tsv: peptide = get_peptide_from_table(t, row)
        else: peptide = row
        peptide.charge = mycharge
        if modify_cysteins: peptide.modify_cysteins()
        peptide.create_fragmentation_pattern(residues)
        if insert_db:
            #insert peptide into db
            insert_peptide_in_db(peptide, db, peptide_table,
                                          transition_group=transition_group)


if use_sqlite: db.commit()

