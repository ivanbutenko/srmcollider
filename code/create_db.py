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
 *
 * Set up the database tables containing the peptide precursors.  If you start
 * out from a FASTA file, you should use the trypsinize.py script in the script
 * folder. Afterwards, use SSRCalc to calculate hydrophobicity and you can use
 * the input from the SSRCalc program directly with this script.
 * 
 * The input in a TSV file should be the following

MTGSVDRPDQNR	12	      8.56	      8.6
MTGSVDRPDQNRGER	15	      7.83	      7.8

 * where the first and the third line will be used (sequence and SSRCalc value).
 *
"""

import sys, csv
sys.path.append('external/')
import MySQLdb
import Residues
import DDB 
import collider

print "Script is deactivated, please edit if you want to use it."
# Since this drops tables, we done want to run it by accident.
sys.exit() #remove this if you want to use this script. 

from optparse import OptionParser, OptionGroup
usage = "usage: %prog [options]\n" 
parser = OptionParser(usage=usage)

group = OptionGroup(parser, "Create db tables Options", "") 
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
group.add_option("--mysql_config", dest="mysql_config", default='~/.my.cnf',
                  help="Location of mysql config (.my.cnf) file" )
group.add_option("--mass_cutoff", dest="mass_cutoff", default='5000',
                  help="M/Z cutoff above which precursors will not be included in the database (default 5000)" )
group.add_option("--max_nr_modifications", dest="max_nr_modifications", default='3',
                  help="Maximal number of modifications per peptide (default 3)" )
group.add_option("--oxidize_methionines", dest="oxidize_methionines", default=False, action="store_true",
                  help="Oxidize Methionines")
group.add_option("--deamidate_asparagine", dest="deamidate_asparagine", default=False, action="store_true",
                  help="Deamindate asparagines")
group.add_option("--doN15", dest="doN15", default=False, action="store_true",
                  help="Modify all proteins with heavy Nitrogen (N15).")
parser.add_option_group(group)

options, args = parser.parse_args(sys.argv[1:])
peptide_table = options.peptide_table
transition_table = options.transition_table
dotransitions = options.dotransitions
tsv_file = options.tsv_file
sqlite_database = options.sqlite_database
mass_cutoff = int(options.mass_cutoff)
# Modifications 
modify_cysteins = True
oxidize_methionines = options.oxidize_methionines 
deamidate_asparagine = options.deamidate_asparagine 
doN15 = options.doN15 
max_nr_modifications = int(options.max_nr_modifications)
# Input file and db
if tsv_file != '': use_tsv = True
else: use_tsv = False
if sqlite_database != '': use_sqlite = True
else: use_sqlite = False

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
    isotope_modification = Residues.NOISOTOPEMODIFICATION
    if doN15:
        isotope_modification = Residues.N15_ISOTOPEMODIFICATION
    #insert peptide into db
    vals = "peptide_key, q1_charge, q1, ssrcalc, modified_sequence, modifications, missed_cleavages, isotopically_modified, transition_group"
    q = "insert into %s (%s) VALUES (%s,%s,%s,%s,'%s', %s, %s, %s, %s)" % (
        peptide_table,
        vals, 
        self.id, self.charge, 
        self.charged_mass, self.ssr_calc, 
        self.get_modified_sequence(),
        self.modifications, 
        self.missed_cleavages(),
        isotope_modification,
        transition_group
    )
    c.execute(q)
    self.parent_id = db.insert_id()

def get_all_modifications(this_sequence, to_modify, replace_with, max_nr_modifications):
    import re
    positions = [m.start() for m in re.finditer(to_modify, this_sequence)]
    for i in range(1,min(max_nr_modifications,1+len(positions))):
        to_replace = list(collider.combinations( positions, i))
        for new_peptide in to_replace:
            it = 0
            curr_seq = ''
            for pos in new_peptide:
                curr_seq += this_sequence[it:pos] + replace_with
                it = pos+1
            curr_seq += this_sequence[it:]
            pep = DDB.Peptide()
            pep.set_sequence(curr_seq)
            pep.modifications = i
            yield pep

def get_all_modified_peptides(peptide, oxidize_methionines, deamidate_asparagine, max_nr_modifications):
    modified = []
    if oxidize_methionines: 
        modified = list(get_all_modifications(peptide.get_modified_sequence(), 'M', 'M[147]', max_nr_modifications))
    ##for pep in modified:
    ##    print pep.get_modified_sequence()
    #continue
    if deamidate_asparagine:
        toappend = []
        for pep in modified:
            tmp_modified = get_all_modifications(pep.get_modified_sequence(), 'N', 'N[115]', max_nr_modifications)
            for pep_new in tmp_modified:
                pep_new.modifications += pep.modifications
                toappend.append(pep_new)
        modified.extend(toappend)
        modified.extend( list(get_all_modifications(peptide.get_modified_sequence(), 'N', 'N[115]', max_nr_modifications) ) )
    return modified



residues = Residues.Residues('mono')
if doN15:
    residues.recalculate_monisotopic_data_for_N15()

###################################
# A) store the transitions
###################################

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
        modifications TINYINT UNSIGNED,
        missed_cleavages TINYINT UNSIGNED,
        isotopically_modified TINYINT UNSIGNED,
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
        modifications TINYINT UNSIGNED,
        missed_cleavages TINYINT UNSIGNED,
        isotopically_modified TINYINT UNSIGNED,
        transition_group INT
    );
    """ % { 'table' : peptide_table} )
    c.execute("alter table %(table)s add index(peptide_key);" % { 'table' : peptide_table} )
    c.execute("alter table %(table)s add index(q1);" % { 'table' : peptide_table} )
    c.execute("alter table %(table)s add index(ssrcalc);" % { 'table' : peptide_table} )
    c.execute("alter table %(table)s add index(transition_group);" % { 'table' : peptide_table} )
    #c.execute('truncate table ' + peptide_table)
    if dotransitions: c.execute('truncate table ' + transition_table)

transition_group = 0
# which modifications:
## methionine oxidation
## de-amidation
## phospho? no
## N-terminal acetylation ?
## methylation
charges = [2,3]  #precursor charge 2 and 3

f = open(tsv_file)
reader = csv.reader(f, delimiter='\t')
for id, line in enumerate(reader):
    if len(line[0]) < 2: continue
    if line[0].startswith("<tr class="): continue
    peptide = DDB.Peptide()
    peptide.set_sequence( line[0] )
    peptide.ssr_calc = line[2] 
    peptide.id = id

    transition_group += 1 
    if modify_cysteins: peptide.modify_cysteins()
    #print transition_group, peptide.sequence
    for mycharge in charges:
        peptide.charge = mycharge
        # calculate charged mass and insert peptide into db
        peptide.create_fragmentation_pattern(residues)
        if peptide.charged_mass > mass_cutoff: continue
        peptide.modifications = 0
        insert_peptide_in_db(peptide, db, peptide_table,
                                      transition_group=transition_group)
    #
    # some heuristics to make the whole thing faster, if already the charged
    # mass if twice as high as the cutoff, even with modifications we will
    # never get into the allowed mass range.
    if peptide.charged_mass / 2.0 > mass_cutoff: continue
    #
    modified = get_all_modified_peptides(peptide, oxidize_methionines, deamidate_asparagine, max_nr_modifications) 
    for mod_peptide in modified:
        if mod_peptide.modifications > max_nr_modifications: continue
        transition_group += 1 
        for mycharge in charges: 
            mod_peptide.ssr_calc = peptide.ssr_calc
            mod_peptide.id = peptide.id
            mod_peptide.charge = mycharge
            # calculate charged mass and insert peptide into db
            mod_peptide.create_fragmentation_pattern(residues)
            if mod_peptide.charged_mass > mass_cutoff: continue
            insert_peptide_in_db(mod_peptide, db, peptide_table,
                                          transition_group=transition_group)
    del modified
    del peptide

db.commit()

