#!/usr/bin/python

"""
There are two different kinds of setups:
    A) store all parent ions and their transitions (srmPeptides and srmTransitions)
    B) also store all collisions (in a table like srmCollisions400710)
"""

import os
import MySQLdb
import sys 
sys.path.append( '/home/ghost/hroest/code/' )
sys.path.append( '/home/ghost/software_libs/' )
sys.path.append( '/home/hroest/lib/' )
sys.path.append( '/home/hroest/msa/code/tppGhost' )
sys.path.append( '/home/hroest/msa/code' )
#import mzXMLreader
#import pepXMLReader
import time
from utils_h import utils
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

exp_key = 3061  #human (800 - 5000 Da)
#exp_key = 3130  #human (all)
#exp_key = 3120  #yeast
all_peptide_query  = """
select distinct peptide.sequence, molecular_weight, ssrcalc,
genome_occurence, peptide.id as peptide_key
from peptide 
inner join peptideOrganism on peptide.id = peptideOrganism.peptide_key
inner join compep.ssrcalc_prediction on ssrcalc_prediction.sequence =
peptide.sequence
where experiment_key = %s
""" % exp_key

#human go from parent_id = 1 to parent_id = 1177958

###################################
# A) store the transitions
###################################
#read in the data from DB and put into bins
insert_db = True
read_from_db = False
c.execute( 'use ddb;' )
t = utils.db_table( c2 )
t.read( all_peptide_query )
print "executed the query"
rows = t.c.fetchall()

#1000 entries / 9 s with fast (executemany), 10.5 with index
#1000 entries / 31 s with normal (execute), 34 with index
transition_table = 'hroest.srmTransitions_test'
peptide_table = 'hroest.srmPeptides_test'
mass_bins = [ []  for i in range(0, 10000) ]
rt_bins = [ []  for i in range(-100, 500) ]
tmp_c  = db.cursor()
start = time.time()
for i,row in enumerate(rows):
    if i % 10000 == 0: print i
    for mycharge in [2,3]:  #precursor charge 2 and 3
        peptide = get_peptide_from_table(t, row)
        peptide.charge = mycharge
        S = silver.Spectrum.Spectrum(SEQUEST_mode =1 )
        S.ion_charge = peptide.charge
        #S.peptide_mass = S.calculate_peptide_mass( 
        #    peptide.get_modified_sequence('SEQUEST'), R.residues, 
        #    None, S.ion_charge)
        S.construct_from_peptide( 
            peptide.get_modified_sequence('SEQUEST'), R.residues, 
            R.res_pairs)
        S.ass_peptide = peptide
        S.fragment_ids = {}
        if insert_db:
            #insert peptide into db
            insert_peptide_in_db(S, db, peptide_table)
            #insert fragment charge 1 and 2 into database
            fast_insert_in_db( S, db, 1, transition_table)
            fast_insert_in_db( S, db, 2, transition_table)
        elif read_from_db:
            #instead of inserting, we read the database
            read_fragment_ids(S, tmp_c, peptide, peptide_table, transition_table)
        #insert peptide into mass_bins hash
        peptide.spectrum = S
        mz = S.peptide_mass / S.ion_charge
        bin = int( mz / 0.7) ;
        mass_bins[ bin ].append( peptide )
        rt_bins[ int(peptide.ssr_calc) ].append( peptide )
        end = time.time()


def read_fragment_ids(self, cursor, peptide, peptide_table, transition_table):
    #reads srmPeptide and srmTransitions in order to store the srm_ids of the 
    #fragments. Speed ~ 300 / s
    qq = """select 
    q3_charge, type, srm_id
    from %s 
    inner join %s on parent_id = parent_key
    where peptide_key = %s and q1_charge = %s
    order by q3_charge, type
    """ % (peptide_table, transition_table, peptide.id, peptide.charge)
    cursor.execute( qq ) 
    all_transitions = cursor.fetchall()
    self.fragment_ids = { 1 : [ [], [] ], 2 : [ [], [] ]}
    tmp = 0
    for r in all_transitions:
        if r[1] == 'b' : tmp =0
        else: tmp = 1
        self.fragment_ids[ r[0] ][tmp].append( r[2] )

###################################
# B) store the collisions
###################################


#here we run through all pairs
reset_pairs_unique(mass_bins)
range_small = range(300,1500)
range_all = range(0,100000)
frag_range = range(1,2) #range for fragment ion charge
MS1_bins = 0.7 /2
MS2_bins = 1.0 /2
ssrcalc_binsize = 4/2 #that may be around 2 minutes, depending on coloumn
#coding == ssrcalc_binsize : MS1_binsize : MS2_binsize
table_all = 'srmCollisions400710_all'
table = 'srmCollisions400710'
lab_peptides_2 = []
pairs_dic = {}
clash_distr = [0 for i in range( 1000 ) ]
for i in range( 1, len( mass_bins )-1):
    if i % 100 ==0:
        print i, len( pairs_dic )
    #ch1 = 1
    #ch2 = 1
    #i = 1500
    current = []
    current.extend( mass_bins[i] )
    current.extend( mass_bins[i-1] )
    current.extend( mass_bins[i+1] )
    current_spectra = []
    for peptide in current: current_spectra.append( peptide.spectrum )
    #for i, S in enumerate(current_spectra):
    for S in current_spectra:
        peptide = S.ass_peptide
        mz = S.peptide_mass / S.ion_charge
        bin = int( mz / 0.7) 
        #only use those as S which are in the current bin => take each S only once
        if bin != i: continue
        #compare query spectra S against all other spectra that are within 
        #a certain tolerance of q1 mass
        for S2 in current_spectra:
            #we are not interested in spectras from the same sequence
            if S2.ass_peptide.sequence == S.ass_peptide.sequence: continue
            ms1diff = abs( mz - S2.peptide_mass / S2.ion_charge)
            ssrdiff = abs( peptide.ssr_calc - S2.ass_peptide.ssr_calc)
            #we need to do S vs S2 BUT ALSO S2 vs S (later one)
            if ( ms1diff < MS1_bins and ssrdiff < ssrcalc_binsize):
                sh = all_calculate_clashes_in_series_insert_db( S, S2, pairs_dic,
                     MS2_bins, clash_distr, range_all, c, table, frag_range )



###########################################################################
#
# Functions
#

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

def insert_peptide_in_db(self, db, peptide_table):
    c = db.cursor()
    peptide = self.ass_peptide
    #insert peptide into db
    vals = "peptide_key, q1_charge, q1, ssrcalc"
    q = "insert into %s (%s) VALUES (%s,%s,%s,%s)" % (
        peptide_table,
        vals, 
        peptide.id, peptide.charge, 
        get_actual_mass(self), peptide.ssr_calc )
    c.execute(q)
    peptide.parent_id = db.insert_id()

def get_actual_mass(self):
    return self.peptide_mass / self.ion_charge

def insert_in_db(self, db, fragment_charge, transition_table):
    assert False #use fast now
    c = db.cursor()
    peptide = self.ass_peptide
    parent_id = peptide.parent_id
    vals = "type, parent_key, q3_charge, q3"
    ch = fragment_charge
    charged_y_series =  [ ( pred + (ch -1)*self.mass_H)/ch for pred in self.y_series ]
    charged_b_series =  [ ( pred + (ch -1)*self.mass_H)/ch for pred in self.b_series ]
    self.fragment_ids[ ch ] = [ [], [] ]
    for i, q3 in enumerate(charged_b_series):
        type = 'b'
        q = "insert into %s (%s) VALUES ('%s',%s,%s,%s)" % (
            transition_table,
            vals, 
            type, parent_id, fragment_charge, q3 )
        c.execute(q)
        self.fragment_ids[ch][0].append( db.insert_id()  )
    pepLen = len(peptide)
    for i, q3 in enumerate(charged_y_series):
        type = 'y'
        q = "insert into %s (%s) VALUES ('%s',%s,%s,%s)" % (
            transition_table,
            vals, 
            type, parent_id, fragment_charge, q3 )
        c.execute(q)
        self.fragment_ids[ch][1].append( db.insert_id()  )

def fast_insert_in_db(self, db, fragment_charge, transition_table):
    c = db.cursor()
    peptide = self.ass_peptide
    parent_id = peptide.parent_id
    vals = "type, parent_key, q3_charge, q3"
    q = "insert into %s (%s)" % (transition_table, vals)  + " VALUES (%s,%s,%s,%s)" 
    ch = fragment_charge
    charged_y =  [ ( pred + (ch -1)*self.mass_H)/ch for pred in self.y_series ]
    charged_b =  [ ( pred + (ch -1)*self.mass_H)/ch for pred in self.b_series ]
    many = [ ['y', parent_id, ch, q3] for i, q3 in enumerate(charged_y) ]
    manyb = [ ['b', parent_id, ch, q3] for i, q3 in enumerate(charged_b) ]
    many.extend( manyb )
    c.executemany( q, many)
    first_id = db.insert_id()
    self.fragment_ids[ ch ] = {
    'y' : [i for i in range(first_id, first_id + len(charged_y )) ] , 
    'b' : [i for i in range(first_id + len(charged_y ),  first_id + 2*len(charged_y )) ]
    }

def all_calculate_clashes_in_series_insert_db( S, S2, pairs_dic, 
                      MS2_bins, clash_distr, range_small, c, table, frag_range):
    for ch1 in frag_range:
        for ch2 in frag_range:
            calculate_clashes_in_series_insert_db( S, S2, ch1, ch2, pairs_dic, 
                                 MS2_bins, clash_distr, range_small, c, table )

class NonExecuteableDatabaseCursor:
    def __init__(self): pass
    def execute(self, t): print t

def calculate_clashes_in_series_insert_db(S, S2, charge1, charge2, pairs_dic, 
                        MS2_bins, clash_distr, myrange, c, table):
    """Compares all b/y ions of spectra S against those of peptide S2.
    If it finds fragments that are within the MS2_bins distance, it will 
    increase share (the return value) by one and flag the corresponding 
    fragment in BOTH associated peptides as non-unique.
    """
    #charge1 = ch1
    #charge2 = ch2
    vals = "coll_srm1, coll_srm2"
    #fragment_ids[ 0 ] = b, [1] = y
    share = 0
    pep1 = S.ass_peptide
    pep2 = S2.ass_peptide
    ch = charge1
    charged_y_series1=[ ( pred + (ch -1)*S.mass_H)/ch for pred in S.y_series ]
    charged_b_series1=[ ( pred + (ch -1)*S.mass_H)/ch for pred in S.b_series ]
    fragment_ids1 = S.fragment_ids[ ch ]
    ch = charge2
    charged_y_series2=[ ( pred + (ch -1)*S2.mass_H)/ch for pred in S2.y_series ]
    charged_b_series2=[ ( pred + (ch -1)*S2.mass_H)/ch for pred in S2.b_series ]
    fragment_ids2 = S2.fragment_ids[ ch ]
    #c = NonExecuteableDatabaseCursor() #for testing
    #every time we count share + 1 we cannot use one b or y
    #only exception: if it clashes with BOTH series, but how
    #often does that happen?
    for i, y in enumerate( charged_y_series1 ):
        if int(y) not in myrange: continue
        for kk, yy in enumerate( charged_y_series2 ):
            if abs( y - yy) < MS2_bins: 
                c.execute("insert into hroest.%s (%s) VALUES (%s,%s)" % ( 
                    table, vals, fragment_ids1[1][i], fragment_ids2[1][kk] )  )
        for kk, bb in enumerate( charged_b_series2 ):
            if abs( y - bb) < MS2_bins: 
                c.execute("insert into hroest.%s (%s) VALUES (%s,%s)" % ( 
                    table, vals, fragment_ids1[1][i], fragment_ids2[0][kk] )  )
    for i, b in enumerate( charged_b_series1 ):
        if int(b) not in myrange: continue
        for kk, yy in enumerate( charged_y_series2 ):
            if abs( b - yy) < MS2_bins: 
                c.execute("insert into hroest.%s (%s) VALUES (%s,%s)" % (
                    table, vals, fragment_ids1[0][i], fragment_ids2[1][kk] )  )
        for kk, bb in enumerate( charged_b_series2 ):
            if abs( b - bb) < MS2_bins: 
                c.execute("insert into hroest.%s (%s) VALUES (%s,%s)" % ( 
                    table, vals, fragment_ids1[0][i], fragment_ids2[0][kk] ) )
    if share > 0:
        pep1.shared_trans.append( pep2 )
        pep2.shared_trans.append( pep1 )
        #if mypeptide.spectrum in (S, S2): 
        #print "share", share, S.ass_peptide.sequence, S2.ass_peptide.sequence, sorted( mypeptide.non_unique)
    #if share > len( S.y_series):
    #    S.ass_peptide.pairs.append( S2 )
    #    S2.ass_peptide.pairs.append( S )
    #    if pairs_dic.has_key( S.peptide ): 
    #        pairs_dic[ S.peptide ].append( S2.peptide )
    #    else: pairs_dic[ S.peptide ] = [ S2.peptide ]
    clash_distr[ share ] += 1
    return share

def reset_pairs_unique(mass_bins):
    #reset pairs and non_unique entries to start level
    i = 0
    for b in mass_bins:
        i += len( b )
        for pep in b:
            pep.non_unique = {}
            pep.pairs = []
            pep.shared_trans = []

def calculate_clashes_in_series(S, S2, charge1, charge2, pairs_dic, 
                        MS2_bins, clash_distr, myrange):
    """Compares all b/y ions of spectra S against those of peptide S2.
    If it finds fragments that are within the MS2_bins distance, it will 
    increase share (the return value) by one and flag the corresponding 
    fragment in BOTH associated peptides as non-unique.
    """
    #charge1 = ch1
    #charge2 = ch2
    share = 0
    pep1 = S.ass_peptide
    pep2 = S2.ass_peptide
    ch = charge1
    charged_y_series1=[ ( pred + (ch -1)*S.mass_H)/ch for pred in S.y_series ]
    charged_b_series1=[ ( pred + (ch -1)*S.mass_H)/ch for pred in S.b_series ]
    ch = charge2
    charged_y_series2=[ ( pred + (ch -1)*S2.mass_H)/ch for pred in S2.y_series ]
    charged_b_series2=[ ( pred + (ch -1)*S2.mass_H)/ch for pred in S2.b_series ]
    #every time we count share + 1 we cannot use one b or y
    #only exception: if it clashes with BOTH series, but how
    #often does that happen?
    for i, y in enumerate( charged_y_series1 ):
        if int(y) not in myrange: continue
        for kk, yy in enumerate( charged_y_series2 ):
            if abs( y - yy) < MS2_bins: 
                share += 1
                pep1.non_unique[ 'y%s' % (i) ] = ''
                pep2.non_unique[ 'y%s' % (kk) ] = ''
        for kk, bb in enumerate( charged_b_series2 ):
            if abs( y - bb) < MS2_bins: 
                share += 1
                pep1.non_unique[ 'y%s' % (i) ] = ''
                pep2.non_unique[ 'b%s' % (kk) ] = ''
    for i, b in enumerate( charged_b_series1 ):
        if int(b) not in myrange: continue
        for kk, yy in enumerate( charged_y_series2 ):
            if abs( b - yy) < MS2_bins: 
                share += 1
                pep1.non_unique[ 'b%s' % (i) ] = ''
                pep2.non_unique[ 'y%s' % (kk) ] = ''
        for kk, bb in enumerate( charged_b_series2 ):
            if abs( b - bb) < MS2_bins: 
                share += 1
                pep1.non_unique[ 'b%s' % (i) ] = ''
                pep2.non_unique[ 'b%s' % (kk) ] = ''
    if share > 0:
        pep1.shared_trans.append( pep2 )
        pep2.shared_trans.append( pep1 )
        #if mypeptide.spectrum in (S, S2): 
        #print "share", share, S.ass_peptide.sequence, S2.ass_peptide.sequence, sorted( mypeptide.non_unique)
    #if share > len( S.y_series):
    #    S.ass_peptide.pairs.append( S2 )
    #    S2.ass_peptide.pairs.append( S )
    #    if pairs_dic.has_key( S.peptide ): 
    #        pairs_dic[ S.peptide ].append( S2.peptide )
    #    else: pairs_dic[ S.peptide ] = [ S2.peptide ]
    clash_distr[ share ] += 1
    return share


#make a list of all peptides

#old analysis


all_pep = []
for pp in mass_bins: 
    for p in pp:
        all_pep.append( p )


#calculate distribution of non-unique transitions
shared_distr = []
shared_rel = []
for b in mass_bins:
    for pep in b:
        shared_distr.append(len(pep.non_unique) )
        shared_rel.append( len( pep.non_unique) / (2.0* len(pep.sequence) - 2) )

h, n = numpy.histogram( shared_rel)

#print distribution in Q1 dimension
f = open( 'mzdist.csv', 'w' )
w = csv.writer(f,  delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL) 
i= 0
for l in bin_length:
    i += 1
    w.writerow( [i * 0.7, l] )

#print distribution in RT dimension
f = open( 'rtdist.csv', 'w' )
w = csv.writer(f,  delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL) 
for i in range(50):
    w.writerow( [i - 50, rt_bin_len[ i - 50] ] )

i= 0
for l in rt_bin_len:
    w.writerow( [i, l] )
    i += 1

f.close()

#####################################################
#here we need the clash analysis


#this looks at the number of pairs per bin
f = open( 'mzclashes.csv', 'w' )
#w = csv.writer(f,  delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL) 
for i, bin in enumerate(mass_bins):
    l = [ p for p in bin if len( p.pairs) > 0]
    #f.write(  '%f %f\n' %  (i * 0.7, len(l)*100.0/(len(bin)+1) ) )
    f.write(  '%f %f\n' %  (i * 0.7, len(l) ) )

#this looks at the clash distribution (also pairs)
f = open( 'clashdist.csv', 'w' )
for i, cl in enumerate(clash_distr):
    f.write(  '%f %f\n' %  (i , cl ) )

#length distribution of the pairs (how many pairs/pep)
f = open( 'pairlendist.csv', 'w' )
sum( h )
pair_distr = [ len( p.pairs) for p in all_pep]
h, n = numpy.histogram( pair_distr, 50, (0,50)  )
for hh, nn in zip( h, n):
    f.write( '%d %f\n' % (nn, hh*100.0/sum(h) ))

f = open( 'lab_pairlendist.csv', 'w' )
pair_distr = [ len( p.pairs) for p in lab_peptides_2]
h, n = numpy.histogram( pair_distr, 50, (0,50)  )
for hh, nn in zip( h, n):
    f.write( '%d %f\n' % (nn, hh*100.0/sum(h) ))









#calculate how many unique transitions are left per peptide
unique_dist = [0 for i in range(1000)]
unique_ratios = []
unique_dist_window = [0 for i in range(1000)]
unique_ratios_window = []
for i, p in enumerate(all_pep):
    unique = len( p.sequence ) - len( p.non_unique)/2 - 1
    ratio = unique * 1.0 / ( len(p.sequence) - 1)
    unique_dist[ unique ] += 1
    unique_ratios.append( ratio)
    if i % 10000 == 0: print i

f = open( 'unique_transitions.csv', 'w' )
for i,u in enumerate(unique_dist):
    f.write( '%d %d\n' % (i,u) )

#calculate transitions per peptide as ratios
f = open( 'unique_ratios.csv', 'w' )
h, n = numpy.histogram( unique_ratios, 10)
for hh, nn in zip( h, n):
    f.write( '%f %f\n' % (nn*100.0 + 5 , hh*100.0/sum(h) ))

f.close()

#All but the last (righthand-most) bin is half-open. In other words, if bins is:
# [1, 2, 3, 4]
#then the first bin is [1, 2) (including 1, but excluding 2) and the second 
#[2, 3). The last bin, however, is [3, 4], which includes 4.]]

[sum( clash_distr[:i] ) for i in range(50)]

lab_peptides_3 = utils.unique( lab_peptides_2)

lab_peptides_3 = []
for pep in lab_peptides_2:
    if not pep in lab_peptides_3: 
        lab_peptides_3.append( pep )

len( lab_peptides_3 )
len( lab_peptides_2 )



tots = {}
seqs = {}
for p in lab_peptides_2: 
    tots[ p.sequence ] = ''
    if len( p.pairs) > 0: 
        print p.sequence
        seqs[ p.sequence ] = ''

f.close()

for p in lab_peptides_2: 
    if p.sequence == 'DIHFMPCSGLTGANIK': break

h, n = numpy.histogram( pair_distr )
[p for p in pair_distr if p >0 ]

len( lab_peptides_2)




