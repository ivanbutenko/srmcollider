#!/usr/bin/python

"""
There are two different kinds of setups:
    A) store all parent ions and their transitions (srmPeptides and srmTransitions)
    B) also store all collisions (in a table like srmCollisions400710)
"""

import MySQLdb
import sys 
sys.path.append( '/home/hroest/lib/' )
sys.path.append( '/home/hroest/srm_clashes/code' ) #Collider
sys.path.append( '/home/hroest/msa/code/tppGhost' ) #DDB
sys.path.append( '/home/hroest/lib' ) #utils
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
#exp_key = 3130  #human (all)
exp_key = 3120  #yeast
exp_key = 3131  #yeast (all)
all_peptide_query  = """
select distinct peptide.sequence, molecular_weight, ssrcalc,
genome_occurence, peptide.id as peptide_key
from peptide 
inner join peptideOrganism on peptide.id = peptideOrganism.peptide_key
inner join compep.ssrcalc_prediction on ssrcalc_prediction.sequence =
peptide.sequence
where experiment_key = %s
and length( peptide.sequence ) > 1
""" % exp_key

#human go from parent_id = 1 to parent_id = 1177958

###################################
# A) store the transitions
###################################
#read in the data from DB and put into bins
# 75 peptides/sec
insert_db = True
modify_cysteins = True
read_from_db = False
c.execute( 'use ddb;' )
t = utils.db_table( c2 )
t.read( all_peptide_query )
print "executed the query"
rows = t.c.fetchall()
#1000 entries / 9 s with fast (executemany), 10.5 with index
#1000 entries / 31 s with normal (execute), 34 with index
transition_table = 'hroest.srmTransitions_test_ionseries'
peptide_table = 'hroest.srmPeptides_test_ionseries'
mass_bins = [ []  for i in range(0, 10000) ]
rt_bins = [ []  for i in range(-100, 500) ]
tmp_c  = db.cursor()
progressm = progress.ProgressMeter(total=len(rows), unit='peptides')
start = time.time()
for i,row in enumerate(rows):
    progressm.update(1)
    #if i >= 1000: break #####FOR TESTING ONLY
    for mycharge in [2,3]:  #precursor charge 2 and 3
        peptide = collider.get_peptide_from_table(t, row)
        peptide.charge = mycharge
        if modify_cysteins: peptide.modify_cysteins()
        S = silver.Spectrum.Spectrum(SEQUEST_mode =1 )
        S.ion_charge = peptide.charge
        S.construct_from_peptide( 
            peptide.get_modified_sequence('SEQUEST'), R.residues, 
            R.res_pairs)
        S.ass_peptide = peptide
        S.fragment_ids = {}
        if insert_db:
            #insert peptide into db
            collider.insert_peptide_in_db(S, db, peptide_table, transition_group=i)
        elif read_from_db:
            #instead of inserting, we read the database
            assert False
            #not implemented with new table layout
            collider.read_fragment_ids(S, tmp_c, peptide, peptide_table, transition_table)
        #insert peptide into mass_bins hash
        peptide.spectrum = S
        mz = S.peptide_mass / S.ion_charge
        bin = int( mz / 0.7) ;
        mass_bins[ bin ].append( peptide )
        rt_bins[ int(peptide.ssr_calc) ].append( peptide )
        end = time.time()
    #we want to insert the fragments only once per peptide
    if insert_db:
        #insert fragment charge 1 and 2 into database
        collider.fast_insert_in_db( S, db, 1, transition_table, transition_group=i)
        collider.fast_insert_in_db( S, db, 2, transition_table, transition_group=i)

#rr = [r for r in rows if r[-1] == 9201171]
#rr = [r for r in rows if r[-1] == 9255505]



#A2 create the additional parent ions for the isotope patterns
cursor = c
vals = "peptide_key, q1_charge, q1, modified_sequence, ssrcalc, isotope_nr"
query = "SELECT parent_id, %s FROM %s" % (vals, peptide_table)
cursor.execute( query )
allpeptides =  cursor.fetchall()
prepared = []
for p in allpeptides:
    q1_charge = p[2]
    for i in range(1,4):
        prepared.append( [p[1], p[2], p[3] + (R.mass_diffC13 * i* 1.0) / q1_charge,
                              p[4], p[5], i] )

q = "INSERT INTO %s (%s)" % (peptide_table, vals)  + \
     " VALUES (" + "%s," *5 + "%s)"
c.executemany( q, prepared)


###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################


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
                sh = collider.all_calculate_clashes_in_series_insert_db( S, S2, pairs_dic,
                     MS2_bins, clash_distr, range_all, c, table, frag_range )

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




