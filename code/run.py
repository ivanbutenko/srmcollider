import os
import MySQLdb
import sys 
sys.path.append( '/home/ghost/hroest/code/' )
sys.path.append( '/home/ghost/software_libs/' )
sys.path.append( '/home/hroest/lib/' )
sys.path.append( '/home/hroest/msa/code/tppGhost' )
sys.path.append( '/home/hroest/msa/code' )
import mzXMLreader
import pepXMLReader
import time
from utils_h import utils
import pipeline 
#db = MySQLdb.connect(read_default_file="~/hroest/.my.cnf")
db = MySQLdb.connect(read_default_file="~/.my.cnf")
c = db.cursor()
c2 = db.cursor()
import silver
import DDB 
import csv 
R = silver.Residues.Residues('mono')

exp_key = 3061  #human
#exp_key = 3120  #yeast
q  = """
select distinct gene.id as gene_id, peptide.sequence, molecular_weight, ssrcalc,
genome_occurence, 
peptide.id as peptide_key
from gene 
inner join experimentOrganism on gene.experiment_key = experimentOrganism.experiment_key
inner join geneProtLink on gene.id = geneProtLink.gene_key
inner join protein on protein.id = geneProtLink.protein_key 
inner join protPepLink on protein.id = protPepLink.protein_key
inner join peptide on protPepLink.peptide_key = peptide.id
inner join peptideOrganism on peptide.id = peptideOrganism.peptide_key
inner join compep.ssrcalc_prediction on ssrcalc_prediction.sequence =
peptide.sequence
where gene.experiment_key = %s
and taxonomy_id = 9606 #human
#and taxonomy_id = 4932 #yeast
#and genome_occurence = 1
""" % exp_key

c.execute( 'select * from hroest.srmClashes where experiment_key = 3061' )


#read in the data from DB
insert_db = False
c.execute( 'use ddb;' )
t = utils.db_table( c2 )
t.read( q )
mass_bins = [ []  for i in range(0, 10000) ]
print "executed the query"
i= 0
while True:
    row = t.fetchone()
    if row == None: break
    if i % 10000 == 0: print i
    peptide = get_peptide_from_table(t, row)
    for mycharge in [2,3]:  #precursor charge 2 and 3
        peptide.charge = mycharge
        S = silver.Spectrum.Spectrum(SEQUEST_mode =1 )
        S.ion_charge = peptide.charge
        S.construct_from_peptide( peptide.get_modified_sequence('SEQUEST'), 
                                 R.residues, R.res_pairs)
        S.peptide = peptide
        if insert_db:
            #insert fragment charge 1 and 2 into database
            insert_in_db( S, c, 1)
            insert_in_db( S, c, 2)
        #insert peptide into mass_bins hash
        S.ass_peptide = peptide
        peptide.spectrum = S
        mz = peptide.mw
        bin = int( mz / 0.7) ;
        mass_bins[ bin ].append( peptide )
    i += 1

def get_peptide_from_table(t, row):
    peptide = DDB.Peptide()
    peptide.set_sequence( t.row(row, 'sequence')  )
    peptide.genome_occurence = t.row( row, 'genome_occurence' )
    peptide.ssr_calc         = t.row( row, 'ssrcalc' )
    peptide.gene_id          = t.row( row, 'gene_id' )
    peptide.mw               = t.row( row, 'molecular_weight' )
    peptide.id               = t.row( row, 'peptide_key' )
    peptide.pairs            = []
    peptide.non_unique = {}
    return peptide

def insert_in_db(self, c, fragment_charge):
    assert False
    peptide = self.peptide
    vals = "type, peptide_key, experiment_key, q1_charge, q3_charge, q1, q3, ssrcalc"
    ch = fragment_charge
    charged_y_series =  [ ( pred + (ch -1)*S.mass_H)/ch for pred in S.y_series ] 
    charged_b_series =  [ ( pred + (ch -1)*S.mass_H)/ch for pred in S.b_series ] 
    for i, q3 in enumerate(charged_b_series):
        type = 'b%s%s' % (i+1, '+' * fragment_charge )
        q = "insert into hroest.srmClashes (%s) VALUES ('%s',%s,%s,%s,%s,%s,%s,%s);" % (
            vals, 
            type, peptide.id, exp_key, peptide.charge, fragment_charge, 
            get_actual_mass(self),q3, peptide.ssr_calc )
        c.execute(q)
    pepLen = len(peptide)
    for i, q3 in enumerate(charged_y_series):
        type = 'y%s%s' % (pepLen - i - 1, '+' * fragment_charge )
        q = "insert into hroest.srmClashes (%s) VALUES ('%s',%s,%s,%s,%s,%s,%s,%s);" % (
            vals, 
            type, peptide.id, exp_key, peptide.charge, fragment_charge, 
            get_actual_mass(self),q3, peptide.ssr_calc )
        c.execute(q)

def get_actual_mass(self):
    return S.peptide_mass / S.ion_charge

def get_all_peptides( mass_bins):
    #make a list of all peptides
    all_pep = []
    for pp in mass_bins: 
        for p in pp:
            all_pep.append( p )
    return all_pep

def reset_pairs_unique(mass_bins):
    #reset pairs and non_unique entries to start level
    i = 0
    for b in mass_bins:
        i += len( b )
        for pep in b:
            pep.non_unique = {}
            pep.pairs = []
            pep.shared_trans = []

###########################################################################

#here we run through all pairs
reset_pairs_unique(mass_bins)
MS1_bins = 0.7 /2
MS2_bins = 1.0 /2
ssrcalc_binsize = 4/2 #that may be around 2 minutes, depending on coloumn
lab_peptides_2 = []
pairs_dic = {}
clash_distr = [0 for i in range( 1000 ) ]
for i in range( 1, len( mass_bins )-1):
    if i % 100 ==0:
        print i, len( pairs_dic )
    current = []
    current.extend( mass_bins[i] )
    current.extend( mass_bins[i-1] )
    current.extend( mass_bins[i+1] )
    current_spectra = []
    for peptide in current: current_spectra.append( peptide.spectrum )
    for S in current_spectra:
        peptide = S.ass_peptide
        if peptide.sequence in lab_peptides: lab_peptides_2.append( peptide )
        mz = peptide.mw
        bin = int( mz / 0.7) 
        #only use those as S which are in the current bin
        if bin != i: continue
        mz = S.peptide_mass / S.ion_charge
        #compare query spectra S against all other spectra that are within 
        #a certain tolerance of q1 mass
        for S2 in current_spectra:
            #we are not interested in spectras from the same sequence
            if S2.ass_peptide.sequence == S.ass_peptide.sequence: continue
            #we measure it only once, always S is the one with the higher mass
            if S2.peptide_mass < S.peptide_mass : continue
            if (abs( mz - S2.peptide_mass / S.ion_charge) < MS1_bins and
               abs(peptide.ssr_calc - S2.ass_peptide.ssr_calc)<ssrcalc_binsize):
                calculate_shared_transitions( S, S2, pairs_dic, 
                                             MS2_bins, clash_distr )

#calculate distribution of non-unique transitions
shared_distr = []
shared_rel = []
for b in mass_bins:
    for pep in b:
        shared_distr.append(len(pep.non_unique) )
        shared_rel.append( len( pep.non_unique) / (2.0* len(pep.sequence) - 2) )

h, n = numpy.histogram( shared_rel)

#read in lab-peptides, bin length etc
bin_length = [ len( b ) for b in mass_bins]
q = "select sequence from hroest.lab_peptides"
t = utils.db_table( c )
t.read( q )
lab_peptides = []
for row in t.rows():
    lab_peptides.append( t.row( row, 'sequence' ))

f = open( 'mzdist.csv', 'w' )
w = csv.writer(f,  delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL) 
i= 0
for l in bin_length:
    i += 1
    w.writerow( [i * 0.7, l] )

f = open( 'mzclashes.csv', 'w' )
#w = csv.writer(f,  delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL) 
for i, bin in enumerate(mass_bins):
    l = [ p for p in bin if len( p.pairs) > 0]
    #f.write(  '%f %f\n' %  (i * 0.7, len(l)*100.0/(len(bin)+1) ) )
    f.write(  '%f %f\n' %  (i * 0.7, len(l) ) )

f = open( 'clashdist.csv', 'w' )
for i, cl in enumerate(clash_distr):
    f.write(  '%f %f\n' %  (i , cl ) )

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
for p in all_pep:
    unique = len( p.sequence ) - len( p.non_unique)/2 - 1
    ratio = unique * 1.0 / ( len(p.sequence) - 1)
    unique_dist[ unique ] += 1
    unique_ratios.append( ratio)

f = open( 'unique_transitions.csv', 'w' )
for i,u in enumerate(unique_dist):
    f.write( '%d %d\n' % (i,u) )

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


def calculate_shared_transitions(S, S2, pairs_dic, MS2_bins, clash_distr):
    """Compares all b/y ions of spectra S against those of peptide S2.
    If it finds fragments that are within the MS2_bins distance, it will 
    increase share (the return value) by one and flag the corresponding 
    fragment in BOTH associated peptides as non-unique.
    """
    share = 0
    pep1 = S.ass_peptide
    pep2 = S2.ass_peptide
    #every time we count share + 1 we cannot use one b or y
    #only exception: if it clashes with BOTH series, but how
    #often does that happen?
    for i, y in enumerate(S.y_series):
        for kk, yy in enumerate(S2.y_series):
            if abs( y - yy) < MS2_bins: 
                share += 1
                pep1.non_unique[ 'y%s' % (i) ] = ''
                pep2.non_unique[ 'y%s' % (kk) ] = ''
        for kk, bb in enumerate(S2.b_series):
            if abs( y - bb) < MS2_bins: 
                share += 1
                pep1.non_unique[ 'y%s' % (i) ] = ''
                pep2.non_unique[ 'b%s' % (kk) ] = ''
    for i, b in enumerate(S.b_series):
        for yy in S2.y_series:
            if abs( b - yy) < MS2_bins: 
                share += 1
                pep1.non_unique[ 'b%s' % (i) ] = ''
                pep2.non_unique[ 'y%s' % (kk) ] = ''
        for kk, bb in enumerate(S2.b_series):
            if abs( b - bb) < MS2_bins: 
                share += 1
                pep1.non_unique[ 'b%s' % (i) ] = ''
                pep2.non_unique[ 'b%s' % (kk) ] = ''
    if share > 0:
        pep1.shared_trans.append( pep2 )
        pep2.shared_trans.append( pep1 )
        #if mypeptide.spectrum in (S, S2): 
        #print "share", share, S.ass_peptide.sequence, S2.ass_peptide.sequence, sorted( mypeptide.non_unique)
    if share > len( S.y_series):
        S.ass_peptide.pairs.append( S2 )
        S2.ass_peptide.pairs.append( S )
        if pairs_dic.has_key( S.peptide ): 
            pairs_dic[ S.peptide ].append( S2.peptide )
        else: pairs_dic[ S.peptide ] = [ S2.peptide ]
    clash_distr[ share ] += 1
    return share




