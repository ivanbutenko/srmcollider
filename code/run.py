#!/usr/bin/python
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

#with this method, I can do 
##method 2
# 22.2 per second, thus do 500k in 6 hours ==> using per peptide methods
##method 1
# 0.5 per second, thus do 500k in 240 hours ==> using per transition methods
do_1vs = True #check only one charge state?
do_vs1 = False #background only one charge state?
q3_range = [300, 1500]
#q3_range = [0, 10000 ]
ssrcalc_window = 2.0 / 2
q1_window = 0.2 / 2
#q1_window = 0.7 / 2
q3_window = 1.0 / 2
do_1_only = "and q1_charge = 2 and q3_charge = 1"
experiment_type = """check all four charge states [%s] vs all four charge states [%s]
with thresholds of SSRCalc %s Q1 %s Q3 %s and a range of %s - %s Da for the q3 transitions.
""" % ( not do_1vs, not do_vs1, ssrcalc_window*2,  q1_window*2, q3_window*2, q3_range[0], q3_range[1])
print experiment_type
if True:
    if do_1vs : query_add = "where q1_charge = 2"; query1_add = do_1_only
    else: query_add = ""
    if do_vs1 : query2_add = do_1_only
    else: query2_add = ""

#
###
#
#here we start
cursor = db.cursor()
query = """
select parent_id
 from hroest.srmPeptide
 %s
""" % query_add
cursor.execute( query )
pepids = cursor.fetchall()
allpeps = [ 0 for i in range(2 * 10**6)]
start = time.time()
for i, p_id in enumerate(pepids):
    if i % 1000 ==0: print i
    query1 = """
    select q1, ssrcalc, q3, srm_id
    from hroest.srmPeptide
    inner join hroest.srmTransitions
      on parent_id = parent_key
    where parent_id = %(parent_id)s
    and q3 > %(q3_low)s and q3 < %(q3_high)s         %(query_add)s
    """ % { 'parent_id' : p_id[0], 'q3_low' : q3_range[0],
           'q3_high' : q3_range[1], 'query_add' : query1_add }
    nr_transitions = cursor.execute( query1 )
    transitions = cursor.fetchall()
    q1 = transitions[0][0]
    ssrcalc = transitions[0][1]
    non_unique = {}
    transitions = [ [t[2], t[3]] for t in transitions ]
    ####method 2 (fast) per protein
    query2 = """
    select q3 from hroest.srmPeptide
    inner join hroest.srmTransitions
      on parent_id = parent_key
    where   ssrcalc > %(ssrcalc)s - %(ssr_window)s 
        and ssrcalc < %(ssrcalc)s + %(ssr_window)s
    and q1 > %(q1)s - %(q1_window)s and q1 < %(q1)s + %(q1_window)s
    and parent_id != %(parent_id)d
    and q3 > %(q3_low)s and q3 < %(q3_high)s         %(query_add)s
    """ % { 'q1' : q1, 'ssrcalc' : ssrcalc, 'parent_id' : p_id[0] ,
           'q3_low':q3_range[0],'q3_high':q3_range[1], 'q1_window' : q1_window,
           'query_add' : query2_add, 'ssr_window' : ssrcalc_window }
    tmp = cursor.execute( query2 )
    collisions = cursor.fetchall()
    #here we loop through all possible combinations of transitions and 
    #potential collisions and check whether we really have a collision
    #TODO it might even be faster to switch the loops
    #TODO since the outer loop can be aborted as soon as a collision is found
    for c in collisions:
        for t in transitions:
            #here, its enough to know that one transition is colliding
            if abs( t[0] - c[0] ) < q3_window: non_unique[ t[1] ] = 0; break
    #####method 1 (slow) per transition
    ##for t in transitions:
    ##    query2 = """
    ##    select q3 from hroest.srmPeptide
    ##    inner join hroest.srmTransitions
    ##      on parent_id = parent_key
    ##    where ssrcalc > %(ssrcalc)s - 2 and ssrcalc < %(ssrcalc)s + 2
    ##    and q1 > %(q1)s - 0.35 and q1 < %(q1)s + 0.35
    ##    and parent_id != %(parent_id)d
    ##    and q3 > %(q3)s - 0.5 and q3 < %(q3)s + 0.5
    ##    ;
    ##    """ % { 'q1' : q1, 'ssrcalc' : ssrcalc, 'parent_id' : p_id[0], 'q3' : t[0] }
    ##    tmp = cursor.execute( query2 )
    ##    if tmp > 0: non_unique[ t[1] ] = 0
    allpeps[ p_id[0] ] = 1.0 - len( non_unique ) * 1.0  / nr_transitions
    end = time.time()


###########################################################################
#
# Analysis and printing
if True:
    #new method, clash distribution
    mydist = []
    for i, p_id in enumerate(pepids):
        #this is the number of UNIQUE transitions
        mydist.append( allpeps[p_id[0] ] ) 
    len(mydist)
    h, n = numpy.histogram( mydist , 10)
    n = [nn * 100.0 + 5 for nn in n]
    h = [ hh *100.0 / len(mydist) for hh in h]
    import gnuplot
    reload( gnuplot )
    filename = '%s_%s_%d%d%d_range%sto%s' % (do_1vs, do_vs1, 
        ssrcalc_window*20,  q1_window*20, q3_window*20, q3_range[0], q3_range[1])
    gnu = gnuplot.Gnuplot.draw_boxes_from_data( [h,n], filename + '.eps',
      'Unique transitions per peptide / %' , 'Occurence / %', keep_data = True,
                                         tmp_csv = filename + '.csv' )
    gnu.add_to_body( "set yrange[0:30]" )
    gnu.draw_boxes()
    #title = '(2+,1+) transitions, Background of 4 combinations [Range 300 to 1500 Da]' )


#get a peptide length distribution
if True:
    cursor = db.cursor()
    query = """
    select parent_id, LENGTH( sequence )
     from hroest.srmPeptide
     inner join ddb.peptide on peptide.id = peptide_key
     %s
    """ % query_add
    cursor.execute( query )
    pepids = cursor.fetchall()
    len_distr  = [ [] for i in range(255) ]
    for zz in pepids:
        p_id, peplen = zz
        len_distr[ peplen ].append( allpeps[ p_id ]  )
    avg_unique_per_length = [ 0 for i in range(255) ]
    for i, l in enumerate(len_distr):
        if len(l) > 0:
            avg_unique_per_length[i] =    sum(l) / len(l) 
    h = avg_unique_per_length[6:26]
    h = [ hh *100.0 for hh in h]
    n = range(6,26)
    reload( gnuplot )
    filename = '%s_%s_%d%d%d_range%sto%s_lendist' % (do_1vs, do_vs1, 
        ssrcalc_window*20,  q1_window*20, q3_window*20, q3_range[0], q3_range[1])
    gnu  = gnuplot.Gnuplot.draw_boxes_from_data( [h,n], filename + '.eps',
      'Peptide length', 'Average number of unique transitions per peptide / %', keep_data = True )
    gnu.add_to_body( "set xrange[%s:%s]"  % (4, 26) )
    gnu.add_to_body( "set yrange[0:100]" )
    gnu.draw_boxes()


cursor = db.cursor()
query = """
select parent_id, LENGTH( sequence )
 from hroest.srmPeptide
 inner join ddb.peptide on peptide.id = peptide_key
 %s
""" % query_add
cursor.execute( query )
pepids = cursor.fetchall()


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




