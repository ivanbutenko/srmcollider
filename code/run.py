#!/usr/bin/python
import os
import MySQLdb
import sys 
sys.path.append( '/home/ghost/hroest/code/' )
sys.path.append( '/home/ghost/software_libs/' )
sys.path.append( '/home/hroest/lib/' )
sys.path.append( '/home/hroest/msa/code/tppGhost' )
sys.path.append( '/home/hroest/msa/code' )
sys.path.append( '/home/hroest/srm_clashes/code' )
#import mzXMLreader
#import pepXMLReader
import time
from utils_h import utils
#import pipeline 
#db = MySQLdb.connect(read_default_file="~/hroest/.my.cnf")
db = MySQLdb.connect(read_default_file="~/.my.cnf")
c = db.cursor()
c2 = db.cursor()
cursor = db.cursor()
import silver
import DDB 
import csv 
R = silver.Residues.Residues('mono')
import numpy
import collider

##testcase
###########################################################################
par  = testcase()
collider = SRMcollider()
collider.find_clashes(cursor, par)

print "I ran the testcase in %ss" % collider.total_time
assert sum( collider.allpeps.values() ) - 975.6326447245566 < 10**(-3)
assert collider.non_unique_count == 26
assert collider.total_count == 12502

#Run the collider
###########################################################################
par = SRM_parameters()
par.q1_window = 0.7 / 2  
par.q3_window = 1.0 / 2  
par.ppm = False
par.eval()
print par.experiment_type
print par.get_common_filename()

collider = SRMcollider()
collider.find_clashes(cursor, par)

collider = SRMcollider()
directory = '/home/hroest/srm_clashes/results/pedro/'
collider.load_from_file( par, directory)
collider.print_unique_histogram(par)
collider.print_cumm_unique(par)
collider.print_q3min(par)
collider.print_q3min_ppm(par)

print "I ran the collider for %s min" % ((end - start)/60 )
print "I ran the collider for %s h" % ((end - start)/3600 )

print_stats(collider)

def print_stats(self):
    #print "I ran the collider for %ss" % self.total_time
    print "Nonunique / Total transitions : %s / %s = %s" % (self.non_unique_count, self.total_count, self.non_unique_count * 1.0 /self.total_count)

total_time = end - start
def print_stats():
    print "I ran the collider for %ss" % total_time
    print "Nonunique / Total transitions : %s / %s = %s" % (non_unique_count, total_count, non_unique_count * 1.0 /total_count)


print_stats()

some 

###########################################################################
#
# Storage of results
self = par

restable = 'hroest.srm_results_' + common_filename
cursor.execute("drop table %s" % restable)
cursor.execute("create table %s (parent_key int, unique_transitions double) " % restable)
cursor.executemany( "insert into %s values " % restable + "(%s,%s)", 
              allpeps.items() )


#p_id = 1
#{1L: 3830996L, 11L: 13837110L, 9L: 3881908L, 12L: 13837111L, 7L: 3881910L}
#{1L: 0.0024949999999535066, 11L: 0.0, 9L: -0.00089500000001407898, 12L: 0.0, 7L: -0.00089500000001407898}

#select * from hroest.srmTransitions_yeast inner join hroest.srmPeptides_yeast
#on parent_key = parent_id inner join  ddb.peptide on peptide_key = peptide.id
#where srm_id in (1, 3830996);



###########################################################################
#
# Analysis and printing

def get_cum_dist(original):
    cumulative = []
    cum  = 0
    for i, val in enumerate(original):
        cum += val
        cumulative.append(  cum )
    return cumulative


import gnuplot
import numpy



#print('Percentage of collisions below 1 ppm: %02.2f~\%%' % (sum( h[40:60] )* 100.0 / sum( h )))
#exact_zero = len( [m for m in q3min_distr if m == 0.0 ]) * 1.0 
#perc_exact_zero = exact_zero / len( q3min_distr ) 
#perc_total_interferences = non_unique_count * 1.0 / total_count
#print "q3min distr: exact_zero, percentage, percentage interferences, total_coun"
#print int(exact_zero), perc_exact_zero, perc_total_interferences, total_count


###########################################################################
## STOP
###########################################################################

#get a peptide length distribution
if True:
    cursor = db.cursor()
    query = """
    select parent_id, LENGTH( sequence )
     from %s
     inner join ddb.peptide on peptide.id = peptide_key
     %s
    """ % (peptide_table, query_add )
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
    filename = 'yeast_%s_%s_%d%d%d_range%sto%s_lendist' % (do_1vs, do_vs1, 
        ssrcalc_window*20,  q1_window*20, q3_window*20, q3_range[0], q3_range[1])
    gnu  = gnuplot.Gnuplot.draw_boxes_from_data( [h,n], filename + '.eps',
      'Peptide length', 'Average number of unique transitions per peptide / %', keep_data = True )
    gnu.add_to_body( "set xrange[%s:%s]"  % (4, 26) )
    gnu.add_to_body( "set yrange[0:100]" )
    gnu.draw_boxes()



reload(gnuplot)
c2.execute( "select floor(q1), count(*) from %s group by floor(q1)" % peptide_table)
bins = c2.fetchall()
filename = 'yeast_mzdist'
n = [nn[0] for nn in bins]
h = [hh[1] for hh in bins]
gnu  = gnuplot.Gnuplot.draw_points_from_data( [h,n], filename + '.eps',
  'Th', 'Occurence', keep_data = True)
gnu.add_to_body( "set xrange[400:]"   )
gnu.draw_points()


reload(gnuplot)
c2.execute( "select floor(ssrcalc), count(*) from %s group by floor(ssrcalc)" % peptide_table)
bins = c2.fetchall()
filename = 'yeast_rtdist'
n = [nn[0] for nn in bins]
h = [hh[1] for hh in bins]
gnu  = gnuplot.Gnuplot.draw_points_from_data( [h,n], filename + '.eps',
  'Retention time', 'Occurence', keep_data = True)
gnu.add_to_body( "set xrange[0:100]"   )
gnu.draw_points(2, True)



#print distribution in RT dimension
f = open( 'rtdist.csv', 'w' )
w = csv.writer(f,  delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL) 
for i in range(50):
    w.writerow( [i - 50, rt_bin_len[ i - 50] ] )



#print distribution in Q1 dimension
f = open( 'mzdist.csv', 'w' )
w = csv.writer(f,  delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL) 
for l in bins:
    w.writerow( [l[0], l[1] ] )

f.close()
