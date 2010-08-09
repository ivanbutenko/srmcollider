#!/usr/bin/python

# This script will plot the distribution of the difference to the next point in
# space for Q1, Q3 and SSRCalc.
# For each dimension, it will thus show, how close on average the next
# transition is. To execute it, you need the tables
# hroest.tmp_srm_get_all_collisions and hroest.tmp_srm_all_collisions_filtered
# which can be created using the SQL script.
import MySQLdb, sys, numpy
db = MySQLdb.connect(read_default_file="~/.my.cnf")
c = db.cursor()
sys.path.append( '/home/hroest/msa/code/tppGhost' )
import gnuplot

coll_table = 'hroest.tmp_srm_all_collisions_filtered'
query_up = """
-- create temporary tables for the plotting
drop table if exists hroest.tmp_analyse_q1diff;
create table hroest.tmp_analyse_q1diff as
select q1diff, q1diff_abs, coll_srm1 from %(coll_table)s
where q1diff < 0.1

order by coll_srm1, q1diff_abs 
;
drop table if exists hroest.tmp_analyse_q3diff;
create table hroest.tmp_analyse_q3diff as
select q3diff, q3diff_abs, coll_srm1 from %(coll_table)s
#where q1diff < 0.1
order by coll_srm1, q3diff_abs 
;
drop table if exists hroest.tmp_analyse_ssrcalcdiff;
create table hroest.tmp_analyse_ssrcalcdiff as
select ssrcalcdiff, ssrcalcdiff_abs, coll_srm1 from %(coll_table)s
#where q1diff < 0.1
order by coll_srm1, ssrcalcdiff_abs 
;

alter table hroest.tmp_analyse_ssrcalcdiff add index(coll_srm1);
alter table hroest.tmp_analyse_q1diff add index(coll_srm1);
alter table hroest.tmp_analyse_q3diff add index(coll_srm1);
select * from hroest.tmp_analyse_q1diff
""" % { 'coll_table' : coll_table}


query_down = """
--#UNDO all of the above
drop table hroest.tmp_analyse_q1diff;
drop table hroest.tmp_analyse_q3diff;
drop table hroest.tmp_analyse_ssrcalcdiff;
""" 

##########################################################################
c.execute( query_up)
print "Query up executed"

c.close()
c = db.cursor()
c.execute("""
select ssrcalcdiff from hroest.tmp_analyse_ssrcalcdiff
group by coll_srm1
"""
)
my_ssrcalcs = [all[0] for all in c.fetchall()]

h, n = numpy.histogram( my_ssrcalcs, 100, (-2.0, 2.0)  )
reload( gnuplot )
filename = 'ssrcalc_diff_distribution_1vs1' 
gnu  = gnuplot.Gnuplot.draw_boxes_from_data( [h,n], filename + '.eps',
  'SSRCalc difference / arbitrary units', 'Occurence', keep_data = True )

print "Query 1 executed"
c.close()
c = db.cursor()
c.execute("""
select q1diff from hroest.tmp_analyse_q1diff 
group by coll_srm1
"""
)
my_q1s = [all[0] for all in c.fetchall()]


h, n = numpy.histogram( my_q1s, 100, (-0.35, 0.35)  )
h = [ q1 * 1.0 / 10**6 for q1 in h]
reload( gnuplot )
filename = 'q1_diff_distribution_1vs1' 
gnu  = gnuplot.Gnuplot.draw_boxes_from_data( [h,n], filename + '.eps',
  'Q1 difference / Th', 'Occurence / 10^6', keep_data = True )
gnu.add_to_body( "set xrange[%s:%s]"  % (-0.35, 0.35) )
gnu.add_to_body( "set xtics -0.3, 0.1" )
gnu.draw_boxes()

c.close()
c = db.cursor()
my_q3s = c.execute("""
select q3diff from hroest.tmp_analyse_q3diff 
group by coll_srm1
"""
)
my_q3s = [all[0] for all in c.fetchall()]

h, n = numpy.histogram( my_q3s, 100, (-0.5, 0.5)  )
h = [q3 * 1.0 / 10**6 for q3 in h]
reload( gnuplot )
filename = 'q3_diff_distribution_1vs1' 
gnu  = gnuplot.Gnuplot.draw_boxes_from_data( [h,n], filename + '.eps',
  'Q3 difference / Th', 'Occurence / 10^6', keep_data = True )



reload( gnuplot )
filename = 'q1_diff_distribution_1vs1' 
gnu  = gnuplot.Gnuplot.draw_boxes_from_data( [h,n], filename + '.eps',
  'Q1 difference / Th', 'Occurence', keep_data = True )


c.close()
c = db.cursor()
c.execute( query_down )
