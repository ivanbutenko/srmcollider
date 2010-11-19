#!/usr/bin/python
# -*- coding: utf-8  -*-
import MySQLdb
import time
import sys 
sys.path.append( '/home/hroest/srm_clashes/code' )
sys.path.append( '/home/hroest/lib/' )
db = MySQLdb.connect(read_default_file="~/.my.cnf")
cursor = db.cursor()
import collider
import hlib

##Run the testcase
###########################################################################
reload( collider )
par  = collider.testcase()
mycollider = collider.SRMcollider()
mycollider.find_clashes_small(db, par) 
print "I ran the testcase in %ss" % mycollider.total_time
assert abs(sum( mycollider.allpeps.values() ) - 975.6326447245566) < 10**(-3)
assert mycollider.non_unique_count == 26
assert mycollider.total_count == 12502
assert abs(mycollider.allpeps[1585] - 0.93333) < 10**(-3)

#verify that with toptrans=False we get the same results
par  = collider.testcase()
mycollider = collider.SRMcollider()
mycollider.find_clashes(db, par, toptrans=False) 
print "I ran the testcase in %ss" % mycollider.total_time
assert abs(sum( mycollider.allpeps.values() ) - 975.6326447245566) < 10**(-3)
assert mycollider.non_unique_count == 26
assert mycollider.total_count == 12502
assert abs(mycollider.allpeps[1585] - 0.93333) < 10**(-3)
assert len( mycollider.q3min_distr ) == 26
assert len( mycollider.q1min_distr ) == 26
assert len( mycollider.found3good ) == 978

#now test with isotopes enabled
par.considerIsotopes = True
par.eval()
mycollider = collider.SRMcollider()
#mycollider.find_clashes_small(db, par) 
mycollider.find_clashes(db, par, toptrans=False) 
print "I ran the testcase in %ss" % mycollider.total_time
assert abs(sum( mycollider.allpeps.values() ) - 971.2154365242601) < 10**(-3)
assert mycollider.non_unique_count == 71
assert mycollider.total_count == 12502
assert abs(mycollider.allpeps[1585] - 0.93333) < 10**(-3)
assert len( mycollider.q3min_distr ) == 71
assert len( mycollider.q1min_distr ) == 71
assert len( mycollider.found3good ) == 978
