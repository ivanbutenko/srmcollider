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



try:
    import c_getnonuis
    #collisions
    #q3, q1, srm_id, peptide_key
    #transitions
    #q3, srm_id
    #
    q3window = 1.0
    ppm = False
    transitions = ( (500.0,1), 
                    (600,2), 
                    (700,3), 
                  )
    collisions = (  (500.4,400.0,101,201),
                    (600.6,400.0,102,201), 
                    (500.6,401.0,103,202),
                    (700.6,401.0,104,202), 
                    (500.7,400.0,105,203),
                    (600.7,400.0,106,203),
                    (700.7,400.0,107,203), 
                 )
    result = c_getnonuis.core_non_unique( transitions, collisions, q3window, ppm)

    assert abs(result[1] - 0.4) < 10**(-3)
    assert abs(result[2] - 0.6) < 10**(-3)
    assert abs(result[3] - 0.6) < 10**(-3)
    #
    #
    result = c_getnonuis.getnonuis( transitions, collisions, q3window, ppm)
    assert result[201] == [1,2]
    assert result[202] == [1,3]
    assert result[203] == [1,2,3]
    #Test 2
    transitions = ( (500.0,1), 
                    (600,2), 
                    (700,3), 
                    (800,4), 
                  )
    collisions = (  (500.4,400.0,101,201),
                    (600.6,400.0,102,201), 
                    (700.6,401.0,103,201),
                    (600.6,401.0,104,202), 
                    (700.7,400.0,105,202),
                    (800.7,400.0,106,202),
                 )
    result = c_getnonuis.getnonuis( transitions, collisions, q3window, ppm)
    assert result[201] == [1,2,3]
    assert result[202] == [2,3,4]
except ImportError:
    print """Module c_getnonuis is not available.
    
    Please compile it if you want to use it."""


