#!/usr/bin/python
# -*- coding: utf-8  -*-
"""
This shows how to run the window finder script
"""
import window_finder
import MySQLdb
db = MySQLdb.connect(read_default_file="~/.my.cnf")
cursor = db.cursor()

################
query = """
select q1
from %s
where
q1 between 400 and 1200 
and ssrcalc between 13 and 47
order by q1
"""

#select the N14 and N15 yeast, then combine 
query14 = query % ('hroest.srmPeptides_yeast')
query15 = query % ('hroest.srmPeptides_yeastN15')
cursor.execute( query14)
v14  = cursor.fetchall()
cursor.execute( query15)
v15  = cursor.fetchall()
joined = list( v14[:] )
joined.extend( v15 )
joined = [j[0] for j in joined]
joined.sort()
mylen = len( joined )

#make 32 bins because we want to sample every 3.2 seconds
mybins = 32
breakafter = numpy.round( mylen * 1.0 / mybins)
values = joined
result = []
for i,v in enumerate(values):
    if i % breakafter == 0: print v; result.append(v)


#naive approach: just sort by Q1 and then make bins 
#this is a good starting point
naive = [result[i+1] - result[i] for i in range(mybins-1)]
naive.append( 1200 - result[-1] )

#go with reduced complexity
reduced = [j for i,j in enumerate(joined) if i%10==0 ]


###########################################################################
###########################################################################
#here the optimization starts

#Strategy 1, just do iterations
w = window_finder.Window_finder(reduced, 32, 400, 1200, overlap=2)
w._set_start( naive )
for i in range(10): w.run(max_change_amount=1, 
                           iterations=1000, temperature=0.2); w.bins.bin_sizes



#Strategy 2, simulated annealing
w = window_finder.Window_finder(reduced, 32, 400, 1200, overlap=2)
w._set_start( naive )
myrange = [10-i for i in range(10)]
window_finder.simulated_annealing(start, w, max_change=1, myrange=myrange)


#Strategy 3, heated chains 
w = window_finder.Window_finder(reduced, 32, 400, 1200, overlap=2)
temperatures = [2.0*i + 0.1 for i in range(3)]
window_finder.start_chains(reduced, naive, temperatures=temperatures, parallel_runs=3)



#best result so far
best_result = [12.20920693318614, 11.021187825939242, 11.202295036902511, 10.825109640504003, 11.498587481932736, 11.858596936213758, 11.789169452746645, 12.171484840021296, 12.684925478005049, 13.310364942294392, 12.949397368590041, 13.253997139645222, 13.992100084656162, 14.919672773233698, 15.145799641393985, 15.967002382049406, 16.847879263750116, 17.674242439490914, 18.943198967946945, 20.059853489763441, 21.556045695687427, 23.398426574814479, 25.178317900770843, 28.415603754408107, 29.202005638995185, 32.4868795924992, 37.389983581921719, 42.481834924231492, 48.534920593374935, 59.790852138254806, 73.28470627319804, 99.813184546577773]

