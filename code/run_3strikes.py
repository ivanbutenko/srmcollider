#!/usr/bin/python
# -*- coding: utf-8  -*-
# vim:set fdm=marker:
"""
select count(distinct parent_key) from result_srmuis where exp_key = 104;
select count(*) from srmPeptides_mouse where isotope_nr = 0 and q1_charge = 2 and q1 between 400 and 1400;

python prepare_rangetree.py 104 60000 9999 hroest.srmPeptides_mouse
sh /tmp/tmp.sh

"""

import MySQLdb
import time
import sys 
import c_rangetree, c_getnonuis
sys.path.append( '/home/hroest/lib/hlib/' )
sys.path.append( '/home/hroest/projects' )
sys.path.append( '/home/hroest/projects/hlib' )
#db_l = MySQLdb.connect(read_default_file="~/.my.cnf.local")
db = MySQLdb.connect(read_default_file="~/.my.cnf")
cursor = db.cursor()
#cursor_l = db_l.cursor()
#cursor = cursor_l
import collider
import progress

exp_key = -1234
min_q1 = 450
max_q1 = 460
ssrcalcwin = 2.0

#11234556 1047 1205 0.25 hroest.srmPeptides_yeast



from optparse import OptionParser, OptionGroup
usage = "usage: %prog experiment_key startQ1 endQ1 [options]"
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "Run uis Options",
                    "None yet")
group.add_option("-f", "--outfile", dest="outfile", default='3strikes.out',
                  help="A different outfile " )
group.add_option("--ssr3strike", dest="ssr3strike", default=0.3, type='float',
                  help="SSRCalc for the 3rd strike" )
group.add_option("--order", dest="myorder", default=3, type='int',
                  help="Order to consider" )
group.add_option("--allow_contamination", dest="allow_contamination", default=3, type='int',
                  help="Number of locally contaminated transitions to be allowed" )
parser.add_option_group(group)



#Run the collider
###########################################################################
#Parse options
db = MySQLdb.connect(read_default_file="~/.my.cnf")
cursor = db.cursor()
par = collider.SRM_parameters()

par = collider.SRM_parameters()
par.parse_cmdl_args(parser)
options, args = parser.parse_args(sys.argv[1:])
par.parse_options(options)

#local arguments
exp_key = sys.argv[1]
min_q1 = float(sys.argv[2])
max_q1 = float(sys.argv[3])
outfile = options.outfile
strike3_ssrcalcwindow = options.ssr3strike
myorder =options.myorder
contamination_allow =options.allow_contamination
par.dontdo2p2f = False #do not look at 2+ parent / 2+ fragment ions

"""

par.q1_window = 0.7 / 2.0 #UIS paper = 1.2
par.q3_window = 0.7 / 2.0 #UIS paper = 2.0
par.ssrcalc_window = 4 / 2.0
par.ppm = False
outfile = '/tmp/out'
strike3_ssrcalcwindow = 0.3
myorder = 4

"""

#par.considerIsotopes = True
#par.max_uis = 5 #turn off uis if set to 0
#par.dontdo2p2f = False #also look at 2+ parent / 2+ fragment ions
#par.q3_range = [400, 1400]
par.eval()
#print par.experiment_type
#par.peptide_table = peptide_table
print par.get_common_filename()
mycollider = collider.SRMcollider()

import Residues
R = Residues.Residues('mono')
isotope_correction = par.isotopes_up_to * R.mass_diffC13 / min(par.parent_charges)
q =  """
select modified_sequence, transition_group, parent_id, q1_charge, q1, ssrcalc, modifications, missed_cleavages
from %(peptide_table)s where q1 between %(lowq1)s - %(isotope_correction)s and %(highq1)s
""" % {'peptide_table' : par.peptide_table, 
              'lowq1'  : min_q1 - par.q1_window, 
              'highq1' : max_q1 + par.q1_window,
              'isotope_correction' : isotope_correction
      } 
cursor.execute(q)


#print "finished query execution: ", time.time() - start
alltuples =  list(cursor.fetchall() )
#from random import shuffle
#shuffle(alltuples)

#print "finished query fetch: ", time.time() - start
mypepids = [
            {
                'sequence'  :  r[0],
                'transition_group' :r[1],
                'parent_id' :  r[2],
                'q1_charge' :  r[3],
                'q1' :         r[4],
                'ssrcalc' :    r[5],
            }
            for r in alltuples
    if r[3] == 2 #charge is 2
    and r[6] == 0 # no modification 
    and r[7] == 0 # no missed cleavages 
    and r[4] >= min_q1
    and r[4] < max_q1
]
parentid_lookup = [ [ r[2], (r[4], r[0], r[1]) ] 
            for r in alltuples
    #if r[3] == 2 and r[6] == 0
]
parentid_lookup  = dict(parentid_lookup)
#print "finished python lookups: ", time.time() - start
pepkey_lookup = [ [ r[1], r[5] ] for r in alltuples ]
pepkey_lookup  = dict(pepkey_lookup)

print "building tree with %s Nodes" % len(alltuples)
c_rangetree.create_tree(tuple(alltuples))
#print "finished build tree : ", time.time() - start


"""
The idea is to find UIS combinations that are globally UIS, locally
clean and also there are no cases in which all the UIS coelute when
they are in different peptides:
    * global UIS = whole RT
    * locally clean = no intereferences around the peptide
    * no coelution: find all peptides globally that share transitions
This is a 3 strikes rule to find good UIS combinations.
"""
#make sure we only get unique peptides
VERY_LARGE_SSR_WINDOW = 9999999
myssrcalc = par.ssrcalc_window

def thisthirdstrike( N, ssrcalcvalues, strike3_ssrcalcwindow):
    M = len(ssrcalcvalues)
    index = [0 for i in range(M)]
    myssr = [0 for i in range(M)]
    discarded_indices = []
    all_nonuis = {}
    # check whether there are any empty ssrcalcvalues
    for k in range(M):
        if len(ssrcalcvalues[k]) == 0:
            discarded_indices.append(k)
    if len(discarded_indices) == M: return {}

    # in each iteration we we either advance one or add one to discarded
    # thus we do this at most sum(N) + len(N) times which is bounded by
    # c*k + c
    while True:
        myssr = [-1 for i in range(M)]
        for k in range(M):
            if k in discarded_indices: continue
            myssr[k] = ssrcalcvalues[k][ index[k] ]

        # find the pivot element (the smallest element that is not yet in a discarded group)
        tmin = max(myssr); piv_i = -1
        for k in range(M):
            if k in discarded_indices: continue
            if myssr[k] <= tmin:
                tmin=myssr[k]
                piv_i = k

        # we need to sort by we also need to have a map back to retrieve the original!
        with_in = [ (a,b) for a,b in enumerate(myssr) ]
        with_in.sort( lambda x,y: cmp(x[1],y[1]))
        sort_idx = [x[0] for x in with_in]
        myssr.sort()
        #if index == [0, 0, 1, 0, 1, 19]: print index, piv_i, tmin, "== ", myssr, sort_idx

        # now find all N different combinations that are not UIS. Since they are
        # sorted we only need to consider elements that have a higher index.
        for k in range(M):
          # or myssr[k] == -1 would also work here
          if sort_idx[k] in discarded_indices: continue 
          nonuis = [k]
          for m in range(k+1,M):
              if not m == k and not (abs(myssr[k] - myssr[m]) > strike3_ssrcalcwindow):
                  nonuis.append(m)
          backsorted = [sort_idx[n] for n in nonuis]
          backsorted.sort()
          all_nonuis[ tuple(backsorted) ] = 0
          #if len(nonuis) > 1: print "added tuple ", tuple(backsorted)
                    
        # Advance the pivot element
        index[piv_i] += 1
        if(index[piv_i] >= len(ssrcalcvalues[piv_i])):
            discarded_indices.append(piv_i)
            #print "wanted to advance", piv_i, "had to append to discarded ", discarded_indices
            if len(discarded_indices) == len(ssrcalcvalues): 
                # break out of loop
                break
    return all_nonuis


print par.experiment_type
self = mycollider
self.mysqlnewtime = 0
self.pepids = mypepids
MAX_UIS = par.max_uis
progressm = progress.ProgressMeter(total=len(self.pepids), unit='peptides')
prepare  = []
at_least_one = 0
f = open(outfile, 'a')
for kk, pep in enumerate(self.pepids):

    p_id = pep['parent_id']
    q1 = pep['q1']
    ssrcalc = pep['ssrcalc']
    q3_low, q3_high = par.get_q3range_transitions()
    #new way to calculate the precursors
    transitions = c_getnonuis.calculate_transitions_ch(
        ((q1, pep['sequence'], p_id),), [1], q3_low, q3_high)
    #fake some srm_id for the transitions
    transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
    nr_transitions = len(transitions)
    q1_low = q1 - par.q1_window 
    q1_high = q1 + par.q1_window
    #correct rounding errors, s.t. we get the same results as before!
    ssrcalc_low = ssrcalc - par.ssrcalc_window + 0.001
    ssrcalc_high = ssrcalc + par.ssrcalc_window - 0.001

    ###############################################################
    #strike 1: it has to be global UIS
    # ssrcalc_low = -999
    # ssrcalc_high = 999
    precursor_ids = tuple(c_rangetree.query_tree( q1_low, ssrcalc_low, 
                                                 q1_high,  ssrcalc_high )  )
    globalprecursors = tuple([parentid_lookup[myid[0]] for myid in precursor_ids
                        #dont select myself 
                       if parentid_lookup[myid[0]][2]  != pep['transition_group']])
    collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide( 
        transitions, globalprecursors, q3_low, q3_high, par.q3_window, par.ppm)
    tmpnonlist = c_getnonuis.get_non_uis( collisions_per_peptide, myorder)
    srm_ids = [t[1] for t in transitions]
    tuples_1strike = collider.get_uis(srm_ids, tmpnonlist, myorder)

    ###############################################################
    #strike 2: it has to be locally clean
    ssrcalc_low = ssrcalc - par.ssrcalc_window + 0.001
    ssrcalc_high = ssrcalc + par.ssrcalc_window - 0.001
    precursor_ids = tuple(c_rangetree.query_tree( q1_low, ssrcalc_low, 
                                                 q1_high,  ssrcalc_high )  )
    precursors = tuple([parentid_lookup[myid[0]] for myid in precursor_ids
                        #dont select myself 
                       if parentid_lookup[myid[0]][2]  != pep['transition_group']])
    collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide( 
        transitions, precursors, q3_low, q3_high, par.q3_window, par.ppm)
    local_interferences = [t[0] for t in c_getnonuis.get_non_uis( collisions_per_peptide, 1).keys()]
    tuples_2strike = []
    for mytuple in tuples_1strike:
        contaminated = 0.0
        for dirty_t in local_interferences:
            if dirty_t in mytuple: contaminated += 1.0
        if contaminated <= contamination_allow: tuples_2strike.append(mytuple)
    #tuples_2strike = tuples_1strike

    ###############################################################
    #strike 3: the transitions in the tuple shall not coelute elsewhere
    # Strike 3 takes usually over 90 % of the time
    tuples_3strike = []
    # 1. For each tuple that we found, find all peptides that collide with
    #    each single transition and store their SSRCalc values. For a tuple of
    #    n transitions, we get n vectors with a number of SSRCalc values.
    # 2. Check whether there exists one retention time (SSRcalc value) that is
    #    present in all vectors.
    ssrcalcvalues_dict = {}
    for t in transitions:
        collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide( 
            (t,), globalprecursors, q3_low, q3_high, par.q3_window, par.ppm)
        svals = [pepkey_lookup[ collkey] for collkey in collisions_per_peptide] 
        svals.sort()
        ssrcalcvalues_dict[t[1]] = svals
    for mytuple in tuples_2strike: 
        thistransitions = [ t for t in transitions if t[1] in mytuple]
        ssrcalcvalues = []
        for t in thistransitions:
            ssrcalcvalues.append(  ssrcalcvalues_dict [t[1]] )
        N = [len(v) for v in ssrcalcvalues]
        # if one of the transitions is contamination-free, the whole set is ok
        if min(N) == 0: tuples_3strike.append( mytuple ); continue
        # This usually takes 90 % of the time of strike 3
        #contaminated = c_getnonuis.thirdstrike( N, ssrcalcvalues, strike3_ssrcalcwindow)
        contaminated = c_getnonuis.thirdstrike_sort( N, ssrcalcvalues, strike3_ssrcalcwindow)
        if not contaminated: tuples_3strike.append( mytuple )

    old_len_tuples = len(tuples_3strike)

    ssrcalcvalues = ssrcalcvalues_dict.values()
    N = [len(v) for v in ssrcalcvalues]
    contaminated = thisthirdstrike(N, ssrcalcvalues, strike3_ssrcalcwindow)
    newdic = dict([(i,list(v)) for i,v in enumerate(contaminated.keys())])
    tmpnonlist.update(c_getnonuis.get_non_uis( newdic, myorder))
    tuples_3strike = collider.get_uis(srm_ids, tmpnonlist, myorder)

    if len(tuples_3strike) > 0: at_least_one += 1
    prepare.append( [ len(tuples_3strike), 
        collider.choose(nr_transitions, min(myorder, nr_transitions)), len(tuples_2strike)-len(tuples_3strike) ] )
    progressm.update(1)

print "Analysed:", kk
print "At least one eUIS of order %s :" % myorder, at_least_one, " which is %s %%" % (at_least_one *100.0/kk)
f.write( "At least one eUIS of order %s : %s which is %s %%\n" % (myorder, at_least_one, at_least_one *100.0/kk) )
f.write( "\n\nt.append([%s,%s])\n\n" % (kk, at_least_one) )

sum_all = sum([p[0]*1.0/p[1] for p in prepare]) 
nr_peptides = len([p for p in prepare])
print  "Random probability to choose good", sum_all*1.0/nr_peptides
sum_all = sum([p[2]*1.0/p[1] for p in prepare]) 
nr_peptides = len([p for p in prepare])
print  "Average lost in strike 3", sum_all*1.0/nr_peptides

f.close()



