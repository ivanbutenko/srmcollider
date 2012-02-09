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
from collider import thisthirdstrike

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

# Run the collider
###########################################################################
# Parse options
par = collider.SRM_parameters()
par.parse_cmdl_args(parser)
options, args = parser.parse_args(sys.argv[1:])
par.parse_options(options)

db = MySQLdb.connect(read_default_file=par.mysql_config)
cursor = db.cursor()
# local arguments
exp_key = sys.argv[1]
min_q1 = float(sys.argv[2])
max_q1 = float(sys.argv[3])
outfile = options.outfile
strike3_ssrcalcwindow = options.ssr3strike
myorder =options.myorder
contamination_allow =options.allow_contamination
par.dontdo2p2f = False #do not look at 2+ parent / 2+ fragment ions
par.eval()
print par.get_common_filename()

# Get the precursors
###########################################################################
from precursor import Precursors
myprecursors = Precursors()
myprecursors.getFromDB(par, db.cursor(), min_q1 - par.q1_window, max_q1 + par.q1_window)
testrange = myprecursors.build_rangetree()

precursors_to_evaluate = [p for p in myprecursors.precursors 
                         if p.q1_charge == 2 
                         and p.modifications == 0
                         and p.missed_cleavages == 0 
                         and p.q1 >= min_q1
                         and p.q1 <= max_q1
                         ]

myprecursors.build_parent_id_lookup()
myprecursors.build_transition_group_lookup()

print "Want to evaluate precursors", len(precursors_to_evaluate)

"""
The idea is to find UIS combinations that are globally UIS, locally
clean and also there are no cases in which all the UIS coelute when
they are in different peptides:
    * global UIS = whole RT
    * locally clean = no intereferences around the peptide
    * no coelution: find all peptides globally that share transitions
This is a 3 strikes rule to find good UIS combinations.
"""

print par.experiment_type
progressm = progress.ProgressMeter(total=len(precursors_to_evaluate), unit='peptides')
prepare  = []
at_least_one = 0
f = open(outfile, 'a')
kk = len(precursors_to_evaluate) -1
for precursor in precursors_to_evaluate:
    q3_low, q3_high = par.get_q3range_transitions()
    transitions = precursor.calculate_transitions(q3_low, q3_high)
    nr_transitions = len(transitions)

    ###############################################################
    #strike 1: it has to be global UIS

    computed_collisions = myprecursors.get_collisions_per_peptide_from_rangetree(precursor, transitions, par)
    collisions_per_peptide = computed_collisions 

    non_useable_combinations = c_getnonuis.get_non_uis( collisions_per_peptide, myorder)
    srm_ids = [t[1] for t in transitions]
    tuples_strike1 = 0
    if not nr_transitions < myorder:
      tuples_strike1 = collider.choose(nr_transitions, myorder ) - len(non_useable_combinations)

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
    # 1. For each tuple that we found, find all peptides that collide with
    #    each single transition and store their SSRCalc values. For a tuple of
    #    n transitions, we get n vectors with a number of SSRCalc values.
    # 2. Check whether there exists one retention time (SSRcalc value) that is
    #    present in all vectors.

    # complexity: nr_tr * k 
    # we get the list of ssrcalc-values for each transition 
    ssrcalcvalues = [ [] for t in transitions]
    for k,v in collisions_per_peptide.iteritems():
        #ssrcalc = trgroup_lookup[k]
        ssrcalc = myprecursors.loopup_by_transition_group(k).ssrcalc
        for tr in v:
            ssrcalcvalues[tr].append(ssrcalc)

    # then sort the ssrcalc values 
    # complexity: nr_tr * k * log(k) 
    for i in range(len(ssrcalcvalues)):
        ssrcalcvalues[i] = sorted(ssrcalcvalues[i])

    # get the list of combinations that are not eUIS
    N = [len(v) for v in ssrcalcvalues]
    cont_comb_list = c_getnonuis.calculate_eUIS(N, ssrcalcvalues, strike3_ssrcalcwindow)
    
    if False:
        #compare c++/python code
        compare = thisthirdstrike(N, ssrcalcvalues, strike3_ssrcalcwindow)
        cont3 = {}
        for c in cont_comb_list:
            cont3[ tuple(c)] = 0
        assert compare == cont3

    # the list of combinations has to be expanded for a specific order (e.g.
    # for order 2, a combination of (1,2,3) has to be expanded into
    # (1,2),(1,3),(2,3) to get all subcombinations.
    newdic = dict([(i,list(v)) for i,v in enumerate(cont_comb_list)])
    non_useable_combinations.update(c_getnonuis.get_non_uis( newdic, myorder))

    tuples_strike3 = 0
    if not nr_transitions < myorder:
      # We are mostly interested in how many tuples are left after strike 3
      tuples_strike3 = collider.choose(nr_transitions, myorder ) - len(non_useable_combinations)
      prepare.append( [ tuples_strike3, collider.choose(nr_transitions, 
        min(myorder, nr_transitions)), tuples_strike1-tuples_strike3 ] )
    # If we have at least one tuple left
    if tuples_strike3 > 0: at_least_one += 1
    progressm.update(1)

print "Analysed:", kk + 1
print "At least one eUIS of order %s :" % myorder, at_least_one, " which is %s %%" % (at_least_one *100.0/(kk+1))
f.write( "At least one eUIS of order %s : %s which is %s %%\n" % (myorder, at_least_one, at_least_one *100.0/(kk+1)) )
f.write( "\n\nt.append([%s,%s])\n\n" % ((kk+1), at_least_one) )

sum_all = sum([p[0]*1.0/p[1] for p in prepare]) 
nr_peptides = len([p for p in prepare])
print  "Random probability to choose good", sum_all*1.0/nr_peptides
sum_all = sum([p[2]*1.0/p[1] for p in prepare]) 
print  "Average lost in strike 3", sum_all*1.0/nr_peptides
sum_all = sum([(p[0]+p[2])*1.0/p[1] for p in prepare]) 
print  "Average without strike 3", 1-sum_all*1.0/nr_peptides



