#!/usr/bin/python
# -*- coding: utf-8  -*-
# vim:set fdm=marker:
# {{{ Main: Run the collider
import MySQLdb
import time
import sys 
sys.path.append( '/home/hroest/srm_clashes/code' )
sys.path.append( '/home/hroest/lib/' )
db = MySQLdb.connect(read_default_file="~/.my.cnf")
cursor = db.cursor()
import collider
import hlib
import copy
#Run the collider
###########################################################################
par = collider.SRM_parameters()
par.q1_window = 0.7 / 2.0 #UIS paper = 1.2
par.q3_window = 0.7 / 2.0 #UIS paper = 2.0
#12 = 90%TPR, 8 = 80%TPR 
par.ssrcalc_window = 4 / 2.0 #for paola do 9999 and 4
par.ssrcalc_window = 9999 / 2.0 #for paola do 9999 and 4
par.ppm = False
par.considerIsotopes = True
par.do_1vs = False
par.transition_table = 'hroest.srmTransitions_yeast_top3'
par.peptide_table = 'hroest.srmPeptides_yeast_top3'
#contains only the actual mrm transitions (and not all b/y ions)
par.transition_table = 'hroest.srmTransitions_yeast_mrmatlas'
par.peptide_table = 'hroest.srmPeptides_yeast_mrmatlas'
par.transition_table = 'hroest.srmTransitions_yeast_mrmatlasall'
par.peptide_table = 'hroest.srmPeptides_yeast_mrmatlasall'
par.transition_table = 'hroest.srmTransitions_yeast'
par.peptide_table = 'hroest.srmPeptides_yeast'
par.max_uis = 5 #turn off uis if set to 0
par.q3_range = [400, 1400]
par.eval()
print par.experiment_type
print par.get_common_filename()
par.query1_add += ' order by Intensity DESC' #only for the toptrans workflow

#paolas workflow
#we need a different background for the peptide atlas
#because there we have differnt peptids and cannot relate
#those with the MRMlink table
reload( collider )
mycollider = collider.SRMcollider()
mycollider.find_clashes_toptrans_paola(db, par) #, bgpar=parbg)

#if you need to have a different background than fg
#only for peptide atlas since those peptide keys map to experiment 3452 and not 
#to 3131 as we would need to then find the transition in the MRM table
parbg = copy.deepcopy(par)
parbg.transition_table = 'hroest.srmTransitions_yeast_pepatlas'
parbg.peptide_table = 'hroest.srmPeptides_yeast_pepatlas'

#calculate for precursor, peptide and protein the min nr transitions necessary
#to measure them without ambiguity
#create a new hroest.experiment
common_filename = par.get_common_filename()
query = """
insert into hroest.experiment  (name, short_description,
description, comment1, comment2, super_experiment_key, ddb_experiment_key)
VALUES (
    'paola 0.7Da/0.7Da', '%s', '%s', '%s', '%s', 3, 0
)
""" %( common_filename + '_' + par.peptide_table.split('.')[1], 
      par.experiment_type, par.peptide_table, par.transition_table)
cursor.execute(query)
myid = db.insert_id()
prepare = [ [p[0], p[1], myid]  for p in mycollider.min_transitions]

#save our result ("the min nr transitions per precursor") linked with 
#our new experiment key
cursor.executemany(
""" insert into hroest.result_srmpaola (parent_key, min_transitions, exp_key) 
    VALUES (%s,%s,%s) """, prepare
)

}}}
#
#this is how do the readout
select id, name, short_description from experiment;
select @exp_key := 18;
select min_transitions, count(*) from hroest.result_srmpaola r 
inner join srmPeptides_yeast srm on srm.parent_id = r.parent_key 
where exp_key = @exp_key 
group by min_transitions ; 

select mt, count(*) from
(
select min(min_transitions) as mt, peptide_key from hroest.result_srmpaola r 
inner join srmPeptides_yeast srm on srm.parent_id = r.parent_key 
where exp_key = @exp_key 
group by peptide_key
) tmp
group by mt ; 

select mt, count(*) from
(
select min(min_transitions) as mt, protein_key from hroest.result_srmpaola r 
inner join srmPeptides_yeast srm on srm.parent_id = r.parent_key 
inner join ddb.protPepLink l on l.peptide_key = srm.peptide_key
where exp_key = @exp_key 
group by protein_key
) tmp
group by mt ; 


#load from file
mycollider = collider.SRMcollider()
directory = '/home/hroest/srm_clashes/results/pedro/'
mycollider.load_from_file( par, directory)

#print results
mycollider.print_unique_histogram(par)
mycollider.print_cumm_unique(par)
mycollider.print_q3min(par)
mycollider.print_q3min_ppm(par)
mycollider.print_stats()
mycollider.print_cumm_unique_all(par, db.cursor() )

#
#
#
reload( collider )
mycollider = collider.SRMcollider()
#pepids = mycollider._get_unique_pepids(par, cursor, ignore_genomeoccurence=True)
pepids = mycollider._get_unique_pepids_toptransitions(par, cursor)

mycollider.find_clashes(db, par, pepids=pepids, toptrans=True, use_per_transition=False )
mycollider.find_clashes_small(db, par, pepids=pepids, UIS_only=False)
mycollider.find_clashes_toptrans_3strike(db, par, pepids=pepids, 
                                         use_per_transition=False )
#mycollider.store_in_file(par)





##Run the testcase => see test_run.py
###########################################################################
# Storage of results {{{
self = par

common_filename = par.get_common_filename()
restable = 'hroest.srm_results_' + common_filename
cursor.execute("drop table %s" % restable)
cursor.execute("create table %s (parent_key int, unique_transitions double) " % restable)
cursor.executemany( "insert into %s values " % restable + "(%s,%s)", 
              mycollider.allpeps.items() )

common_filename = par.get_common_filename()
restable = 'hroest.srmCollisions_' + common_filename
cursor.execute("drop table %s" % restable)
cursor.execute("create table %s (coll_srm1 int, coll_srm2 int) " % restable)
cursor.executemany( "insert into %s values " % restable + "(%s,%s)", 
              mycollider.allcollisions )

}}}
###########################################################################
# find 2 colliding peptides in the 1200 peptides by Pedro {{{

query = """
select parent_id, q1, q1_charge, ssrcalc
 from %s
 inner join
 ddb.peptide on peptide.id = %s.peptide_key
 %s
""" % (par.peptide_table, par.peptide_table, par.query_add )
cursor.execute( query )
mypepids =  cursor.fetchall()

mycollider = collider.SRMcollider()
mycollider.find_clashes(db, par, toptrans=False, pepids=mypepids)

[pep for pep, c in mycollider.allpeps.iteritems() if c < 0.7].__len__()

#all peptides with more than 3 collisions
interest2 = [pep for pep, c in mycollider.colliding_peptides.iteritems() if c > 3]
interest_pep2 = [pep for pep in pepids if pep['parent_id'] in interest2]

#all peptides with less than 70% unique transitions
interest = [pep for pep, c in mycollider.allpeps.iteritems() if c < 0.7]
interest = [pep for pep in interest if pep not in interest2]
interest_pep = [pep for pep in pepids if pep['parent_id'] in interest]


cursor.execute( """select peptide.id, Tr_original from ddb.peptide inner
               join hroest.yeast_dp_data_original on peptide.sequence =
               stripped_sequence where experiment_key = 3352
               group by sequence"""
)

cursor.execute( """select peptide.id, RT from ddb.peptide inner
               join hroest.ludovic_compiled2_mascot20101001 on peptide.sequence =
               pep_seq
               where experiment_key = 3352
               group by sequence"""
)
real_rts = cursor.fetchall()
real_rts = dict(real_rts )
interest = pepids
cursor.execute( """select sequence, RT from ddb.peptide inner
               join hroest.ludovic_compiled2_mascot20101001 on peptide.sequence =
               pep_seq
               where experiment_key = 3352
               group by sequence"""
)
realrt_seq = cursor.fetchall()
realrt_seq = dict(realrt_seq )



MAX_UIS = 10
self.UIS_redprob = [0 for i in range(MAX_UIS+1)]


pep = interest_pep[0]

manually_peps = """
DGPWDVMLK
DVYLLDLR
FDELIPSLK
FNYAWGLIK
FSEGLLSVIK
IQTYGETTAVDALQK
LGEELTAIAK
LLDTLGFQK
LLENGITFPK
NRPTLQVGDLVYAR
""".split()


out_csv = ""
par.ssrcalc_window = 100000
REALDIFF =  0.5
metatext = ''
count = 0
for pep in interest:
    transitions = self._get_all_transitions(par, pep, cursor)
    nr_transitions = len( transitions )
    collisions = self._get_all_collisions(par, pep, cursor)
    collisions_per_peptide = {}
    q3_window_used = par.q3_window
    for t in transitions:
        for c in collisions:
            if abs( t[0] - c[0] ) <= q3_window_used:
                #gets all collisions
                if collisions_per_peptide.has_key(c[3]):
                    if not t[1] in collisions_per_peptide[c[3]]:
                        collisions_per_peptide[c[3]].append( t[1] )
                else: collisions_per_peptide[c[3]] = [ t[1] ] 
    pairs = []
    for key in collisions_per_peptide.keys():
        for key2 in collisions_per_peptide.keys():
            try:
                realdiff = real_rts[ key ] - real_rts[ key2 ]
                if abs( realdiff ) < REALDIFF and key < key2:
                    print realdiff; count += 1
                    pairs.append( [key, key2])
            except:pass #we may not have rts for all of them
    print "="*75
    if len( pairs ) == 0: continue
    flatpairs = flatten_list( pairs )
    interest_list  =  dict( [ [p,collisions_per_peptide[p]] for p in flatpairs] )
    non_uids = get_all_non_uis(self, pep, par, cursor, 2)
    #collisions_per_peptide.values() 
    cursor.execute( "select sequence from ddb.peptide where id =%s" % pep['peptide_key'])
    seq = cursor.fetchall()[0][0]
    if seq not in manually_peps: continue
    real_rt = -1
    try:
        real_rt = real_rts[ pep['peptide_key']]
    except: pass
    text = """Peptide %s\nQ1: %s\tSSRCalc %s\tReal %s\nTo Check:\n""" % (seq,
        pep['q1'], pep['ssrcalc'], real_rt)
    colliding_seqs = [seq]
    coll_energy = calculate_coll_energy( 2, pep['q1'])
    for collpep, coll in interest_list.iteritems():
        cursor.execute( "select sequence from ddb.peptide where id =%s" % collpep )
        sequence = cursor.fetchall()[0][0]
        colliding_seqs.append( sequence )
        for cc in coll:
            cursor.execute( "select q3 from hroest.srmTransitions_yeast_1200 where srm_id =%s" % cc )
            q3 = cursor.fetchall()[0][0]
            text += "%s\t%s\t %s\t%s\n" % (pep['q1'], q3, sequence, realrt_seq[sequence] )
            out_csv += "%s,%s,%s,coll_%s\n" % (pep['q1'], q3, coll_energy, seq + "_" + sequence )
    text += "We cannot use these combinations at the same time (not UIS):\n" 
    for pair in non_uids:
        for coll in pair:
            cursor.execute( "select q3 from hroest.srmTransitions_yeast_1200 where srm_id =%s" % coll )
            text += "%s\t%s\t\t" % (pep['q1'], cursor.fetchall()[0][0])
        text += "\n"
    for s in colliding_seqs:
        #cursor.execute( """select Q1,Q3, stripped_sequence, isotype,
        #               relative_intensity, rank from
        #               hroest.yeast_dp_data_original where stripped_sequence =
        #               '%s' and isotype='light' order by prec_z, rank;
        #              """ % s)
        cursor.execute( """select Q1,Q3, naked, prec_z,
                       intensities from
                       hroest.ludovic_compiled2_daveDirtyPeP20101001 where naked =
                       '%s' order by prec_z, intensities;
                      """ % s)
        text += "Best transitions for %s\n" % s
        for tr in cursor.fetchall():
            text += "%s\t"*4 % (tr[0], tr[1], tr[2], tr[4])
            text += "\n"
            if abs(pep['q1'] - tr[0]) > 10: continue
            coll_energy = calculate_coll_energy( tr[3], tr[0])
            out_csv += "%s,%s,%s,%s\n" % (tr[0], tr[1], coll_energy, tr[2])
    #print text
    metatext += text + '\n' * 3




open('real_close.txt', 'w').write( metatext)
open('more3coll.txt', 'w').write( metatext)
open('less70pct_unique.txt', 'w').write( metatext)


open('20101004_hannes_srmcollider.csv', 'w').write( out_csv )
open('20101004_hannes_srmcollider_2.csv', 'w').write( out_csv )



collider.permutations


def calculate_coll_energy(charge, mz):
    #Fuer 2+
    # CE (2+)= 0.034*(m/z) +3.314
    #Fuer 3+
    # CE (3+)= 0.044*(m/z) +3.314
    if charge ==2: 
        return 0.034*(mz) +3.314
    elif charge ==3: 
        return 0.044*(mz) +3.314
    else: assert False

def flatten_list(l):
    return [item for sublist in l for item in sublist]






def get_all_non_uis(self, pep, par, cursor, order, ssrcalc_window = 100000):
    old_win = par.ssrcalc_window
    par.ssrcalc_window = ssrcalc_window
    transitions = self._get_all_transitions(par, pep, cursor)
    nr_transitions = len( transitions )
    collisions = self._get_all_collisions(par, pep, cursor)
    q3_window_used = par.q3_window
    collisions_per_peptide = {}
    for t in transitions:
        for c in collisions:
            if abs( t[0] - c[0] ) <= q3_window_used:
                #gets all collisions
                if collisions_per_peptide.has_key(c[3]):
                    if not t[1] in collisions_per_peptide[c[3]]:
                        collisions_per_peptide[c[3]].append( t[1] )
                else: collisions_per_peptide[c[3]] = [ t[1] ] 







#here we calculate the UIS for this peptide with the given RT-range
reslist = set()
for pepc in collisions_per_peptide.values():
    collider.get_non_uis(pepc, reslist, order)

    par.ssrcalc_window = old_win 
    return reslist




newl = []
for p, c in mycollider.count_pair_collisions.iteritems():
    parent, peptide = p.split(':')
    newl.append( [parent, c] )

newl.sort( lambda x,y: -cmp(x[1], y[1]) )
interesting = [n for n in newl if n[1] > 2]

def get_tbl(id):
    trpep = collider.get_trans_collisions(par, db, p_id=id, q3_low = 300, q3_high = 5000)
    latex_res  = ''
    self = trpep
    for tr in self.transitions:
        latex_res += tr.print_latex() + '\n'
    mytable = r"""
    \begin{table}[h]
    \centering
    \caption[%(c1)s]{\textbf{%(c1)s} %(c2)s}
    \label{tab:%(label)s}
    \begin{tabular}{ %(rows)s }
    %%\maketablespace
    %(tabletop)s \\
    %%\toprule
    %(content)s
    \end{tabular}
    \end{table}
    """ % { 'c1' : 'Analysis for peptide number ' + str(id) + '.', 'c2' : """ 
          The peptide with sequence \url{%s}, 
           a q1 mass of at %s and ssrcalc %s was analyzed.
           """ % (trpep.sequence, trpep.q1, trpep.ssrcalc) , 
           'label' : 'l' + str(id),
           'rows' : 'c c |' + 'c'*7,
           'tabletop' : 'q3 & ion & q3 & ion & q1 & SSRcal & dq3 /1000 & sequence & q3charge',
           'content' : latex_res
    }
    return mytable

#latex table
reload( collider )
mystr = ''
for i in interesting:
    mystr += get_tbl(int(i[0]) )
    mystr += '\n'

latex.create_document( mystr )


}}}
###########################################################################
## cumulative distr. of how much each ppm resolution contributes {{{


mycollider = collider.SRMcollider()
directory = '/home/hroest/srm_clashes/results/pedro/'
mycollider.load_from_file( par, directory)
self = mycollider
h, n = numpy.histogram( [abs(q) for q in self.q3min_distr_ppm], 100 )
#h = [gg*100.0 / len( self.q3min_distr_ppm) for gg in h]
cumh = collider.get_cum_dist( h )
filename = par.get_common_filename() + "_q3cum_abs"
plt_cmd = 'with lines lt -1 lw 2'
gnu = gnuplot.Gnuplot.draw_from_data( [cumh,n], plt_cmd, filename + '.eps',
  'Closest difference in Q3 / ppm' , 'Cumulative Occurence / %', 
  keep_data=True, tmp_csv = filename + '.csv')
gnu.add_to_body( "set yrange[0:100]" )
gnu.draw(plt_cmd, keep_data=True)



len( [q3 for q3 in self.q3min_distr_ppm if abs(q3) < 5 ] ) * 1.0 / len( self.q3min_distr_ppm )
}}}
###########################################################################
# number of shared collisions between two peptides {{{
more_one_coll = [ {} for i in range(10)]
for p, c in mycollider.count_pair_collisions.iteritems():
    try:
        parent, peptide = p.split(':')
        if c > 0: more_one_coll[0][ parent ] = ''
        if c > 1: more_one_coll[1][ parent ] = ''
        if c > 2: more_one_coll[2][ parent ] = ''
        if c > 3: more_one_coll[3][ parent ] = ''
        if c > 4: more_one_coll[4][ parent ] = ''
        if c > 5: more_one_coll[5][ parent ] = ''
        if c > 6: more_one_coll[6][ parent ] = ''
    except AttributeError: pass


more_one_coll_len = [ len(more) for more in more_one_coll ]
more_one_coll_len 

for myc in range(10):
    number  = prepare_int(more_one_coll_len[myc], 4 )
    print get_latex_row( [myc + 1, number ] )


total = len(mycollider.count_pair_collisions) - 134029
for myc in range(10):
    mylen = len([p for p, c in mycollider.count_pair_collisions.iteritems() if c == myc] )
    mylen_latex = prepare_int( mylen, 7) 
    prct_latex = prepare_float( mylen * 100.0 / total, 2, 3)
    print get_latex_row( [myc, mylen_latex, prct_latex ] )

}}}


