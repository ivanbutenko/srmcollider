#!/usr/bin/python
# -*- coding: utf-8  -*-
# vim:set fdm=marker:
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

# {{{ Main: Run the collider
###########################################################################
par = collider.SRM_parameters()
par.q1_window = 1600 / 2.0 #UIS paper = 1.2 => we do 0.7 and 1.0
par.q3_window = 10.0 / 2.0 #UIS paper = 2.0
par.q1_window = 1.2 / 2.0 #UIS paper = 1.2 => we do 0.7 and 1.0
par.q3_window = 2.0 / 2.0 #UIS paper = 2.0
##12 = 90%TPR, 8 = 80%TPR 
par.ssrcalc_window = 4 / 2.0 #for paola do 9999 and 4
#par.ssrcalc_window = 9999 / 2.0 #for paola do 9999 and 4
#par.ssrcalc_window = 2 / 2.0 #for pedro
par.ppm = False
par.considerIsotopes = True
par.do_1vs = False
par.dontdo2p2f = False #do not look at 2+ parent / 2+ fragment ions
par.transition_table = 'hroest.srmTransitions_yeast_top3'
par.peptide_table = 'hroest.srmPeptides_yeast_top3'
#contains only the actual mrm transitions (and not all b/y ions)
#par.transition_table = 'hroest.srmTransitions_yeast_mrmatlas'
#par.peptide_table = 'hroest.srmPeptides_yeast_mrmatlas'
par.transition_table = 'hroest.srmTransitions_yeast_mrmatlasall'
par.peptide_table = 'hroest.srmPeptides_yeast_mrmatlasall'
par.transition_table = 'hroest.srmTransitions_yeast'
par.peptide_table = 'hroest.srmPeptides_yeast'
#par.transition_table = 'hroest.srmTransitions_yeast_N14N15'
#par.peptide_table = 'hroest.srmPeptides_yeast_N14N15'
#par.transition_table = 'hroest.srmTransitions_yeast_pepatlas_union'
#par.peptide_table = 'hroest.srmPeptides_yeast_pepatlas_union'
##
#par.transition_table = 'hroest.srmTransitions_human'
#par.peptide_table = 'hroest.srmPeptides_human'
#par.transition_table = 'hroest.srmTransitions_mouse'
#par.peptide_table = 'hroest.srmPeptides_mouse'
par.max_uis = 10 #turn off uis if set to 0
par.q3_range = [400, 1200]
par.eval()
print par.experiment_type
print par.get_common_filename()

#{{{ paolas workflow
#we need a different background for the peptide atlas
#because there we have differnt peptide_keys and cannot relate
#those with the MRMlink table
#only for peptide atlas since those peptide keys map to experiment 3452 and not
#to 3131 as we would need to then find the transition in the MRM table
#we only use the parbg for the calculations of the transitions. We still want 
parbg = copy.deepcopy(par)
parbg.transition_table = 'hroest.srmTransitions_yeast_pepatlas'
parbg.peptide_table = 'hroest.srmPeptides_yeast_pepatlas'

par.query1_add += ' order by Intensity DESC' #only for the toptrans workflow

reload( collider )
mycollider = collider.SRMcollider()

mycollider.find_clashes_toptrans_paola(db, par)#, bgpar=parbg)


#calculate for precursor, peptide and protein the min nr transitions necessary
#to measure them without ambiguity
#create a new hroest.experiment

common_filename = par.get_common_filename()
query = """
insert into hroest.experiment  (name, short_description,
description, comment1, comment2, comment3, super_experiment_key, ddb_experiment_key)
VALUES (
    'paola 1Da/1Da', '%s', '%s', '%s', '%s', '%s', 3, 0
)
""" %( common_filename + '_' + par.peptide_table.split('.')[1], #parbg.peptide_table.split('.')[1], 
      par.experiment_type, par.peptide_table, par.transition_table, 
     'using new MRMAtlas_qtrap_final_no_pyroGlu')
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


{{{ Print plot / result

{{{ sherman result 1
select @total ;


select @total := count(distinct gene_key) 
from  ddb.peptide p 
inner join ddb.peptideOrganism o on p.id = o.peptide_key
inner join ddb.protPepLink l on l.peptide_key = p.id
inner join ddb.protein  on l.protein_key = protein.id
#inner join ddb.geneProtLink ll on ll.protein_key = l.protein_key
where p.experiment_key = 3130
and genome_occurence = 1;



select count(distinct gene_key)   / @total 
from result_srmuis 
inner join srmPeptides_human pep on pep.parent_id = parent_key
inner join ddb.peptide p on p.id = pep.peptide_key
inner join ddb.peptideOrganism o on p.id = o.peptide_key
inner join ddb.protPepLink l on p.id = l.peptide_key
inner join ddb.geneProtLink ll on ll.protein_key = l.protein_key
where exp_key = 35
and uisorder = 2
and genome_occurence = 1
and experiment_key = 3130
and total_UIS - non_useable_UIS  > 1
and q1 between 400 and 2000
;



}}}
    #{{{Sherman fig 1

    cursor.execute(' drop table hroest.srmuis_tmp32 ')
    cursor.execute(
    """
    create temporary table hroest.srmuis_tmp32 as 
    #select sum(non_useable_UIS) sum_non, sum(total_UIS) sum_total , uisorder, exp_key
    select @c := count(*) from hroest.result_srmuis     where exp_key = 35 and uisorder =1;

    select sum(non_useable_UIS) sum_non, sum(total_UIS) sum_total , uisorder, exp_key,
    sum(non_useable_UIS) / sum(total_UIS) sum_total ,
    count(*),
    sum(non_useable_UIS / total_UIS) / count(*) as sumperpep
    from hroest.result_srmuis 
    where exp_key = 132
    group by uisorder, exp_key
    ;

    select exp_key
    sum(non_useable_UIS) / sum(total_UIS) sum_total ,
    count(*),
    sum(non_useable_UIS / total_UIS) / count(*) as sumperpep
    from hroest.result_srmuis 




    where exp_key = 132
    group by uisorder, exp_key
    ;

    """)

    cursor.execute( """
    select sum_non / sum_total as pcnt_non, uisorder, exp_key 
    from hroest.srmuis_tmp32 
    order by exp_key, uisorder
    """ )
    result = cursor.fetchall()


select sequence, molecular_weight, pi from ddb.peptide p
inner join ddb.peptideOrganism o on o.peptide_key = p.id
#our strain H73Rv
where experiment_key = 3474
and genome_occurence = 1
#not in human
and sequence not in
(select sequence from ddb.peptide where experiment_key = 3130)
#but in the other strain
and sequence in
(select sequence from ddb.peptide where experiment_key = 3475)

select count(*)
from ddb.peptide p
inner join ddb.peptideOrganism o on o.peptide_key = p.id
#our strain H73Rv
where experiment_key = 3474
and length(sequence) > 5 
and length(sequence) < 21 
and genome_occurence = 1
#not in human
and sequence not in
(select sequence from ddb.peptide where experiment_key = 3130)
#but in the other strain
and sequence in
(select sequence from ddb.peptide where experiment_key = 3475)


    """
    select sum(non_useable_UIS) sum_non, sum(total_UIS) sum_total , uisorder, exp_key,
    sum(non_useable_UIS) / sum(total_UIS) sum_total ,
    count(*),
    sum(non_useable_UIS / total_UIS) / count(*) as sumperpep
    from hroest.result_srmuis 
    where exp_key = 132
    and uisorder = 1
    group by uisorder, exp_key
    ;
    """

    #this is correct for Figure 1!!
    q = """
    select sum(non_useable_UIS) sum_non, sum(total_UIS) sum_total , uisorder, exp_key,
    sum(non_useable_UIS) / sum(total_UIS) sum_total ,
    count(*),
    sum(non_useable_UIS / total_UIS) / count(*) as sumperpep
    from hroest.result_srmuis 
    where exp_key = 132
    group by uisorder, exp_key
    ;
    """
    cursor.execute(q)

    result = cursor.fetchall()

    data = []
    e = [float(r[0])*100 for r in result if r[2] == 32]
    n = [r[1] for r in result if r[2] == 32]
    data.append( [h,n])
    h = [float(r[0])*100 for r in result if r[2] == 35]
    n = [r[1] for r in result if r[2] == 35]
    data.append( [h,n])
    """
    h = [float(r[0])*100 for r in result if r[2] == 105]
    n = [r[1] for r in result if r[2] == 105]
    data.append( [h,n])
    """
    h = [float(r[0])*100 for r in result if r[2] == 52]
    n = [r[1] for r in result if r[2] == 52]
    data.append( [h,n])
    h = [float(r[0])*100 for r in result if r[2] == 51]
    n = [r[1] for r in result if r[2] == 51]
    data.append( [h,n])
    """
    h = [float(r[0])*100 for r in result if r[2] == 104]
    n = [r[1] for r in result if r[2] == 51]
    data.append( [h,n])
    """

    data = []
    for ekey in range(132,137):
        h = [float(r[0])*100 for r in result if r[2] == ekey]
        n = [r[1] for r in result if r[2] == ekey]
        data.append( [h,n])

    import gnuplot
    titles = [
        'Human (no RT)',
        'Human (16 units)',
        'Human (12 units)',
        'Human (8 units)',
        'Human (4 units)']
    output = 'out.eps'
    xlabel = 'Number of transitions'
    ylabel = 'Probability of non-unique transition set / %'
    title = 'Probability to select a redundant SRM assay (transitions between 400 and 1400)'
    cmd = 'with lines lw 6'
    ba = """
        set yrange[0:100]
        set xtics 1
        set mxtics 1
        set key top right
        """
    gnu = gnuplot.Gnuplot.draw_from_multiple_data(data, cmd,
            '/tmp/fig1.eps', xlabel, ylabel, '/tmp/htmp', title, 
            keep_data = True, datatitles=titles, body_add = ba, nocolor=True, 
                                                      addls = False)
    #}}}


#{{{ Roest paper figure 2 (sherman fig 1)

    # there are two ways to calculate it:
    # 1. sum up all the nonuseable and total and then make the average
    # 2. create the average per peptide and then sum up the averages
    #
    # Here we chose the second way
    cursor.execute(' drop table hroest.srmuis_tmp32 ')
    cursor.execute(
    """
    create temporary table hroest.srmuis_tmp32 as 
    #select sum(non_useable_UIS) sum_non, sum(total_UIS) sum_total , uisorder, exp_key
    select sum(non_useable_UIS/total_UIS) / count(*)  sum_non, uisorder, exp_key
    from hroest.result_srmuis 
    where exp_key in (32, 132, 110)
    group by uisorder, exp_key;
    """)
    cursor.execute( """
    #select sum_non / sum_total as pcnt_non, uisorder, exp_key 
    select sum_non as pcnt_non, uisorder, exp_key 
    from hroest.srmuis_tmp32 
    order by exp_key, uisorder
    """ )
    result = cursor.fetchall()
    # 35 was the first human attempt but has no limit of 400-1400 init
    # 132 has that limit for human
    ## for yeast we dont have that data

    data = []
    ## h = [float(r[0])*100 for r in result if r[2] == 32]
    ## n = [r[1] for r in result if r[2] == 32]
    h = [float(r[0])*100 for r in result if r[2] == 132]
    n = [r[1] for r in result if r[2] == 132]
    data.append( [h,n])
    h = [ 0.99999598 , 0.90668559 , 0.36609717 , 0.10413078 , 0.03534957 ]
    h = [hh*100 for hh in h]
    n = [1,2,3,4,5]
    data.append( [h,n])
    # h = [float(r[0])*100 for r in result if r[2] == 110]
    # n = [r[1] for r in result if r[2] == 110]
    # data.append( [h,n])

    import gnuplot
    titles = [
        'Human', 
        'Yeast',
        #'Mouse'
]
    output = '/tmp/fig1.eps'
    xlabel = 'Number of transitions'
    ylabel = 'Probability of non-unique transition set / %'
    title = ''
    cmd = 'with lines lw 6'
    ba = """
        set yrange[0:100]
        set xtics 1
        set mxtics 1
        set key top right
        """
    gnu = gnuplot.Gnuplot.draw_from_multiple_data(data, cmd,
            output, xlabel, ylabel, '/tmp/htmp', title, 
            keep_data = True, datatitles=titles, body_add = ba, nocolor=True, 
                                                      addls = False)
    #}}}

#{{{ Roest paper figure 3 (rt influence)

    # there are two ways to calculate it:
    # 1. sum up all the nonuseable and total and then make the average
    # 2. create the average per peptide and then sum up the averages
    #
    # Here we chose the second way
    cursor.execute(' drop table hroest.srmuis_tmp32 ')
    cursor.execute(
    """
    create temporary table hroest.srmuis_tmp32 as 
    #select sum(non_useable_UIS) sum_non, sum(total_UIS) sum_total , uisorder, exp_key
    select sum(non_useable_UIS/total_UIS) / count(*)  sum_non, uisorder, exp_key
    from hroest.result_srmuis 
    where exp_key between 132 and 137
    group by uisorder, exp_key;
    """)
    cursor.execute( """
    #select sum_non / sum_total as pcnt_non, uisorder, exp_key 
    select sum_non as pcnt_non, uisorder, exp_key 
    from hroest.srmuis_tmp32 
    order by exp_key, uisorder
    """ )
    result = cursor.fetchall()

    data = []
    for ekey in range(132,137):
        h = [float(r[0])*100 for r in result if r[2] == ekey]
        n = [r[1] for r in result if r[2] == ekey]
        data.append( [h,n])

    import gnuplot
    titles = [
        'Human (no RT)',
        'Human (16 units)',
        'Human (12 units)',
        'Human (8 units)',
        'Human (4 units)']
    output = 'out.eps'
    xlabel = 'Number of transitions'
    ylabel = 'Probability of non-unique transition set / %'
    #title = 'Probability to select a redundant SRM assay (transitions between 400 and 1400)'
    title = ''
    cmd = 'with lines lw 6'
    ba = """
        set yrange[0:100]
        set xtics 1
        set mxtics 1
        set key top right
        """
    gnu = gnuplot.Gnuplot.draw_from_multiple_data(data, cmd,
            '/tmp/fig1.eps', xlabel, ylabel, '/tmp/htmp', title, 
            keep_data = True, datatitles=titles, body_add = ba, nocolor=True, 
                                                      addls = False)





    #}}}


    {{{hroest style figure (cumulative plot )


#regular 
select count(*) / 121009
, exp_key from hroest.result_srmuis
where exp_key between 238 and 261
and exp_key in (238, 239, 241, 245, 249, 250, 251, 252, 253, 254)
and total_UIS - non_useable_UIS >= 5
group by exp_key
order by exp_key;


#peptide atlas
select count(*) / 53002
, exp_key from hroest.result_srmuis
where exp_key between 238 and 261
and exp_key not in (238, 239, 240, 241, 245, 249, 250, 251, 252, 253, 254)
and total_UIS - non_useable_UIS >= 3
group by exp_key;



    cursor.execute(
    """
    select total_UIS - non_useable_UIS, exp_key 
    from hroest.result_srmuis
    where exp_key between 238 and 261
    """ )
    result = cursor.fetchall()

    s  = ''
    for i in range(25): s += '%s,' % i

    s  += '\n'
    for i in range(238,262):
        h, n = numpy.histogram( [r[0] for r in result if r[1] == i], 25, (0,25) )
        for aa in h: s += '%s,' % aa
        cursor.execute(" select short_description, comment3 from hroest.experiment where id = %s" % i)
        try: s += '%s, %s\n'  % cursor.fetchall()[0]
        except IndexError: s += '\n'


===========================================================================




    cursor.execute(
    """
    create temporary table hroest.nr_transitions as 
    select parent_id, count(*) as transitions from 
    hroest.srmPeptides_yeast srm 
    inner join hroest.srmTransitions_yeast tr on srm.transition_group = tr.group_id
    and q1 between 400 and 1400
    and q3 between 400 and 1400
    and q1_charge = 2 and q3_charge = 1
    group by parent_id
    ;
    """ )
    cursor.execute(" alter table hroest.nr_transitions add index(parent_id);" )

    cursor.execute(
    """
    select parent_key, free_transitions, exp_key, n.transitions, modified_sequence ,
    free_transitions * n.transitions free
    from hroest.result_srmcompare r
    inner join hroest.srmPeptides_yeast srm on srm.parent_id = r.parent_key
    inner join hroest.nr_transitions n on n.parent_id = r.parent_key
    #inner join ddb.peptide on peptide.id = srm.peptide_key
    #   and length(sequence) between 7 and 21
       and q1 between 400 and 1400
    """ )
    result = cursor.fetchall()

    data = []
    h, n = numpy.histogram( [r[1]*100 for r in result if r[2] == 50], 100 )
    cumh = collider.get_cum_dist( [hq*100.0/numpy.sum(h) for hq in h] )
    data.append( [cumh,n])
    h, n = numpy.histogram( [r[1]*100 for r in result if r[2] == 53], 100 )
    cumh = collider.get_cum_dist( [hq*100.0/numpy.sum(h) for hq in h] )
    data.append( [cumh,n])

    data = []
    h, n = numpy.histogram( [r[5] for r in result if r[2] == 50], 20, (0,20) )
    cumh = collider.get_cum_dist( [hq*100.0/numpy.sum(h) for hq in h] )
    data.append( [cumh,n])
    h, n = numpy.histogram( [r[5] for r in result if r[2] == 53], 20, (0,20) )
    cumh = collider.get_cum_dist( [hq*100.0/numpy.sum(h) for hq in h] )
    data.append( [cumh,n])




    import gnuplot
    titles = ['29 Da, 50 ppm', '0.7 Da, 0.7 Da']
    output = '/tmp/test.eps'
    xlabel = 'Unique transitions per peptide / %'
    ylabel = 'Cummulative Occurence / %'
    title = 'Comparison at 2 SSRCalc unit with yeast N^{14} + N^{15} background'
    cmd = 'with lines '
    ba = """
        #set xrange[0:20]
        set key top left
        """
    gnu = hlib.gnuplot.Gnuplot.draw_from_multiple_data(data, cmd,
            output, xlabel, ylabel, '/tmp/htmp', title, 
            keep_data = True, datatitles=titles, body_add = ba, nocolor=True)

}}}

    {{{hroest style figure (from result_srmuis ) Figure 4 Hroest 

    cursor.execute(
    """
    select total_UIS - non_useable_UIS as useable, 1 as dummy, exp_key from 
    hroest.result_srmuis
        where (exp_key between 122 and 130 or exp_key between 238 and 271 
        or exp_key between 580 and 600 )
        and uisorder = 1
    """ )
    result = cursor.fetchall()

    data = []
    # for ek in range(124,130):
    #     h, n = numpy.histogram( [r[0] for r in result if r[2] == ek], 30, (0,30) )
    #     cumh = collider.get_cum_dist( [hq*100.0/numpy.sum(h) for hq in h] )
    #     data.append( [cumh,n])
    #     print ek, len(n)


    titles = [ 
        '800 Da, 10ppm, 2SSRCalc',
        '800 Da, 10ppm',
        '800 Da, 1ppm',
        '800 Da, 0.1ppm',
        '1.2 Da, 2.0 Da',
        '1.0 Da, 1.0 Da',
        '25 Da 10ppm', 
        '0.7 Da, 0.7 Da',
        'test', 
        ]
    data_ranges = [239,262,263,264, 249, 250, 245]

    titles = [ 
        '800 Da, 10ppm 2SSR',
        '800 Da, 1ppm 2SSR',
        '800 Da, 0.1ppm 2SSR',
  #      '800 Da, 10ppm ',
        '1.2 Da, 2.0 Da',
        '1.0 Da, 1.0 Da',
        '25 Da 10ppm', 
        '0.7 Da, 0.7 Da',
        'test', 
        ]
    data_ranges = [ 239,240,241, 249,250 ]

    data = []
    #for ek in [ 245 , 250 , 252 , 268 , 269 , 583 , 586 , 587 , 589 ]:
    #for ek in [ 268, 269, 250 , 245, 252 ]:
    for ek in data_ranges:
        h, n = numpy.histogram( [r[0] for r in result if r[2] == ek], 30, (0,30) )
        cumh = collider.get_cum_dist( [hq*100.0/numpy.sum(h) for hq in h] )
        data.append( [cumh,n])

    import gnuplot
    output = '/tmp/fig1.eps'
    xlabel = 'Unique transitions per peptide'
    ylabel = 'Cummulative Occurence / %'
    title = 'Comparison at 2 SSRCalc unit with yeast background'
    cmd = 'with lines lw 1'
    ba = """
        #set xrange[0:20]
        set key bottom right 
        """
    gnu = hlib.gnuplot.Gnuplot.draw_from_multiple_data(data, cmd,
            output, xlabel, ylabel, '/tmp/htmp', title, 
            keep_data = True, datatitles=titles, body_add = ba, nocolor=False)


    cursor.execute( """
    select sum_non / sum_total as pcnt_non, uisorder, exp_key 
    from hroest.srmuis_tmp32 
    where exp_key between 122 and 130
    order by exp_key, uisorder
    """ )
    result = cursor.fetchall()
    data = []
    for ek in range(123,130):
        h = [float(r[0])*100 for r in result if r[2] == ek]
        n = [r[1] for r in result if r[2] == ek]
        data.append( [h,n])

    import gnuplot
    xlabel = 'Number of transitions'
    ylabel = 'Probability of non-unique transition set / %'
    output = '/tmp/test.eps'
    title = 'Comparison at 2 SSRCalc unit with yeast background'
    cmd = 'with lines'
    ba = """
        set xrange[1:3]
        set yrange[0:60]
        set xtics 1
        set mxtics 1
        set key top right
        """
    gnu = hlib.gnuplot.Gnuplot.draw_from_multiple_data(data, cmd,
            output, xlabel, ylabel, '/tmp/htmp', title, 
            keep_data = True, datatitles=titles, body_add = ba, nocolor=False)

}}}

    {{{ hroest complete Q1 Q3 map




    pointmap = {
    #all 9999 SSRCalc
    138 : (1 , 0.01),
    139 : (1 , 0.05),
    140 : (1 , 0.1),
    141 : (1 , 0.5),
    142 : (1 , 1),
    143 : (10 , 0.01),
    144 : (10 , 0.05),
    145 : (10 , 0.1),
    146 : (10 , 0.5),
    147 : (10 , 1),
    148 : (20 , 0.01),
    149 : (20 , 0.05),
    150 : (20 , 0.1),
    151 : (20 , 0.5),
    152 : (20 , 1),
    153 : (50 , 0.01),
    154 : (50 , 0.05),
    155 : (50 , 0.1),
    156 : (50 , 0.5),
    157 : (50 , 1),
    158 : (150 , 0.01),
    159 : (150 , 0.05),
    160 : (150 , 0.1),
    161 : (150 , 0.5),
    162 : (150 , 1),
    163 : (150 , 0.001),
 190 : (1.0 , 0.5     ),
 191 : (1.0 , 0.55    ),
 192 : (1.0 , 0.6     ),
 193 : (1.0 , 0.65    ),
 194 : (1.0 , 0.7     ),
 195 : (1.0 , 0.75    ),
 196 : (1.0 , 0.8     ),
 197 : (1.0 , 0.85    ),
 198 : (1.0 , 0.9     ),
 199 : (1.0 , 0.95    ),
 200 : (10.0 , 0.125  ),
 201 : (10.0 , 0.15   ),
 202 : (10.0 , 0.175  ),
 203 : (10.0 , 0.2    ),
 204 : (10.0 , 0.25   ),
 205 : (10.0 , 0.4    ),
 206 : (20.0 , 0.05   ), 
 207 : (20.0 , 0.055  ), 
 208 : (20.0 , 0.06   ), 
 209 : (20.0 , 0.065  ), 
 210 : (20.0 , 0.07   ), 
 211 : (20.0 , 0.075  ), 
 212 : (20.0 , 0.08   ), 
 213 : (20.0 , 0.085  ), 
 214 : (20.0 , 0.09   ), 
 215 : (20.0 , 0.095  ), 
 216 : (50.0 , 0.005  ), 
 217 : (50.0 , 0.01   ), 
 218 : (50.0 , 0.02   ), 
 219 : (50.0 , 0.03   ), 
 220 : (50.0 , 0.04   ), 
 221 : (50.0 , 0.045  ), 
 222 : (50.0 , 0.055  ), 
 223 : (50.0 , 0.06   ), 
 224 : (150.0 , 0.02  ), 
 225 : (150.0 , 0.03  ), 
 226 : (150.0 , 0.04  ), 
 227 : (150.0 , 0.045 ), 
 228 : (150.0 , 0.055 ), 
 229 : (150.0 , 0.06  ), 
 230 : (250.0 , 0.02  ), 
 231 : (250.0 , 0.03  ), 
 232 : (250.0 , 0.04  ), 
 233 : (250.0 , 0.045 ), 
 234 : (250.0 , 0.055 ), 
 235 : (250.0 , 0.06  ), 
 591 : (500.0 , 0.01  ), 
 592 : (500.0 , 0.02  ), 
 593 : (500.0 , 0.03  ), 
 594 : (500.0 , 0.04  ), 
 595 : (500.0 , 0.045 ), 
 596 : (500.0 , 0.055 ), 
 597 : (500.0 , 0.06  ), 
 598 : (500.0 , 0.1   ), 
    }




    select avg((total_UIS - non_useable_UIS) / total_UIS )
    from hroest.result_completegraph
    where uisorder = 2
    and exp_key = 212



select count(*)
, exp_key from hroest.result_completegraph
#where exp_key between 200 and 205
#where exp_key between 206 and 215
where exp_key between 216 and 223
and uisorder = 1
and total_UIS - non_useable_UIS >= 2
group by exp_key;


create table hroest.result_completegraph_analysis 
select avg( (total_UIS - non_useable_UIS) / total_UIS) as useable_pcnt,
 uisorder, exp_key, comment3 
from hroest.result_completegraph      
inner join hroest.experiment on experiment.id = exp_key  
group by exp_key, uisorder;





update hroest.experiment set comment3 = 'q1: 1.0; q3: 0.01' where id = 138;
update hroest.experiment set comment3 = 'q1: 1.0; q3: 0.05' where id = 139;
update hroest.experiment set comment3 = 'q1: 1.0; q3: 0.10' where id = 140;
update hroest.experiment set comment3 = 'q1: 1.0; q3: 0.50' where id = 141;
update hroest.experiment set comment3 = 'q1: 1.0; q3: 1.00' where id = 142;
update hroest.experiment set comment3 = 'q1: 10.0; q3: 0.01' where id = 143;
update hroest.experiment set comment3 = 'q1: 10.0; q3: 0.05' where id = 144;
update hroest.experiment set comment3 = 'q1: 10.0; q3: 0.10' where id = 145;
update hroest.experiment set comment3 = 'q1: 10.0; q3: 0.50' where id = 146;
update hroest.experiment set comment3 = 'q1: 10.0; q3: 1.00' where id = 147;
update hroest.experiment set comment3 = 'q1: 20.0; q3: 0.01' where id = 148;
update hroest.experiment set comment3 = 'q1: 20.0; q3: 0.05' where id = 149;
update hroest.experiment set comment3 = 'q1: 20.0; q3: 0.10' where id = 150;
update hroest.experiment set comment3 = 'q1: 20.0; q3: 0.50' where id = 151;
update hroest.experiment set comment3 = 'q1: 20.0; q3: 1.00' where id = 152;
update hroest.experiment set comment3 = 'q1: 50.0; q3: 0.01' where id = 153;
update hroest.experiment set comment3 = 'q1: 50.0; q3: 0.05' where id = 154;
update hroest.experiment set comment3 = 'q1: 50.0; q3: 0.10' where id = 155;
update hroest.experiment set comment3 = 'q1: 50.0; q3: 0.50' where id = 156;
update hroest.experiment set comment3 = 'q1: 50.0; q3: 1.00' where id = 157;
update hroest.experiment set comment3 = 'q1: 150.0; q3: 0.01' where id = 158;
update hroest.experiment set comment3 = 'q1: 150.0; q3: 0.05' where id = 159;
update hroest.experiment set comment3 = 'q1: 150.0; q3: 0.10' where id = 160;
update hroest.experiment set comment3 = 'q1: 150.0; q3: 0.50' where id = 161;
update hroest.experiment set comment3 = 'q1: 150.0; q3: 1.00' where id = 162;







select useable_pcnt, uisorder, exp_key, experiment.comment3 
from hroest.result_completegraph_analysis 
inner join  hroest.experiment on experiment.id = exp_key
where exp_key in (142, 149, 150, 154, 159) order by uisorder, exp_key









0
1
2
3
4
5
6
7
8
9

    cursor.execute(
    """
    select (total_UIS - non_useable_UIS) / total_UIS as useable_pcnt, uisorder, exp_key from 
    hroest.result_completegraph
    where uisorder = 2
    #and exp_key between  138 and 163
    """ )
    result = cursor.fetchall()


    data1 = {}
    data2 = {}
    data3 = {}
    for ek in pointmap.keys():
        #currdata = [float(r[0]) for r in result if r[2] == ek and r[1] == 1]
        #data1[ek] = numpy.average(currdata)
        currdata = [float(r[0]) for r in result if r[2] == ek and r[1] == 2]
        data2[ek] = numpy.average(currdata)
        #currdata = [float(r[0]) for r in result if r[2] == ek and r[1] == 3]
        #data3[ek] = numpy.average(currdata)

        h, n = numpy.histogram( , 30, (0,30) )
        cumh = collider.get_cum_dist( [hq*100.0/numpy.sum(h) for hq in h] )
        data.append( [cumh,n])



for k in data2:
    print pointmap[k][0], pointmap[k][1],data2[k] 


for k in data1:
    if data1[k] < 0.0020 and data1[k] > 0.0:
        print k,pointmap[k],data1[k] *100

for k in data2:
    if data2[k] < 0.6 and data2[k] > 0.2:
        print k,pointmap[k],data2[k] 


for k in data2:
    if data2[k] < 0.85 and data2[k] > 0.6:
        print k,pointmap[k],data2[k] 


for k in data3:
    if data3[k] < 0.98 and data3[k] > 0.8:
        print k,pointmap[k],data3[k] 




############higher == better
############

#1data
#1Da/1Da
0.0324046971713
#20Da / 50mDa 
0.708949416986
#20Da / 100mDa 
0.15957606459
#50Da / 50mDa 
0.248670925303
#150Da / 50mDa 
0.110951958159
#150Da / 10mDa 
0.606860008601


#2data
#1Da/1Da
0.422533811535
#10Da / 100mDa 
0.467749669859
#20Da / 50mDa 
0.596599187664
#20Da / 100mDa 
0.32510815146
#50Da / 50mDa 
0.427
#150Da / 50mDa 
0.298477131917
#150Da / 10mDa 
0.759989313708




#3data
#1Da/1Da
0.900423977555
#10Da / 100mDa 
0.928531104298
#20Da / 50mDa 
0.955863807651
#20Da / 100mDa 
0.896544923105
#50Da / 50mDa 
0.928519673743
#150Da / 50mDa 
0.928715289319
#150Da / 10mDa 
0.984139073013


f = open('/tmp/data' , 'w')
for k,v in data2.iteritems():
    f.write('%s %s %s\n'  % (pointmap[k][0], pointmap[k][1], v))



pointrowhash  = []
for k,v in data2.iteritems():
    pointrowhash.append( (pointmap[k][0], pointmap[k][1], v) )


pointrowhash.sort( lambda x,y: cmp( x[0]+x[1], y[0]+y[1]) )
for k,v in zip(pointrowhash[:-1], pointrowhash[1:]):
    if k[0] == v[0]:
        print "set arrow from %s,%s,%s to %s,%s,%s nohead" % (k[0], k[1], k[2], v[0], v[1], v[2])

for k,v in zip(pointrowhash[:-5], pointrowhash[6:]):
    if True:
        print "set arrow from %s,%s,%s to %s,%s,%s nohead" % (k[0], k[1], k[2], v[0], v[1], v[2])

for k,v in pointrowhash.iteritems():
    print "set arrow from 1,0.01,0.94 to 1,0.05,0.85 nohead"




    q1vals = [1, 10, 20, 50, 150]
    q3vals = [ 1000, 125, 70, 50, 40 ]
    q3vals = [q/1000.0 for q in q3vals]
    #q3vals = numpy.log(q3vals ) 
    data = [ [q3vals,  q1vals ] ]
    q1vals = [10, 20, 50, 150]
    q3vals = [ 500, 200, 100, 80 ]
    q3vals = [q/1000.0 for q in q3vals]
    data.append( [q3vals,  q1vals ] )
    q1vals = [1, 10, 20, 50, 150]
    q3vals = [ 100, 50, 25, 10 , 7]
    q3vals = [q/1000.0 for q in q3vals]
    data.append( [q3vals,  q1vals ] )
    q1vals = [1, 10, 20, 50, 150]
    q3vals = [ 500, 50, 25, 10 , 7]
    q3vals = [q/1000.0 for q in q3vals]
    data.append( [q3vals,  q1vals ] )
    #from data3
    q1vals = [10, 20, 50, 150]
    q3vals = [ 1000, 500, 200, 100 ]
    q3vals = [q/1000.0 for q in q3vals]
    data.append( [q3vals,  q1vals ] )
    newd = []
    for d in data:
        newd.append([ numpy.log(d[0])/numpy.log(10), numpy.log(d[1]) /numpy.log(10)])

    #data = newd
    import gnuplot
    titles = [ 
        'Contour line 0.42',
        'Contour line 0.7',
        'Contour line 0.17',
        'Contour line data3 0.8',
        ]
    output = '/tmp/test.eps'
    xlabel = 'Q1 Value / Da'
    ylabel = 'Q3 Value / Da'
    title = 'Comparison at 2 SSRCalc unit with yeast background'
    cmd = 'with linespoints '
    ba = """
        set xrange[1:200]
        set key top right
        set log xy
        """
    gnu = hlib.gnuplot.Gnuplot.draw_from_multiple_data(data, cmd,
            output, xlabel, ylabel, '/tmp/htmp', title, 
            keep_data = True, datatitles=titles, body_add = ba, nocolor=False)


    cursor.execute( """
    select sum_non / sum_total as pcnt_non, uisorder, exp_key 
    from hroest.srmuis_tmp32 
    where exp_key between 122 and 130
    order by exp_key, uisorder
    """ )
    result = cursor.fetchall()
    data = []
    for ek in range(123,130):
        h = [float(r[0])*100 for r in result if r[2] == ek]
        n = [r[1] for r in result if r[2] == ek]
        data.append( [h,n])

    import gnuplot
    xlabel = 'Number of transitions'
    ylabel = 'Probability of non-unique transition set / %'
    output = '/tmp/test.eps'
    title = 'Comparison at 2 SSRCalc unit with yeast background'
    cmd = 'with lines'
    ba = """
        set xrange[1:3]
        set yrange[0:60]
        set xtics 1
        set mxtics 1
        set key top right
        """
    gnu = hlib.gnuplot.Gnuplot.draw_from_multiple_data(data, cmd,
            output, xlabel, ylabel, '/tmp/htmp', title, 
            keep_data = True, datatitles=titles, body_add = ba, nocolor=False)

}}}

}}}

}}}

#
#this is how do the readout
select id, name, short_description from experiment
where name like 'paola%';

#{{{ Paola readout
create temporary table allowed_sequences as
select modified_sequence from hroest.srmPeptides_yeast ;
alter table allowed_sequences add index(modified_sequence);

#in result_srmpaola we store the parent_id of the precursor
#this is usually based on the one found in srmPeptides_yeast but not all cases
#e.g. if we use the srmPeptides_yeast_pepatlas_union table, we have different 
#ids and thus need to join with a different table!

#43 - 45/46

#new paola : 66 - 73
select @exp_key := 72; 

select @exp_key := 177; 
select @exp_key := 178; 
select @exp_key := 182; 

select @exp_key := 619; 
select exp_key, min_transitions, count(*) from hroest.result_srmpaola r 
#inner join hroest.srmPeptides_yeast_pepatlas_union srm on srm.parent_id = r.parent_key 
#inner join allowed_sequences seqs on seqs.modified_sequence = srm.modified_sequence
where exp_key = @exp_key 
#where exp_key > 74 and exp_key < 91
#where exp_key > 90 and exp_key < 102
#and srm.modified_sequence in (select modified_sequence from  allowed_sequences )
#comment out to see total
group by min_transitions
; 

drop table tmp34 ;
create temporary table tmp34 as 
select exp_key, min_transitions, count(*) as occurence from hroest.result_randomtransitions r 
where exp_key > 0
group by min_transitions, exp_key;


select * from tmp34;
select min_transitions, std(occurence), AVG(occurence), occurence from tmp34
group by min_transitions;

select std(occurence) from tmp34
group by min_transitions;

# SRM Atlas selection
+---------+-----------------+----------+
| exp_key | min_transitions | count(*) |
+---------+-----------------+----------+
|     599 |               1 |        4 |
|     599 |               2 |     7425 |
|     599 |               3 |    22030 |
|     599 |               4 |     9395 |
|     599 |               5 |     2252 |
|     599 |               6 |      529 |
|     599 |               7 |      141 |
|     599 |               8 |       49 |
|     599 |               9 |       19 |
|     599 |              10 |        9 |
|     599 |              11 |       37 |


# y ions above precursor 
+---------+-----------------+----------+
|     601 |               2 |     4243 | 
|     601 |               3 |    23871 | 
|     601 |               4 |    11327 | 
|     601 |               5 |     2038 | 
|     601 |               6 |      411 | 
+---------+-----------------+----------+


#random selection
+---------+-----------------+----------+
| exp_key | min_transitions | count(*) |
+---------+-----------------+----------+
|     603 |               1 |       72 | 
|     603 |               2 |     9040 | 
|     603 |               3 |    18692 | 
|     603 |               4 |     9495 | 
|     603 |               5 |     2981 | 
|     603 |               6 |      902 | 
|     603 |               7 |      289 | 
|     603 |               8 |      132 | 
|     603 |               9 |       59 | 
|     603 |              10 |       25 | 
|     603 |              11 |      203 | 
+---------+-----------------+----------+



{{{ Python code (automated code execution)

cursor.execute('use hroest')
cursor.execute(
""" create temporary table allowed_sequences as
select modified_sequence from hroest.srmPeptides_yeast ;
 """)
cursor.execute('alter table allowed_sequences add index(modified_sequence);')


s = ''
#new paola : 66 - 73
#pepatlas are  burned up to and including 184
#compare 187 to 186
#for key in [186, 185, 179, 180, 236, 237]:
for key in [599, 601, 603]:
    cursor.execute( 'select @exp_key := %s' % key)
    cursor.execute(
    """
    select min_transitions, count(*) from hroest.result_srmpaola r 
    inner join hroest.srmPeptides_yeast srm on srm.parent_id = r.parent_key 
    #inner join allowed_sequences seqs on seqs.modified_sequence = srm.modified_sequence
    where exp_key = @exp_key 
    #and srm.modified_sequence in (select modified_sequence from hroest.allowed_sequences )
    #comment out to see total
    group by min_transitions
    ; 
    """
    )
    r = cursor.fetchall()
    for rr in r: s += ( '%s;' % rr[1])
    s += '\n'

#new paola : 66 - 73
for key in [236, 237]:
    cursor.execute( 'select @exp_key := %s' % key)
    cursor.execute(
    """
    select min_transitions, count(*) from hroest.result_srmpaola r 
    inner join hroest.srmPeptides_yeast_pepatlas_union srm on srm.parent_id = r.parent_key 
    #inner join allowed_sequences seqs on seqs.modified_sequence = srm.modified_sequence
    where exp_key = @exp_key 
    and srm.modified_sequence in (select modified_sequence from hroest.allowed_sequences )
    #comment out to see total
    group by min_transitions
    ; 
    """
    )
    r = cursor.fetchall()
    for rr in r: s += ( '%s;' % rr[1])
    s += '\n'

print s


#for key in [72,73,71,70,68,69]:
for key in [67,66]:
    cursor.execute( 'select @exp_key := %s' % key)
    cursor.execute(
    """
#peptides
select mt, count(*) from
(
select min(min_transitions) as mt, peptide_key from hroest.result_srmpaola r 
inner join srmPeptides_yeast_pepatlas_union srm on srm.parent_id = r.parent_key 
#inner join srmPeptides_yeast srm on srm.parent_id = r.parent_key 
where exp_key = @exp_key 
and srm.modified_sequence in (select modified_sequence from  allowed_sequences )
group by peptide_key
) tmp
group by mt ; 
"""
    )
    r = cursor.fetchall()
    for rr in r: s += ( '%s;' % rr[1])
    s += '\n'

#for key in [67,66]:
for key in [72,73,71,70,68,69]:
    cursor.execute( 'select @exp_key := %s' % key)
    cursor.execute(
    """
    #proteins
select mt, count(*) from
(
select min(min_transitions) as mt, protein_key from hroest.result_srmpaola r 
inner join srmPeptides_yeast srm on srm.parent_id = r.parent_key 
#inner join srmPeptides_yeast_pepatlas_union srm on srm.parent_id = r.parent_key 
inner join ddb.protPepLink l on l.peptide_key = srm.peptide_key
where exp_key = @exp_key 
#and srm.modified_sequence in (select modified_sequence from  allowed_sequences )
group by protein_key
) tmp
group by mt ; 
"""
    )
    r = cursor.fetchall()
    for rr in r: s += ( '%s;' % rr[1])
    s += '\n'




}}}

#peptides
select mt, count(*) from
(
select min(min_transitions) as mt, peptide_key from hroest.result_srmpaola r 
inner join srmPeptides_yeast srm on srm.parent_id = r.parent_key 
where exp_key = @exp_key 
and srm.modified_sequence in (select modified_sequence from  allowed_sequences )
group by peptide_key
) tmp
group by mt ; 

#proteins
select mt, count(*) from
(
select min(min_transitions) as mt, protein_key from hroest.result_srmpaola r 
#inner join srmPeptides_yeast srm on srm.parent_id = r.parent_key 
inner join srmPeptides_yeast_pepatlas_union srm on srm.parent_id = r.parent_key 
inner join ddb.protPepLink l on l.peptide_key = srm.peptide_key
where exp_key = @exp_key 
and srm.modified_sequence in (select modified_sequence from  allowed_sequences )
group by protein_key
) tmp
group by mt ; 



#create the table with min_nr transitions for each precursor
select @exp_key := 13;
select 
RIGHT(ac, LENGTH(ac) - INSTR(ac, '.')) as protein, 
sequence as stripped_sequence, 
CONCAT(modified_sequence, '/', q1_charge) as modified_sequence,
min_transitions 
from hroest.result_srmpaola r 
inner join srmPeptides_yeast srm on srm.parent_id = r.parent_key 
inner join ddb.peptide p on srm.peptide_key = p.id
inner join ddb.protPepLink l on l.peptide_key = p.id
inner join ddb.protein prot on prot.id = l.protein_key
inner join ddb.isbAc ac on prot.sequence_key = ac.sequence_key
where exp_key = @exp_key 
and db='eggnog'
into outfile '/home/hroest/data/allbackground_nossrcalc_1Da1Da.csv'
; 

}}}

{{{ SWATH workflow


reload( collider )
mycollider = collider.SRMcollider()
mycollider.find_clashes(db, par, exp_key = None)


db = MySQLdb.connect(read_default_file="~/.my.cnf")
cursor = db.cursor()

common_filename = par.get_common_filename()
query = """
insert into hroest.experiment  (name, short_description,
description, comment1, comment2, super_experiment_key, ddb_experiment_key)
VALUES (
    'SWATH yeast', '%s', '%s', '%s', '%s', 47, 0
)
""" %( common_filename + '_' + par.peptide_table.split('.')[1], 
      par.experiment_type, par.peptide_table, par.transition_table)
cursor.execute(query)
myid = db.insert_id()

prepare = [ [k,v, myid]  for k,v in mycollider.allpeps.iteritems() ]
#save our result ("the min nr transitions per precursor") linked with 
#our new experiment key
cursor.executemany(
""" insert into hroest.result_srmcompare (parent_key, free_transitions, exp_key) 
    VALUES (%s,%s,%s) """, prepare
)

create temporary table nr_transitions as 
select parent_id, count(*) as transitions from 
hroest.srmPeptides_yeast srm 
inner join hroest.srmTransitions_yeast tr on srm.transition_group = tr.group_id
and q1 between 400 and 1400
and q3 between 400 and 1400
and q1_charge = 2 and q3_charge = 1
group by parent_id
;
alter table nr_transitions add index(parent_id);


select parent_key, free_transitions, exp_key, n.transitions, modified_sequence ,
free_transitions * n.transitions free
from hroest.result_srmcompare r
inner join hroest.srmPeptides_yeast srm on srm.parent_id = r.parent_key
inner join hroest.nr_transitions n on n.parent_id = r.parent_key
#inner join ddb.peptide on peptide.id = srm.peptide_key
where exp_key = 49
#   and length(sequence) between 7 and 21
   and q1 between 400 and 1400
limit 10;


   create table hroest.result_srmcompare (
   parent_key int, 
   free_transitions double, 
   exp_key int
   )
   alter table hroest.result_srmcompare add index(exp_key, parent_key);



#
#this is how do the readout in SWATH for ludovic
select id, name, short_description from experiment
where id between 166 and 173
;




select count(*)
, exp_key from hroest.result_srmuis
where exp_key between 170 and 173
#and total_UIS - non_useable_UIS >= 3
group by exp_key;

select count(*)
, exp_key from hroest.result_srmuis
where exp_key between 166 and 169
#and total_UIS - non_useable_UIS >= 3
group by exp_key;






select id, short_description, name, comment3 from experiment where name like 'ludo%' order by id;

select count(*)
, exp_key from hroest.result_srmuis
where exp_key between 238 and 261
#and total_UIS - non_useable_UIS >= 3
group by exp_key;


#regular 
select count(*) / 121009
, exp_key from hroest.result_srmuis
where exp_key between 238 and 261
and exp_key in (238, 239, 241, 245, 249, 250, 251, 252, 253, 254)
and total_UIS - non_useable_UIS >= 5
group by exp_key
order by exp_key;


#peptide atlas
select count(*) / 53002
, exp_key from hroest.result_srmuis
where exp_key between 238 and 261
and exp_key not in (238, 239, 240, 241, 245, 249, 250, 251, 252, 253, 254)
and total_UIS - non_useable_UIS >= 3
group by exp_key;


order by exp_key



# regular
 ### new data with smaller RT, e.g. 0.33 instead of 2
select count(*) as occurence #/ 121009 * 100 as prcnt
, exp_key , total_UIS - non_useable_UIS as useable_transitions from hroest.result_srmuis
where exp_key between 262 and 267
#and total_UIS - non_useable_UIS >= 3
group by total_UIS - non_useable_UIS, exp_key
order by exp_key, total_UIS - non_useable_UIS;


s = ''
cursor.execute(
"""
select  short_description
, id from hroest.experiment
where id between 166 and 173
order by id
"""
)
for r in cursor.fetchall(): s += '%s;' % r[0]

s = ''
cursor.execute(
"""
#peptide atlas
select count(*) / 53002
, exp_key from hroest.result_srmuis
where exp_key between 238 and 261
and exp_key not in (238, 239, 240, 241, 245, 249, 250, 251, 252, 253, 254)
and total_UIS - non_useable_UIS >= 5
group by exp_key;
"""
)
for r in cursor.fetchall(): s += '%s;' % r[0]

print s


{{{ #use only peptides between 400 and 1200 



s = ''
cursor.execute(
"""
#
select count(*) / 53002
, exp_key from hroest.result_srmuis
where exp_key between 238 and 261
and exp_key not in (238, 239, 240, 241, 245, 249, 250, 251, 252, 253, 254)
and total_UIS - non_useable_UIS >= 5
group by exp_key;
"""
)
for r in cursor.fetchall(): s += '%s;' % r[0]


#for key in [246, 256, 258, 270, 271]:
#for key in [245, 250, 252, 268, 269]:
#for key in [589]:
#for key in [607, 610, 611]:
for key in [620,622,624,628,630]:
    s = ''
    tmp = cursor.execute(
    """
    select count(*) /  111880 
    , total_UIS - non_useable_UIS , exp_key
    from hroest.result_srmuis
    inner join hroest.srmPeptides_yeast
    pep on  pep.parent_id = parent_key
    where exp_key in (%s)
    and isotope_nr = 0 and q1 between 400.0 and 1200.0
    group by exp_key, total_UIS - non_useable_UIS 
    """ % key
    )
    all = cursor.fetchall()
    for i,r in enumerate(all): s += '%s;' % sum([ rr[0] for rr in all[:i]] )
    print key, s



#for key in [246, 256, 258, 270, 271]:
#for key in [245, 250, 252, 268, 269]:
#for key in [590]:
#for key in [ 608, 609, 612]:
for key in [621,623,625,629,631]:
    s = ''
    tmp = cursor.execute(
    """
    #peptide atlas
    select count(*) /  48087 #111880 #
    , total_UIS - non_useable_UIS , exp_key
    from hroest.result_srmuis
    inner join hroest.srmPeptides_yeast_pepatlas
    pep on  pep.parent_id = parent_key
    where exp_key in (%s)
    and isotope_nr = 0 and q1 between 400.0 and 1200.0
    group by exp_key, total_UIS - non_useable_UIS 
    """ % key
    )
    all = cursor.fetchall()
    for i,r in enumerate(all): s += '%s;' % sum([ rr[0] for rr in all[:i]] )
    print key, s




select count(*), exp_key from hroest.result_srmuis 
inner join hroest.srmPeptides_yeast_pepatlas pep on  pep.parent_id = parent_key
where exp_key = 258 
and isotope_nr = 0 and q1 between 400.0 and 1200.0
group by exp_key;

select count(*), exp_key from hroest.result_srmuis 
inner join hroest.srmPeptides_yeast_pepatlas pep on  pep.parent_id = parent_key
where exp_key = 258 
and isotope_nr = 0 and q1 between 400.0 and 1200.0
group by exp_key;


|   |     268 |
|   111880 |     269 |
|    48087 |     270 |


    
}}}

}}}

create table tmp_notinrangetree as 
select distinct parent_key from result_srmuis
where exp_key = 103 and parent_key not in
(select parent_key from tmp123)
limit 1


create temporary table tmp123 as
select parent_key from  result_srmuis where exp_key = 102;
alter table tmp123 add index(parent_key)

{{{ srmuis filler workflow

#31 yeast, 34 human
db = MySQLdb.connect(read_default_file="~/.my.cnf")
cursor = db.cursor()
common_filename = par.get_common_filename()
query = """
insert into hroest.experiment  (name, short_description,
description, comment1, comment2, super_experiment_key, ddb_experiment_key)
VALUES (
    'uis_perpeptide', '%s', '%s', '%s', '%s', 31, 0
)
""" %( common_filename + '_' + par.peptide_table.split('.')[1], 
      par.experiment_type, par.peptide_table, par.transition_table)
cursor.execute(query)
myid = db.insert_id()

reload( collider )
mycollider = collider.SRMcollider()
mycollider.find_clashes_small(db, par, UIS_only=True, exp_key = myid,
                              calc_collisions=False)


#select collisions
par.query2_add = 'and q1 between 400 and 1200'

#all peptides to search
par.query_add = ' and q1 between 400 and 1200'

{{{Verify srmuis filler


#see that the two methods produce the same results
#differences because
#old methods takes ALL peptides and excludes the one with genome_occurence > 1
select count(distinct parent_key) from result_srmuis where exp_key = 112;
select count(distinct parent_id) from srmPeptides_yeast where q1_charge = 2 and q1 between 400 and 1400 and isotope_nr = 0;

cursor.execute( " drop table hroest.tmp145;")
cursor.execute( " drop table hroest.tmp115;")
cursor.execute(
"""
create temporary table hroest.tmp145 as
    select distinct parent_id from hroest.srmPeptides_yeast where q1_charge = 2 and q1
    between 400 and 1400 and isotope_nr = 0;
"""
)
cursor.execute(
"""
alter table hroest.tmp145 add index (parent_id)
"""
)

cursor.execute(
"""
create temporary table hroest.tmp115 as
    select distinct parent_key as parent_id from hroest.result_srmuis where exp_key = 102
"""
)
cursor.execute(
"""
alter table hroest.tmp115 add index (parent_id)
"""
)


cursor.execute(
"""
select  non_useable_UIS, parent_key , total_UIS 
    from hroest.result_srmuis where exp_key in ( 102, 10113) 
    and parent_key in (select parent_id from hroest.tmp145) 
    and parent_key in (select parent_id from hroest.tmp115) 
order by parent_key, uisorder
"""
)
result = cursor.fetchall()
for i in range(1, len(result), 2):
    print result[i][0] ,  result[i-1][0], i
    assert result[i][0] == result[i-1][0]



}}}

}}}

{{{ diffdistribution analysis

common_filename = par.get_common_filename()
query = """
insert into hroest.experiment  (name, short_description,
description, comment1, comment2, super_experiment_key, ddb_experiment_key)
VALUES (
    'diffdistr', '%s', '%s', '%s', '%s', 24, 0
)
""" %( common_filename + '_' + par.peptide_table.split('.')[1], 
      par.experiment_type, par.peptide_table, par.transition_table)
cursor.execute(query)
myid = db.insert_id()



reload( collider )
mycollider = collider.SRMcollider()
pepids = mycollider._get_unique_pepids(par, cursor)

import random
random.shuffle(pepids)
mycollider.find_clashes(db, par, pepids = pepids[:10000], exp_key = -100)





q1r = [-par.q1_window, par.q1_window] 
q3r = [-par.q3_window, par.q3_window] 

import numpy
self = mycollider
#self.q1min_hist += numpy.histogram( self.q1min_distr, 1000, q1r)[0]
bins = numpy.histogram( self.q1min_distr, 1000, q1r)[1]
prepare = [[myid, b, v, 1] for b,v in zip(bins, self.q1min_hist)]
cursor.executemany(
""" insert into hroest.result_diffhist (exp_key, bin, value, type ) 
    VALUES (%s,%s,%s,%s) """, prepare
)


#self.q3min_hist += numpy.histogram( self.q3min_distr, 1000, q3r)[0]
bins = numpy.histogram( self.q3min_distr, 1000, q3r)[1]
prepare = [[myid, b, v, 2] for b,v in zip(bins, self.q3min_hist)]
cursor.executemany(
""" insert into hroest.result_diffhist (exp_key, bin, value, type ) 
    VALUES (%s,%s,%s,%s) """, prepare
)



#self.q1all_hist += numpy.histogram( self.q1all_distr, 1000, q1r)[0]
bins = numpy.histogram( self.q1all_distr, 1000, q1r)[1]
prepare = [[myid, b, v, 3] for b,v in zip(bins, self.q1all_hist)]
cursor.executemany(
""" insert into hroest.result_diffhist (exp_key, bin, value, type ) 
    VALUES (%s,%s,%s,%s) """, prepare
)

#self.q3all_hist += numpy.histogram( self.q3all_distr, 1000, q3r)[0]
bins = numpy.histogram( self.q3all_distr, 1000, q3r)[1]
prepare = [[myid, b, v, 4] for b,v in zip(bins, self.q3all_hist)]
cursor.executemany(
""" insert into hroest.result_diffhist (exp_key, bin, value, type ) 
    VALUES (%s,%s,%s,%s) """, prepare
)


create table hroest.result_diffhist 
(
exp_key int, 
bin double,
value double,
type int(4)
);
alter table hroest.result_diffhist add index(exp_key, type);

}}}

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
###########################################################################

select pr.sequence_key,sequence,rank,label,q1,q3,q1_charge,q3_charge
from protein pr inner join protPepLink ln on ln.protein_key = pr.id
inner join peptide pe on ln.peptide_key = pe.id 
inner join mrmAssay a on a.peptide_key = pe.id 
inner join transition on mrmAssay_key = a.id 
where pr.experiment_key = 3373 
and repl_by = 0 
and q1_charge = 2 
and label = 'none' 
#and rank = 3
and sequence = 'LPSTSGSEGVPFR'  
#limit 1
;



create temporary table all_assays as 
select pr.sequence_key,sequence,rank,label,q1,q3,q1_charge,q3_charge
#select count(distinct sequence_key) 
from protein pr inner join protPepLink ln on ln.protein_key = pr.id
inner join peptide pe on ln.peptide_key = pe.id 
inner join mrmAssay a on a.peptide_key = pe.id 
inner join transition on mrmAssay_key = a.id 
where pr.experiment_key = 3373  # master
# where pr.experiment_key = 1358  # ncbi genome
and repl_by = 0 
and q1_charge = 2 
and label = 'none' 
#and rank = 3
#and sequence = 'LPSTSGSEGVPFR'  
#limit 1
;

alter table all_assays add index(sequence);


# 3446 peptides 
select count(distinct peptide.sequence) from ddb.peptide 
inner join all_assays on all_assays.sequence = peptide.sequence
and rank = 3
where experiment_key = 1389 ;

# 1578 proteins out of 1905
select count(distinct ln.protein_key) from ddb.peptide 
inner join all_assays on all_assays.sequence = peptide.sequence
inner join protPepLink ln on ln.peptide_key = peptide.id
and rank = 3
where experiment_key = 1389 ;


# 18026 coordinates 
select count(*) from 
(
select count(*) from ddb.peptide 
inner join all_assays on all_assays.sequence = peptide.sequence
where experiment_key = 1389 
and rank < 100
#and peptide.sequence = 'MNPLIQSLTEGQLR'
group by peptide.sequence, rank
) t
;



