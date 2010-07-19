##########################################################################
##########################################################################
#get the full set of collisions for all transitions
#calculate fraction of collisions for certain parent_ids

#for each transition, see whether there are collisions
drop table hroest.tmp_srm_get_all_collisions;
create table hroest.tmp_srm_get_all_collisions as
select
parent_key as p_parent_key, q3_charge as p_q3_charge,
q3 as p_q3, type as p_type , coll_srm1 , coll_srm2,
ssrcalc as p_ssrcalc, q1 as p_q1
from srmTransitions
inner join hroest.srmCollisions400710 on coll_srm1 = srm_id
inner join hroest.srmPeptide on parent_id = parent_key
where 
    #parent_key < 1001
    q3_charge = 1
    and q1_charge = 2
#group by parent_key, q3
;
alter table hroest.tmp_srm_get_all_collisions add index( coll_srm2 );

#12 minutes
#we need to join over srmTransitisions again to get the 
#charge states of the collisions and filter for them
drop table hroest.tmp_srm_all_collisions_filtered;
create table hroest.tmp_srm_all_collisions_filtered as
select 
p_parent_key , p_q3_charge , p_q3        , p_type , coll_srm1 , coll_srm2 , 
srm_id   , parent_key , q3_charge , q3          , type , sequence, q1, ssrcalc,
p_ssrcalc- ssrcalc as ssrcalcdiff,
p_q3- q3 as q3diff,
p_q1- q1 as q1diff,
ABS( p_ssrcalc- ssrcalc) as ssrcalcdiff_abs,
ABS( p_q3- q3) as q3diff_abs,
ABS( p_q1- q1) as q1diff_abs
from hroest.tmp_srm_get_all_collisions
inner join srmTransitions on srmTransitions.srm_id = coll_srm2
inner join srmPeptide on parent_key = parent_id
inner join ddb.peptide on id = peptide_key
where
    q3_charge = 1
    and q1_charge = 2
;
alter table hroest.tmp_srm_all_collisions_filtered add index(coll_srm1);



##%%f = open( 'mydist26.py', 'w' )
##%%f.write( 'mydist26 = [' )
##%%for d in mydist:
##%%    f.write( '%f, ' % d)
##%%
##%%f.write( ']' )
##%%f.close()




###
#QUESTION: does the q1 accuracy make a big difference?
#10866
#8057
#20% less collisions ==> when 1vs4 is searched
#then the q1diff < 0.1 makes a big difference
######
#3671
#3645
# 1 % less collisions ==> when only 1vs1 is searched
#20098260
#20040839
# 0.3 % fewer collisions
#each of those takes four minutes
#they contain for each transition "coll_srm1"
#the best/closest collision in dimension q1, q3, ssrcalc
#this can be used for plotting
select count(*) from hroest.tmp_srm_all_collisions_filtered
where q1diff < 0.1
;

drop table hroest.tmp_analyse_q1diff;
create table hroest.tmp_analyse_q1diff as
select q1diff, q1diff_abs, coll_srm1 from hroest.tmp_srm_all_collisions_filtered
where q1diff < 0.1

order by coll_srm1, q1diff_abs 
;
drop table hroest.tmp_analyse_q3diff;
create table hroest.tmp_analyse_q3diff as
select q3diff, q3diff_abs, coll_srm1 from hroest.tmp_srm_all_collisions_filtered
#where q1diff < 0.1
order by coll_srm1, q3diff_abs 
;
drop table hroest.tmp_analyse_ssrcalcdiff;
create table hroest.tmp_analyse_ssrcalcdiff as
select ssrcalcdiff, ssrcalcdiff_abs, coll_srm1 from hroest.tmp_srm_all_collisions_filtered
#where q1diff < 0.1
order by coll_srm1, ssrcalcdiff_abs 
;

alter table hroest.tmp_analyse_ssrcalcdiff add index(coll_srm1);
alter table hroest.tmp_analyse_q1diff add index(coll_srm1);
alter table hroest.tmp_analyse_q3diff add index(coll_srm1);
select * from hroest.tmp_analyse_q1diff
group by coll_srm1
limit 4;

--##########################################################################
--                          PYTHON CODE
--##########################################################################
#-- make a picture of this
reload( gnuplot )

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


--##########################################################################
--                          PYTHON CODE
--##########################################################################



drop table hroest.not_in_lower_cutoff;
create temporary table hroest.not_in_lower_cutoff as
select coll_srm1 from hroest.tmp_srm_all_collisions_filtered
where coll_srm1 not in
(select coll_srm1 from hroest.tmp_q1cutoff)
;
select * from hroest.not_in_lower_cutoff limit 10;

select * from hroest.srmTransitions where srm_id = 276;
select * from hroest.srmPeptide where parent_id = 5;

select * from hroest.srmCollisions400710
inner join srmTransitions on srmTransitions.srm_id = coll_srm2
inner join srmPeptide on parent_key = parent_id
inner join ddb.peptide on id = peptide_key
 where coll_srm1 = 276;


#how many collisions did we detect?
select distinct coll_srm1  from hroest.tmp_srm_all_collisions_filtered 
where q1
#order by q3diff
;

#which collide most
select count(*), sequence from hroest.tmp_srm_all_collisions_filtered
group by sequence
order by count(*);

drop table hroest.srm_groupcollisions23;
create temporary table hroest.srm_groupcollisions23 as
select p_parent_key, count(*) as nr_collisions
from hroest.tmp_srm_all_collisions_filtered
group by p_parent_key
;
select count(*) from hroest.srm_groupcollisions23;
select * from hroest.srm_groupcollisions23;




set @mypID = 1001;
set @mypID = 100015;
    select q1, ssrcalc, q3, srm_id, type, sequence
    from hroest.srmPeptide
    inner join hroest.srmTransitions
      on parent_id = parent_key
    inner join ddb.peptide on id = peptide_key
    where parent_id = @mypId
    and q3 > 0 and q3 < 10000         and q1_charge = 2 and q3_charge = 1
order by q3
;



set @q1_accuracy = 0.01;
    select q3, peptide_key, type, sequence from hroest.srmPeptide
    inner join hroest.srmTransitions
      on parent_id = parent_key
    inner join ddb.peptide on id = peptide_key
    where   ssrcalc > 6.74 - 1.0
        and ssrcalc < 6.74 + 1.0
    and q1 > 406.699 - @q1_accuracy and q1 < 406.699 + @q1_accuracy
    and parent_id != 1001
    and q3 > 0 and q3 < 10000         and q1_charge = 2 and q3_charge = 1
group by sequence
order by q3
;



#my peptide is EHLEER
#EEHLER
#candidates are 
##GDTGLGHR, GPAPCDR, SHQGLDR , SNLEGHR , CLPSQHK , NPHLGTSS, NDVDPPR , 
##QMYNTR  , EEHLER  , MSNVPHK , HEELER  
#
| q1          | ssrcalc | q3         | srm_id |
+-------------+---------+------------+--------+
| 406.6990325 |    6.74 | 130.050415 |  62097 |   b   COLL  : EEHLER
| 406.6990325 |    6.74 |   175.1195 |  62106 |   y   COLL  : all
| 406.6990325 |    6.74 | 267.109325 |  62098 |   b   COLL  : HEELER
| 406.6990325 |    6.74 |  304.16209 |  62105 |   y   COLL  : EEHLER, HEELER
| 406.6990325 |    6.74 | 380.193385 |  62099 |   b    unique 
| 406.6990325 |    6.74 |  433.20468 |  62104 |   y    unique
| 406.6990325 |    6.74 | 509.235975 |  62100 |   b    COLL : EEHLER, HEELER
| 406.6990325 |    6.74 |  546.28874 |  62103 |   y   COLL  : HEELER 
| 406.6990325 |    6.74 | 638.278565 |  62101 |   b    COLL : all
| 406.6990325 |    6.74 |  683.34765 |  62102 |   y    COLL : EEHLER




#my peptides is EPISVSSEQVLK
#collisions occur with
##VVSGGLCPVLESR, TQYVHSPYDRPGWNPR , STAAEVQQVLNR, LLAGQQVWDASK     ,
##ENPLNGASLSWK,  TTHLMFHTVTK      , VQAVQLCQSALR, EEYGHSEVVEYYCNPR ,
##ALSSEWKPEIR,   LEALLDECANPK     , NIATSLHEICSK,
  | p_q3        | p_type |   | q3          | type | sequence         |
--+-------------+--------+---+-------------+------+------------------+
1 |  340.187235 | b      | 1 |  340.162095 | b    | TTHLMFHTVTK      |
1 |   359.26582 | y      | 1 |  359.193055 | b    | ALSSEWKPEIR      |
1 |  427.219265 | b      | 1 |  427.255645 | b    | LEALLDECANPK     |
1 |    487.3244 | y      | 1 |  487.251635 | b    | NIATSLHEICSK     |
1 |  526.287675 | b      | 1 |  526.298915 | b    | VQAVQLCQSALR     |
1 |   616.36699 | y      | 1 |  616.236705 | b    | EEYGHSEVVEYYCNPR |
1 |  700.351735 | b      | 1 |   700.39935 | y    | VVSGGLCPVLESR    |
1 |   703.39902 | y      | 1 |  703.268735 | b    | EEYGHSEVVEYYCNPR |
1 |   790.43105 | y      | 1 |   790.42453 | y    | VQAVQLCQSALR     |
1 |  829.394325 | b      | 1 |   829.42419 | y    | NIATSLHEICSK     |
1 |   889.49946 | y      | 1 |   889.40893 | y    | LEALLDECANPK     |
1 |  957.452905 | b      | 1 |   957.51577 | y    | ALSSEWKPEIR      |
1 |   976.53149 | y      | 1 |  976.452855 | b    | TQYVHSPYDRPGWNPR |
1 | 1056.521315 | b      | 1 |  1056.58017 | y    | STAAEVQQVLNR     |
1 |  1089.61555 | y      | 1 |  1089.53288 | y    | LLAGQQVWDASK     |
1 | 1169.605375 | b      | 1 | 1169.562475 | b    | NIATSLHEICSK     |
1 |  1186.66831 | y      | 1 |  1186.62203 | y    | ENPLNGASLSWK     |





    select *
    from hroest.srmPeptide
    inner join hroest.srmTransitions
      on parent_id = parent_key
    inner join hroest.srmCollisions400710
      on coll_srm1 = srm_id
    where parent_id = 1001
    and q3 > 0 and q3 < 10000         and q1_charge = 2 and q3_charge = 1
group by q3


select * from 



