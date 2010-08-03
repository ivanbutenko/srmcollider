# --
# -- This script creates the two tables tmp_srm_get_all_collisions and
# -- tmp_srm_all_collisions_filtered which use srmCollisions400710 to find all
# -- possible collisions for each transition.
# --
# --

#undo this script:
#drop table hroest.tmp_srm_get_all_collisions;
#drop table hroest.tmp_srm_all_collisions_filtered;


#get the full set of collisions for all transitions
#calculate fraction of collisions for certain parent_ids
#

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

#this takes about 12 minutes
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
# -- only compare to 1 out of 4 charge state combinations
    q3_charge = 1
    and q1_charge = 2
;
alter table hroest.tmp_srm_all_collisions_filtered add index(coll_srm1);

#how many collisions did we detect?
select distinct coll_srm1  from hroest.tmp_srm_all_collisions_filtered 
#where q1
#order by q3diff
;

#which collide most
select count(*), sequence from hroest.tmp_srm_all_collisions_filtered
group by sequence
order by count(*);

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
select count(*) from hroest.tmp_srm_all_collisions_filtered
where q1diff < 0.1
;

--get all transitions for a specific peptide
set @mypID = 1001;
set @mypID = 100015;
    select q1, ssrcalc, q3, srm_id, type, sequence
    from hroest.srmPeptide
    inner join hroest.srmTransitions
      on parent_id = parent_key
    inner join ddb.peptide on id = peptide_key
    where parent_id = @mypId
    and q3 > 300 and q3 < 10000         and q1_charge = 2 and q3_charge = 1
order by q3
;

--get all peptides that are close
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

--# get all transitions of a peptide together with the best collisions
--# unique transitions will get a NULL for all except the first coloumn
create temporary table my_pep_collisions as 
select p_q3, p_type, q3, type, sequence
from hroest.tmp_srm_get_all_collisions 
inner join srmTransitions on coll_srm2 = srm_id
inner join srmPeptide on parent_key = parent_id
inner join ddb.peptide on peptide_key = peptide.id
where p_parent_key =  @mypId  group by coll_srm1 
order by p_q3;
create temporary table my_pep_transitions as
select q3 as all_q3 
from hroest.srmPeptide
inner join hroest.srmTransitions
  on parent_id = parent_key
inner join ddb.peptide on id = peptide_key
where parent_id = @mypId
and q3 > 300 and q3 < 10000         and q1_charge = 2 and q3_charge = 1
order by q3 ;
select all_q3 as p_q3, p_type, q3, type, sequence
from my_pep_transitions
left join  my_pep_collisions
on all_q3 = p_q3;


