#-- srmPeptide contains peptide information 
#-- srmTransitions contains transitions
#-- they are joined over parent_id = parent_key

#-- peptide_key links the peptides to the ddb.peptide table
#-- type is either b or y ion

#-- old table == srmClashes (non-normalized; contains ALL transitions)

drop table hroest.srmPeptide;
truncate table hroest.srmPeptide;
create table hroest.srmPeptide(
    parent_id INT PRIMARY KEY AUTO_INCREMENT,
    peptide_key INT,
    q1_charge INT,
    q1 DOUBLE,
    ssrcalc DOUBLE
);
alter table hroest.srmPeptide add index(peptide_key);
alter table hroest.srmPeptide add index(q1);
alter table hroest.srmPeptide add index(ssrcalc);

drop table hroest.srmTransitions;
truncate table hroest.srmTransitions;
create table hroest.srmTransitions(
    srm_id INT PRIMARY KEY AUTO_INCREMENT,
    parent_key INT,
    q3_charge INT,
    q3 DOUBLE,
    type VARCHAR(3)
);
alter table hroest.srmTransitions add index(parent_key);
alter table hroest.srmTransitions add index(q3);

#-- create a table that has bins in MS1 dimenstion of size 0.7 Da and in 
#-- MS2 dimensions of size 1 Da and in SSRCalc dimension 4 units
drop table hroest.srmCollisions400710;
truncate table hroest.srmCollisions400710;
create table hroest.srmCollisions400710(
    coll_srm1 INT,
    coll_srm2 INT
);
alter table hroest.srmCollisions400710 add unique index(coll_srm1, coll_srm2);
drop table hroest.srmCollisions400710_all;
truncate table hroest.srmCollisions400710_all;
create table hroest.srmCollisions400710_all(
    coll_srm1 INT,
    coll_srm2 INT
);

drop table hroest.srm_counts;
create temporary table hroest.srm_counts as
select round(q1 / 0.7 ) * 0.7, round(q3), count(*) as occurence from hroest.srmClashes
where q1 between 300 and 1500 and q3 between 300 and 1500 
and experiment_key = 3061
group by round(q1 / 0.7), round(q3) #order by count(*)
;

