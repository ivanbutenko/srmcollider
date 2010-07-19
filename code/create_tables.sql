#-- srmPeptide contains peptide information 
#-- srmTransitions contains transitions
#-- they are joined over parent_id = parent_key

#-- peptide_key links the peptides to the ddb.peptide table
#-- type is either b or y ion


select * from hroest.srmClashes 
inner join ddb.peptide on peptide.id = srmClashes.peptide_key
where q1 between 518 and 519 
and q3 between 740 and 760
and ssrcalc between 20 and 25
#limit 100
;

#drop table hroest.srmClashes;
#truncate table hroest.srmClashes;
#create table hroest.srmClashes(
#    id INT PRIMARY KEY AUTO_INCREMENT,
#    type varchar(100),
#    peptide_key INT,
#    experiment_key INT,
#    q1_charge INT,
#    q3_charge INT,
#    q1 DOUBLE,
#    q3 DOUBLE, 
#    ssrcalc DOUBLE
#);
#alter table hroest.srmClashes add index(experiment_key);
#alter table hroest.srmClashes add index(peptide_key);
#alter table hroest.srmClashes add index(q1, q3, ssrcalc);

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

#--create a table that has bins in MS1 dimenstion of size 0.7 Da and in 
#--MS2 dimensions of size 1 Da
drop table hroest.srm_counts;
create temporary table hroest.srm_counts as
select round(q1 / 0.7 ) * 0.7, round(q3), count(*) as occurence from hroest.srmClashes
where q1 between 300 and 1500 and q3 between 300 and 1500 
and experiment_key = 3061
group by round(q1 / 0.7), round(q3) #order by count(*)
;

alter table hroest.srmTest add index(experiment_key);
alter table hroest.srmTest add index(peptide_key);
alter table hroest.srmTest add index(q1, q3, ssrcalc);


