#-- srmPeptides_test contains peptide information 
#-- srmTransitions_test contains transitions
#-- they are joined over parent_id = parent_key

#-- peptide_key links the peptides to the ddb.peptide table
#-- type is either b or y ion

#-- old table == srmClashes (non-normalized; contains ALL transitions)

#drop table hroest.srmPeptides_test;
truncate table hroest.srmPeptides_test;
create table hroest.srmPeptides_test(
    parent_id INT PRIMARY KEY AUTO_INCREMENT,
    peptide_key INT,
    q1_charge INT,
    q1 DOUBLE,
    ssrcalc DOUBLE
);
alter table hroest.srmPeptides_test add index(peptide_key);
alter table hroest.srmPeptides_test add index(q0);
alter table hroest.srmPeptides_test add index(ssrcalc);

#drop table hroest.srmTransitions_test;
truncate table hroest.srmTransitions_test;
create table hroest.srmTransitions_test(
    srm_id INT PRIMARY KEY AUTO_INCREMENT,
    parent_key INT,
    q3_charge INT,
    q3 DOUBLE,
    type VARCHAR(3)
);
alter table hroest.srmTransitions_test add index(parent_key);
alter table hroest.srmTransitions_test add index(q3);

#-- create a table that has bins in MS1 dimenstion of size 0.7 Da and in 
#-- MS2 dimensions of size 1 Da and in SSRCalc dimension 4 units
#drop table hroest.srmCollisions4010250_test;
#truncate table hroest.srmCollisions4010250_test;
create table hroest.srmCollisions4010250_test(
    coll_srm1 INT,
    coll_srm2 INT
);
alter table hroest.srmCollisions4010250_test add unique index(coll_srm1, coll_srm2);

drop table hroest.srm_counts;
create temporary table hroest.srm_counts as
select round(q1 / 0.7 ) * 0.7, round(q3), count(*) as occurence from hroest.srmClashes
where q1 between 300 and 1500 and q3 between 300 and 1500 
and experiment_key = 3061
group by round(q1 / 0.7), round(q3) #order by count(*)
;

