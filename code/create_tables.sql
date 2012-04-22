#-- srmPeptides_test contains peptide information 
#-- srmTransitions_test contains transitions
#-- they are joined over parent_id = parent_key

#-- peptide_key links the peptides to the ddb.peptide table
#-- type is either b or y ion

#-- old table == srmClashes (non-normalized; contains ALL transitions)

drop table hroest.srmPeptides_test;
truncate table hroest.srmPeptides_test;
create table hroest.srmPeptides_test(
    parent_id INT PRIMARY KEY AUTO_INCREMENT,
    peptide_key INT,
    modified_sequence VARCHAR(255),
    q1_charge TINYINT,
    q1 DOUBLE,
    ssrcalc DOUBLE,
    isotope_nr TINYINT,
    transition_group INT
);
alter table hroest.srmPeptides_test add index(peptide_key);
alter table hroest.srmPeptides_test add index(q1);
alter table hroest.srmPeptides_test add index(ssrcalc);
alter table hroest.srmPeptides_test add index(transition_group);


#-- delete all isotopes
#-- delete from hroest.srmPeptides_test where isotope_nr > 0;

drop table hroest.srmTransitions_test;
truncate table hroest.srmTransitions_test;
create table hroest.srmTransitions_test(
    srm_id INT PRIMARY KEY AUTO_INCREMENT,
    group_id INT,
    q3_charge TINYINT ,
    q3 DOUBLE,
    type VARCHAR(8),
    fragment_number TINYINT
);
alter table hroest.srmTransitions_test add index(group_id);
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


create table result_srmuis (
exp_key int, 
parent_key int, 
non_useable_UIS int,
total_UIS int, 
uisorder int(4)
)


ALTER TABLE hroest.srmPeptides_human ADD isotopically_modified tinyint(3) unsigned;
UPDATE hroest.srmPeptides_human SET isotopically_modified = 0;

ALTER TABLE hroest.srmPeptides_human ADD modifications tinyint(3) unsigned;
UPDATE hroest.srmPeptides_human SET modifications = 0;
ALTER TABLE hroest.srmPeptides_human ADD missed_cleavages tinyint(3) unsigned;
UPDATE hroest.srmPeptides_human SET missed_cleavages = 0;


ALTER TABLE hroest.srmPeptides_yeast ADD isotopically_modified tinyint(3) unsigned;
UPDATE hroest.srmPeptides_yeast SET isotopically_modified = 0;
ALTER TABLE hroest.srmPeptides_yeast ADD modifications tinyint(3) unsigned;
UPDATE hroest.srmPeptides_yeast SET modifications = 0;
ALTER TABLE hroest.srmPeptides_yeast ADD missed_cleavages tinyint(3) unsigned;
UPDATE hroest.srmPeptides_yeast SET missed_cleavages = 0;

ALTER TABLE hroest.srmPeptides_test ADD isotopically_modified tinyint(3) unsigned;
UPDATE hroest.srmPeptides_test SET isotopically_modified = 0;
ALTER TABLE hroest.srmPeptides_test ADD modifications tinyint(3) unsigned;
UPDATE hroest.srmPeptides_test SET modifications = 0;
ALTER TABLE hroest.srmPeptides_test ADD missed_cleavages tinyint(3) unsigned;
UPDATE hroest.srmPeptides_test SET missed_cleavages = 0;

ALTER TABLE srmPeptides_yeast_oxMetDeamid_miss1 ADD isotopically_modified tinyint(3) unsigned;
UPDATE srmPeptides_yeast_oxMetDeamid_miss1 SET isotopically_modified = 0;

ALTER TABLE srmPeptides_human_oxMetDeamid_miss1 ADD isotopically_modified tinyint(3) unsigned;
UPDATE srmPeptides_human_oxMetDeamid_miss1 SET isotopically_modified = 0;

ALTER TABLE srmPeptides_humanpepatlas_oxMetDeamid_miss0 ADD isotopically_modified tinyint(3) unsigned;
UPDATE srmPeptides_humanpepatlas_oxMetDeamid_miss0 SET isotopically_modified = 0;

ALTER TABLE srmPeptides_yeastpepatlas_oxMetDeamid_miss0 ADD isotopically_modified tinyint(3) unsigned;
UPDATE srmPeptides_humanpepatlas_oxMetDeamid_miss1 SET isotopically_modified = 0;
