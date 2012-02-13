DROP TABLE hroest.yeast_dp_data_original;
CREATE TABLE hroest.yeast_dp_data_original (
id INT AUTO_INCREMENT PRIMARY KEY, 
Q1 DOUBLE, 
Q3 DOUBLE,
Tr_recalibrated DOUBLE, 
Tr_original DOUBLE, 
transition_name TEXT, 
protein_name    VARCHAR(255), 
isotype VARCHAR(255), 
CE DOUBLE, 
relative_intensity  INT,
rank    INT, 
mysql_assay_id  INT, 
mysql_transitions_id    INT,
mysql_peptides_id INT,
stripped_sequence VARCHAR(255),
modification VARCHAR(255), 
PI      DOUBLE, 
spectrum_score  VARCHAR(255),
score_type      VARCHAR(255), 
machine VARCHAR(255), 
JPT_plate_id    VARCHAR(255),
row     VARCHAR(255),
col     VARCHAR(255),
redundancy      INT,
original_spec_name      VARCHAR(255),
dummy   INT,
prec_z  INT,
frg_type VARCHAR(1),
loss_type       VARCHAR(255),
frg_z   INT,
frg_nr  INT,
transition_group_id     VARCHAR(255),
decoy   INT,
mymod VARCHAR(255)
);


LOAD DATA LOCAL INFILE
'/home/hroest/srm_clashes/data/YEAST_spectrast_DP_pep.xls'
#'/home/hroest/srm_clashes/data/test.xls'
INTO TABLE hroest.yeast_dp_data_original 
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
IGNORE 1 LINES 
( Q1 ,
Q3 ,
Tr_recalibrated ,
Tr_original ,
transition_name ,
protein_name    ,
isotype ,
CE ,
relative_intensity  ,
rank    ,
mysql_assay_id  ,
mysql_transitions_id    ,
mysql_peptides_id ,
stripped_sequence ,
modification ,
PI      ,
spectrum_score  ,
score_type      ,
machine ,
JPT_plate_id    ,
row     ,
col     ,
redundancy      ,
original_spec_name      ,
dummy   ,
prec_z  ,
frg_type ,
loss_type       ,
frg_z   ,
frg_nr  ,
transition_group_id     ,
decoy   ,
mymod )
;


SHOW WARNINGS limit 200;


DROP TABLE hroest.yeast_dp_data_light;
CREATE TABLE hroest.yeast_dp_data_light as
SELECT
id,
Q1,
Q3,
Tr_original,
CE,
relative_intensity  ,
rank,
stripped_sequence,
modification as modified_sequence,
prec_z  ,
frg_type ,
frg_z   ,
frg_nr 
from hroest.yeast_dp_data_original
where isotype = 'light';
alter table hroest.yeast_dp_data_light add primary key (id);
alter table hroest.yeast_dp_data_light modify id int auto_increment ;
alter table hroest.yeast_dp_data_light add index(stripped_sequence);
alter table hroest.yeast_dp_data_light add index(modified_sequence);


#how many occur in both datasets? 17893 sequences out of a total of 18240
select count(*) 
from hroest.yeast_dp_data_light
inner join hroest.srmPeptides_yeast_2
on srmPeptides_yeast_2.modified_sequence  = yeast_dp_data_light.modified_sequence
group by srmPeptides_yeast_2.modified_sequence;

create table tmp_all_seqs as
select distinct modified_sequence
from hroest.yeast_dp_data_light;
alter table tmp_all_seqs add index(modified_sequence);

create table tmp_all_seqs2 as
select distinct modified_sequence
from hroest.srmPeptides_yeast_2;
alter table tmp_all_seqs2 add index(modified_sequence);

create table tmp_not_found_seqs as
select distinct modified_sequence
from hroest.yeast_dp_data_light
where modified_sequence not in (select modified_sequence 
from hroest.srmPeptides_yeast_2
);



drop table tmp_not_found_seqs;

create table tmp_not_found_seqs as
select tmp_all_seqs2.modified_sequence as mod2, hroest.tmp_all_seqs.modified_sequence as mod1
from hroest.tmp_all_seqs
left join hroest.tmp_all_seqs2 
on tmp_all_seqs2.modified_sequence  = tmp_all_seqs.modified_sequence
;
delete from tmp_not_found_seqs where mod2 is NOT NULL;

select count(*) from tmp_not_found_seqs;
select mod1 from tmp_not_found_seqs;




#how many peptides do we have 192794
select count(*) from ddb.peptide
where experiment_key = 3131
and length( peptide.sequence ) > 1
;

#how many peptides have ssrcalc predictions?
#192792
select distinct peptide.sequence
from ddb.peptide 
inner join compep.ssrcalc_prediction on ssrcalc_prediction.sequence =
peptide.sequence
where experiment_key = 3131
and length( peptide.sequence ) > 1









"""

select  prec_z, parent_id, srmPeptides_yeast_allCAM.q1 as q1, q1_charge,
ssrcalc          
#select *
from hroest.srmPeptides_yeast_allCAM          inner join
ddb.peptide on peptide.id = srmPeptides_yeast_allCAM.peptide_key      inner
join ddb.peptideOrganism on peptide.id = peptideOrganism.peptide_key
inner join hroest.yeast_dp_data_light  on
yeast_dp_data_light.modified_sequence =
srmPeptides_yeast_allCAM.modified_sequence          
where genome_occurence = 1
and q1_charge = prec_z  
group by parent_id, q1_charge
limit 10
;


select * from srmPeptides_yeast_allCAM where parent_id = 122;
select * from yeast_dp_data_light where modified_sequence = 'ENLMYQASVFNQDVGK';



and q1_charge = 2         
and prec_z = 3 limit 250;



"""





###################################
# A minus B
SELECT DISTINCT a.member_id, a.name
FROM a LEFT JOIN b USING (member_id, name)
WHERE b.member_id IS NULL

#ddb.peptide minus 
SELECT DISTINCT peptide.sequence 
FROM ddb.peptide 
LEFT JOIN compep.ssrcalc_prediction ON ssrcalc_prediction.sequence = peptide.sequence
WHERE ssrcalc_prediction.sequence IS NULL
and experiment_key = 3130
;


select TABLE_SCHEMA,TABLE_NAME, 
concat(round(data_length/(1024*1024),2),'M') DATA,
concat(round(INDEX_LENGTH/(1024*1024),2),'M') IDX
from information_schema.TABLES 
where TABLE_SCHEMA = 'hroest'
#where TABLE_SCHEMA = 'saccharomyces_cerevisiae_core_57_1j'
order by data_length ;




#big join over srm stuff and dirty peptides
select parent_id, peptide_key, q1_charge, srmPeptides_yeast_allCAM.q1, 
yeast_dp_data_light.Q1, srm_id, q3_charge ,
srmTransitions_yeast_allCAM.q3, yeast_dp_data_light.Q3, relative_intensity, rank,
frg_type, type, frg_z, frg_nr, fragment_number as frag_nr,
srmPeptides_yeast_allCAM.modified_sequence as sequence, Tr_original, ssrcalc
from hroest.srmPeptides_yeast_allCAM
inner join hroest.srmTransitions_yeast_allCAM
ON parent_key = parent_id
inner join hroest.yeast_dp_data_light ON 
yeast_dp_data_light.modified_sequence = srmPeptides_yeast_allCAM.modified_sequence
where 
    prec_z = q1_charge
and frg_type = type
and frg_z = q3_charge
and frg_nr = fragment_number
limit 40;


#to correleate ssrcalc with RT
select distinct
Tr_original, ssrcalc
from hroest.srmPeptides_yeast_allCAM
inner join hroest.srmTransitions_yeast_allCAM
ON parent_key = parent_id
inner join hroest.yeast_dp_data_light ON 
yeast_dp_data_light.modified_sequence = srmPeptides_yeast_allCAM.modified_sequence
where 
    prec_z = q1_charge
and frg_type = type
and frg_z = q3_charge
and frg_nr = fragment_number
;




















LOAD DATA LOCAL INFILE
'/home/hroest/Found_mixAll.csv'
INTO TABLE ddb.peptide 
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
( sequence )
SET experiment_key = 3352,
peptide_type = 'norm';



LOAD DATA LOCAL INFILE
'/tmp/tmp.out'
INTO TABLE hroest.srmPeptides_yeast 
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
( sequence )
SET experiment_key = 3352,
peptide_type = 'norm';




DROP TABLE hroest.openmsModel;
CREATE TABLE hroest.openmsModel (
sequence    VARCHAR(255), 
openmsrt DOUBLE
);


LOAD DATA LOCAL INFILE
'/tmp/modeltestout'
INTO TABLE hroest.openmsModel 
FIELDS TERMINATED BY ' '
LINES TERMINATED BY '\n'


create table hroest.tmp_seqs (sequence varchar(255) );
LOAD DATA LOCAL INFILE
'/tmp/sequences3456'
INTO TABLE hroest.tmp_seqs 
FIELDS TERMINATED BY ' '
LINES TERMINATED BY '\n'






truncate TABLE hroest.manndata_lfq_filemapping ;
LOAD DATA LOCAL INFILE
'/tmp/filemapping.csv'
INTO TABLE hroest.manndata_lfq_filemapping 
FIELDS TERMINATED BY ','
LINES TERMINATED BY '\n'
IGNORE 1 LINES
;



create TABLE hroest.tppsearchresults
(id int primary key auto_increment,
 exp_key int, 
 otherid int,
 name text, 
 rt double, 
 someint int, 
 someotherint int,
 sequence varchar(255),
 mybefore varchar(1), 
 myafter varchar(1), 
 protein varchar(255), 
 morert double, 
 d1 double, 
 d2 double, 
 d3 double, 
 d4 double, 
 d5 double, 
 d6 double
)
alter table hroest.tppsearchresults add index(exp_key);
alter table hroest.tppsearchresults add index(peptide);


#50201 - 50205
'/home/hroest/data/searchresults/2010Nov/CS_PUB_TPP_v05_hroest_1289397066450.csv'
'/home/hroest/data/searchresults/2010Nov/CS_PUB_TPP_v05_hroest_1289412551195.csv'
'/home/hroest/data/searchresults/2010Nov/CS_PUB_TPP_v05_hroest_1289423320075.csv'
'/home/hroest/data/searchresults/2010Nov/CS_PUB_TPP_v05_hroest_1289423680212.csv'
'/home/hroest/data/searchresults/2010Nov/CS_PUB_TPP_v05_hroest_1289424056378.csv'


LOAD DATA LOCAL INFILE
'/home/hroest/data/searchresults/2010Nov/CS_PUB_TPP_v05_hroest_1289424056378.csv'
INTO TABLE hroest.tppsearchresults
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
(
 otherid ,
 name ,
 rt ,
 someint ,
 someotherint ,
 sequence ,
 mybefore ,
 myafter ,
 protein ,
 morert ,
 d1 ,
 d2 ,
 d3 ,
 d4 ,
 d5 ,
 d6 
)
SET exp_key = 50205
;



drop table tmp_tppsearch_analysis ; 
create table tmp_tppsearch_analysis as 
select name, sequence,
@common_name := LEFT(name,  INSTR( name, '.') - 1) as common_name,
# length(@common_name) - locate("_", reverse(@common_name)) as acs,
@contfrac := LEFT(@common_name, length(@common_name) - locate("_", reverse(@common_name)) ) as contains_frag,
RIGHT(@contfrac, locate("_", reverse(@contfrac))-1  ) as fraction , 
#
LEFT( 
       RIGHT(@common_name, locate("_", reverse(@common_name))-1 ) , 
INSTR( RIGHT(@common_name, locate("_", reverse(@common_name))-1 ), '-')-1 )
as sds,
RIGHT(@common_name,  1) as replicate
from hroest.tppsearchresults 
# where sequence = 'AAAAAAAAK'
#limit 10;
group by fraction, sds, replicate, sequence
;

alter table tmp_tppsearch_analysis add index(fraction);
alter table tmp_tppsearch_analysis add index(sds);
alter table tmp_tppsearch_analysis add index(replicate);
alter table tmp_tppsearch_analysis add index(sequence);

select * from tmp_tppsearch_analysis  a
inner join manndata_lfq_filemapping map on  map.file = a.common_name
limit 1;

drop table tmp_tppsearch_analysis_2 ;
create table tmp_tppsearch_analysis_2 as
select count(*) as obs, sequence, fraction, sds, map.id as file
from tmp_tppsearch_analysis a
inner join manndata_lfq_filemapping map on map.file = CONCAT(CONCAT(a.contains_frag,'_'), sds)
group by sequence, fraction, sds;

alter table tmp_tppsearch_analysis_2 add index(fraction);
alter table tmp_tppsearch_analysis_2 add index(file);
alter table tmp_tppsearch_analysis_2 add index(sds);
alter table tmp_tppsearch_analysis_2 add index(sequence);

select * from  tmp_tppsearch_analysis_2 limit 10;






'NDLLANIVLNSTAFENR',
'YVFGLEFLR',
'TGTLTENVMTVVR',
'LSLGLQPGELQYLR',
'TETALLSLAR',
'LFGYESNSLFK',
'GLILDGLLGIQDPLR',
'LLVETLK',
'SFLQLVWAAFNDK'




create temporary table hroest.ruthspeps (sequence varchar(255) );
LOAD DATA LOCAL INFILE
'/home/hroest/data/ruthspeptides'
INTO TABLE hroest.ruthspeps 
FIELDS TERMINATED BY ' '
LINES TERMINATED BY '\n';

select count(*) from hroest.ruthspeps;
select count(distinct gene_key) 
from hroest.ruthspeps
inner join ddb.peptide on peptide.sequence = ruthspeps.sequence 
inner join ddb.protPepLink l on l.peptide_key = peptide.id
inner join ddb.geneProtLink ll on l.protein_key = ll.protein_key
where experiment_key = 3130;

# peptides
#3995 in human
#1558 in mouse

#proteins 
# 1312 in human
# 842 in mouse




drop table hroest.tca_peps ;
create temporary table hroest.tca_peps (sequence varchar(255) );
LOAD DATA LOCAL INFILE
'/home/hroest/data/tca_peptides.csv'
INTO TABLE hroest.tca_peps 
FIELDS TERMINATED BY ' '
LINES TERMINATED BY '\n';


drop table hroest.tca_peps_u ;
create temporary table hroest.tca_peps_u as 
select distinct sequence from hroest.tca_peps;
select * from hroest.tca_peps_u ;



select u.sequence, o.genome_occurence from hroest.tca_peps_u  u 
inner join ddb.peptide p on u.sequence = p.sequence 
inner join ddb.peptideOrganism o on o.peptide_key = p.id
where experiment_key = 3131
and genome_occurence = 1
;



select p.sequence, o.genome_occurence from ddb.peptide p 
inner join ddb.peptideOrganism o on o.peptide_key = p.id
where experiment_key = 3131
and sequence = 'SLVPNIPFQMLLR'
;
