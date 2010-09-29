f = open('/home/hroest/tmp.sql', 'w')
f.write( """
DROP TABLE hroest.orinner_scores;
CREATE TABLE hroest.orinner_scores (
id INT AUTO_INCREMENT PRIMARY KEY, 
file_id INT, 
run_id TEXT,
assays_id INT,
pg_rank INT, 
Tr DOUBLE,
decoy INT,
group_id TEXT,
pepseq VARCHAR(255) );
        """)

text = """
LOAD DATA LOCAL INFILE
'/home/hroest/data/A_D100819_SPLASMAVER_TRID1521_%s_0%s_scores.xls'
INTO TABLE hroest.orinner_scores 
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
IGNORE 1 LINES 
( 
run_id,
assays_id,
group_id, 
pepseq,
@dummy,
decoy,
@dummy,
@dummy,
@dummy,
pg_rank,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
    Tr)
SET file_id = %s;
"""

for i in range(1,28):
    if i != 12: f.write( text % (i, 1, i) )
    else: f.write( text % (i, 2, i) )


f.close()




bigtable = \
"""

DROP TABLE hroest.orinner_peakgroups;
CREATE TABLE hroest.orinner_peakgroups (
id INT AUTO_INCREMENT PRIMARY KEY, 
run_id TEXT,
assays_id INT,
pg_rank INT, 
Tr_min DOUBLE,
Tr DOUBLE,
decoy INT,
group_id TEXT,
pepseq VARCHAR(255) ,
discrimination_score double,
m_score double
);

#decoy = 0 means true
LOAD DATA LOCAL INFILE
'/home/hroest/data/mInteract_peakgroups_mod.xls'
INTO TABLE hroest.orinner_peakgroups 
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
IGNORE 1 LINES 
( 
@dummy,
decoy,
Tr_min,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
run_id,
assays_id,
group_id, 
pepseq,
@dummy,
@dummy,
@dummy,
@dummy,
pg_rank,
@dummy,
@dummy,
@dummy,
    Tr,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
discrimination_score,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
@dummy,
m_score
);

"""






select sequence from ddb.peptide 
inner join hroest.orinner_peakgroups on peptide.sequence = orinner_peakgroups.pepseq
where experiment_key = 3394




select sequence, Tr from ddb.peptide 
inner join hroest.orinner_scores on peptide.sequence = orinner_scores.pepseq
where experiment_key = 3394
order by sequence, Tr

