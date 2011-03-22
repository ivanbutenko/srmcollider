import sqlite
conn = sqlite.connect('/tmp/example')
c = conn.cursor()

c.execute(
"""
create table srmPeptides_yeast(
    parent_id INT PRIMARY KEY ,
    peptide_key INT,
    modified_sequence VARCHAR(255),
    q1_charge TINYINT,
    q1 DOUBLE,
    ssrcalc DOUBLE,
    isotope_nr TINYINT,
    transition_group INT
);
create index pepkey   on srmPeptides_yeast (peptide_key);
create index q1       on srmPeptides_yeast (q1);
create index ssrcalc  on srmPeptides_yeast (ssrcalc);
create index trgroup  on srmPeptides_yeast (transition_group);
"""
)


import csv
#parent_id   peptide_key modified_sequence   q1_charge   q1  ssrcalc isotope_nr  transition_group
creader = csv.reader( open('/tmp/yeast_table.out'), delimiter='\t')
creader.next()
for rr in creader:
    c.execute(
        """
insert into srmPeptides_yeast 
(modified_sequence, peptide_key, parent_id, q1_charge, q1, ssrcalc, isotope_nr)
values ('%s', %s, %s, %s, %s, %s, %s)
        """ % (rr[0], rr[1], rr[2], rr[3], rr[4], rr[5], rr[6]       )
    )

conn.commit()

c.execute("select count(*) from srmPeptides_yeast")
c.fetchall()

c.close()
