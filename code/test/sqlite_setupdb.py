import sqlite
conn = sqlite.connect('/tmp/testdb')
c = conn.cursor()

table = 'srmPeptides_test'
ttable = 'srmTransitions_test'

c.execute(
"""
create table %(table)s(
    parent_id INT PRIMARY KEY ,
    peptide_key INT,
    modified_sequence VARCHAR(255),
    q1_charge TINYINT,
    q1 DOUBLE,
    ssrcalc DOUBLE,
    isotope_nr TINYINT,
    transition_group INT,
    isotopically_modified TINYINT UNSIGNED
);
create index testpepkey   on %(table)s (peptide_key);
create index testq1       on %(table)s (q1);
create index testssrcalc  on %(table)s (ssrcalc);
create index testtrgroup  on %(table)s (transition_group);
    """ % {'table' : table}
)


import csv
creader = csv.reader( open('sqltestp.out'), delimiter='\t')
#creader.next()
for rr in creader:
    query = 'insert into %s' % table + """
( parent_id , peptide_key, modified_sequence, q1_charge, q1, ssrcalc, isotope_nr, transition_group, isotopically_modified)
values (%s, %s, '%s', %s, %s, %s, %s, %s, 0)
        """ % (rr[0], rr[1], rr[2], rr[3], rr[4], rr[5], rr[6] , rr[7] )
    c.execute( query
    )

conn.commit()

c.execute(
"""
create table %(table)s (
    srm_id INT PRIMARY KEY,
    group_id INT,
    q3_charge TINYINT ,
    q3 DOUBLE,
    type VARCHAR(8),
    fragment_number TINYINT
);
create index testgi  on %(table)s (group_id);
create index testq3  on %(table)s (q3);
    """ % {'table' : ttable}
)
creader = csv.reader( open('sqltestt.out'), delimiter='\t')
#creader.next()
for rr in creader:
    query = 'insert into %s' % ttable + """
(srm_id, group_id, q3_charge, q3, type, fragment_number)
values (%s, %s, %s, %s, '%s', %s)
        """ % (rr[0], rr[1], rr[2], rr[3], rr[4], rr[5])
    c.execute( query
    )

conn.commit()



c.close()












exit()


