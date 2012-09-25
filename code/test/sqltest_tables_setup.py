import MySQLdb
import sqlite
import test_shared
import test_db, sys

if sys.argv[1] == "mysql":
    try: 
      conn = MySQLdb.connect(read_default_file=test_db.mysql_conf_file)
      c = conn.cursor()
      c.connection.autocommit(True)
    except MySQLdb.OperationalError as e:
        print "Could not connect to database: Please check the configuration in test/test_db.py!\n", e
        sys.exit()
elif sys.argv[1] == "sqlite":
    conn = sqlite.connect(test_shared.SQLITE_DATABASE_LOCATION)
    c = conn.cursor()
else:
    print "First argument must be either 'mysql' or 'sqlite'"
    sys.exit()

table = 'srmPeptides_test'

try:
    c.execute("DROP TABLE %s" % table)
except MySQLdb.OperationalError as e:
    # If the table doesnt exist, mysql throws an error. We cannot use IF EXISTS
    # in the SQL query because Sqlite doesnt understand that!
    pass
except sqlite.DatabaseError as e:
    pass

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

    modifications TINYINT UNSIGNED,
    missed_cleavages TINYINT UNSIGNED,
    isotopically_modified TINYINT UNSIGNED
); """ % {'table' : table})
c.execute("create index testpepkey   on %(table)s (peptide_key);" % {'table': table})
c.execute("create index testq1       on %(table)s (q1)" % {'table': table})
c.execute("create index testssrcalc  on %(table)s (ssrcalc)" % {'table': table})
c.execute("create index testtrgroup  on %(table)s (transition_group)" % {'table': table})

import csv
creader = csv.reader( open('sqltestp.out'), delimiter='\t')
#creader.next()
for rr in creader:
    # skip all that have isotope_nr other than zero (we used to store higher
    # isotopes explicitely but dont do this any more)
    if int(rr[6]) != 0: continue
    query = 'insert into %s' % table + """
( parent_id , peptide_key, modified_sequence, q1_charge, q1, ssrcalc, isotope_nr, transition_group, isotopically_modified, modifications, missed_cleavages)
values (%s, %s, '%s', %s, %s, %s, %s, %s, 0, 0, 0)
        """ % (rr[0], rr[1], rr[2], rr[3], rr[4], rr[5], rr[6] , rr[7] )
    c.execute(query)
conn.commit()

c.close()
exit()


