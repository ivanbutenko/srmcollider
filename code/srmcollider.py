#!/usr/bin/python

import MySQLdb
import sys 
sys.path.append( '/home/hroest/msa/code' )
from utils_h import utils
db = MySQLdb.connect(read_default_file="~/.my.cnf")
c = db.cursor()
c2 = db.cursor()

def get_collisions( q1, q3, ssrcalc, q1_window, q3_window, ssrcalc_window, 
                   exp_key, db):
    select_temp = """
    select * from hroest.srmClashes 
    inner join ddb.peptide on peptide.id = srmClashes.peptide_key
    where q1 between %s and %s
    and q3 between %s and %s
    and ssrcalc between %s and %s
    and srmClashes.experiment_key = %s
    ;
    """
    q1_low      =  q1 - q1_window /2.0
    q1_high     =  q1 + q1_window /2.0
    q3_low      =  q3 - q3_window /2.0   
    q3_high     =  q3 + q3_window /2.0
    ssrcalc_low =  ssrcalc - ssrcalc_window /2.0 
    ssrcalc_high=  ssrcalc + ssrcalc_window /2.0   
    ##q1_low      =   518
    ##q1_high     =   519
    ##q3_low      =   740
    ##q3_high     =   760
    ##ssrcalc_low =   20
    ##ssrcalc_high=   25
    c2 = db.cursor()
    t = utils.db_table( c2 )
    t.read( select_temp, [q1_low, q1_high, q3_low, q3_high, 
                            ssrcalc_low, ssrcalc_high, exp_key])
    #TODO implement
    #res = t.fetchall()
    res = []
    while True:
        r = t.fetchone()
        if r == None: break
        res.append( r)
    return t, res


def main(input, q1_w, q3_w, ssr_w, exp_key, db, high, low):
    #sanitize input
    import re
    seqs = "'"
    for inp in input.split():
        #only alphanumeric
        seqs += filter(str.isalnum, inp)+ "','"
    seqs = seqs[:-2]

    #get input transitions
    input_q  = """
    select * from ddb.peptide
    inner join ddb.peptideOrganism on peptideOrganism.peptide_key = peptide.id
    inner join hroest.srmClashes on srmClashes.peptide_key = peptide.id
    where peptide.experiment_key = %s
    and peptide.sequence in 
    (%s)
    """ % (exp_key, seqs)
    c = db.cursor()
    c.execute( 'use ddb;')
    t = utils.db_table( c )
    t.read( input_q )
    result = t.fetchall_groupBy('sequence')

    #print some links to csv file and input/output validation
    print "<a href ='/hroest/srmcollider.csv'>Download csv file</a>"
    print "<br/>"
    print "input: %s peptides" % (len( input.split() )) 
    print "<br/>"
    print "found: %s peptides" % (len( result )) 

    #also prepare a csv
    import csv
    f = open('/nas/www/html/hroest/srmcollider.csv', 'w')
    w = csv.writer(f, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    w.writerow( ['sequence', 'q1', 'q3', 'type'] )
    
    for sequence in result.values():
        useable = []
        collisions= []
        for row in sequence:
            q1  = t.row( row, 'q1')
            q3  = t.row( row, 'q3')
            #print q3, low
            if q3 < low or q3 > high: continue
            ssr = t.row( row, 'ssrcalc')
            occurence = t.row( row, 'genome_occurence')
            mysequence = t.row( row, 'sequence')
            tt, rr = get_collisions( q1, q3, ssr, q1_w, q3_w, ssr_w, exp_key, db)
            assert len( rr) > 0  #this would mean the sequence is not in the DB
            if len( rr) == 1: useable.append( row )
            else: 
                collisions.append( [row, [rrr for rrr in rr 
                      if tt.row(rrr, 'sequence') != mysequence]  ] )
        print "<br/>"*2, mysequence, q1, round(ssr, 2)
        print "<br/>"
        print "useable: <br/>"
        for u in useable:
            mytype = t.row(u, 'type') 
            q3 = t.row(u, 'q3') 
            print mytype + " " + str( round(q3, 2 ) ) + "<br/>"
            w.writerow( [ mysequence, q1, q3, mytype ] )

        print "<br/>"*1
        print "collisions:<br/>"
        for col in collisions:
            myrow = col[0]
            print t.row( myrow, 'type'), '::::'
            for c in col[1]:
                print tt.row(c, 'sequence') 
                print tt.row(c, 'type') 
                print round( tt.row(c, 'q3'), 2)
                print tt.row(c, 'ssrcalc') 
                print ";;"
            print "<br/>"*1 

        print "\n\n\n\n"
        print "%% useable: %d" % (len(useable)*100.0 / ( len(useable) 
            + len(collisions)) )
    f.close()

input = """
DDGSGVDIIDRPSMCLEYTTSK   
DLEILPAGDLTEIGEK         
AVGIGFIAVGIIGYAIK        
DQTSNLR                  
DQLIENCSK                
NFMPLEAK                 
DLTETCR                  
CNATPETFLQDIDK           
SLGQELIR                 
DGISQEGTIENAPANPSK       
SFLVILPAGLK              
TQMSHTLCVTTQLPR          
DIFTTLVSPELLSTCK         
GILQYVWFKPFYCFGTLICSAWK  
ENGFSMFK                 
HLQGPSDLVK               
QTHPVIIPSFASWFDISK       
IEDDLLFR                 
SGLVPAQFIEPVR            
AISSKPLSVR               
"""
q1_w = 1
q3_w = 1
ssr_w = 4
exp_key = 3120 #yeast
high = 1500
low = 300

warm_welcome = """
###########################################################################
###########################################################################
#                              SRM Collider                               #
###########################################################################
#                                  alpha                                  #
###########################################################################
"""
warm_welcome = """
###########################################################################<br/>
###########################################################################<br/>
#
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;
                              SRM Collider
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
#<br/>
###########################################################################<br/>
###########################################################################<br/>
"""

# 2+, 3+ precursor
# 1+, 2+ fragmente
# interessante faelle raussuchen => in dp liste
sample_peptides = """
AFGIPVNTFSSEVVTLWYR
AIPAPHEILTSNVVTR
VTDISTGIYK
GYSENPVENSQFLTEYVATR
ETLVGFMTEYVATR
IQDPQMTGYVSTR
ATMVGTPYWMAPEIVNQK
TNSFVGTEEYLAPEVIR
TNSFVGTEEYIAPEVIR
LINSIADTFVGTSTYMSPER
"""
sample_peptides_html = ''
for s in sample_peptides.split():
    sample_peptides_html += s + '<br/>'



print 'Content-type: text/html\n\n'
print warm_welcome

import cgi
form = cgi.FieldStorage()   # FieldStorage object to
if form.has_key('peptides'):
    peptides = form.getvalue('peptides')
    #print peptides
    #peptides = input
    main( peptides, q1_w, q3_w, ssr_w, exp_key, db, high, low)
else:
  print """
<FORM action="/cgi-bin/hroest/srmcollider.py" method="post">
    <P>
    <label for="peptides">Peptides</label><br />
    <textarea cols="60" name="peptides" rows="20"></textarea>
    <br/>
    <br/>
    <INPUT type="submit" value="Send"> 
    </P>
 </FORM>
<br/>
The following parameters are set:
<br/>
q1_window = 1 Da <br/>
q3_window = 1 Da <br/>
ssr_window = 4 units <br/>
genome = yeast<br/>
high_mass = 1500 Da <br/>
low_mass = 300 Da <br/>
peptides are fully tryptic only <br/>
<br/><br/>
To try this tool, you could use the following sample peptides:
<br/>%s    
""" % sample_peptides_html

