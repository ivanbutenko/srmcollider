#!/usr/bin/python

import cgitb; cgitb.enable()

import MySQLdb, time
import sys 
sys.path.append( '/home/hroest/projects/msa/code' )
sys.path.append( '/home/hroest/projects/srm_clashes/code' )
from utils_h import utils
#db = MySQLdb.connect(read_default_file=".my.cnf")
#db = MySQLdb.connect("orl.ethz.ch", "hroest", "", db="hroest")
db = MySQLdb.connect("orl.ethz.ch", "srmcollider", "srmcollider", db="hroest")
c = db.cursor()
c2 = db.cursor()
import collider

#myCSVFile = '/nas/www/html/hroest/srmcollider.csv'
myCSVFile = '/var/www/documents/srmcollider.csv'
myUIS_CSVFile = '/var/www/documents/uis_srmcollider.csv'
myUIS_CSVFile_rel = '/../documents/uis_srmcollider.csv'
myCSVFile_rel = '/../documents/srmcollider.csv'

#def get_collisions( q1, q3, ssrcalc, q1_window, q3_window, ssrcalc_window, 
#                   exp_key, db):
#    select_temp = """
#    select * from hroest.srmClashes 
#    inner join ddb.peptide on peptide.id = srmClashes.peptide_key
#    where q1 between %s and %s
#    and q3 between %s and %s
#    and ssrcalc between %s and %s
#    and srmClashes.experiment_key = %s
#    ;
#    """
#    q1_low      =  q1 - q1_window /2.0
#    q1_high     =  q1 + q1_window /2.0
#    q3_low      =  q3 - q3_window /2.0   
#    q3_high     =  q3 + q3_window /2.0
#    ssrcalc_low =  ssrcalc - ssrcalc_window /2.0 
#    ssrcalc_high=  ssrcalc + ssrcalc_window /2.0   
#    ##q1_low      =   518
#    ##q2_high     =   519
#    ##q3_low      =   740
#    ##q3_high     =   760
#    ##ssrcalc_low =   20
#    ##ssrcalc_high=   25
#    c2 = db.cursor()
#    t = utils.db_table( c2 )
#    t.read( select_temp, [q1_low, q1_high, q3_low, q3_high, 
#                            ssrcalc_low, ssrcalc_high, exp_key])
#    #TODO implement
#    #res = t.fetchall()
#    res = []
#    while True:
#        r = t.fetchone()
#        if r == None: break
#        res.append( r)
#    return t, res


def main(input, q1_w, q3_w, ssr_w, exp_key, db, high, low, genome, isotope, uis):

    #sanitize input
    import re
    seqs = "'"
    for inp in input.split():
        #only alphanumeric
        seqs += filter(str.isalnum, inp)+ "','"

    seqs = seqs[:-2]

    #figure out which db to use 
    db_used = 'hroest'
    if genome == 'yeast':
        table_used =  'yeast'
    elif genome == 'yeastN15':
        table_used =  'yeastN15'
    elif genome == 'human':
        table_used =  'human'
    else: print "Genome not recognized"; exit()
    #create the parameter object
    par = collider.SRM_parameters()
    par.q1_window = q1_w / 2.0
    par.q3_window = q3_w / 2.0
    par.ssrcalc_window = ssr_w / 2.0 
    par.ppm = False
    par.considerIsotopes = True
    par.q3_range = [low, high]
    par.peptide_table = db_used + '.srmPeptides_' + table_used
    par.transition_table = db_used + '.srmTransitions_' + table_used
    par.eval()
    extra_table = 'ddb.peptide'
    mycollider = collider.SRMcollider()

    if uis > 0:
        print "<p>Calculated UIS for peptide %s" % input.split()[0]
        print "<a href ='%s'>Download csv file with UIS</a></p>" % myUIS_CSVFile_rel
        if uis >  5 or len( input.split() ) > 1:
            print "Can only calculate up to order 5 and only 1 peptide"
            exit()

    

    #get input transitions
    input_q  = """
    select parent_id, q1, q1_charge, ssrcalc, b.id, b.sequence
    from %(ptable)s a
    inner join %(et)s b on a.peptide_key = b.id
    where b.sequence in 
    (%(sequences)s)
    %(query_add)s
    """ % { 'ptable' : par.peptide_table, 'sequences' : seqs, 'et' : extra_table, 
           'query_add' : par.query_add  }
    cursor = db.cursor()
    cursor.execute(input_q)
    res = cursor.fetchall()
    result = [
        {
            'parent_id' :  r[0],
            'q1' :         r[1],
            'q1_charge' :  r[2],
            'ssrcalc' :    r[3],
            'peptide_key' :r[4],
            'sequence'    :r[5]
        }
        for r in res
    ]

    #lets see whether we have any results:

    if len( result ) == 0:
        print "Unfortunately, none of the peptides were found in the database."
        print "<a href ='srmcollider.py'>Return to main page.</a>" 
        exit()

    #print some links to csv file and input/output validation
    #print "<a href ='/../documents/srmcollider.csv'>Download csv file</a>"
    unique = utils.unique(input.split() )
    print "<a href ='%s'>Download csv file</a>" % myCSVFile_rel
    print "<br/>"
    print "input: %s peptides" % (len( input.split() )) 
    print "<br/>"
    print "unique: %s peptides" % (len( unique )) 
    print "<br/>"
    print "found: %s peptides" % (len( result )) 
    print "<div class='toc'><ul>"
    for u in result: print '<li><a href="#%s">%s</a></li>' % (u['sequence'],u['sequence'])
    print "</ul></div>"


    #also prepare a csv
    import csv
    f = open( myCSVFile, 'w')
    w = csv.writer(f, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    w.writerow( ['sequence', 'q1', 'q3', 'type'] )

    #do the part for UIS
    if uis > 0:
        #prepare another csv
        import csv
        fuis = open( myUIS_CSVFile, 'w')
        wuis = csv.writer(fuis, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        MAX_UIS = uis
        pep = result[0]
        transitions = mycollider._get_all_transitions(par, pep, cursor, values="q3, srm_id, type, fragment_number")
        collisions = mycollider._get_all_collisions(par, pep, cursor, 
            values="q3, q1, srm_id, peptide_key, type, fragment_number, modified_sequence, ssrcalc, isotope_nr")
        non_uis_list = [set() for i in range(MAX_UIS+1)]
        collisions_per_peptide = {}
        q3_window_used = par.q3_window
        for t in transitions:
            #if par.ppm: q3_window_used = par.q3_window * 10**(-6) * t[0]
            this_min = q3_window_used
            for c in collisions:
                if abs( t[0] - c[0] ) <= q3_window_used:
                    if collisions_per_peptide.has_key(c[3]):
                        if not t[1] in collisions_per_peptide[c[3]]:
                            collisions_per_peptide[c[3]].append( t[1] )
                    else: collisions_per_peptide[c[3]] = [ t[1] ] 
        #here we calculate the UIS for this peptide with the given RT-range
        for pepc in collisions_per_peptide.values():
            for i in range(1,MAX_UIS+1):
                collider.get_non_uis(pepc, non_uis_list[i], i)
        #we could actually calculate the UIS combinations 
        #from the non-UIS combinations
        srm_ids = [t[1] for t in transitions]
        srm_lookup = [ (t[1] , t) for t in transitions]
        srm_lookup = dict(srm_lookup) 
        #start = time.time()
        for i in range(1,MAX_UIS+1):
            uis_list = collider.get_uis(srm_ids, non_uis_list[i], i)
            #end = time.time()
            for comb in uis_list:
                tmp = [ srm_lookup[elem] for elem in comb]  
                myrow = []
                for tt in tmp:
                    myrow.extend( [ tt[0], tt[2] + str(tt[3] )]  )
                wuis.writerow(myrow)

        fuis.close()

    for pep in result:
        useable = []
        collisions= []
        non_unique = {}
        mysequence = pep['sequence']
        transitions = mycollider._get_all_transitions(par, pep, cursor, values="q3, srm_id, type, fragment_number")
        collisions = mycollider._get_all_collisions(par, pep, cursor, 
            values="q3, q1, srm_id, peptide_key, type, fragment_number, modified_sequence, ssrcalc, isotope_nr")
        q3_window_used = par.q3_window
        for t in transitions:
            #if par.ppm: q3_window_used = par.q3_window * 10**(-6) * t[0]
            this_min = q3_window_used
            for c in collisions:
                if abs( t[0] - c[0] ) <= this_min:
                    try: 
                        non_unique[ t[1] ][1].append( c )
                    except:
                        non_unique[ t[1] ] = [t, [c]]
        useable = [t for t in transitions if not t[1] in non_unique]
        #
        #
        print "<h2 id='%s'>Peptide %s</h2>" % (mysequence, mysequence)
        print "<div class='pep_property'>"
        print "<p>Q1: %s</p>" % pep['q1']
        print "<p>SSRCalc: %s</p>" % round(pep['ssrcalc'], 2)
        print "<p>Percent Useable: %d</p>" % (len(useable)*100.0 / len(transitions)  )
        print "</div>"

        if len( useable ) > 0:
            print "<h3>Useable transitions</h3>"
            print "<table class='tr_table'>"
            print "<tr><td>Type</td> <td>Q3</td> </tr>"
            for u in useable:
                print "<tr><td>"
                mytype = u[2] + str(u[3])
                print mytype
                print "</td><td>"
                print str( round( u[0], 2 ) ) 
                #q1_charge = t.row(u, 'q1_charge') 
                #q1 = t.row(u, 'q1') 
                #q3 = t.row(u, 'q3') 
                #print "(" + str( q1_charge ) + ")" + mytype + " " + " " +\
                #                     str( round(q3, 2 ) ) + "<br/>"
                w.writerow( [ mysequence, pep['q1'], u[0], mytype ] )
                print "</td></tr>"
            print "</table>"
        else: 
            print "<p>No useable transitions for this peptide!</p>"
            w.writerow( ['"Sorry: no useable transitions for this peptide. Try to calculate UIS."' ] )

        print "<h3>Unuseable transitions (Collisions)</h3>"
        print "<table class='col_table'>"
        for u, coll in non_unique.values():
            print "<tr><td>"
            mytype = u[2] + str(u[3])
            print mytype
            print str( round( u[0], 2 ) ) 
            print "</td><td>"
            #myrow = col[0]
            #print "(" + str(t.row( myrow, 'q1_charge')) + ")", \
            #      t.row( myrow, 'type'), \
            #        '::::'
            for c in coll:
                print '(', round(c[1], 2), round(c[0], 2),  ')'
                print c[7]
                print c[4] + str(c[5])
                if c[8] != 0:
                    print c[6], 'isotope', c[8], ' </br>'
                else: print c[6], ' </br>'
                #print tt.row(c, 'sequence') 
                #print tt.row(c, 'type') 
                #print round( tt.row(c, 'q1'), 2)
                #print round( tt.row(c, 'q3'), 2)
                #print tt.row(c, 'ssrcalc') 
                #print ";;"
            print "</td></tr>"
        print "</table>"


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

warm_welcome = """
<div class="top">
    <div class="header">
        SRM Collider
    </div>
    <div class="version">
        version 0.1
        </br>
        alpha
        </br>
        Hannes Roest 2010
    </div>
    <div class="navigation">
        <span class="active-nav"> 
            <a href="srmcollider.py">Collider </a>
        </span>
        <span class="inactive-nav"> <a href="download.html">Download </a> </span>
        <span class="inactive-nav"> <a href="about.html">About </a> </span>
    </div>
</div>
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




###########################################################################
###########################################################################
# START OF HTML


print 'Content-type: text/html\n\n'
print """
<!DOCTYPE html PUBtdC "-//W3C//DTD XHTML 1.0 Transitional//EN"
       "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"> 
 
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en"> 
<head> 
  <meta http-equiv="content-type" content="text/html;charset=UTF-8" /> 
  <title>SRM Collider</title> 
  <link href="/stylesheets/srmcollider.css" media="screen" rel="stylesheet" type="text/css" /> 
  <link href="/stylesheets/%s.css" media="screen" rel="stylesheet" type="text/css" /> 
</head> 

<body>
<div class="whole">
""" % 'brown'
#options for the css are
## first
## second
## brown (quite nice)
## ocker 
## green
print warm_welcome
print "<div class='main'>"

import cgi
form = cgi.FieldStorage()   # FieldStorage object to
if form.has_key('peptides'):
    peptides = form.getvalue('peptides')
    q1_w = float(form.getvalue('q1_window') )
    q3_w = float(form.getvalue('q3_window') )
    ssr_w = float(form.getvalue('ssr_window') )
    high = float(form.getvalue('high_mass') )
    low = float(form.getvalue('low_mass') )
    genome = form.getvalue('genome') 
    #TODO 
    #isotope = int(form.getvalue('isotope') )
    isotope = 3
    uis = int(form.getvalue('uis') )
    #print peptides
    #peptides = input
    start = time.time()
    main( peptides, q1_w, q3_w, ssr_w, exp_key, db, high, low, genome, isotope, uis)
    print "<hr> <br/>This query took: %s s" % (time.time() - start)
else:
  print """
<form action="/cgi-bin/srmcollider.py" method="post">
    <p class='input_field'>
        <label for="peptides">Please enter the peptide sequences here:</label><br />
        <textarea id="pep_input" name="peptides" rows="20"></textarea>
    </p>

    <p class='input_field'>
        <label class="mylabel" for="ssr_window">SSRCalc window</label>
        <input class="number_input" type="text" name="ssr_window" value="4"> arbitrary units
    </p>

    <p class='input_field'>
        <label class="mylabel" for="q1_window">Q1 mass window</label>
        <input class="number_input" type="text" name="q1_window" value="0.7"> Th
    </p>

    <p class='input_field'>
        <label class="mylabel" for="q3_window">Q3 mass window</label>
        <input class="number_input" type="text" name="q3_window" value="1.0"> Th
    </p>

    <p class='input_field'>
        <label class="mylabel" for="low_mass">Low mass threshold for transitions</label>
        <input class="number_input" type="text" name="low_mass" value="300"> Th
    </p>

    <p class='input_field'>
        <label class="mylabel" for="high_mass">High mass threshold for transitions</label>
        <input class="number_input" type="text" name="high_mass" value="1500"> Th
    </p>

    <p class='input_field'>
        <label class="mylabel" for="genome">Genome</label>
        <select name="genome" class="number_input">
          <option value="yeast">Yeast (tryptic)</option>
          <option value="yeastN15">Yeast N15 (tryptic)</option>
          <option value="human">Human (tryptic)</option>
        </select>
    </p>

    <p class='input_field'>
        <label class="mylabel" for="isotope">Consider isotopes up to </label>
        <!-- >
        <input class="number_input" type="text" name="isotope" value="3"> amu
        </!-->
        <input class="number_input" type="text" disabled="False" name="isotope" value="3"> amu
    </p>


    <p class='input_field'>
        <label class="mylabel" for="uis">Find UIS up to order* </label>
        <input class="number_input" type="text" name="uis" value="0"> 
    </p>

    <INPUT type="submit" value="Send"> 

 </form>
* due to computational constraints, only one peptide at a time can be searched in UIS mode

<!-- 
q1_window = 1 Da <br/>
q3_window = 1 Da <br/>
ssr_window = 4 units <br/>
high_mass = 1500 Da <br/>
low_mass = 300 Da <br/>
genome = yeast<br/>
peptides are fully tryptic only <br/>
<br/><br/>
-->
<!-- >
To try this tool, you could use the following sample peptides:
<br/>%s    
</!-->
""" % sample_peptides_html


print "</div>"
print """
</div>
</body>
</html>"""
