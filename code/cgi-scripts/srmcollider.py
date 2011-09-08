#!/usr/bin/python
# -*- coding: utf-8  -*-

"""
 *
 * Program       : SRMCollider
 * Author        : Hannes Roest <roest@imsb.biol.ethz.ch>
 * Date          : 05.02.2011 
 *
 *
 * Copyright (C) 2011 Hannes Roest
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307, USA
 *
"""

import MySQLdb, time, sys 
import csv, re
import cgitb; cgitb.enable()
import cgi

sys.path.append( '/home/hroest/projects/msa/code' )
sys.path.append( '/home/hroest/projects/srm_clashes/code' )
import DDB; from Residues import Residues
import collider
import c_getnonuis
import sharedhtml as shared

# some options that can be changed locally for your convenience
default_mysql = "/IMSB/users/hroest/.srm.cnf"
db_used = 'hroest'
default_org_prefix = '.srmPeptides_'
genomes_that_require_N15_data = ['yeastN15']
default_ssrcalc = 'hroest.ssrcalc_pr_copy'
ssrcalc_path = '/tmp/' 

#myCSVFile = '/nas/www/html/hroest/srmcollider.csv'
myCSVFile = '/var/www/documents/srmcollider.csv'
myUIS_CSVFile = '/var/www/documents/uis_srmcollider.csv'
myUIS_CSVFile_rel = '/../documents/uis_srmcollider.csv'
myCSVFile_rel = '/../documents/srmcollider.csv'

# If you want to add additional genomes, edit the genome_select HTML and the
# map_db_tables function
genome_select = """
    <option value="yeast">Yeast (tryptic)</option>
    <option value="yeastN15">Yeast N15 (tryptic)</option>
    <option value="mouse">Mouse (tryptic)</option>
    <option value="human">Human (tryptic)</option>
"""

def map_db_tables(genome):
    # figure out which db to use 
    # you will have to change that if you have a different setup
    if genome == 'yeast':
        table_used =  'yeast'
    elif genome == 'yeastN15':
        table_used =  'yeastN15'
    elif genome == 'mouse':
        table_used =  'mouse'
    elif genome == 'human':
        table_used =  'human'
    else: 
        print "Genome not recognized";
        exit()
    return table_used


###########################################################################
# No changes after here

db = MySQLdb.connect(read_default_file=default_mysql)
c = db.cursor()
cursor = c

class Peak():
    def __init__(self, t=None): 
        if not t is None: 
            self.transition = t
    @property
    def pQ3(self):
        return round(self.transition[0], 2)


def get_ssrcalc_values(seqs, input_sequences, default_ssrcalc):
    import os, csv
    ssr_query = """
    select sequence, ssrcalc
    from %(ssrcalc_table)s
    where sequence in (%(seqs)s)
    """ % { 'seqs' : seqs, 'ssrcalc_table' : default_ssrcalc}
    cursor.execute( ssr_query )
    pepmap = dict( cursor.fetchall() )
    #pepmap = {}

    # TODO: the used version in the TPP is 3.0 which is old and cannot be used
    # any more!? Is there a new pl script?

    not_found = []
    for ii,s in enumerate(input_sequences):
        try: ssrcalc = pepmap[filter(str.isalpha,s)]
        except KeyError: not_found.append(s)

    # SSRCalc finds the parameter file with ENV
    shellfile = '/tmp/ssrfile%s.sh' % os.getpid()
    outfile = '/tmp/ssrout%s.out' % os.getpid()
    env = {'SSRCalc' : ssrcalc_path } 
    cmd = """/SSRCalc3.pl --alg 3.0 --seq %s --output tsv --B 1 --A 0 > """ % " / ".join(not_found)
    cmd = ssrcalc_path + cmd + outfile

    f = open(shellfile, 'w')
    f.write(cmd)
    f.close()

    os.spawnlpe(os.P_WAIT, "/bin/bash", "bash", shellfile, env)
    r = csv.reader( open(outfile), delimiter='\t')
    for line in r:
        pepmap[line[0]] = float(line[2])

    os.system("rm %s" % shellfile)
    os.system("rm %s" % outfile)
    return pepmap

def unique_values(seq): 
    # order preserving
    checked = []
    for e in seq:
        if e not in checked:
            checked.append(e)
    return checked

def print_collding_peptides(collisions_per_peptide, precdic, ii, peaks):
    coll = [(k,v) for k,v in collisions_per_peptide.iteritems()]

    print """
    <a title="Show Tables" href="javascript:toggleDisplay('col_peptides_%s')">
        <h3>Colliding peptides</h3>
    </a>""" % ii

    #helper functions
    def sortpep_by_ssr(x,y):
        if len(x[1]) != len(y[1]): return -cmp( len(x[1]), len(y[1]) )
        #by SSRCalc
        return cmp( precdic[x[0]][3], precdic[y[0]][3]) 
    def sortpep_by_q1(x,y):
        #by difference in q1
        return cmp( abs(precdic[x[0]][0] - pep['q1']), abs(precdic[y[0]][0] - pep['q1'])) 

    coll.sort( sortpep_by_ssr )
    print "<table class='col_table' id='col_peptides_%i'>" % ii
    print "<tr> <td>Q1</td> <td>RT</td> <td>Sequence</td>  <td>Transitions</td> </tr>"
    for c in coll:
        print "<tr><td>"
        print round(precdic[c[0]][0], 2)
        print "</td><td>"
        print precdic[c[0]][3]
        print "</td><td>"
        print precdic[c[0]][1]
        print "</td><td>"
        for t in c[1]:
            print peaks[t].annotation
        print ' <br />'
        print "</td></tr>"
    print "</table>"

def print_unuseable(unuseable, nonunique, ii):

    print """
    <a title="Show Tables" href="javascript:toggleDisplay('col_transitions_%s')"> 
        <h3>Unuseable transitions (Collisions)</h3>
    </a>""" % ii
    print "<table class='col_table' id='col_transitions_%s'>" % ii
    print "<tr> <td>Type / Q3</td> <td>Colliding (Q1, Q3) / SSRCalc / Type / Sequence / Isotope </td> </tr>"
    unuseable.sort( lambda x,y: cmp(len(nonunique[x.transition[1]]), len(nonunique[y.transition[1]])))
    for peak in unuseable: 
        print "<tr><td>"
        print peak.annotation
        print peak.pQ3
        print "</td><td>"
        for c in nonunique[ peak.transition[1] ]:
            print '(%s,%s)' % (round(c[1], 2), round(c[0], 2))
            print c[7]
            print c[4] + str(c[5])
            if c[8] != 0:
                print c[6], 'isotope', c[8], ' </br>'
            else: print c[6], ' </br>'

        print "</td></tr>"
    print "</table>"

def main(myinput, q1_w, q3_w, ssr_w, db, high, low, genome, isotope, uis, ions):

    q3_low = low
    q3_high = high
    cursor = db.cursor()

    # sanitize input: all input is already sanitized except myinput and genome
    seqs = "'"
    input_sequences = []
    for inp in myinput.split():
        #only alphanumeric and [ ]
        sanitized = "".join( [i for i in inp if (str.isalnum(i) or i in [ '[', ']']  )] )
        #to look ssrcalc up in the db, we need no modifications
        seqs += filter(str.isalpha, inp) + "','"
        input_sequences.append(sanitized)
    seqs = seqs[:-2]

    table_used = map_db_tables(genome)

    #Now all input should be sane

    R = Residues('mono')
    if genome in genomes_that_require_N15_data: 
        R.recalculate_monisotopic_data_for_N15()

    #prepare a csv
    f = open( myCSVFile, 'w')
    w = csv.writer(f, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    w.writerow( ['sequence', 'q1', 'q3', 'type', 'charge'] )
    fuis = open( myUIS_CSVFile, 'w')
    wuis = csv.writer(fuis, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    MAX_UIS = uis

    #create the parameter object
    par = collider.SRM_parameters()
    par.dontdo2p2f = False #do not look at 2+ parent / 2+ fragment ions
    par.q1_window = q1_w / 2.0
    par.q3_window = q3_w / 2.0
    par.ssrcalc_window = ssr_w / 2.0 
    par.ppm = False
    par.considerIsotopes = True
    par.isotopes_up_to = isotope
    par.q3_range = [low, high]
    par.peptide_table = db_used + default_org_prefix + table_used
    par.transition_table = db_used + '.srmTransitions_' + table_used
    par.__dict__.update( ions )
    par.eval()
    mycollider = collider.SRMcollider()

    if uis > 0:
        if uis > 5:
            print "I will only calculate UIS up to order 5 \
                    (otherwise your Excel might break)." #and we dont want that.
            exit()
        print "<a href ='%s'>Download csv file with UIS</a>" % myUIS_CSVFile_rel
        print "<br/>"

    # Print headers, link to csv file and table of content
    unique = unique_values(myinput.split() )
    print "<a href ='%s'>Download csv file</a>" % myCSVFile_rel
    print "<br/>"
    print "input: %s peptides" % (len( myinput.split() )) 
    print "<br/>"
    print "unique: %s peptides" % (len( unique )) 
    print "<div class='toc'><ul>"
    for u in input_sequences: print '<li><a href="#%s">%s</a></li>' % (u,u)
    print "</ul></div>"
    print shared.toggleDisplay # Javascript function to toggle a div
    

    ###########################################################################
    # Start analysis
    # 1. Find SSRCalc values for all peptides
    # 2. Iterate through all sequences and calculate the b / y series
    # 3. Find all (potentially) interferring precursors
    # 4. For each precursors, find the list of transitions that interfers
    # 5. For each transition, find the precursors that interfere 
    ###########################################################################

    pepmap = get_ssrcalc_values(seqs, input_sequences, default_ssrcalc)
 
    for ii,s in enumerate(input_sequences):
        try: ssrcalc = pepmap[filter(str.isalpha,s)]
        except KeyError: ssrcalc = 25

        peptide = DDB.Peptide()
        peptide.set_sequence(s)
        peptide.charge = 2
        peptide.create_fragmentation_pattern(R)
        b_series = peptide.b_series
        y_series = peptide.y_series
        peaks = []

        #
        # Step 2 : calculate b/y ion series of the target
        #
        pcount = 0
        for ch in [1]:
            scount = 0
            for pred in y_series:
                scount += 1
                q3 = ( pred + (ch -1)*R.mass_H)/ch
                if q3 < q3_low or q3 > q3_high: continue
                p = Peak( [q3, pcount] ) 
                p.annotation = 'y' + str(len(y_series) + 1 - scount) 
                p.charge = ch
                pcount += 1
                peaks.append( p )
            scount = 0
            for pred in b_series:
                scount += 1
                q3 = ( pred + (ch -1)*R.mass_H)/ch
                if q3 < q3_low or q3 > q3_high: continue
                p = Peak( [q3, pcount] ) 
                p.annotation = 'b' + str(scount)
                p.charge = ch
                pcount += 1
                peaks.append( p )

        q1 = peptide.charged_mass 
        transitions = c_getnonuis.calculate_transitions_ch(
            ((q1, s, 1),), [1], q3_low, q3_high)
        transitions = [ (p[0], i) for i,p in enumerate(transitions)]
        nr_transitions = len( transitions )
        if nr_transitions == 0: continue #no transitions in this window

        #
        # Step 3 : find all potentially interfering precursors
        #
        pep = {
            'mod_sequence' : s,
            'parent_id' :  -1,
            'q1' :         q1,
            'q1_charge' :  2,
            'ssrcalc' :    ssrcalc,
            'peptide_key' :-1,
        }
        precursors = mycollider._get_all_precursors(par, pep, cursor, 
            values="q1, modified_sequence, peptide_key, ssrcalc, isotope_nr, q1_charge")
        precursors = [p for p in precursors if p[1] != pep['mod_sequence'] ]
        precdic = dict( [ (p[2], p) for p in precursors] )

        # 
        # Step 4 and 5
        #
        # 
        # Find interferences per precursor, then find interferences per
        # transition (see the two readouts in the html)
        if par.do_b_y_only():
            collisions_per_peptide = \
            c_getnonuis.calculate_collisions_per_peptide( tuple(transitions),
                tuple(precursors), q3_low, q3_high, par.q3_window, par.ppm)
            nonunique = c_getnonuis._find_clashes_forall( tuple(transitions),
                tuple(precursors), q3_low, q3_high, par.q3_window, par.ppm)
        else:
            collisions_per_peptide = \
            c_getnonuis.calculate_collisions_per_peptide_other_ion_series(
                tuple(transitions), tuple(precursors), q3_low, q3_high,
                par.q3_window, par.ppm, par) 
            nonunique = c_getnonuis._find_clashes_forall_other_series( 
                tuple(transitions), tuple(precursors), q3_low, q3_high, 
                par.q3_window, par.ppm, par)

        non_uis = []
        if uis > 0:
            wuis.writerow([s])

            srm_ids = [peak.transition[1] for peak in peaks]
            srm_lookup = [ (peak.transition[1] , peak) for peak in peaks]
            srm_lookup = dict(srm_lookup) 
            for order in range(1,uis+1): 
                non_uis = c_getnonuis.get_non_uis(collisions_per_peptide, order)
                really_calculate_uis = True
                if not really_calculate_uis:
                    # here we just output the non-UIS combinations. Usually
                    # these are more informative and are preferable to a list
                    # of UIS combinations.
                    for comb in non_uis:
                        tmp = [ srm_lookup[elem] for elem in comb]  
                        myrow = []
                        for tt in tmp:
                            myrow.extend( [ tt.transition[0], tt.annotation ])
                        wuis.writerow(myrow)
                if really_calculate_uis:
                    # if you want the real deal, go ahead. 
                    uis_list = collider.get_uis(srm_ids, non_uis, order)
                    for comb in uis_list:
                        tmp = [ srm_lookup[elem] for elem in comb]  
                        myrow = []
                        for tt in tmp:
                            myrow.extend( [ tt.transition[0], tt.annotation ])
                        wuis.writerow(myrow)

        useable = []
        unuseable = []
        for peak in peaks: 
            if peak.transition[1] in nonunique:
                unuseable.append(peak)
            else: useable.append(peak)

        ########################################################################
        # Printing
        ########################################################################
        print "<h2 id='%s'>Peptide %s</h2>" % (s, s)
        print "<div class='pep_property'>"
        print "<p>Q1: %s</p>" % pep['q1']
        print "<p>SSRCalc: %s</p>" % round(pep['ssrcalc'], 2)
        print "<p>Percent Useable: %d %%</p>" % (len(useable)*100.0 / len(transitions)  )
        print "</div>"

        if len( useable ) > 0:
            print "<h3>Useable transitions</h3>"
            print "<table class='tr_table'>"
            print "<tr><td>Type</td> <td>Q3</td> </tr>"
            for peak in useable:
                print "<tr><td>"
                print peak.annotation
                print "</td><td>"
                print peak.pQ3
                w.writerow( [ s, pep['q1'], round(peak.pQ3), peak.annotation, 
                             peak.charge ] )
                print "</td><td>"
            print "</table>"
        else: 
            print "<p>No useable transitions for this peptide!</p>"
            w.writerow( ['"Sorry: no useable transitions for this peptide. Try to calculate UIS."' ] )

        print_collding_peptides(collisions_per_peptide, precdic, ii, peaks)
        print_unuseable(unuseable, nonunique, ii)

    fuis.close()
    f.close()

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
print shared.header
print shared.warm_welcome
print "<div class='main'>"

form = cgi.FieldStorage()   # FieldStorage object to
if form.has_key('peptides'):
    peptides = form.getvalue('peptides')
    q1_w = float(form.getvalue('q1_window') )
    q3_w = float(form.getvalue('q3_window') )
    ssr_w = float(form.getvalue('ssr_window') )
    high = float(form.getvalue('high_mass') )
    low = float(form.getvalue('low_mass') )
    genome = form.getvalue('genome') 
    isotope = int(form.getvalue('isotope') )
    uis = int(form.getvalue('uis') )
    ions = {'aions'    : bool(form.getvalue('aions'      )),
            'aMinusNH3': bool(form.getvalue('aMinusNH3'  )),
            'bions'    : bool(form.getvalue('bions'      )),
            'bMinusH2O': bool(form.getvalue('bMinusH2O'  )),
            'bMinusNH3': bool(form.getvalue('bMinusNH3'  )),
            'bPlusH2O' : bool(form.getvalue('bPlusH2O'   )),
            'cions'    : bool(form.getvalue('cions'      )),
            'xions'    : bool(form.getvalue('xions'      )),
            'yions'    : bool(form.getvalue('yions'      )),
            'yMinusH2O': bool(form.getvalue('yMinusH2O'  )),
            'yMinusNH3': bool(form.getvalue('yMinusNH3'  )),
            'zions'    : bool(form.getvalue('zions'      ))}
    start = time.time()
    main( peptides, q1_w, q3_w, ssr_w, db, high, low, genome, isotope, uis, ions)
    print "<hr> <br/>This query took: %s s" % (time.time() - start)
else:
  ions = ['aions'      ,
         'aMinusNH3'  ,
         'bions'      ,
         'bMinusH2O'  ,
         'bMinusNH3'  ,
         'bPlusH2O'   ,
         'cions'      ,
         'xions'      ,
         'yions'      ,
         'yMinusH2O'  ,
         'yMinusNH3'  ,
         'zions'      ]
  html_ions = ''
  for ion in ions:
      #<label class="mylabel" for="%s">%s</label>
      check = ''
      if ion in ['bions', 'yions']: check = 'checked="yes" '
      html_ions += """
      <input %(check)s type="checkbox" name="%(ion)s" value="s(ion)s"> %(ion)s<br>
      """ %  {'ion' : ion, 'check' : check}
  print shared.toggleDisplay # Javascript function to toggle a div
  print """
<form action="/srmcollider/srmcollider.py" method="post">
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
            %(genome_select)s
        </select>
    </p>

    <p class='input_field'>
        <label class="mylabel" for="isotope">Consider isotopes up to </label>
        <input class="number_input" type="text" name="isotope" value="3"> amu
        <!-- >
        <input class="number_input" type="text" disabled="False" name="isotope" value="3"> amu
        </!-->
    </p>


    <p class='input_field'>
        <label class="mylabel" for="uis">Find UIS up to order* </label>
        <input class="number_input" type="text" name="uis" value="0"> 
    </p>

        
    <p>
    <a title="Show Tables" href="javascript:toggleDisplay('bg_ion')">
        Background Ion Series
    </a>

    <p id="bg_ion" style="display:none;">
      %(ion_series)s

    </p>

    </p>

    <INPUT type="submit" value="Send"> 

 </form>
<p>
* This will print out a list of all UIS combinations of transitions for all
peptides. Preferably, only combinations of less than 5 transitions are
considered, otherwise the result might get quite large.
</p>

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
<br/>%(sample_peptides_html)s    
</!-->
""" % {'sample_peptides_html' : sample_peptides_html, 'genome_select':
       genome_select, 'ion_series' : html_ions} 

print "</div>"
print """
</div>
</body>
</html>"""
