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

# All changes that should be done by the user are in the collider_config.py
from collider_config import *

###########################################################################
# No changes after here

import MySQLdb, time, sys 
import csv, re
import random, string
import cgi
import os

sys.path.append(SRMCOLLIDER_HOME)
import DDB
from Residues import Residues
import collider
import c_getnonuis
import sharedhtml as shared
from precursor import Precursor

db = MySQLdb.connect(read_default_file=default_mysql)
c = db.cursor()
cursor = c
R = Residues('mono')

def get_ssrcalc_values(seqs, input_sequences, default_ssrcalc):
    if default_ssrcalc != '':
        ssr_query = """
        select sequence, ssrcalc
        from %(ssrcalc_table)s
        where sequence in (%(seqs)s)
        """ % { 'seqs' : seqs, 'ssrcalc_table' : default_ssrcalc}
        cursor.execute( ssr_query )
        pepmap = dict( cursor.fetchall() )
    else: pepmap = {}

    not_found = []
    for ii,s in enumerate(input_sequences):
        try: ssrcalc = pepmap[filter(str.isalpha,s)]
        except KeyError: not_found.append(s)

    # TODO: the used version in the TPP is 3.0 which is old and cannot be used
    # any more online. It makes it hard to compare. Is there a new pl script?

    # SSRCalc finds the parameter file with ENV
    shellfile = '/tmp/ssrfile%s.sh' % os.getpid()
    outfile = '/tmp/ssrout%s.out' % os.getpid()
    env = {'SSRCalc' : ssrcalc_path } 
    cmd = """/SSRCalc3.pl --alg 3.0 --seq "%s" --output tsv --B 1 --A 0 > """ % " / ".join(not_found)
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

def print_collding_peptides(collisions_per_peptide, precdic, ii, fragments):
    coll = [(k,v) for k,v in collisions_per_peptide.iteritems()]

    print """
    <a title="Show Tables" href="javascript:toggleDisplay('col_peptides_%s')">
        <h3>Interfering peptides <small><small>(Click to fold in.)</small></small> </h3>

    </a>""" % ii

    #helper functions
    def sortpep_by_ssr(x,y):
        if len(x[1]) != len(y[1]): return -cmp( len(x[1]), len(y[1]) )
        #by SSRCalc
        return cmp( precdic[x[0]].ssrcalc, precdic[y[0]].ssrcalc) 
    def sortpep_by_q1(x,y):
        #by difference in q1
        return cmp( abs(precdic[x[0]].q1 - pep['q1']), abs(precdic[y[0]].q1 - pep['q1'])) 

    coll.sort( sortpep_by_ssr )
    print "<table class='col_table' id='col_peptides_%i'>" % ii
    print "<tr> <td>Q1</td> <td>RT</td> <td>Sequence</td>  <td>Transitions</td> </tr>"
    for c in coll:
        # Q1, RT, Sequence, Transitions
        precursor = precdic[c[0]]
        print "<tr><td>"
        print round(precursor.q1, 2)
        print "</td><td>"
        print precursor.ssrcalc
        print "</td><td>"
        print precursor.modified_sequence
        print "</td><td>"
        for t in c[1]:
            print fragments[t].annotation
        print ' <br />'
        print "</td></tr>"
    print "</table>"

def print_unuseable(unuseable, nonunique, ii):
    print """
    <a title="Show Tables" href="javascript:toggleDisplay('col_transitions_%s')"> 
        <h3>All transitions with Collisions <small><small>(Click to fold in.)</small></small> </h3>

    </a>""" % ii
    print "<table class='col_table' id='col_transitions_%s'>" % ii
    print "<tr> <td>Type / Q3</td> <td>Colliding (Q1, Q3) / SSRCalc / Type / Sequence / Isotope </td> </tr>"
    unuseable.sort( lambda x,y: cmp(len(nonunique[x.fragment_count]), len(nonunique[y.fragment_count])))
    for peak in unuseable: 
        print "<tr><td>"
        print peak.annotation
        print peak.pQ3
        print "</td><td>"
        this_interference = nonunique[ peak.fragment_count ]
        this_interference.sort( lambda x,y: cmp(x[1], y[1]))
        for c in this_interference:
            print '(%s,%s)' % (round(c[1], 2), round(c[0], 2)) # print Q1/Q3
            print c[7] # print ssrcalc
            print c[4] +'_' + str(c[5]) # print type + fragment nr
            if c[8] != 0:
                print c[6], 'isotope', c[8], ' </br>'
            else: print c[6], ' </br>'

        print "</td></tr>"
    print "</table>"

def main(myinput, q1_w, q3_w, ssr_w, db, high, low, genome, isotope, uis, ions, 
         missed, oxMet, Deamid, chargeCheck):

    q3_low = low
    q3_high = high
    cursor = db.cursor()

    # create unique files
    thisrandom = "".join( [string.ascii_letters[ int(random.random() * len(string.ascii_letters) )] for i in range(10)])
    myCSVFile         = myCSVFile_         + thisrandom + '.csv'
    myCSVFile_rel     = myCSVFile_rel_     + thisrandom + '.csv'
    myUIS_CSVFile     = myUIS_CSVFile_     + thisrandom + '.csv'
    myUIS_CSVFile_rel = myUIS_CSVFile_rel_ + thisrandom + '.csv'

    # sanitize input: all input is already sanitized except myinput and genome
    seqs = "'"
    input_sequences = []
    for inp in myinput.split():
        #only alphanumeric and [ ]
        sanitized = "".join( [i for i in inp if (str.isalnum(i) or i in [ '[', ']']  )] )
        #to look ssrcalc up in the db, we need no modifications
        if len(sanitized) == 0: continue
        seqs += filter(str.isalpha, inp) + "','"
        input_sequences.append(sanitized.upper())
    seqs = seqs[:-2]

    table_used = map_db_tables(genome)

    #Now all input should be sane
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
    par.query2_add = ''
    if not oxMet and not Deamid:
        par.query2_add = ' and modifications = 0 '
    elif oxMet and not Deamid:
        par.query2_add = " and modified_sequence not like '%N[115]%' "
    elif not oxMet and Deamid:
        par.query2_add = " and modified_sequence not like '%M[147]%' "

    if missed == 0:
        par.query2_add += ' and missed_cleavages = 0 '
    if missed == 1:
        par.query2_add += ' and missed_cleavages <= 1 '

    print shared.resultInterpretation
    if uis > 0:
        if uis > 5:
            print "I will only calculate UIS up to order 5 \
                    (otherwise your Excel might break)." #and we dont want that.
            exit()
        print "<a href ='%s'>Download csv file with UIS</a>" % myUIS_CSVFile_rel
        print "<br/>"

    # Print headers, link to csv file and table of content
    unique = unique_values(myinput.split() )
    # print "<a href ='%s'>Download csv file</a>" % myCSVFile_rel
    # print "<br/>"
    print "input: %s peptides" % (len( myinput.split() )) 
    print "<br/>"
    print "unique: %s peptides" % (len( unique )) 
    print "<div class='toc'><ul>"
    for u in input_sequences: print '<li><a href="#%s">%s</a></li>' % (u,u)
    print "</ul></div>"
    print """<a title="Toggle all" href="javascript:toggleAll()">
    Toggle all <small>Click to fold/unfold all</small>
    </a>"""

    print shared.toggleDisplay # Javascript function to toggle a div

    do_analysis(input_sequences, seqs, q3_low, q3_high, par, wuis, w)
    fuis.close()
    f.close()

def do_analysis(input_sequences, seqs, q3_low, q3_high, par, wuis, w):
 
    ###########################################################################
    # Do analysis
    # 1. Find SSRCalc values for all peptides
    # 2. Iterate through all sequences and calculate the b / y series
    # 3. Find all (potentially) interfering precursors
    # 4. For each precursors, find the list of transitions that interfers
    # 5. For each transition, find the precursors that interfere 
    ###########################################################################

    pepmap = get_ssrcalc_values(seqs, input_sequences, default_ssrcalc)
    toggle_all_str = '<script language="javascript"> function toggleAll(){ '
    mycollider = collider.SRMcollider()

    for ii,s in enumerate(input_sequences):
        try: ssrcalc = pepmap[filter(str.isalpha,s)]
        except KeyError: ssrcalc = 25

        peptide = DDB.Peptide()
        peptide.set_sequence(s)
        peptide.charge = 2
        peptide.create_fragmentation_pattern(R)

        #
        # Step 2 : calculate b/y ion series of the target and get the transitions
        #
        fragments = list(peptide.get_fragment_objects(reversed(peptide.y_series),
            'y', 1, R, q3_low, q3_high))
        fragments.reverse()
        fragments.extend(list( peptide.get_fragment_objects(peptide.b_series, 
            'b', 1, R, q3_low, q3_high)))
        for fcount, f in enumerate(fragments): 
            f.fragment_count = fcount

        transitions = [ (f.q3, f.fragment_count) for f in fragments]
        if len( transitions ) == 0: continue # no transitions in this window

        #
        # Step 3 : find all potentially interfering precursors
        #  Create precursor and use db to find all interfering precursors
        #
        precursor = Precursor(
            modified_sequence = s, parent_id = -1,
            q1 = peptide.charged_mass, q1_charge = 2, 
            ssrcalc = ssrcalc, transition_group = -1)
        precursors_obj = mycollider. _get_all_precursors(
            par, precursor, cursor, bysequence=True)
        precursor_obj_dic = dict( [ (p.transition_group, p) for p in precursors_obj] )

        # 
        # Step 4 and 5
        #
        # 
        # Find interferences per precursor, then find interferences per
        # transition (see the two readouts in the html)

        collisions_per_peptide = \
        c_getnonuis.calculate_collisions_per_peptide_other_ion_series(
            tuple(transitions), precursors_obj, par, q3_low, q3_high,
            par.q3_window, par.ppm, chargeCheck) 

        nonunique = c_getnonuis._find_clashes_forall_other_series( 
            tuple(transitions), precursors_obj, par, q3_low, q3_high,
            par.q3_window, par.ppm, peptide.charged_mass - par.q1_window, chargeCheck)

        non_uis = []
        if uis > 0:
            wuis.writerow([s])

            srm_ids = [f.fragment_count for f in fragments]
            srm_lookup = [ (fragment.fragment_count, fragment) for fragment in fragments]
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
                            myrow.extend( [ tt.q3, tt.annotation ])
                        wuis.writerow(myrow)
                if really_calculate_uis:
                    # if you want the real deal, go ahead. 
                    uis_list = collider.get_uis(srm_ids, non_uis, order)
                    if(len(uis_list) == 0): wuis.writerow([ 'Sorry, no UIS found for order %s' % order ])
                    for comb in uis_list:
                        tmp = [ srm_lookup[elem] for elem in comb]  
                        myrow = []
                        for tt in tmp:
                            myrow.extend( [ tt.q3, tt.annotation ])
                        wuis.writerow(myrow)

        useable = []
        unuseable = []
        for fragment in fragments: 
            if fragment.fragment_count in nonunique:
                unuseable.append(fragment)
            else: useable.append(fragment)

        ########################################################################
        # Printing
        ########################################################################
        print "<h2 id='%s'>Peptide %s</h2>" % (s, s)
        print "<div class='pep_property'>"
        print "<p>Q1: %s</p>" % precursor.q1
        print "<p>SSRCalc: %s</p>" % round(precursor.ssrcalc, 2)
        # print "<p>Percent Useable: %d %%</p>" % (len(useable)*100.0 / len(transitions)  )
        print "</div>"

        if len(useable ) > 0:
            print "<h3>Transitions without any interferences</h3>"
            print "<table class='tr_table'>"
            print "<tr><td>Type</td> <td>Q3</td> </tr>"
            for peak in useable:
                print "<tr><td>"
                print peak.annotation
                print "</td><td>"
                print peak.pQ3
                w.writerow( [s, precursor.q1, round(peak.pQ3), peak.annotation, 
                             peak.charge ] )
                print "</td><td>"
            print "</table>"
        else: 
            print "<p>No transitions that have no interference at all!</p>"
            w.writerow( ['"Sorry: no useable transitions for this peptide. Try to calculate UIS."' ] )

        print_collding_peptides(collisions_per_peptide, precursor_obj_dic,
                                ii, fragments)
        print_unuseable(unuseable, nonunique, ii)
        toggle_all_str += "toggleDisplay('col_peptides_%s'); toggleDisplay('col_transitions_%s');\n" % (ii,ii)

    toggle_all_str += "};</script>"
    print toggle_all_str


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
            'zions'    : bool(form.getvalue('zions'      )),
            'MMinusH2O': bool(form.getvalue('MMinusH2O'  )),
            'MMinusNH3': bool(form.getvalue('MMinusNH3'  )),
           }
    missed = int(form.getvalue('missed') )
    oxMet =  bool(form.getvalue('oxMet') )
    Deamid = bool(form.getvalue('Deamid') )
    chargeCheck = bool(form.getvalue('chargeCheck') )
    start = time.time()
    main( peptides, q1_w, q3_w, ssr_w, db, high, low, genome, isotope, uis, ions,
        missed, oxMet, Deamid, chargeCheck)
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
         'zions'      ,
         'MMinusH2O'  ,
         'MMinusNH3'  ,
         ]
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
        <input class="number_input" type="text" name="ssr_window" value="10"> arbitrary units
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
        <label class="mylabel" for="missed">Missed Cleavages</label>
        <select name="missed" class="number_input">
            <option value="0">0</option>
            <option value="1">1</option>
        </select>
    </p>


    <p class='input_field'>
        <label class="mylabel" for="uis">Find UIS up to order* </label>
        <input class="number_input" type="text" name="uis" value="2"> 
    </p>

    <p>
    
    <label> Charge check: </label>
      <input type="checkbox" name="chargeCheck" value="chargeCheck"> Check that 
      interfering signal can actually hold charge (e.g. 2+ charge) </br>
    </p>

    <label> Modifications: </label>
      <input type="checkbox" name="oxMet" value="oxMet"> oxidized Methionines
      <input type="checkbox" name="Deamid" value="Deamid"> deamidated Asparagines<br>
    </p>
        
    <p>
    <a title="Show Tables" href="javascript:toggleDisplay('bg_ion')">
        Background Ion Series
    </a>

    <p id="bg_ion" style="display:block;">
      %(ion_series)s

    </p>
    Please note that this server may take some time to respond to the
    <b>first</b> query, this is due to the time needed to populate the MySQL
    cache with the corresponding indices. After running a couple of peptides
    (10-20), the computation time per peptide should be well below 100 ms
    for yeast (other proteomes are bigger and might take longer, also
    non-default options like higher order UIS might increase the computing
    time).
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
