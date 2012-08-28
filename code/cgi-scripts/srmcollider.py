#!/usr/bin/python
# -*- coding: utf-8  -*-

"""
 *
 * Program       : SRMCollider
 * Author        : Hannes Roest <roest@imsb.biol.ethz.ch>
 * Date          : 05.02.2011 
 *
 *
 * Copyright (C) 2011 - 2012 Hannes Roest
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
import collider
import c_getnonuis
import sharedhtml as shared
from precursor import Precursor
from srmcollider_website_helper import getSRMParameter, get_ssrcalc_values
from srmcollider_website_helper import unique_values, write_csv_row
from srmcollider_website_helper import SRMColliderController

backw_compatible = False
rounding_precision = 4

db = MySQLdb.connect(read_default_file=default_mysql)

controller = SRMColliderController()
controller.initialize(db_used=db_used, default_org_prefix=default_org_prefix,
    db_tables_map=db_tables_map)

def get_html_ions():
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
  return html_ions

def print_header(input_sequences):
    """
    Print headers, link to csv file and table of content
    """
    unique = unique_values(input_sequences)
    # print "<a href ='%s'>Download csv file</a>" % myCSVFile_rel
    # print "<br/>"
    print "input: %s peptides" % (len(input_sequences)) 
    print "<br/>"
    print "unique: %s peptides" % (len( unique )) 
    print "<div class='toc'><ul>"
    for u in input_sequences: print '<li><a href="#%s">%s</a></li>' % (u,u)
    print "</ul></div>"
    print """<a title="Toggle all" href="javascript:toggleAll()">
    Toggle all <small>Click to fold/unfold all</small>
    </a>"""

def prepare_extern_graph(data, q1, ssrcalc, name):
    thisrandom = "".join( [string.ascii_letters[ int(random.random() * len(string.ascii_letters) )] for i in range(10)])
    filename_part = 'srmcollider_data_%s_peptide' % (thisrandom)
    filename = '/tmp/' + filename_part
    f = open(filename, 'w')
    for d in data:
        f.write('%s\t%s\n' % (d[0], d[1]) )
    f.close()
    args = "data=%s&mz=%s&rt=%s&pepname=%s" % (
        filename_part, q1, ssrcalc, name)
    return args

def print_peptide_header(current_sequence, precursor, par, precursors_obj):
    nr_interferences = len(precursors_obj)
    print "<h2 id='%s'>Peptide %s</h2>" % (current_sequence, current_sequence)
    print "<div class='pep_property'>"

    print "<table border='1' class='col_table' id='header_peptide'>" 
    print """<tr> 
      <td>Sequence</td> 
      <td>Q1</td> 
      <td>Q1 window</td> 
      <td><a href='http://hs2.proteome.ca/SSRCalc/SSRCalcX.html'>SSRCalc</a></td>  
      <td>SSRCalc window</td> 
      <td>Interfering precursors</td> 
      <td>Background</td> 
      <td>Graph</td> 
      </tr>"""
    #<td>Interferences</td> 
    print "<tr><td>"
    print precursor.modified_sequence
    print "</td><td>"
    print precursor.q1
    print "</td><td> &#x00B1;"
    print par.q1_window
    print "</td><td>"
    print round(precursor.ssrcalc, 2)
    print "</td><td> &#x00B1;"
    print par.ssrcalc_window
    print "</td><td>"
    print nr_interferences
    print "</td><td>"
    print par.genome
    print "</td><td>"
    args = prepare_extern_graph([ (c.q1, c.ssrcalc) for c in precursors_obj],
      precursor.q1, precursor.ssrcalc, precursor.modified_sequence)
    print "<a rel='lightbox' href='plot_graph_dynamic.py?%s' >Graph</a>" % args
    print "</td><tr>" 
    print "</table>"


    print "</div>"

def print_transition_overview(fragments, precursor, nonunique_obj):
    ii = precursor.seq_id
    print "<h3>Transition Overview</h3>" 
    ## #print "<a title="Show Tables" href="javascript:toggleDisplay('col_transitions2_%s')"> " % ii
    ## #print "<h3>Transition Overview<small><small> (Click to fold)</small></small> </h3>" 
    ## #print '</a>' 
    print "<table border='1' class='col_table' id='col_transitions2_%s'>" % ii
    print "<tr> <td>Transition </td> <td>Q3</td>  <td>Interferences</td> <td>Graph</td> </tr>"
    ## #print "<tr> <td>Transition </td> <td>Q3</td>  <td>Interferences</td> </tr>"
    fragments.sort( lambda x,y: cmp(len(nonunique_obj[x.fragment_count]), len(nonunique_obj[y.fragment_count])))
    for peak in fragments: 
        interferences = nonunique_obj[peak.fragment_count]
        print "<tr><td>"
        print peak.annotation
        if False: # debug
            print "</td><td>"
            print peak.fragment_count
        print "</td><td>"
        print peak.pQ3
        print "</td><td>%s" % len(nonunique_obj[ peak.fragment_count ])
        # print "</td><td>"
        # print [ (i.q3, i.ssrcalc) for i in interferences]
        print "</td><td>"
        this_interference = nonunique_obj[ peak.fragment_count ]
        args = prepare_extern_graph( [ (c.q3, c.ssrcalc) for c in this_interference], 
          peak.q3, precursor.ssrcalc, "%s_%s" % (precursor.modified_sequence, peak.annotation) )
        print "<a rel='lightbox' href='plot_graph_dynamic.py?%s' >Graph</a>" % args
        print "</td></tr>"
    print "</table>"

def print_collding_peptides(collisions_per_peptide, precursors_obj, ii, fragments):
    precdic = dict( [ (p.transition_group, p) for p in precursors_obj] )
    coll = [(k,v) for k,v in collisions_per_peptide.iteritems()]

    print """
    <a title="Show Tables" href="javascript:toggleDisplay('col_peptides_%s')">
        <h3>Interfering peptides <small><small>(Click to fold)</small></small> </h3>

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

def print_uniqueness_analysis(collisions_per_peptide, peptide):
    if len([0 for f in peptide.fragments if f.library_intensity > 0]) == 0: return

    nonunique_set = set([ tuple(v) for k,v in collisions_per_peptide.iteritems()])
    transitions_order = [ (f.fragment_count, f.library_intensity, f) for f in peptide.fragments]
    transitions_order.sort( lambda x,y: -cmp(x[1], y[1]))
    ii = 0
    print "<h3>Uniqueness of Top n transitions</h3>" 
    print "<table border='1' class='col_table' id='x_col_transitions2_%s'>" % ii
    print "<tr> <td>Transitions</td><td>Nr transitions</td> <td>Combined Uniqueness</td> </tr>"
    for i in range(len(transitions_order)):
      #print [f for f in peptide.fragments if f.fragment_count ==  transitions_order[i][0]]
      mytransitions = tuple(sorted( [t[0] for t in transitions_order[:i+1]]))
      tr = [t[2] for t in transitions_order[:i+1]]


      # mytransitions = tuple(sorted([t[1] for t in transitions[:j]]))
      # unuseable = False
      # for k,v in collisions_per_peptide.iteritems():
      #     if tuple(sorted(v)[:j]) == mytransitions: unuseable=True
      # if not unuseable: min_needed = j

      unuseable = False
      for k,v in collisions_per_peptide.iteritems():
          if tuple(sorted(v)[:i+1]) == mytransitions: unuseable=True

      print "<tr><td>"
      print " + ".join([t.annotation for t in tr])
      print "</td><td>%s</td><td>" % (i+1)
      if unuseable:
        print "No"
      else:
        print "Yes"
      print "</td></tr>"
    print "</table>"

def print_transition_detail(unuseable, nonunique_obj, ii):
    print """
    <a title="Show Tables" href="javascript:toggleDisplay('col_transitions_%s')"> 
        <h3>All transitions with Collisions <small><small>(Click to fold)</small></small> </h3>

    </a>""" % ii
    print "<table class='col_table' id='col_transitions_%s'>" % ii
    print "<tr> <td>Type / Q3</td> <td>Colliding (Q1, Q3) / SSRCalc / Type / Sequence / Isotope </td> </tr>"

    unuseable.sort( lambda x,y: cmp(len(nonunique_obj[x.fragment_count]), len(nonunique_obj[y.fragment_count])))
    for peak in unuseable: 
        if backw_compatible and len(nonunique_obj[ peak.fragment_count ]) == 0: continue
        print "<tr><td>"
        print peak.annotation
        print round(peak.q3, rounding_precision)
        print "</td><td>"

        this_interference = nonunique_obj[ peak.fragment_count ]
        if backw_compatible: this_interference.sort( lambda x,y: cmp(x.q1, y.q1))
        else: this_interference.sort( lambda x,y: cmp(x.q3, y.q3))
        for c in this_interference:
            print '(%s,%s)' % (round(c.q1, rounding_precision), round(c.q3, rounding_precision)) # print Q1/Q3
            print c.ssrcalc # print ssrcalc
            print c.ion_type +'_' + str(c.ion_number) # print type + fragment nr
            if c.isotope_nr != 0:
                print c.sequence, 'isotope', c.isotope_nr, ' </br>'
            else: print c.sequence, ' </br>'

        print "</td></tr>"
    print "</table>"

def main(par):
    local_cursor = db.cursor()

    # create unique files and prepare a csv
    thisrandom = "".join( [string.ascii_letters[ int(random.random() * len(string.ascii_letters) )] for i in range(10)])
    myUIS_CSVFile     = myUIS_CSVFile_     + thisrandom + '.csv'
    myUIS_CSVFile_rel = myUIS_CSVFile_rel_ + thisrandom + '.csv'
    fuis = open( myUIS_CSVFile, 'w')
    writer_uis = csv.writer(fuis, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    writer_uis.writerow(['Sequence', 'UIS Order', 'Q3', 'Annotation', '(For higher order UIS, the pattern Q3/Annotation repeats)'])

    print shared.resultInterpretation
    if par.uis > 0:
        if par.uis > 5:
            print "I will only calculate UIS up to order 5 \
                    (otherwise your Excel might break)." #and we dont want that.
            exit()
        print "<a href ='%s'>Download csv file with UIS</a>" % myUIS_CSVFile_rel
        print "<br/>"

    print_header(par.input_sequences)
    print shared.toggleDisplay # Javascript function to toggle a div

    do_analysis(par.input_sequences, par.seqs, par, writer_uis, local_cursor)
    fuis.close()

def do_analysis(input_sequences, seqs, par, wuis, local_cursor):
    """
    ###########################################################################
    # Do analysis
    # 1. Find SSRCalc values for all peptides
    # 2. Iterate through all sequences and calculate the b / y series
    # 3. Find all (potentially) interfering precursors
    # 4. For each precursors, find the list of transitions that interfers
    # 5. For each transition, find the precursors that interfere 
    # 6. Print information
    ###########################################################################
    """

    q3_low, q3_high = par.q3_range
    uis = par.uis
    pepmap = get_ssrcalc_values(seqs, input_sequences, default_ssrcalc, local_cursor, ssrcalc_path)
    toggle_all_str = '<script language="javascript"> function toggleAll(){ '
    mycollider = collider.SRMcollider()

    for seq_id, peptide in enumerate(controller.peptides):
        #
        # Step 2 : find the SSRCalc values for this sequence
        #
        try: ssrcalc = pepmap[filter(str.isalpha, peptide.sequence)]
        except KeyError: ssrcalc = 25
        transitions = [ (f.q3, f.fragment_count) for f in peptide.fragments]
        if len( transitions ) == 0: continue # no transitions in this window

        #
        # Step 3 : find all potentially interfering precursors
        #  Create precursor and use db to find all interfering precursors
        #
        precursor = Precursor(modified_sequence = peptide.get_modified_sequence(), parent_id = -1,
            q1 = peptide.charged_mass, q1_charge = 2, ssrcalc = ssrcalc, transition_group = -1)
        precursor.seq_id = seq_id
        precursors_obj = mycollider._get_all_precursors(par, precursor, local_cursor)

        # 
        # Step 4 and 5: Find interferences per precursor, then find
        # interferences per transition (see the two readouts in the html)
        #
        collisions_per_peptide = \
        c_getnonuis.calculate_collisions_per_peptide_other_ion_series(
            tuple(transitions), precursors_obj, par, q3_low, q3_high,
            par.q3_window, par.ppm, par.chargeCheck) 

        nonunique = c_getnonuis._find_clashes_forall_other_series( 
            tuple(transitions), precursors_obj, par, q3_low, q3_high,
            par.q3_window, par.ppm, peptide.charged_mass - par.q1_window, par.chargeCheck)

        # also add those that have no interference
        for fragment in peptide.fragments: 
            if not fragment.fragment_count in nonunique:
                nonunique[fragment.fragment_count] = []
        nonunique_obj = controller.getNonuniqueObjects(nonunique)

        # 
        # Step 6: printing
        #
        do_all_print(peptide, collisions_per_peptide, 
                 wuis, precursor, par, precursors_obj, nonunique_obj)
        toggle_all_str += "toggleDisplay('col_peptides_%s'); toggleDisplay('col_transitions_%s');\n" % (seq_id,seq_id)

    toggle_all_str += "};</script>"
    print toggle_all_str
    print """
    <script language="javascript">
        window.onload = toggleAll();
    </script>"""

def do_all_print(peptide, collisions_per_peptide, 
    wuis, precursor, par, precursors_obj, nonunique_obj):
    fragments = peptide.fragments
    uis = par.uis
    current_sequence = precursor.modified_sequence
    if uis > 0: write_csv_row(fragments, collisions_per_peptide, current_sequence, uis, wuis)
    print_peptide_header(current_sequence, precursor, par, precursors_obj ) 
    print_transition_overview(fragments[:], precursor, nonunique_obj)
    print_uniqueness_analysis(collisions_per_peptide, peptide)
    print_collding_peptides(collisions_per_peptide, precursors_obj, precursor.seq_id, fragments)
    print_transition_detail(fragments[:], nonunique_obj, precursor.seq_id)

###########################################################################
###########################################################################
# START OF HTML

print 'Content-type: text/html\n\n'
print shared.header
print shared.warm_welcome
print "<div class='main'>"



sample_peptides_html = controller.get_sample_peptides_html()
form = cgi.FieldStorage()   # FieldStorage object to
if form.has_key('peptides') and not form.has_key("from_api"):
    par = controller.parse_srmcollider_form(form, genomes_that_require_N15_data)
    start = time.time()
    main(par)
    print "<hr> <br/>This query took: %s s" % (time.time() - start)
else:

  html_ions = get_html_ions()
  textfield_peptides = ""
  if form.has_key('peptides'):
    textfield_peptides = cgi.escape(form.getvalue("peptides"))

  print shared.toggleDisplay # Javascript function to toggle a div
  print """
<form action="/srmcollider/srmcollider.py" method="post">
    <p class='input_field'>
        <label for="peptides">Please enter the peptide sequences here:</label><br />
        <textarea id="pep_input" name="peptides" rows="20">%(textfield_peptides)s</textarea>
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
       genome_select, 'ion_series' : html_ions, 'textfield_peptides' : textfield_peptides} 

print "</div>"
print """
</div>
</body>
</html>"""
