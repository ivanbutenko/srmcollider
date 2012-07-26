#!/usr/bin/python
# -*- coding: utf-8  -*-
# vim:set fdm=marker:

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

""" {{{
This program allows to input a list of peptides with relative transition
intensity information and will output the minimal number of transitions needed
to create a unique assay.

When using mProphet files, it is also possible to use experimental intensities
to check whether the measured transitions are still sufficient to form an UIS.
Accepted inputs are srmAtlas tsv files, mProphet csv files and SPECTRAST
spectral library files.

Example Workflow: 
    1. go to https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetTransitions
    2. search for some protein, e.g. YLR304C  
    3. download file "best_peptides..." and open it in some editor 
    4. run "python runcollider.py filename background --srmatlas_tsv --max_uis 10"
}}}
"""

import MySQLdb, sys, csv, time
sys.path.append('code')
import collider, progress, precursor; from Residues import Residues
#from code import collider, Residues, progress

# some options that can be changed locally for your convenience
default_mysql = "~/.srm.cnf"
default_org_prefix = 'hroest.srmPeptides_'
default_ssrcalc = 'srmcollider.ssrcalc_pr_copy'

# Cmd line option parsing
# {{{
from optparse import OptionParser, OptionGroup
usage = "usage: %prog spectrallibrary [options]\n" +\
        " *      spectrallibrary: path to library WITHOUT the splib ending\n" +\
        "\n" 

parser = OptionParser(usage=usage)
group = OptionGroup(parser, "SRMCollider Options", "")
group.add_option("--safety", dest="safetytransitions", default=3, type="float",
    help="Number of transitions to add above the absolute minimum. " + 
    "Defaults to 3." , metavar='3')
group.add_option("-f", "--file", dest="outfile", default='outfile', metavar='out',
    help="Output file"   )
group.add_option("-l", "--log", dest="logfile", default='', metavar='LOG',
    help="Log file"   )
group.add_option("--db_prefix", dest="organism_prefix", default=default_org_prefix,
    help="The DB table will be prefix+organism"   )
group.add_option("--ssrcalc_table", dest="ssrcalc_table", default=default_ssrcalc,
    help="The DB table that holds ssrcalc values" )
group.add_option("--mapfile", dest="pepmapfile", default='', metavar='File',
    help="Text file with one sequence and SSRCalc value per line"   )
group.add_option("--csv", dest="csv", default=False, action='store_true',
                  help="""The file passed as the first option is a mprophet style csv file
 The format is a header row which MUST contain the following fields:
 'ModifiedSequence',
 'PrecursorCharge',
 'PeptideSequence',
 'PrecursorMz',
 'LibraryIntensity',
 'ProductMz'.
All other fields except 'protein_name' are ignored.
"""
                )
group.add_option("--srmatlas_tsv", dest="srmatlas_tsv", default=False, action='store_true',
    help="""The file passed as the first option is a srmatlas style tsv file
    The format is one header line and then the following fields: Protein, Pre
    AA, Sequence, Fol AA, Adj SS, Source, q1_mz, q1_chg, q3_mz,
    q3_chg, Label, Rank, RI, SSRT, n_obs
"""
                )
group.add_option("--exp", dest="exp_resultfile", default='', #metavar='out',
    help="""Experimental result file (mProphet output) that contains information,
                 which peptides could be measured with how many transitions. 
 The format is a header row which MUST contain the following fields:
'target_height'
'transition_group_id'
'transition_group_pepseq'
'target_trs'. Optionally the fields can be filtered with a threshold for target_height.
                 """ )
group.add_option("--thr", dest="threshold", default=0, 
    help="Threshold for target_height (e.g. the noise level)")
parser.add_option_group(group)

# }}}
###########################################################################

# parameter evaluation
# {{{
parameters = collider.SRM_parameters()
parameters.parse_cmdl_args(parser, default_mysql=default_mysql)
options, args = parser.parse_args(sys.argv[1:])
parameters.parse_options(options)

#local arguments
safetytransitions = options.safetytransitions
outfile = options.outfile
logfile = options.logfile
if len(options.logfile) != 0:
    sys.stdout = open(options.logfile, 'w')
pepmapfile = options.pepmapfile
libfile = args[0]
use_experimental_height = False
if options.exp_resultfile != '': use_experimental_height = True

par = parameters
parameters.dontdo2p2f = False
parameters.print_query = False
parameters.eval()

if par.max_uis == 0:
    err = "Error: --max_uis needs to be larger than 0\n"
    print err
    sys.stderr.write(err) 
    sys.exit()

# determine whether we have the C++ modules available or not
try:
    import c_getnonuis
    use_cpp = True
    print "Imported the C++ libraries with success"
except ImportError: 
    use_cpp = False
    print "Could not import the C++ libraries, please compile them."

db = MySQLdb.connect(read_default_file=par.mysql_config)
cursor = db.cursor()
try:
    cursor.execute("desc %s" % parameters.peptide_tables[0])
except Exception:
    err = "Could not access database '%s'.\n" % parameters.peptide_tables[0] +\
    "Make sure that your --db_prefix and configuration file are correct.\n" 
    print err
    sys.stderr.write(err) 
    sys.exit()

###########################################################################
# }}}
###########################################################################

# File parsing routines
from Fileparser import parse_srmatlas_file, parse_mprophet_resultfile, parse_mprophet_methodfile

# here we parse the files
# {{{
if not options.csv and not options.srmatlas_tsv:
    import speclib_db_lib
    sptxt = libfile + ".splib"
    pepidx = libfile + ".pepidx"
    file_used = sptxt

    print "Experiment Type"
    print ' '*5, 'Check transitions from a spectral library with'
    #peak into file, see how many peptides there are 
    for i,l in enumerate(open(pepidx)): pass
    if i > 10000: 
        print '\n'*5, "Do you really want to create assays for %s peptides?" % i
        print "Please filter the spectral library first."
        print '\n'*1, "Sorry, I give  up."
        sys.exit(9)

    library_key = 4242 
    library = speclib_db_lib.Library(library_key)
    library.read_sptxt_pepidx( sptxt, pepidx, library_key )
elif options.csv: 
    print "Experiment Type"
    print ' '*5, 'Check transitions from a csv transitions list with'
    file_used = libfile
    library = parse_mprophet_methodfile(libfile)
elif options.srmatlas_tsv: 
    print "Experiment Type"
    print ' '*5, 'Check transitions from a csv transitions list with'
    file_used = libfile
    library = parse_srmatlas_file(libfile)

print ' '*5, par.thresh
print ' '*5, 'Da for the q3 transitions.'
print ' '*5, 'Using file : %s' % file_used
print ' '*5, 'Using background organism(s): %s' % par.peptide_tbl_identifier
print ' '*5,  "with %s modifications and %s missed cleavages" % (par.max_mods, par.max_MC)



if use_experimental_height:
    infile = options.exp_resultfile
    threshold = float(options.threshold)
    parse_mprophet_resultfile(library, infile, threshold)

# }}}

###########################################################################
# Start program
###########################################################################

##
## First create the mapping of the SSRcalc values to the peptide sequences 
## {{{
mycollider = collider.SRMcollider()
pepmap = {}
seqs = ''
seqdic = {}
for icount, spectrum in enumerate(library):
    seqs += "'%s'," % spectrum.sequence
    seqdic[ spectrum.sequence ] = 0

# try to find SSRCalc values
try:
    ssr_query = """
    select sequence, ssrcalc
    from %(ssrcalc_table)s
    where sequence in (%(seqs)s)
    """ % { 'seqs' : seqs[:-1], 'ssrcalc_table' : options.ssrcalc_table}
    cursor.execute( ssr_query )
    pepmap = dict( cursor.fetchall() )
except Exception: pepmap = {}

# add more SSRCalc values from the pepmap file
if pepmapfile != '':
    try:
        f = open(pepmapfile)
        for l in f:
            if l != '': pepmap[ l.split()[0] ] = l.split()[1]
    except Exception:
        print "Could not read pepmapfile. Please provide a file that has one "+\
              "sequence and SSRCalc value per line separated by a space."

print "%s peptide precursors provided" % icount
print "%s unique peptide sequences provided" % len(seqdic)
print "%s unique peptide sequences have SSRCalc values " % len(pepmap)
if len(seqdic) > len(pepmap):
    print "If you want to use SSRCalc predictions, please provide a file with SSRCalc values."
    print "Otherwise your results will be wrong because you are missing some SSRcalc values."

#}}}

## 
## START the main loop
## {{{
progressm = progress.ProgressMeter(total=icount+1, unit='peptides')
get_precursor_time = 0
get_min_tr_time = 0
for counter,spectrum in enumerate(library):
    tmp_time = time.time()
    spectrum.score = -99
    spectrum.min_needed = -1
    # Get the spectrum peaks, sort by intensity
    # Assign a ssrcalc value to (if not present, use default value).
    # Then create a peptide hash and get the transitions
    peaks = spectrum.get_peaks()
    peaks.sort( lambda x,y: -cmp(x.intensity, y.intensity))
    try: ssrcalc = pepmap[spectrum.sequence]
    except KeyError: ssrcalc = 25
    #
    #
    peptide_obj = precursor.Precursor(
        modified_sequence=spectrum.modified_sequence,
        parent_id=-1,
        transition_group=-1,
        q1=spectrum.precursorMZ,
        q1_charge=spectrum.name.split('/')[1],
        ssrcalc=ssrcalc
    )
    transitions = [ (p.peak, i) for i,p in enumerate(peaks)]
    if use_experimental_height: 
        transitions = [ (p.peak, i) for i,p in enumerate(peaks) if p.experimental_height > threshold]
    nr_transitions = len( transitions )
    if nr_transitions == 0: continue #no transitions in this window
    if (nr_transitions > 31):
        print "Skipping", spectrum.modified_sequence, "too many transitions (%s)" % nr_transitions
        continue

    peptide_obj.fix_mprophet_sequence_bug() # fix mProphet-bug (modified cysteins are [C160] instead of C[160]

    #
    # Get all interfering precursors, (w/o the current peptide)
    precursors = mycollider._get_all_precursors(par, peptide_obj, cursor)

    R = Residues('mono')
    q3_low, q3_high = par.get_q3range_collisions()
    #
    # check for q3_low and q3_high values, if the transitions are outside this
    # range, we will not find collisions for those == incorrect results
    #
    for t in transitions: 
        if t[0] < q3_low or t[0] > q3_high:
            err = "\n"+ "Error: Found transition %s out of range\n" % t[0] + \
            "Please adjust the --q3_low and --q3_high parameters to get correct results\n"
            print err
            sys.stderr.write(err) 
            sys.exit()

    get_precursor_time += time.time() - tmp_time; tmp_time = time.time()

    if not use_experimental_height and use_cpp:
        # We dont have experimental height data and use C++ code
        import c_integrated
        old_prec = [(0,p.modified_sequence) for p in precursors]
        min_needed = c_integrated.getMinNeededTransitions(tuple(transitions), tuple(old_prec), 
            par.max_uis, par.q3_window, par.ppm, par)
    elif not use_experimental_height:
        # We dont have experimental height data and cannot use C++ code
        collisions_per_peptide = collider.get_coll_per_peptide(mycollider, 
            transitions, par, peptide_obj, cursor)
        min_needed = mycollider._sub_getMinNeededTransitions(par, transitions, collisions_per_peptide)
        #min_needed = mycollider.getMinNeededTransitions_direct(par, transitions, precursors)
    else:
        # here we consider the case that we have measured a number of
        # transitions experimentally and want to know how many of them are
        # sufficient to establish uniqueness. For this, all we need is
        # that one tuple of transitions establishes uniqueness since we
        # were able to measure it above the background noise.
        collisions_per_peptide = collider.get_coll_per_peptide(mycollider, 
            transitions, par, pep, cursor)
        for order in range(1,nr_transitions+1): 
            mymax = collider.choose(nr_transitions, order)
            if use_cpp: non_uis = c_getnonuis.get_non_uis(collisions_per_peptide, order)
            else: 
                non_uis = set()
                for pepc in collisions_per_peptide.values():
                    get_non_uis(pepc, non_uis, i)
            if len(non_uis) < mymax: break
        if len(non_uis) < mymax: min_needed  = order
        else: min_needed = -1
    spectrum.score = min_needed * nr_transitions
    spectrum.min_needed = min_needed
    if min_needed != -1: spectrum.score = nr_transitions - min_needed
    if not par.quiet: progressm.update(1)
    get_min_tr_time += time.time() - tmp_time; tmp_time = time.time()


# }}}

totaltime = get_precursor_time + get_min_tr_time
print "Time for getting the precursors %0.2fs (%.0f%%) vs the time to calculate the minimal nr transitions %0.2fs (%.0f%%)" % (
    get_precursor_time, 100*get_precursor_time/totaltime, get_min_tr_time, 100*get_min_tr_time/totaltime,)
##
## PRINT some statistics of our results
## For each precursor, peptide and protein print the necessary transitions
## {{{
if not use_experimental_height:
    # default case
    mycollider.perprecursor = [0 for i in range(par.max_uis+1)]
    mycollider.perpeptide = [0 for i in range(par.max_uis+1)]
    mycollider.perprotein = [0 for i in range(par.max_uis+1)]

    precursors = {}
    peptides = {}
    proteins = {}
    for spectrum in library:
        peptides[spectrum.sequence] = []
        proteins[spectrum.protein] = []

    for spectrum in library:
        if spectrum.min_needed == -1: continue
        peptides[spectrum.sequence].append( spectrum.min_needed )
        proteins[spectrum.protein].append( spectrum.min_needed )
        try:
          mycollider.perprecursor[ spectrum.min_needed ] += 1
        except IndexError:
            sys.exit( "Error, increase --max_uis to at least %s" % spectrum.min_needed )

    for k,v in proteins.iteritems():
        if len(v) == 0: continue
        mycollider.perprotein[ min(v) ] += 1
    for k,v in peptides.iteritems():
        if len(v) == 0: continue
        mycollider.perpeptide[ min(v) ] += 1

    #write report
    sumtrans = 0
    for i in range(par.max_uis+1):
        print "Precursors needing %s transitions: " % i, mycollider.perprecursor[i], \
                '\t(Peptides %s)' % mycollider.perpeptide[i], \
                '\t(Proteins %s)' % mycollider.perprotein[i]
        sumtrans += i*mycollider.perprecursor[i]

    measured_precursors = len(library)
    measured_peptides = len(dict( [(s.sequence,0) for s in library]))
    measured_proteins = len(dict( [(s.protein,0) for s in library ]))
    print "Precursors not possible to observe (too few transitions):  %s / %s" %( 
        measured_precursors - sum( mycollider.perprecursor), measured_precursors)
    print "Peptides not possible to observe (too few transitions):  %s / %s" % (
        measured_peptides - sum( mycollider.perpeptide) , measured_peptides)
    print "Proteins not possible to observe (too few transitions):  %s / %s" % (
        measured_proteins - sum( mycollider.perprotein) , measured_proteins)

    print "Minimal number of transitions needed to measure these %s peptides:" % (
        icount - measured_precursors + sum(mycollider.perprecursor) ), sumtrans
else:
    # we have experimental data
    mycollider.perprecursor = [0 for i in range(par.max_uis+1)]
    mycollider.perpeptide = [0 for i in range(par.max_uis+1)]
    mycollider.perprotein = [0 for i in range(par.max_uis+1)]

    precursors = {}
    peptides = {}
    proteins = {}
    for spectrum in library:
        peptides[spectrum.sequence] = []
        proteins[spectrum.protein] = []

    for spectrum in library:
        if spectrum.min_needed == -1: continue
        peptides[spectrum.sequence].append( spectrum.score )
        proteins[spectrum.protein].append( spectrum.score )
        mycollider.perprecursor[ spectrum.score ] += 1

    for k,v in proteins.iteritems():
        if len(v) == 0: continue
        mycollider.perprotein[ max(v) ] += 1
    for k,v in peptides.iteritems():
        if len(v) == 0: continue
        mycollider.perpeptide[ max(v) ] += 1

    #write report
    sumtrans = 0
    for i in range(par.max_uis+1):
        print "Precursors having %s extra transitions: " % i, mycollider.perprecursor[i], \
                '\t(Peptides %s)' % mycollider.perpeptide[i], \
                '\t(Proteins %s)' % mycollider.perprotein[i]
        sumtrans += i*mycollider.perprecursor[i]

    measured_precursors = len([s for s in library 
        if s.peaks[0].__dict__.has_key('measured')])
    measured_peptides = len(dict( [(s.sequence,0) for s in library 
        if s.peaks[0].__dict__.has_key('measured')]))
    measured_proteins = len(dict( [(s.protein,0) for s in library 
        if s.peaks[0].__dict__.has_key('measured')]))

    print "Precursors not possible to observe (too few transitions):  %s / %s" %( 
        measured_precursors - sum( mycollider.perprecursor), measured_precursors)
    print "Peptides not possible to observe (too few transitions):  %s / %s" % (
        measured_peptides - sum( mycollider.perpeptide) , measured_peptides)
    print "Proteins not possible to observe (too few transitions):  %s / %s" % (
        measured_proteins - sum( mycollider.perprotein) , measured_proteins)

# }}}

# here we write to the outfile
# {{{
f = open(outfile + '_transitions.csv', 'w')
w = csv.writer(f)
for spectrum in library:
    peaks = spectrum.get_peaks()
    peaks.sort( lambda x,y: -cmp(x.intensity, y.intensity))
    for i,p in enumerate(peaks):
        if i > int(spectrum.min_needed-1 + safetytransitions): break
        if not options.csv and not options.srmatlas_tsv:
            w.writerow( [spectrum.precursorMZ, p.peak, p.peak_annotation, spectrum.name])
        else: w.writerow(p.mprophetline)

f.close()
f = open(outfile + '_peptides.csv', 'w')
w = csv.writer(f)
for spectrum in library:
    w.writerow( [spectrum.min_needed, spectrum.name, spectrum.protein])

f.close()
print "Wrote transition list into file ", outfile

#}}}

