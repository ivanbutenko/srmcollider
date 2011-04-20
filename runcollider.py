#!/usr/bin/python
import MySQLdb, time, sys
from string import Template
import numpy
sys.path.append( '/IMSB/users/hroest/projects')
sys.path.append( '/IMSB/users/hroest/projects/hlib')
sys.path.append( '/IMSB/users/hroest/projects/msa/code/')
sys.path.append( '/IMSB/users/hroest/projects/srm_clashes/code/')
db = MySQLdb.connect(read_default_file="/IMSB/users/hroest/.srm.cnf")
cursor = db.cursor()

import collider
import speclib_db_lib
import progress
import silver
import DDB 
import csv

from optparse import OptionParser, OptionGroup
usage = "usage: %prog spectrallibrary backgroundorganism [options]\n" +\
        " *      spectrallibrary: path to library WITHOUT the splib ending\n" +\
        " *      backgroundorganism: e.g. 'yeast' or 'human', others need to be requested " +\
        "\n" + \
        "This program overrides the options --q3_low and --q3_high default values" + \
        "with values of 0 to "

    
parser = OptionParser(usage=usage)

group = OptionGroup(parser, "Spectral library Options", "")
group.add_option("--safety", dest="safetytransitions", default=3, type="float",
                  help="Number of transitions to add above the absolute minimum. " + 
                  "Defaults to 3" , metavar='3')
group.add_option("-f", "--file", dest="outfile", default='outfile', metavar='out',
                  help="Output file"   )
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
All other fields are ignored.
"""
                )
group.add_option("--exp", dest="exp_resultfile", default='', #metavar='out',
                  help="""Experimental result file (mProphet output) that contains information, which peptides could be measured with how many transitions. 
 The format is a header row which MUST contain the following fields:
'target_height'
'transition_group_id'
'transition_group_pepseq'
'target_trs'. Optionally the fields can be filtered with a threshold for target_height.
                 """ )
group.add_option("--thr", dest="threshold", default=0, 
                 help="Threshold for target_height (e.g. the noise level)")
parser.add_option_group(group)




#parameters evaluation
parameters = collider.SRM_parameters()
parameters.parse_cmdl_args(parser)
options, args = parser.parse_args(sys.argv[1:])

#local arguments
safetytransitions = options.safetytransitions
outfile = options.outfile
pepmapfile = options.pepmapfile
libfile = args[0]
peptide_table = 'hroest.srmPeptides_'  + args[1]

par = parameters
par.__dict__.update( options.__dict__ )
par.q3_range = [options.q3_low, options.q3_high]
par.q1_window /= 2.0
par.q3_window /= 2.0
par.ssrcalc_window /= 2.0
if par.ppm == 'True': par.ppm = True
elif par.ppm == 'False': par.ppm = False

parameters.peptide_table = peptide_table
parameters.dontdo2p2f = False
parameters.eval()
sptxt = libfile + ".splib"
pepidx = libfile + ".pepidx"
print par.get_common_filename()
#print par.experiment_type

if not options.csv:
    print "Experiment Type"
    print ' '*5, 'Check transitions from a spectral library with'
    print ' '*5, par.thresh
    print ' '*5, 'Da for the q3 transitions.'
    print ' '*5, 'Using spectral library: %s' % sptxt
    print ' '*5, 'Using background organism: %s' % args[1]

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


else: 
    csvfile = libfile
    print "Experiment Type"
    print ' '*5, 'Check transitions from a csv transitions list with'
    print ' '*5, par.thresh
    print ' '*5, 'Da for the q3 transitions.'
    print ' '*5, 'Using csv file library: %s' % csvfile
    print ' '*5, 'Using background organism: %s' % args[1]

    # some dummy classes that can do the same thing as the spectral library
    # classes
    class Peak():
        def __init__(self):
            self.experimental_height  = 0

    class Spectrum():
        def __init__(self):
            self.peaks = []
            self.sequence = None
            self.name = None
            self.precursorMZ = None
            self.protein = ''
        def get_peaks(self):
            return self.peaks

    r = csv.reader( open(csvfile, 'rU') , delimiter='\t') 
    header = r.next()
    headerdict = dict([(t,i) for i,t in enumerate(header)])
    specdict = {}
    for line in r:
        key = line[ headerdict['ModifiedSequence']] + '/' + line[headerdict['PrecursorCharge']]
        if specdict.has_key( key ): specdict[key].append(line)
        else: specdict[key] = [ line ]

    library = []
    for key in specdict:
        lines = specdict[key]
        s = Spectrum()
        s.name        = key
        s.sequence    = lines[0][headerdict['PeptideSequence']]
        s.precursorMZ = float(lines[0][headerdict['PrecursorMz']])
        try: s.protein = lines[0][headerdict['protein_name']]
        except Exception: pass
        for peak in lines:
            p = Peak()
            try:
                p.peak             = float(peak[headerdict['ProductMz']])
                p.intensity        = float(peak[headerdict['LibraryIntensity']])
            except Exception: p.intensity = 0
            #
            try: p.type        = peak[headerdict['frg_type']] + peak[headerdict['frg_nr']] + '_' +  peak[headerdict['frg_z']]
            except Exception: pass
            p.mprophetline = peak
            s.peaks.append(p)
        library.append(s)


###########################################################################
###########################################################################
use_experimental_height = False
if options.exp_resultfile != '': use_experimental_height = True
if use_experimental_height:
    #parse the experimental infile 
    # here we try to map back the outfile of the mProphet to the
    # originally used transition list. Of course this fails sometimes.
    # libdict cotains the PEPTIDE.2 identifier for sequence.charge_state
    # the column "target_trs" should contain y7_1 typeNr_charge
    exp_peptides = {}
    libdict = dict([ (s.name.replace('/', '.'), s) for s in library ] )
    # infile = '/IMSB/users/hroest/detected_assays_depl_plasma_trs_selected_columns.csv'
    # infile = '/IMSB/users/hroest/detected_assays_depl_urine_trs_selected_columns.csv'
    # outfile = '/tmp/outtest'
    # threshold = 200
    infile = options.exp_resultfile
    threshold = float(options.threshold)
    #
    experimental_infile = csv.reader( open(infile, 'rU') , delimiter='\t') 
    print infile
    header = experimental_infile.next()
    headerdict = dict([(t,i) for i,t in enumerate(header)])
    headerdict['target_height']
    ## #here we write to the outfile
    ## f = open(outfile, 'w')
    ## w = csv.writer(f)
    ## w.writerow(header)
    sequence_n_found = []
    tr_n_found = []
    for l in experimental_infile:
        if float(l[ headerdict['target_height'] ] ) > threshold:
            gr_id = l[ headerdict['transition_group_id'] ]
            charge = gr_id.split('.')[1]
            gr_id = l[ headerdict['transition_group_pepseq'] ] + '.' + charge
            exp_peptides[ l[ headerdict['transition_group_pepseq'] ] ] = 0
            try:
                pp = [p for p in libdict[ gr_id ].peaks 
                      if p.type == l[ headerdict['target_trs'] ] ]
                assert len(pp) == 1
                pp[0].experimental_height = l[headerdict['target_height']]
                pp[0].measured = True
            except KeyError: sequence_n_found.append( gr_id)
            except AssertionError: tr_n_found.append( [
                libdict[ gr_id ],
                l[ headerdict['target_trs'] ], 
                [p.type for p in libdict[ gr_id ].peaks]
            ]  )
            # w.writerow(l)
    print "found %s peptides experimentally" % len( exp_peptides)
    print "found %s peptides not in the library input file" % len( sequence_n_found)
    print "found %s transitions not in the library input file" % len(tr_n_found)
    f = open('error.log', 'w')
    f.write('Could not find the following peptides:\n')
    for s in sequence_n_found:
        f.write(s + '\n')
    f.write('Could not find the following transitions:\n')
    for t in tr_n_found:
        f.write( t[0].sequence + '\t'  + t[1] + '\n')
    f.close()

###########################################################################

mycollider = collider.SRMcollider()
mycollider.min_transitions = []

pepmap = {}
seqs = ''
seqdic = {}
for icount, spectrum in enumerate(library):
    seqs += "'%s'," % spectrum.sequence
    seqdic[ spectrum.sequence ] = 0

ssr_query = """
select sequence, ssrcalc
from hroest.ssrcalc_pr_copy
where sequence in (%(seqs)s)
""" % { 'seqs' : seqs[:-1], }
cursor.execute( ssr_query )
pepmap = dict( cursor.fetchall() )

#add more SSRCalc values from the pepmap file
if pepmapfile != '':
    try:
        f = open(pepmapfile)
        for l in f:
            if l != '': pepmap[ l.split()[0] ] = l.split()[1]
    except Exception:
        print "Could not read pepmapfile. Please provide a file that has one "+\
              "sequence and SSRCalc value per line separated by a space."

start = time.time()
print "%s peptide precursors provided" % icount
print "%s unique peptide sequences provided" % len(seqdic)
print "%s unique peptide sequences have SSRCalc values " % len(pepmap)
if len(seqdic) > len(pepmap):
    print "If you want to use SSRCalc predictions, please provide a file with SSRCalc values."
    print "Otherwise your results will be wrong because you are missing some SSRcalc values."

if par.max_uis == 0:
    err = "Error: --max_uis needs to be larger than 0"
    print err
    sys.stderr.write(err) 
    sys.exit()

parameters.aions      =  True
parameters.aMinusNH3  =  True
parameters.bMinusH2O  =  True
parameters.bMinusNH3  =  True
parameters.bPlusH2O   =  True
parameters.yMinusH2O  =  True
parameters.yMinusNH3  =  True
parameters.cions      =  True
parameters.xions      =  True
parameters.zions      =  True










progressm = progress.ProgressMeter(total=icount+1, unit='peptides')
#mycollider.minstats = [0 for i in range(par.max_uis+1)]
#mycollider.persequence = dict([ (s, []) for s in seqdic])
#mycollider.peptides = dict([ (s, []) for s in seqdic])
for counter,spectrum in enumerate(library):
    spectrum.score = -99
    spectrum.min_needed = -1
    peaks = spectrum.get_peaks()
    peaks.sort( lambda x,y: -cmp(x.intensity, y.intensity))
    bgpar = parameters
    par = parameters
    msequence = spectrum.name.split('/')[0]
    try: ssrcalc = pepmap[spectrum.sequence]
    except KeyError: ssrcalc = 25
    pep = {
        'mod_sequence' :   spectrum.name.split('/')[0],
        'parent_id' :  -1,
        'q1' :         spectrum.precursorMZ,
        'q1_charge' :  spectrum.name.split('/')[1],
        'ssrcalc' :    ssrcalc,
        'peptide_key' :-1,
    }
    transitions = [ (p.peak, i) for i,p in enumerate(peaks)]
    if use_experimental_height: transitions = [ (p.peak, i) for i,p in enumerate(peaks) if p.experimental_height > threshold]
    nr_transitions = len( transitions )
    if nr_transitions == 0: continue #no transitions in this window
    pep['mod_sequence'] = pep['mod_sequence'].replace('[C160]', 'C[160]')
    precursors = mycollider._get_all_precursors(par, pep, cursor)
    precursors = [p for p in precursors if p[1] != pep['mod_sequence'] ]
    R = silver.Residues.Residues('mono')
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
    #
    if False:
        collisions = list( mycollider._get_all_collisions_calculate_sub(precursors, par, R, q3_low, q3_high) )
        min_needed = mycollider._getMinNeededTransitions(par, transitions, collisions)
    if True:
        import c_getnonuis
        collisions_per_peptide = c_getnonuis.calculate_collisions_per_peptide(
            tuple(transitions), tuple(precursors), q3_low, q3_high, 
            par.q3_window, par.ppm)
        if not use_experimental_height:
            min_needed = mycollider._sub_getMinNeededTransitions(par, transitions, collisions_per_peptide)
        else:
            # here we consider the case that we have measured a number of
            # transitions experimentally and want to know how many of them are
            # sufficient to establish uniqueness. For this, all we need is
            # that one tuple of transitions establishes uniqueness since we
            # were able to measure it above the background noise.
            for order in range(1,nr_transitions+1): 
                mymax = collider.choose(nr_transitions, order)
                non_uis = c_getnonuis.get_non_uis(collisions_per_peptide, order)
                if len(non_uis) < mymax: break
            if len(non_uis) < mymax: min_needed  = order
            else: min_needed = -1
    #mycollider.min_transitions.append( [spectrum.name, min_needed, nr_transitions] )
    spectrum.score = min_needed * nr_transitions
    spectrum.min_needed = min_needed
    if min_needed != -1: spectrum.score = nr_transitions - min_needed
    #    #mycollider.minstats[ min_needed ] += 1
    #    #mycollider.persequence[ spectrum.sequence ].append( spectrum.score)
    #mycollider.peptides[ spectrum.sequence ].append( spectrum )
    end = time.time()
    if not par.quiet: progressm.update(1)


sys.exit()

if False:
    print '\n'
    print numpy.histogram([s.score for s in library if s.score > -99], 20, (-10,10) )
    print mycollider.persequence
    mycollider.persequenceminstats = [0 for i in range(par.max_uis+1)]
    for k,v in mycollider.persequence.iteritems():
        if len(v) == 0: continue
        #print k, min(v)
        mycollider.persequenceminstats[ min(v) ] += 1

    mycollider.perprotein = [0 for i in range(par.max_uis+1)]
    try: 
        proteins = {}
        for spectrum in library:
            #print spectrum.protein
            if spectrum.min_needed == -1: continue
            if proteins.has_key( spectrum.protein ): proteins[spectrum.protein].append( spectrum.score )
            else: proteins[ spectrum.protein ]  = [ spectrum.score]
        for k,v in proteins.iteritems():
            if len(v) == 0: continue
            #print k, min(v)
            mycollider.perprotein[ min(v) ] += 1

    except Exception: pass
    # print mycollider.persequenceminstats
    # print sum(mycollider.persequenceminstats)
    # print len(seqdic) - sum(mycollider.persequenceminstats)

    #write report
    sumtrans = 0
    for i,v in enumerate(mycollider.minstats):
        print "Precursors needing %s transitions: " % i, v, \
                '\t(Peptides %s)' % mycollider.persequenceminstats[i], \
                '\t(Proteins %s)' % mycollider.perprotein[i]
        sumtrans += i*v

    not_observeable = len( [0 for a,b,c in mycollider.min_transitions if b == -1])
    print "Precursors not possible to observe (too few transitions): ",  not_observeable
    print "Proteins not possible to observe (too few transitions): ",  \
            len(dict( [(s.protein,0) for s in library])) - sum( mycollider.perprotein)
    print "Peptides not possible to observe (too few transitions): ",  len(seqdic) - sum(mycollider.persequenceminstats)


if not use_experimental_height:
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
        mycollider.perprecursor[ spectrum.min_needed ] += 1

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

    ## a, b = numpy.histogram([c for a,b,c in mycollider.min_transitions if b == -1], 10, (0,10) )
    ## print "Transition distribution of not observeable precursors"
    ## print a, b
    print "Minimal number of transitions needed to measure these %s peptides:" % (
        icount - measured_precursors + sum(mycollider.perprecursor) ), sumtrans
else:
    #use experimental height
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

    # not_observeable = len( [0 for a,b,c in mycollider.min_transitions if b == -1])
    # print "Precursors not possible to observe (too few transitions): %s / %s" % (
    #   len([s for s in library if s.score < 0]) , len(library))
    # print "Peptides not possible to observe (too few transitions): %s / %s" % (
    #         len(peptides)- sum( mycollider.perpeptide), len(peptides) )
    # print "Proteins not possible to observe (too few transitions): %s / %s" % (
    #         len(proteins)- sum( mycollider.perprotein), len(proteins) )

    # #print "Peptides not possible to observe (too few transitions): ",  len(seqdic) - sum(mycollider.persequenceminstats)

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


#here we write to the outfile
f = open(outfile, 'w')
w = csv.writer(f)
for spectrum in library:
    peaks = spectrum.get_peaks()
    peaks.sort( lambda x,y: -cmp(x.intensity, y.intensity))
    for i,p in enumerate(peaks):
        if i > int(spectrum.min_needed-1 + safetytransitions): break
        if not options.csv: w.writerow( [spectrum.precursorMZ, p.peak, p.peak_annotation, spectrum.name])
        else: w.writerow(p.mprophetline)

f.close()
print "Wrote transition list into file ", outfile



"""

drop table hroest.ssrcalc_pr_copy ;
create table 
hroest.ssrcalc_pr_copy as
select * from compep.ssrcalc_prediction;
alter table hroest.ssrcalc_pr_copy add index(sequence);

"""
