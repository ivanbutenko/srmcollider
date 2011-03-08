#!/usr/bin/python
import MySQLdb, time, sys
from string import Template
sys.path.append( '/IMSB/users/hroest/projects')
sys.path.append( '/IMSB/users/hroest/projects/hlib')
sys.path.append( '/IMSB/users/hroest/projects/msa/code/')
sys.path.append( '/IMSB/users/hroest/projects/srm_clashes/code/')
db = MySQLdb.connect(read_default_file="/IMSB/users/hroest/.srm.cnf")
cursor = db.cursor()

import collider
import utils
import speclib_db_lib
import progress
import silver
import DDB 

from optparse import OptionParser, OptionGroup
usage = "usage: %prog spectrallibrary backgroundorganism [options]\n" +\
        " *      spectrallibrary: path to library WITHOUT the splib ending\n" +\
        " *      backgroundorganism: e.g. 'yeast' or 'human', others need to be requested "
parser = OptionParser(usage=usage)

group = OptionGroup(parser, "Spectral library Options", "")
group.add_option("--safety", dest="safetytransitions", default=3, type="float",
                  help="Number of transitions to add above the absolute minimum. " + 
                  "Defaults to 3" , metavar='3')
group.add_option("-f", "--file", dest="outfile", default='outfile',
                  help="Output file"   )
group.add_option("--mapfile", dest="pepmapfile", default='',
                  help="Text file with one sequence and SSRCalc value per line"   )
group.add_option("--csv", dest="csv", default=False, action='store_true',
                  help="File given by first option is a mprophet style csv file")
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
    import csv
    print "Experiment Type"
    print ' '*5, 'Check transitions from a csv transitions list with'
    print ' '*5, par.thresh
    print ' '*5, 'Da for the q3 transitions.'
    print ' '*5, 'Using csv file library: %s' % csvfile
    print ' '*5, 'Using background organism: %s' % args[1]

    #another way to do it, parsing csv files
    class Peak():
        def __init__(self): pass

    class Spectrum():
        def __init__(self):
            self.peaks = []
            self.sequence = None
            self.name = None
            self.precursorMZ = None
        def get_peaks(self):
            return self.peaks

    #r = csv.reader( open('/IMSB/users/hroest/mprophet.txt', 'rU') , delimiter='\t') 
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
        for peak in lines:
            p = Peak()
            p.intensity        = float(peak[headerdict['LibraryIntensity']])
            p.peak             = float(peak[headerdict['ProductMz']])
            p.mprophetline = peak
            s.peaks.append(p)
        library.append(s)



mycollider = collider.SRMcollider()
mycollider.min_transitions = []

pepmap = {}
seqs = ''
for icount, spectrum in enumerate(library):
    seqs += "'%s'," % spectrum.sequence

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
print "%s Peptides provided" % icount
print "%s Peptides have SSRCalc values " % len(pepmap)
if icount > len(pepmap):
    print "If you want to use SSRCalc predictions, please provide a file with SSRCalc values."
    print "Otherwise your results will be wrong because you are missing some SSRcalc values."

if par.max_uis == 0:
    print "max uis needs to be larger than 0"
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
mycollider.minstats = [0 for i in range(par.max_uis+1)]
for spectrum in library:
    peaks = spectrum.get_peaks()
    peaks.sort( lambda x,y: -cmp(x.intensity, y.intensity))
    bgpar = parameters
    par = parameters
    msequence = spectrum.name.split('/')[0]
    try: ssrcalc = pepmap[spectrum.sequence]
    except KeyError: ssrcalc = 25
    pep = {
        'sequence' :   spectrum.name.split('/')[0],
        'parent_id' :  -1,
        'q1' :         spectrum.precursorMZ,
        'q1_charge' :  spectrum.name.split('/')[1],
        'ssrcalc' :    ssrcalc,
        'peptide_key' :-1,
    }
    #
    transitions = [ (p.peak, i) for i,p in enumerate(peaks)]
    nr_transitions = len( transitions )
    if nr_transitions == 0: continue #no transitions in this window
    precursors = mycollider._get_all_precursors(par, pep, cursor)
    precursors = [p for p in precursors if p[1] != pep['sequence'] ]
    R = silver.Residues.Residues('mono')
    q3_low, q3_high = par.get_q3range_collisions()
    collisions = list( mycollider._get_all_collisions_calculate_sub(precursors, par, R, q3_low, q3_high) )
    #print len(collisions)
    nr_used_tr = min(par.max_uis+1, nr_transitions)
    mytransitions = transitions[:nr_used_tr]
    collisions_per_peptide = {}
    q3_window_used = par.q3_window
    for t in mytransitions:
        if par.ppm: q3_window_used = par.q3_window * 10**(-6) * t[0]
        for c in collisions:
            if abs( t[0] - c[0] ) <= q3_window_used:
                #gets all collisions
                if collisions_per_peptide.has_key(c[3]):
                    if not t[1] in collisions_per_peptide[c[3]]:
                        collisions_per_peptide[c[3]].append( t[1] )
                else: collisions_per_peptide[c[3]] = [ t[1] ] 
    #take the top j transitions and see whether they, as a tuple, are
    #shared
    min_needed = par.max_uis
    for j in range(par.max_uis,0,-1): 
        #take top j transitions
        mytransitions = tuple(sorted([t[1] for t in transitions[:j]]))
        unuseable = False
        for k,v in collisions_per_peptide.iteritems():
            if tuple(sorted(v)[:j]) == mytransitions: unuseable=True
            #if len(v) > 1: print 'v', v
        if not unuseable: min_needed = j
        #print j, unuseable, min_needed
    #print 'min needed', min_needed
    mycollider.min_transitions.append( [spectrum.name, min_needed] )
    mycollider.minstats[ min_needed ] += 1
    spectrum.min_needed = min_needed
    #if min_needed  > 1: print spectrum.name
    #if spectrum.name == 'FAGGPSDALAALSK/2': pass
    #if spectrum.name == 'LHATHLAAEATK/3':
    #    print ""
    #    print "asdfasd", min_needed
    #    print len(collisions)
    #    for k,v in collisions_per_peptide.iteritems():
    #        if len(v) > 2: print v
    end = time.time()
    if not par.quiet: progressm.update(1)


#write report
sumtrans = 0
for i,v in enumerate(mycollider.minstats):
    print "Peptides needing %s transitions: " % i, v
    sumtrans += i*v

print "Minimal number of transitions needed to measure these %s peptides:" % icount, sumtrans

#here we write to the outfile
import csv
f = open(outfile, 'w')
w = csv.writer(f)
for spectrum in library:
    peaks = spectrum.get_peaks()
    peaks.sort( lambda x,y: -cmp(x.intensity, y.intensity))
    if spectrum.name == 'LHATHLAAEATK/3':
        #print "asdfasd", spectrum.min_needed
        pass
    for i,p in enumerate(peaks):
        #i starts at 0
        if spectrum.name == 'LHATHLAAEATK/3':
            #print i , int(spectrum.min_needed-1 + safetytransitions)
            #print p.mprophetline
            pass
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
