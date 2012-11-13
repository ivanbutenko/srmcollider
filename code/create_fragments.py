#!/usr/bin/python
# -*- coding: utf-8  -*-
"""
 *
 * Program       : SRMCollider
 * Author        : Hannes Roest <roest@imsb.biol.ethz.ch>
 * Date          : 13.11.2012
"""
import srmcollider_website_helper
import DDB, collider, sys

from Residues import Residues
R = Residues('mono')
peptideparser = srmcollider_website_helper.PeptideParser(R)

###########################################################################
# Parse options
from optparse import OptionParser, OptionGroup
usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "Run options")
group.add_option("--in", dest="inputfile", type="str", help="Input file with one peptide sequence per line")
group.add_option("--out", dest="outputfile", type="str", help="Output file in mProphet format")
parser.add_option_group(group)

mycollider = collider.SRMcollider()
par = collider.SRM_parameters()
par.parse_cmdl_args(parser)
options, args = parser.parse_args(sys.argv[1:])
par.parse_options(options)
par.eval()

assumed_precursor_charge = 2
assumed_library_intensity = 0

f = open(options.inputfile)
peptides_raw = f.read()
seqs, input_sequences = peptideparser.sanitize_peptide_input(peptides_raw)
peptides = []
for s in input_sequences: 
    peptide = DDB.Peptide(); 
    peptide.set_sequence(s)
    peptide.charge = assumed_precursor_charge
    peptides.append(peptide)

peptideparser.calculate_default_fragmenation(peptides, par)


import csv
outstream = csv.writer(open(options.outputfile, "w"))
outstream.writerow( ['ModifiedSequence',
      'PrecursorCharge',  'PeptideSequence',  'PrecursorMz',
      'LibraryIntensity',  'ProductMz'] )
for seq_id, peptide in enumerate(peptides):
        transitions = [ (f.q3, f.fragment_count) for f in peptide.fragments]
        for fragment in peptide.fragments:
          outstream.writerow( [peptide.sequence, assumed_precursor_charge, peptide.sequence, peptide.charged_mass, assumed_library_intensity, fragment.q3] )

