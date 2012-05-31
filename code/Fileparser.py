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

import csv

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


# File parsing routines
def parse_srmatlas_file(csvfile):
    myfile = open(csvfile, 'rU') 
    dialect = csv.Sniffer().sniff(myfile.read(1024))
    myfile.seek(0)
    r = csv.reader(myfile, dialect)
    header = r.next()
    # We know the header
    headerdict = {}
    headerdict['ModifiedSequence'] = 2
    headerdict['PrecursorCharge'] = 7
    headerdict['PrecursorMz'] = 6
    headerdict['ProductMz'] = 8
    headerdict['frg_z'] = 9
    headerdict['frg_type'] = 10
    headerdict['LibraryIntensity'] = 11
    headerdict['protein_name'] = 0
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
        s.sequence    = lines[0][headerdict['ModifiedSequence']].replace( 'C[160]', 'C').replace('C[+57]', 'C')
        s.modified_sequence    = lines[0][headerdict['ModifiedSequence']].replace('C[+57]', 'C[160]')
        s.precursorMZ = float(lines[0][headerdict['PrecursorMz']])
        s.protein = lines[0][headerdict['protein_name']]
        for peak in lines:
            p = Peak()
            p.peak             = float(peak[headerdict['ProductMz']])
            p.intensity        = 10000.0 * 1.0/float(peak[headerdict['LibraryIntensity']]) 
            p.type             = peak[headerdict['frg_type']] + '_' +  peak[headerdict['frg_z']]
            p.mprophetline     = peak
            s.peaks.append(p)
        library.append(s)
    return library
 
def parse_mprophet_methodfile(csvfile):

    myfile = open(csvfile, 'rU') 
    dialect = csv.Sniffer().sniff(myfile.read(1024))
    myfile.seek(0)
    r = csv.reader( myfile, dialect)
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
        s.modified_sequence    = lines[0][headerdict['ModifiedSequence']].replace('C[+57]', 'C[160]')
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
    return library

def parse_mprophet_resultfile(library, infile, threshold):
    # parse the experimental infile 
    # here we try to map back the outfile of the mProphet to the
    # originally used transition list. Of course this fails sometimes.
    # libdict cotains the PEPTIDE.2 identifier for sequence.charge_state
    # the column "target_trs" should contain y7_1 typeNr_charge
    exp_peptides = {}
    libdict = dict([ (s.name.replace('/', '.'), s) for s in library ] )
    myfile = open(infile, 'rU') 
    dialect = csv.Sniffer().sniff(myfile.read(1024))
    myfile.seek(0)
    experimental_infile = csv.reader(myfile, dialect)

    header = experimental_infile.next()
    headerdict = dict([(t,i) for i,t in enumerate(header)])
    headerdict['target_height']
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
    print "found %s peptides experimentally" % len( exp_peptides)
    print "found %s peptides not in the library input file" % len( sequence_n_found)
    print "found %s transitions not in the library input file" % len(tr_n_found)
    #
    err = 'Could not find the following peptides:\n'
    for s in sequence_n_found: err += s + '\n'
    err += 'Could not find the following transitions:\n'
    for t in tr_n_found: err += t[0].sequence + '\t'  + t[1] + '\n'
    sys.stderr.write(err) 


