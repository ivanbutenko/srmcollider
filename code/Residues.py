
import string

class Residues:

    # http://www.sisweb.com/referenc/source/exactmaa.htm
    # http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
    average_elements = {
        'H' : 1.007825   * 99.99/100 + 2.014102  * 0.015/100,
        'N' : 14.003074  * 99.63/100 + 15.000109 * 0.37/100,
        'O' : 15.994915  * 99.76/100 + 16.999131 * 0.038/100  + 17.999159 * 0.20/100,
        'C' : 12.000000  * 98.90/100 + 13.003355 * 1.10,
        'P' : 30.973763 
    }

    monoisotopic_elements = {
        'H'   : 1.007825032,
        'H2'  : 2.01410178,

        'C'   : 12.000000,
        'C13' : 13.00335484,

        'N'   : 14.003074005,
        'N15' : 15.000108898,

        'O'   : 15.994914620,
        'O17' : 16.999132,
        'O18' : 17.999161,

        'P'   : 30.973762,
        'S'   : 31.972071
    }

    aa_codes = {
        'A' : 'Ala',
        'R' : 'Arg',
        'N' : 'Asn',
        'D' : 'Asp',
        'C' : 'Cys',
        'E' : 'Glu',
        'Q' : 'Gln',
        'G' : 'Gly',
        'H' : 'His',
        'I' : 'Ile',
        'L' : 'Leu',
        'K' : 'Lys',
        'M' : 'Met',
        'F' : 'Phe',
        'P' : 'Pro',
        'S' : 'Ser',
        'T' : 'Thr',
        'W' : 'Trp',
        'Y' : 'Tyr',
        'V' : 'Val',
        'C[160]' : 'Cys+CAM',
        'M[147]' : 'Met+Ox',
    }

    aa_names = {
        'A': 'Alanine',
        'B': 'Aspartic Acid or Asparagine',
        'C': 'Cysteine',
        'c': 'Modified cysteine' ,
        'D': 'Aspartate',
        'E': 'Glutamate',
        'F': 'Phenylalanine',
        'G': 'Glycine',
        'H': 'Histidine',
        'I': 'Isoleucine',
        'K': 'Lysine',
        'k': 'Lys->Cys substitution and carbamidomethylation (903)',
        'L': 'Leucine',
        'M': 'Methionine',
        'm': 'Modified methionine' ,
        'N': 'Asparagine',
        'P': 'Proline',
        'Q': 'Glutamine',
        'R': 'Arginine',
        'S': 'Serine',
        'T': 'Threonine',
        'V': 'Valine',
        'W': 'Tryptophan',
        'X': 'Leucine/Isoleucine',
        'Y': 'Tyrosine',
        'Z': 'Glutamic acid'
        }

    aa_sum_formulas_text = {
        'A' : 'C3H5ON',
        'R' : 'C6H12ON4',
        'N' : 'C4H6O2N2',
        'D' : 'C4H5O3N',
        'C' : 'C3H5ONS',
        'E' : 'C5H7O3N',
        'Q' : 'C5H8O2N2',
        'G' : 'C2H3ON',
        'H' : 'C6H7ON3',
        'I' : 'C6H11ON',
        'L' : 'C6H11ON',
        'K' : 'C6H12ON2',
        'M' : 'C5H9ONS',
        'F' : 'C9H9ON',
        'P' : 'C5H7ON',
        'S' : 'C3H5O2N',
        'T' : 'C4H7O2N',
        'W' : 'C11H10ON2',
        'Y' : 'C9H9O2N',
        'V' : 'C5H9ON'
    }

    #from http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html
    aa_sum_formulas = {
        'A' : { 'C' : 3,  'H' : 5  , 'O' : 1, 'N' : 1  },
        'R' : { 'C' : 6,  'H' : 12 , 'O' : 1, 'N' : 4  },
        'N' : { 'C' : 4,  'H' : 6  , 'O' : 2, 'N' : 2  },
        'D' : { 'C' : 4,  'H' : 5  , 'O' : 3, 'N' : 1  },
        'C' : { 'C' : 3,  'H' : 5  , 'O' : 1, 'N' : 1, 'S' : 1  },
        'E' : { 'C' : 5,  'H' : 7  , 'O' : 3, 'N' : 1  },
        'Q' : { 'C' : 5,  'H' : 8  , 'O' : 2, 'N' : 2  },
        'G' : { 'C' : 2,  'H' : 3  , 'O' : 1, 'N' : 1  },
        'H' : { 'C' : 6,  'H' : 7  , 'O' : 1, 'N' : 3  },
        'I' : { 'C' : 6,  'H' : 11 , 'O' : 1, 'N' : 1  },
        'L' : { 'C' : 6,  'H' : 11 , 'O' : 1, 'N' : 1  },
        'K' : { 'C' : 6,  'H' : 12 , 'O' : 1, 'N' : 2  },
        'M' : { 'C' : 5,  'H' : 9  , 'O' : 1, 'N' : 1, 'S' : 1   },
        'F' : { 'C' : 9,  'H' : 9  , 'O' : 1, 'N' : 1  },
        'P' : { 'C' : 5,  'H' : 7  , 'O' : 1, 'N' : 1  },
        'S' : { 'C' : 3,  'H' : 5  , 'O' : 2, 'N' : 1  },
        'T' : { 'C' : 4,  'H' : 7  , 'O' : 2, 'N' : 1  },
        'W' : { 'C' : 11, 'H' : 10 , 'O' : 1, 'N' : 2  },
        'Y' : { 'C' : 9,  'H' : 9  , 'O' : 2, 'N' : 1  },
        'V' : { 'C' : 5,  'H' : 9  , 'O' : 1, 'N' : 1  },
        'C[160]' : { 'C' : 3+2,  'H' : 5+3  , 'O' : 1+1, 'N' : 1+1, 'S' : 1  }, # + CAM = H(3) C(2) N O 
        'M[147]' : { 'C' : 5,  'H' : 9  , 'O' : 1+1, 'N' : 1, 'S' : 1   },
    }

    mass_H = monoisotopic_elements['H']
    mass_N = monoisotopic_elements['N']
    mass_O = monoisotopic_elements['O']
    mass_C = monoisotopic_elements['C']
    mass_S = monoisotopic_elements['S']
    mass_P = monoisotopic_elements['P']
    mass_NH2 = mass_N + 2*mass_H
    mass_NH3 = mass_N + 3*mass_H
    mass_CO =  mass_C +   mass_O
    mass_H2O = mass_O + 2*mass_H
    mass_OH =  mass_O +   mass_H
    mass_H3PO4 = mass_P + mass_O * 4 + mass_H * 3
    mass_H1PO4 = mass_P + mass_O * 4 + mass_H * 1
    mass_H1PO3 = mass_P + mass_O * 3 + mass_H * 1
    mass_CAM = 2* mass_C + 4*mass_H + mass_O + mass_N #CH2-CONH2

    mass_C13 = monoisotopic_elements['C13']
    mass_N15 = monoisotopic_elements['N15']
    mass_diffC13 = mass_C13 - mass_C
    mass_diffN15 = mass_N15 - mass_N

    average_data = { 
        # Key on abbreviation, give name, molecular weight (in daltons).
        'A': ('Alanine', 71.0788),
        'B': ('Aspartic Acid or Asparagine', 114.5962), 
        'C': ('Cysteine', 103.1448), 
        'c': ('Modified cysteine' , 160.1448), # Add 57
        'D': ('Aspartate', 115.0886),
        'E': ('Glutamate', 129.1155),
        'F': ('Phenylalanine', 147.1766),
        'G': ('Glycine', 57.0519),
        'H': ('Histidine', 137.1411),
        'I': ('Isoleucine', 113.1594),
        'K': ('Lysine', 128.1741),
        'k': ('Lys->Cys substitution and carbamidomethylation (903)', 128.09496 + 32.0219),
        'L': ('Leucine', 113.1594),
        'M': ('Methionine', 131.1986),
        'm': ('Modified methionine' , 147.1986), # add 16
        'N': ('Asparagine', 114.1038),
        'P': ('Proline', 97.1167),
        'Q': ('Glutamine', 128.1307),
        'R': ('Arginine', 156.1875),
        'S': ('Serine', 87.0782),
        'T': ('Threonine', 101.1051),
        'V': ('Valine', 99.1326),
        'W': ('Tryptophan', 186.2132),
        'X': ('Leucine/Isoleucine', 113.1594), # Can't distinguish leucine/isoleucine.
        'Y': ('Tyrosine', 163.176),
        'Z': ('Glutamic acid, or glutamine', 128),
        }
    
    #e.g. from http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html
    monoisotopic_data = { 
        # Key on abbreviation, give name, molecular weight (in daltons).
        'A': ('Alanine', 71.03711),
        'B': ('Aspartic Acid or Asparagine', 114.04293), 
        'C': ('Cysteine', 103.00919), 
        #'c': ('Modified cysteine' , 160.00919), # Add 57
        'c': ('Modified cysteine' , 103.00919 + mass_CAM - mass_H ), # CAM replaces H
        'C[160]': ('Modified cysteine' , 103.00919 + mass_CAM - mass_H ), # CAM replaces H
        'D': ('Aspartate', 115.02694),
        'E': ('Glutamate', 129.04259),
        'F': ('Phenylalanine', 147.06841),
        'G': ('Glycine', 57.02146),
        'H': ('Histidine', 137.05891),
        'I': ('Isoleucine', 113.08406),
        'K': ('Lysine', 128.09496),
        'k': ('Lys->Cys substitution and carbamidomethylation (903)', 128.09496 + 31.935685),
        'L': ('Leucine', 113.08406),
        'M': ('Methionine', 131.04049),
        #'m': ('Modified methionine', 147.04049), # add 16
        'm': ('Modified methionine', 131.04049 + mass_O), # oxygen
        'M[147]': ('Modified methionine', 131.04049 + mass_O), # oxygen
        'N': ('Asparagine', 114.04293),
        'N[115]': ('Asparagine', 114.04293 - mass_N - mass_H + mass_O),
        'P': ('Proline', 97.05276),
        'Q': ('Glutamine', 128.05858),
        'R': ('Arginine', 156.10111),
        'S': ('Serine', 87.03203),
        'T': ('Threonine', 101.04768),
        'V': ('Valine', 99.06841),
        'W': ('Tryptophan', 186.07931),
        'X': ('Leucine/Isoleucine', 113.08406), # Can't distinguish leucine/isoleucine
        'Y': ('Tyrosine', 163.06333),
        'Z': ('Glutamic acid, or glutamine', 128.05858),
        }

#E[111] 351      => pyroGlu
#C[169] 58       => ? 
#R[166] 12599    => SILAC?
#C[143] 319      => Pyro-carbamidomethyl 
#C[152] 2        => ?
#C[553] 4        => ICAT?
#K[136] 15131    => SILAC?
#Q[111] 1730     => pyroGlu
#C[330] 37       => ICAT?
#C[339] 20       => ICAT?
#C[545] 4        => ICAT?
#W[202] 23       => Oxidation?


    hydrophobicity = {
        'F': 5.00,
        'W': 4.88,
        'L': 4.76,
        'X': 4.59,
        'I': 4.41,
        'M': 3.23,
        'V': 3.02,
        'C': 2.50,
        'Y': 2.00,
        'A': 0.16,
        'T': -1.08,
        'E': -1.50,
        'Z': -2.13,
        'D': -2.49,
        'Q': -2.76,
        'R': -2.77,
        'S': -2.85,
        'B': -3.14,
        'G': -3.31,
        'N': -3.79,
        'H': -4.63,
        'P': -4.92,
        'K': -5.00
        }

    basicity = {
        'G': 202.7,
        'C': 206.2,
        'A': 206.4,
        'S': 207.6,
        'D': 208.6,
        'V': 208.7,
        'L': 209.6,
        'X': 210.2,
        'B': 210.7,
        'I': 210.8,
        'T': 211.7,
        'F': 212.1,
        'N': 212.8,
        'Y': 213.1,
        'M': 213.3,
        'Q': 214.2,
        'P': 214.4,
        'Z': 214.9,
        'E': 215.6,
        'W': 216.1,
        'K': 221.8,
        'H': 223.7,
        'R': 237.0
        }

    helicity = {
        'F': 1.26,
        'W': 1.07,
        'L': 1.28,
        'X': 1.29, #avg L,I
        'I': 1.29,
        'M': 1.22,
        'V': 1.27,
        'C': 0.79,
        'Y': 1.11,
        'A': 1.24,
        'T': 1.09,
        'E': 0.85,
        'D': 0.89,
        'Z': 0.91, #avg Q,E
        'B': 0.92, #avg N,D
        'Q': 0.96,
        'R': 0.95,
        'S': 1.00,
        'G': 1.15,
        'N': 0.94,
        'H': 0.97,
        'P': 0.57,
        'K': 0.88,
        }
    
    pI = {
        'G': 6.0,
        'A': 6.0,
        'V': 6.0,
        'L': 6.0,
        'X': 6.0, #L or I
        'I': 6.0,
        'F': 5.5,
        'P': 6.3,
        'S': 5.7,
        'T': 5.6,
        'Y': 5.7,
        'C': 5.0,
        'M': 5.7,
        'N': 5.4,
        'B': 4.1, #avg N and D
        'Q': 5.7,
        'Z': 4.5, #avg Q,E
        'W': 5.9,
        'D': 2.8,
        'E': 3.2,
        'K': 9.7,
        'R': 10.8,
        'H': 7.6
    }

    def __init__(self, type="mono"):
        """Set up the residue data structure."""
        #add the phosphorylations
        self.monoisotopic_data[ 's' ] = ('Phospho-S', 
        self.monoisotopic_data[ 'S' ][1] + self.mass_H1PO3)
        self.monoisotopic_data[ 't' ] = ('Phospho-T', 
        self.monoisotopic_data[ 'T' ][1] + self.mass_H1PO3)
        self.monoisotopic_data[ 'y' ] = ('Phospho-Y',
        self.monoisotopic_data[ 'Y' ][1] + self.mass_H1PO3)

        self.average_data[ 's' ] = ('Phospho-S', 
        self.average_data[ 'S' ][1] + self.mass_H1PO3)
        self.average_data[ 't' ] = ('Phospho-T', 
        self.average_data[ 'T' ][1] + self.mass_H1PO3)
        self.average_data[ 'y' ] = ('Phospho-Y',
        self.average_data[ 'Y' ][1] + self.mass_H1PO3)

        if not type:
            self.residues = self.average_data
        elif type.startswith("mono"):
            self.residues = self.monoisotopic_data
        elif type.startswith("av"):
            self.residues = self.average_data
        else:
            raise ValueError("Type of residue must be one of: mono[isotopic], av[erage] (characters within [] are optional.")
        keys = self.residues.keys()
        self.res_pairs = [ string.join((r, s), '') for r in keys for s in keys ]

    def recalculate_monisotopic_data(self):

        self.monoisotopic_data = {}
        for abbrev, formula in self.aa_sum_formulas.iteritems(): 
            mysum = 0.0
            for key, value in formula.iteritems():
                mysum += self.monoisotopic_elements[ key ] * value
            self.monoisotopic_data[ abbrev ] = ( self.aa_codes[abbrev] , mysum )

        #
        self.monoisotopic_data['c'] = self.monoisotopic_data['C'] + self.mass_CAM - self.mass_H
        self.monoisotopic_data['c'] = ( 'Modified cystein', 
             self.monoisotopic_data['C'][1] + self.mass_CAM - self.mass_H)
        self.monoisotopic_data['k'] = ( 'Lys->Cys substitution and carbamidomethylation (903)',
            self.monoisotopic_data['K'][1] + 31.935685)
        self.monoisotopic_data['m'] = ( 'Modified methionine', 
             self.monoisotopic_data['M'][1] + self.mass_O)
        self.monoisotopic_data[ 's' ] = ('Phospho-S', 
            self.monoisotopic_data[ 'S' ][1] + self.mass_H1PO3)
        self.monoisotopic_data[ 't' ] = ('Phospho-T', 
            self.monoisotopic_data[ 'T' ][1] + self.mass_H1PO3)
        self.monoisotopic_data[ 'y' ] = ('Phospho-Y',
            self.monoisotopic_data[ 'Y' ][1] + self.mass_H1PO3)
        self.residues = self.monoisotopic_data

    def recalculate_monisotopic_data_for_N15(self):

        self.monoisotopic_data = {}
        for abbrev, formula in self.aa_sum_formulas.iteritems(): 
            mysum = 0.0
            for key, value in formula.iteritems():
                #replace N with N15
                if key == 'N': key = 'N15'
                mysum += self.monoisotopic_elements[ key ] * value
            self.monoisotopic_data[ abbrev ] = ( self.aa_codes[abbrev] , mysum )

        #IMPORTANT: CAM is added afterwards and is NOT heavy
        #
        self.monoisotopic_data['C[160]'] = ( 'Modified cystein', 
             self.monoisotopic_data['C'][1] + self.mass_CAM - self.mass_H)
        #
        self.monoisotopic_data['c'] = ( 'Modified cystein', 
             self.monoisotopic_data['C'][1] + self.mass_CAM - self.mass_H)
        self.monoisotopic_data['k'] = ( 'Lys->Cys substitution and carbamidomethylation (903)',
            self.monoisotopic_data['K'][1] + 31.935685)
        self.monoisotopic_data['m'] = ( 'Modified methionine', 
             self.monoisotopic_data['M'][1] + self.mass_O)
        self.monoisotopic_data[ 's' ] = ('Phospho-S', 
            self.monoisotopic_data[ 'S' ][1] + self.mass_H1PO3)
        self.monoisotopic_data[ 't' ] = ('Phospho-T', 
            self.monoisotopic_data[ 'T' ][1] + self.mass_H1PO3)
        self.monoisotopic_data[ 'y' ] = ('Phospho-Y',
            self.monoisotopic_data[ 'Y' ][1] + self.mass_H1PO3)
        self.residues = self.monoisotopic_data

