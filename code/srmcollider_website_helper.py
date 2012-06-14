import os, csv
import collider
import c_getnonuis

import DDB
from Residues import Residues

"""
Expected input format


Format A)
[Sequence]
[Sequence]
[...]

ex:
MTMDRK
MMDRK
MTMRK

Format B)

Name: [Sequence]
Charge: [Charge]
[Q1] [LibraryIntensity] [Annotation] [FragmentCharge]

ex:

Name: MTMDRK
Charge: 2
701 20 b6 2
801 100 y7 4

"""

class NonUnique():

    def __init__(self):
            self.q3            = None
            self.q1            = None
            self.dummy         = None
            self.peptide_key   = None
            self.ion_type      = None
            self.ion_number    = None
            self.sequence      = None
            self.ssrcalc       = None
            self.isotope_nr    = None
            self.charge        = None

class PeptideParser():

    def __init__(self, R):
        self.R  = R

    def parse_stack(self, stack):
        assert len(stack) > 2
        peptide = DDB.Peptide()
        s = stack[0][5:].strip() 
        peptide.charge = int(stack[1][7:].strip())
        sanitized = "".join( [i for i in s if (str.isalnum(i) or i in [ '[', ']']  )] )
        peptide.set_sequence(sanitized)
        peptide.fragments = []
        cnt = 0
        for tr in stack[2:]:
            if len(tr.strip()) == 0: continue
            newtr = tr.split()
            ann = "".join( [i for i in newtr[2] if (str.isalnum(i))] )
            fr = DDB.Fragment(float(newtr[0]), ann, int(newtr[3]) )
            fr.library_intensity = float(newtr[1])
            fr.fragment_count = cnt
            cnt += 1
            peptide.fragments.append(fr)
        # to make sure that it has a charged mass!
        peptide.create_fragmentation_pattern(self.R)
        return peptide

    def parse_transition_list(self, data):
        peptides = []
        stack = []
        for line in data.splitlines(True):
            if len(line.strip() ) == 0: continue 
            if line[:5] == 'Name:':
              if len(stack) > 0:
                peptide = self.parse_stack(stack)
                peptides.append(peptide)
              #
              stack = [line]
            else: stack.append(line)

        if len(stack) > 0:
          peptide = self.parse_stack(stack)
          peptides.append(peptide)

        return peptides

    def sanitize_peptide_input(self, myinput):

        seqs = "'"
        input_sequences = []
        peptides = []
        for inp in myinput.split():
            #only alphanumeric and [ ]
            sanitized = "".join( [i for i in inp if (str.isalnum(i) or i in [ '[', ']']  )] )
            #to look ssrcalc up in the db, we need no modifications
            if len(sanitized) == 0: continue
            seqs += filter(str.isalpha, inp) + "','"
            input_sequences.append(sanitized.upper())

            peptide = DDB.Peptide()
            peptide.set_sequence(sanitized)
            peptide.fragments = []

            peptides.append(peptide)

        seqs = seqs[:-2]
        return seqs, input_sequences

    def calculate_default_fragmenation(self, peptides, par):

        q3_low, q3_high = par.q3_range
        for peptide in peptides:
            peptide.charge = 2
            peptide.create_fragmentation_pattern(self.R)
            fragments = list(peptide.get_fragment_objects(reversed(peptide.y_series),
                'y', 1, self.R, q3_low, q3_high))
            fragments.reverse()
            fragments.extend(list( peptide.get_fragment_objects(peptide.b_series, 
                'b', 1, self.R, q3_low, q3_high)))
            for fcount, f in enumerate(fragments): f.fragment_count = fcount
            peptide.fragments = fragments

    def get_seqs(self, peptides):
        seqs = "'"
        input_sequences = []
        for p in peptides:
            seqs += filter(str.isalpha, p.sequence) + "','"
            input_sequences.append(p.sequence.upper())
        seqs = seqs[:-2]
        return seqs, input_sequences

    def get_sql_sequences(self, peptides):
        seqs = "'"
        for p in peptides:
            seqs += filter(str.isalpha, p.sequence) + "','"
        seqs = seqs[:-2]
        return seqs

    def get_upper_modified_sequences(self, peptides):
        input_sequences = []
        for p in peptides:
            input_sequences.append(p.sequence.upper())
        return input_sequences

class SRMColliderController():

    def __init__(self):
        self.peptides = []
        self.R = Residues('mono')

    def initialize(self, db_used, default_org_prefix, db_tables_map):
        self.db_used = db_used
        self.default_org_prefix = default_org_prefix
        self.db_tables_map = db_tables_map

    def get_sample_peptides_html(self):
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
        return sample_peptides_html 

    def parse_srmcollider_form(self, form, genomes_that_require_N15_data):
        peptides_raw = form.getvalue('peptides')
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

        try:
          table_used = self.db_tables_map[ genome ]
        except KeyError:
            print "Genome not recognized";
            exit()

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
        par.peptide_table = self.db_used + self.default_org_prefix + table_used
        par.transition_table = self.db_used + '.srmTransitions_' + table_used
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

        par.chargeCheck = chargeCheck
        par.genome = genome
        par.uis = uis

        # check whether we need to recalculate the residues for N15
        if genome in genomes_that_require_N15_data: 
            self.R.recalculate_monisotopic_data_for_N15()

        parser = PeptideParser(self.R)

        # try to parse transition list. If it doesnt work, try to parse the
        # input as single peptide sequences.
        try:
          self.peptides = parser.parse_transition_list(peptides_raw)
          seqs, input_sequences = parser.get_seqs(self.peptides)
        except AssertionError:
          self.peptides = []

        if len(self.peptides) == 0:
            # sanitize input: all input is already sanitized except myinput and genome
            seqs, input_sequences = parser.sanitize_peptide_input(peptides_raw)
            parser.calculate_default_fragmenation(self.peptides, par)

        par.seqs = seqs
        par.input_sequences = input_sequences

        return par

    def getNonuniqueObjects(self, nonunique):

        # Syntax
        ## tmplist.append( python::make_tuple(series[k],
        ## q1_used, 0, peptide_key, curr_ion, snumber, precursor.attr("modified_sequence"), ssrcalc, isotope_nr, ch));

        res = {}
        for key in nonunique:
            nonunique_obj = []
            for n in nonunique[key]:
                o = NonUnique()
                o.q3            = n[0]
                o.q1            = n[1]
                o.dummy         = n[2]
                o.peptide_key   = n[3]
                o.ion_type      = n[4]
                o.ion_number    = n[5]
                o.sequence      = n[6]
                o.ssrcalc       = n[7]
                o.isotope_nr    = n[8]
                o.charge        = n[9]
                nonunique_obj.append(o)
            res[key] = nonunique_obj
        return res

def getSRMParameter(q1_w, q3_w, ssr_w, high, low, isotope, ions,
         missed, oxMet, Deamid, db_used, default_org_prefix, table_used):

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

    return par

def get_ssrcalc_values(seqs, input_sequences, default_ssrcalc, cursor, ssrcalc_path):
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

def write_csv_row(fragments, collisions_per_peptide, current_sequence, uis, wuis):
    srm_ids = [f.fragment_count for f in fragments]
    srm_lookup = [ (fragment.fragment_count, fragment) for fragment in fragments]
    srm_lookup = dict(srm_lookup) 
    for order in range(1,uis+1): 
        non_uis = c_getnonuis.get_non_uis(collisions_per_peptide, order)
        if False:
            # here we just output the non-UIS combinations. Usually
            # these are more informative and are preferable to a list
            # of UIS combinations.
            for comb in non_uis:
                tmp = [ srm_lookup[elem] for elem in comb]  
                myrow = []
                for tt in tmp:
                    myrow.extend( [ tt.q3, tt.annotation ])
                wuis.writerow(myrow)
        else:
            # if you want the real deal, go ahead. 
            uis_list = collider.get_uis(srm_ids, non_uis, order)
            #if(len(uis_list) == 0): wuis.writerow([ 'Sorry, no UIS found for order %s' % order ])
            for comb in uis_list:
                tmp = [ srm_lookup[elem] for elem in comb]  
                myrow = [current_sequence, order]
                for tt in tmp:
                    myrow.extend( [ tt.q3, tt.annotation ])
                wuis.writerow(myrow)

