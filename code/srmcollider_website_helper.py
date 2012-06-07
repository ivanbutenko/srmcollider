import os, csv
import collider
import c_getnonuis

def sanitize_peptide_input(myinput):

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

    return seqs, input_sequences

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
