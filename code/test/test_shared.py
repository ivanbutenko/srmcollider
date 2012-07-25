#
# vim:set fdm=marker:

import sys
sys.path.extend( ['.', '../', '../external'])
SQLITE_DATABASE_LOCATION = '/tmp/srmcollider_testdb'

def ignoreImportError_rangetree(f):
    def new(*args):
        try:
            import c_rangetree
            return f(*args)
        except ImportError: pass
    return new

def ignoreImportError_cget(f):
    def new(*args):
        try:
            import c_getnonuis
            return f(*args)
        except ImportError: pass
    return new

def dummy(): 
    pass
def check_cgetnonuis_availability(function):
    try:
        import c_getnonuis
        return function
    except ImportError:
        return dummy

def check_crangetree_availability(function):
    try:
        import c_rangetree
        return function
    except ImportError:
        return dummy

def check_cintegrated_availability(function):
    try:
        import c_integrated
        return function
    except ImportError:
        return dummy



from precursor import Precursor
class ThreePeptideExample():

    precursor = Precursor(modified_sequence='YYLLDYR', q1=503.256187374, q1_charge=2, transition_group = 34, parent_id = 69, isotopically_modified=0, ssrcalc = 25)
    interfering_precursors = [
      Precursor(modified_sequence='GGLIVELGDK', q1=500.787837374, q1_charge=2, transition_group = 665, parent_id = 1331, isotopically_modified=0, ssrcalc = 25),
      Precursor(modified_sequence='NGTDGGLQVAIDAMR', q1=506.58461326, q1_charge=3, transition_group = 618, parent_id = 1238, isotopically_modified=0, ssrcalc = 25)
    ]

    expected_transitions = (
          (842.44121971600021, 0), 
          (679.37788971600014, 1), 
          (566.29382971600012, 2), (453.2097697160001, 3), (440.21854503200007, 4), 
          (553.30260503200009, 5), (668.32954503200006, 6), (831.39287503200012, 7) )

#{{{
#collisions
#q3, q1, srm_id, peptide_key
#transitions
#q3, srm_id
#


#{{{
transitions_def1 = ( (500.0,1), 
                (600,2), 
                (700,3), 
              )
collisions_def1 = (  (500.4,400.0,101,201),
                (600.6,400.0,102,201), 
                (500.6,401.0,103,202),
                (700.6,401.0,104,202), 
                (500.7,400.0,105,203),
                (600.7,400.0,106,203),
                (700.7,400.0,107,203), 
             )

refcollperpep1 = {
    201 : [1,2],
    202 : [1,3],
    203 : [1,2,3],
}

lennonuis1 = [3,3,1,0,0]

refnonuis1 = [
    set([]),
    set([(1,), (2,), (3,)]),
    set([(1, 2), (1, 3), (2, 3)]),
    set([(1, 2, 3)]),
    set([]),
    set([]),
]

#}}}
#{{{
transitions_def2 = ( (500.0,1), 
                (600,2), 
                (700,3), 
                (800,4), 
              )
transitions_def2_unsorted = ( (500.0,1), 
                (600,2), 
                (800,4), 
                (700,3), 
              )
#peptide 201 shares transitions 1-3 and 
#peptide 202 shares transitions 2-4
collisions_def2 = (  (500.4,400.0,101,201),
                (600.6,400.0,102,201), 
                (700.6,401.0,103,201),
                (600.6,401.0,104,202), 
                (700.7,400.0,105,202),
                (800.7,400.0,106,202),
             )

refcollperpep2 = {
    201 : [1,2,3],
    202 : [2,3,4],
}

lennonuis2 = [4,5,2,0,0]

refnonuis2_sorted = [
    set([]),
    set([(1,), (2,), (3,), (4,)]),
    set([(1, 2), (2, 3),  (1, 3), (2, 4), (3, 4)]),
    set([(1, 2, 3), (2, 3, 4), ]),
    set([]),
    set([]),
]

refnonuis2_unsorted = [
    set([]),
    set([(1,), (2,), (3,), (4,)]),
    set([(1, 2), (2, 3),  (1, 3), (2, 4), (4, 3)]),
    set([(1, 2, 3), (2, 4, 3), ]),
    set([]),
    set([]),
]





#}}}
#{{{
transitions_def3 = ( (500.0,1), 
                (600,2), 
                (700,3), 
                (800,4), 
                (900,5), 
                (1000,6), 
              )
#peptide 201 shares transitions 1-3 and 
#peptide 202 shares transitions 2-4
#peptide 203 hares transitions 1-6
collisions_def3 = (  (500.4,400.0,101,201),
                (600.6,400.0,102,201), 
                (700.6,401.0,103,201),
                (600.6,401.0,104,202), 
                (700.7,400.0,105,202),
                (800.7,400.0,106,202),
                (500.6,401.0,104,203), 
                (600.6,401.0,104,203), 
                (700.7,400.0,105,203),
                (800.7,400.0,106,203),
                (900.7,400.0,106,203),
                (1000.7,400.0,106,203),
             )



lennonuis3 = [6, 15, 20, 15, 6]

refnonuis3 = [
    set([]),
    set([(1,), (2,), (3,), (4,), (5,), (6,)]),
    set([
                (5,6),
                (4,6),
                (3,6),
                (2,6),
                (1,6),
                (4,5),
                (3,5),
                (2,5),
                (1,5),
                (3,4),
                (2,4),
                (1,4),
                (2,3),
                (1,3),
                (1,2),
                ]),
    set([       (3, 4, 6), (2, 3, 5), (1, 2, 6),
                (2, 5, 6), (4, 5, 6), (2, 3, 6), (1, 3, 6), (2, 4, 6), (1, 4,
                5), (1, 2, 5), (1, 2, 3), (1, 3, 5), (3, 5, 6), (2, 4, 5), (1,
                3, 4), (3, 4, 5), (1, 4, 6), (1, 5, 6), (1, 2, 4), (2, 3, 4)]),
    set([
                # present     #absent
                (1, 2, 3, 4), # 5,6
                (1, 2, 3, 5), # 4,6
                (1, 2, 4, 5), # 3,6
                (1, 3, 4, 5), # 2,6
                (2, 3, 4, 5), # 1,6
                (1, 2, 3, 6), # 4,5 
                (1, 2, 4, 6), # 3,5
                (1, 3, 4, 6), # 2,5
                (2, 3, 4, 6), # 1,5
                (1, 2, 5, 6), # 3,4
                (1, 3, 5, 6), # 2,4
                (2, 3, 5, 6), # 1,4
                (1, 4, 5, 6), # 2,3
                (2, 4, 5, 6), # 1,3
                (3, 4, 5, 6), # 1,2
                ]),

    set([(1, 2, 3, 4, 5),
                (1, 2, 3, 4, 6),
                (1, 2, 3, 5, 6),
                (1, 2, 4, 5, 6),
                (1, 3, 4, 5, 6),
                (2, 3, 4, 5, 6),
                ])

]
#}}}
#{{{

transitions_def4 = ( (500.0,1), 
                (600,2), 
                (700,3), 
                (800,4), 
                (900,5), 
                (1000,6), 
              )
#peptide 201 shares transitions 1-3 and 
#peptide 202 shares transitions 2-4
#peptide 203 shares transitions 4-6
collisions_def4 = (  (500.4,400.0,101,201),
                (600.6,400.0,102,201), 
                (700.6,401.0,103,201),
                (600.6,401.0,104,202), 
                (700.7,400.0,105,202),
                (800.7,400.0,106,202),
                (800.7,400.0,106,203),
                (900.7,400.0,106,203),
                (1000.7,400.0,106,203),
             )


refnonuis4 = [
    set([]),
    #1
    set([(1,), (2,), (3,), (4,), (5,), (6,)]),
    #2
    set([
        #203
        (5,6),
        (4,6),
        (4,5),
        #202
        (3,4),
        (2,4),
        #201
        (2,3),
        (1,3),
        (1,2),
        ]) ,
    #3
    set([(1, 2, 3), (2, 3, 4), (4, 5, 6) ]),
    set([]),
    set([])
]

lennonuis4 = [6 ,8 ,3 ,0 ,0]

#}}}


transitions_def5 = ( (500.0,1), 
                (600,2), 
                (700,3), 
                (800,4), 
                (900,5), 
                (1000,6), 
              )
#peptide 201 is shares transitions 1-3,6 and 
#peptide 202 is shares transitions 2-4
#peptide 203 is shares transitions 4-6
collisions_def5 = (  (500.4,400.0,101,201),
                (600.6,400.0,102,201), 
                (700.6,401.0,103,201),
                (1000.7,400.0,106,201),
                (600.6,401.0,104,202), 
                (700.7,400.0,105,202),
                (800.7,400.0,106,202),
                (800.7,400.0,106,203),
                (900.7,400.0,106,203),
                (1000.7,400.0,106,203),
             )
#}}}

#{{{ peptide 1
peptide1 = (500, 'PEPTIDE', 1, 1)
#the singly and doubly charged b and y series between 300-1500
transitions_12_between300_1500 = [703.31500971599996,
                                       574.27241971600006,
                                       477.21965971600002,
                                       376.17197971600001,
                                       324.155935032,
                                       425.20361503200002,
                                       538.28767503200004,
                                       653.31461503200001,
                                       352.161417374,
                                       327.16122003200002]
#this is the double charged b and y series (peptide 1)
pep1_yseries = [352.161417374, 287.64012237400004,
                     239.11374237400003, 188.58990237400002,
                     132.04787237400001, 74.53440237400001]
pep1_bseries = [49.534205032000003, 114.055500032,
                     162.58188003200002, 213.10572003200002,
                     269.64775003200003, 327.16122003200002]

#}}}

peptide2 = (400, 'CEPC[160]IDM[147]E',2,2)

from test_shared_large import *

runprecursors_obj1 = []
for p in runprecursors1:
    runprecursors_obj1.append(Precursor(
    modified_sequence = p[1], transition_group = p[2], q1_charge = p[3], isotopically_modified = p[4]))

runprecursors_obj2 = []
for p in runprecursors2:
    runprecursors_obj2.append(Precursor(
    modified_sequence = p[1], transition_group = p[2], q1_charge = p[3], isotopically_modified = p[4]))

runpep_obj1 = Precursor( q1 = runpep1['q1'], modified_sequence = 
  runpep1['mod_sequence'], transition_group = runpep1['transition_group'], 
  isotopically_modified = 0, ssrcalc = runpep1['ssrcalc']) 

runpep_obj2 = Precursor( q1 = runpep2['q1'], modified_sequence = 
  runpep2['mod_sequence'], transition_group = runpep2['transition_group'], 
  isotopically_modified = 0, ssrcalc = runpep2['ssrcalc']) 

from SRM_parameters import SRM_parameters
def get_default_setup_parameters():
        par = SRM_parameters()
        par.q1_window = 1 / 2.0
        par.q3_window = 1 / 2.0
        par.ssrcalc_window = 10 / 2.0
        par.ppm = False
        par.isotopes_up_to = 3
        par.q3_low = 400
        par.q3_high = 1400
        par.max_uis = 5
        par.peptide_tables = ['srmPeptides_test']
        par.mysql_config = '~/.my.cnf'
        par.sqlite_database = SQLITE_DATABASE_LOCATION
        par.use_sqlite = False
        par.quiet = False

        par.bions      =  True
        par.yions      =  True
        par.aions      =  False
        par.aMinusNH3  =  False
        par.bMinusH2O  =  False
        par.bMinusNH3  =  False
        par.bPlusH2O   =  False
        par.yMinusH2O  =  False
        par.yMinusNH3  =  False
        par.cions      =  False
        par.xions      =  False
        par.zions      =  False
        par.MMinusH2O  =  False
        par.MMinusNH3  =  False
        par.q3_range = [par.q3_low, par.q3_high]
        par.set_default_vars()
        par.eval()
        return par

#########################
#########################
#########################
##### Legacy functions
##### used by speed test
#########################
#########################
def _get_unique_pepids(par, cursor, ignore_genomeoccurence=False):
    query = """
    select parent_id, q1, q1_charge, ssrcalc, peptide.id, modified_sequence, transition_group
     from %s
     inner join
     ddb.peptide on peptide.id = %s.peptide_key
     inner join ddb.peptideOrganism on peptide.id = peptideOrganism.peptide_key 
     where genome_occurence = 1
     %s
    """ % (par.peptide_table, par.peptide_table, par.query_add )
    if ignore_genomeoccurence:
        query = """
        select parent_id, q1, q1_charge, ssrcalc, peptide_key, modified_sequence, transition_group
         from %s
         where 4 = 4
         %s
        """ % (par.peptide_table, par.query_add )
    if par.print_query: print query
    #print query
    cursor.execute( query )
    res = cursor.fetchall()
    return [
        {
            'parent_id' :  r[0],
            'q1' :         r[1],
            'q1_charge' :  r[2],
            'ssrcalc' :    r[3],
            'peptide_key' :r[4],
            'mod_sequence':r[5],
            'transition_group':r[6],
        }
        for r in res
    ]

def _get_all_collisions(self, par, pep, cursor,
            values="q3, q1, srm_id, peptide_key", transitions=None,
            bysequence=False):
    q3_low, q3_high = par.get_q3range_collisions()
    vdict = { 'q1' : pep['q1'], 'ssrcalc' : pep['ssrcalc'],
            'peptide_key' : pep['peptide_key'], 'q1_window' : par.q1_window,
            'query_add' : par.query2_add, 'ssr_window' : par.ssrcalc_window,
            'pep' : par.peptide_table, 'values' : values,
            'pepseq' : pep['mod_sequence'],
            'q3_low':q3_low,'q3_high':q3_high,
            'trans' : par.transition_table,
            }
    #we compare the parent ion against 4 different parent ions
    #thus we need to take the PEPTIDE key here
    if bysequence: selectby = "and %(pep)s.modified_sequence != '%(pepseq)s'" % vdict
    else: selectby = "and %(pep)s.peptide_key != %(peptide_key)d" % vdict
    vdict['selectby'] = selectby
    query2 = """
    select %(values)s
    from %(pep)s
    inner join %(trans)s
      on %(pep)s.transition_group = %(trans)s.group_id
    where ssrcalc > %(ssrcalc)s - %(ssr_window)s 
        and ssrcalc < %(ssrcalc)s + %(ssr_window)s
    and q1 > %(q1)s - %(q1_window)s and q1 < %(q1)s + %(q1_window)s
    %(selectby)s
    and q3 > %(q3_low)s and q3 < %(q3_high)s
    %(query_add)s
    """ % vdict
    if transitions is None:
        #if par.select_floor: query2 += "\n GROUP BY FLOOR(q3)"
        if par.print_query: print query2
        cursor.execute( query2 )
    #Here we add N conditions to the SQL query where N = 2*len(transitions)
    #we only select those transitions that also could possible interact
    #We gain about a factor two to three [the smaller the window, the better]
    else:
        txt = 'and ('
        q3_window_used = par.q3_window
        for q3, id in transitions:
            if par.ppm: q3_window_used = par.q3_window * 10**(-6) * q3
            txt += '(q3 > %s and q3 < %s) or\n' % (q3 - q3_window_used, q3 + q3_window_used)
        txt = txt[:-3] + ')'
        #if par.select_floor: txt += "\n GROUP BY FLOOR(q3)"
        if par.print_query: print query2 + txt
        # print query2
        # print txt
        cursor.execute( query2 + txt )
    return cursor.fetchall()

def getnonuis(transitions, collisions, q3_window, ppm):
        collisions_per_peptide = {}
        q3_window_used = q3_window
        for t in transitions:
            if ppm: q3_window_used = q3_window * 10**(-6) * t[0]
            this_min = q3_window_used
            for c in collisions:
                if abs( t[0] - c[0] ) <= q3_window_used:
                    #gets all collisions
                    if collisions_per_peptide.has_key(c[3]):
                        if not t[1] in collisions_per_peptide[c[3]]:
                            collisions_per_peptide[c[3]].append( t[1] )
                    else: collisions_per_peptide[c[3]] = [ t[1] ] 
        return collisions_per_peptide

def get_non_UIS_from_transitions(transitions, collisions, par, MAX_UIS, 
                                forceset=False):
    """ Get all combinations that are not UIS 
    
    Note that the new version returns a dictionary. To convert it to a set, one 
    needs to force the function to return a set.
    """
    try: 
        #using C++ functions for this == faster
        import c_getnonuis
        non_uis_list = [{} for i in range(MAX_UIS+1)]
        collisions_per_peptide = getnonuis(transitions, collisions, par.q3_window, par.ppm)
        for order in range(1,MAX_UIS+1):
            non_uis_list[order] = c_getnonuis.get_non_uis(
                collisions_per_peptide, order)

        if forceset: return [set(k.keys()) for k in non_uis_list]
        return non_uis_list

    except ImportError:
        #old way of doing it
        return get_non_UIS_from_transitions_old(transitions, collisions, par, MAX_UIS)

import sys
sys.path.extend(['..'])
import uis_functions
def get_non_UIS_from_transitions_old(transitions, collisions, par, MAX_UIS, unsorted=False):
    """ Get all combinations that are not UIS """
    #collisions
    #q3, q1, srm_id, peptide_key
    #transitions
    #q3, srm_id
    collisions_per_peptide = {}
    non_uis_list = [set() for i in range(MAX_UIS+1)]
    q3_window_used = par.q3_window
    for t in transitions:
        if par.ppm: q3_window_used = par.q3_window * 10**(-6) * t[0]
        this_min = q3_window_used
        for c in collisions:
            if abs( t[0] - c[0] ) <= q3_window_used:
                #gets all collisions
                if collisions_per_peptide.has_key(c[3]):
                    if not t[1] in collisions_per_peptide[c[3]]:
                        collisions_per_peptide[c[3]].append( t[1] )
                else: collisions_per_peptide[c[3]] = [ t[1] ] 
    #here we calculate the UIS for this peptide with the given RT-range
    for pepc in collisions_per_peptide.values():
        for i in range(1,MAX_UIS+1):
            if unsorted: get_non_uis_unsorted(pepc, non_uis_list[i], i)
            else: uis_functions.get_non_uis(pepc, non_uis_list[i], i)
    return non_uis_list

