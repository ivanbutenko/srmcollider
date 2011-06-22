#!/usr/bin/python
# -*- coding: utf-8  -*-

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

import sys, os, time
sys.path.append( '/home/hroest/projects/msa/code/tppGhost' )
sys.path.append( '/home/hroest/projects/hlib/' )
sys.path.append( '/IMSB/users/hroest/projects/hlib/' )
sys.path.append( '/IMSB/users/hroest/projects/msa/code/tppGhost' )
sys.path.append( '/IMSB/users/hroest/projects/tppGhost' )
import numpy
import progress
import gnuplot
import DDB

class SRM_parameters(object):

    def __init__(self): 
        self.do_1vs = True #check only one charge state?
        self.do_vs1 = False #background only one charge state?
        self.dontdo2p2f = True #do not look at 2+ parent / 2+ fragment ions
        self.considerIsotopes = False #do not consider the C13 isotopes
        self.isotopes_up_to = 0
        self.ppm = True #measure q3 in ppm
        self.transition_table = 'hroest.srmTransitions_yeast'
        self.peptide_table = 'hroest.srmPeptides_yeast'
        self.q3_range = [400, 1200]
        self.ssrcalc_window = 2.0 / 2
        self.q1_window = 25.0 / 2.0
        self.q3_window = 10.0 / 2.0
        self.max_uis = 5 #maximal order of UIS to calculate (no UIS => set to 0)
        self.do_1_only = "and q1_charge = 2 and q3_charge = 1"
        self.print_query = False
        #self.select_floor = False
        self.quiet = False
        self.aions      =  False
        self.aMinusNH3  =  False
        self.bMinusH2O  =  False
        self.bMinusNH3  =  False
        self.bPlusH2O   =  False
        self.yMinusH2O  =  False
        self.yMinusNH3  =  False
        self.cions      =  False
        self.xions      =  False
        self.zions      =  False

    def parse_cmdl_args(self, parser, default_mysql = "~/.my.cnf"):
        from optparse import OptionGroup
        group = OptionGroup(parser, "General Options",
                            "These are the general options for the SRM  Collider")
        #group.add_option( dest="q1_window", default=1, help="Q1 window (e.g. use 1 for +- 0.5 Da). Defaults to 1", metavar="Q1WIN")
        group.add_option("--q1_window", dest="q1_window", default=1, type="float",
                          help="Q1 window (e.g. use 1 for +- 0.5 Da). " + 
                          "Defaults to 1", metavar="Q1WIN")
        group.add_option("--q3_window", dest="q3_window", default=1, type="float",
                          help="Q3 window (e.g. use 1 for +- 0.5 Da). " + 
                          "Defaults to 1", metavar="Q3WIN")
        group.add_option("--ssrcalc_window", dest="ssrcalc_window", default=9999, type="float",
                          help="RT (retention time) window (e.g. use 1 for +- 0.5 units)." + 
                          " Defaults to 9999 (infinite.)", metavar="RTWIN")
        group.add_option("--ppm", dest="ppm", default=False, 
                          help="Interpret Q3 window as PPM (default False)")
        group.add_option("-i", "--considerIsotopes", dest="considerIsotopes", default=True,
                          help="Consider isotopes of the precursors, up to 3amu " +
                          "(default True)")
        group.add_option("--q3_low", dest="q3_low", default=400, type="float",
                          help="Start of transition range to analyse (default 400)", metavar="Q3LOW")
        group.add_option("--q3_high", dest="q3_high", default=1400, type="float",
                          help="End of transition range to analyse (default 1400)", metavar="Q3HIGH")
        group.add_option("--max_uis", dest="max_uis", default=0, type='int',
                          help="maximal order of UIS to calculate " +
                          "(defaults to 0 == no UIS)" )
        group.add_option("--peptide_table", dest="peptide_table", default='hroest.srmPeptides_yeast',
                          help="MySQL table containing the peptides" )
        group.add_option("--transition_table", dest="transition_table", default='hroest.srmTransitions_yeast',
                          help="MySQL table containing the transitions" )
        group.add_option("--mysql_config", dest="mysql_config", default=default_mysql,
                          help="Location of mysql config file, defaults to %s" % default_mysql )
        group.add_option("-q", "--quiet", dest="quiet", default=False,
                          help="don't print status messages to stdout")
        parser.add_option_group(group)

    def parse_options(self, options):
        self.q3_range = [options.q3_low, options.q3_high]
        self.q1_window /= 2.0
        self.q3_window /= 2.0
        self.ssrcalc_window /= 2.0
        if self.ppm == 'True': par.ppm = True
        elif self.ppm == 'False': par.ppm = False

    def eval(self):
        #query will get all parent ions to consider
        #query1 will get all interesting transitions
        #query2 will get all colliding transitions
        self.query_add = "and isotope_nr = 0 "
        self.query1_add = "and isotope_nr = 0"
        self.query2_add = ""
        self.ppm_string = "Th"
        if self.do_1vs :
            self.query_add += "and q1_charge = 2 "
            self.query1_add = self.do_1_only
        if self.do_vs1 : 
            self.query2_add += self.do_1_only
        elif self.dontdo2p2f:
            self.query2_add += "and not (q1_charge = 2 and q3_charge = 2) "
        if not self.considerIsotopes: self.query2_add += "and isotope_nr = 0"
        if self.ppm: self.ppm_string = "PPM"
        self.experiment_type = """Experiment Type:
        check all four charge states [%s] vs all four charge states [%s] with
        thresholds of SSRCalc %s, Q1 %s (Th), Q3 %s (%s) and a range of %s to %s
        Da for the q3 transitions.  Ignore 2+ parent / 2+ fragment ions %s, 
        selecting from %s and %s.
        Consider Isotopes: %s""" % ( 
            not self.do_1vs, not self.do_vs1,
          self.ssrcalc_window*2,  self.q1_window*2, self.q3_window*2, 
          self.ppm_string, self.q3_range[0], self.q3_range[1], self.dontdo2p2f, 
          self.peptide_table, self.transition_table, self.considerIsotopes)
        self.thresh = \
        "thresholds of SSRCalc %s, Q1 %s (Th), Q3 %s (%s) and a range of %s to %s" % (
          self.ssrcalc_window*2,  self.q1_window*2, self.q3_window*2, 
          self.ppm_string, self.q3_range[0], self.q3_range[1]
        )

    def get_copy(self):
        import copy
        return copy.deepcopy(self)

    def get_q3range_transitions(self):
        return self.q3_range[0], self.q3_range[1]

    def get_q3range_collisions(self):
        if not self.ppm:
            return self.q3_range[0]-self.q3_window,\
                   self.q3_range[1]+self.q3_window
        else:
            return self.q3_range[0]-self.q3_window*10**(-6)*self.q3_range[0],\
                   self.q3_range[1]+self.q3_window*10**(-6)*self.q3_range[1]

    def get_q3_window_transitions(self, q3):
        if self.ppm:
            return [q3 - self.q3_window * 10**(-6) *q3, q3 + self.q3_window * 10**(-6) * q3]
        else: return [q3 - self.q3_window, q3 + self.q3_window]

    def get_common_filename(self):
        background = self.do_vs1 
        if not self.do_vs1 and self.dontdo2p2f: background = 'vs3'
        common_filename = '%s_%s_%s_%d_%d' % (self.peptide_tbl_identifier, 
            self.do_1vs, background,
            self.ssrcalc_window*20,  self.q1_window*20)
        if self.ppm:
            common_filename += '_%dppm' % (self.q3_window*20)
        else:
            common_filename += '_%d' % (self.q3_window*20)
        if self.q3_range[1] < 0:
            common_filename += "_range%stomin%s" % (self.q3_range[0], abs( self.q3_range[1]) )
        else:
            common_filename += "_range%sto%s" % (self.q3_range[0], abs( self.q3_range[1]) )
        return common_filename

    @property
    def transition_db(self): return self.transition_table.split('.')[0]

    @property
    def transition_tbl(self): return self.transition_table.split('.')[1]

    @property
    def peptide_db(self): return self.peptide_table.split('.')[0]

    @property
    def peptide_tbl(self): 
        try:
            return self.peptide_table.split('.')[1]
        except IndexError:
            #not in database.table format
            return self.peptide_table

    @property
    def peptide_tbl_identifier(self): return self.peptide_tbl[12:] #cut off 'srmPeptides'

def testcase():
    par = SRM_parameters()
    par.q1_window = 0.7 / 2
    par.q3_window = 1.0 / 2
    par.ppm = False
    par.transition_table = 'hroest.srmTransitions_test'
    par.peptide_table = 'hroest.srmPeptides_test'
    par.dontdo2p2f = False #do not look at 2+ parent / 2+ fragment ions
    #default 
    par.considerIsotopes = False #do not consider the C13 isotopes
    par.do_1vs = True #check only one charge state?
    par.do_vs1 = False #background only one charge state?
    par.q3_range = [400, 1200]
    par.ssrcalc_window = 2.0 / 2
    par.do_1_only = "and q1_charge = 2 and q3_charge = 1"
    #
    par.eval()
    return par

def get_cum_dist(original):
    cumulative = []
    cum  = 0
    for i, val in enumerate(original):
        cum += val
        cumulative.append(  cum )
    return cumulative

class SRMcollider(object):

    def find_clashes_small(self, db, par, use_per_transition=False,
                           UIS_only=False, exp_key=None, calc_collisions=False, 
                           pepids=None):
        #make sure we only get unique peptides
        cursor = db.cursor()
        if pepids is not None: self.pepids = pepids
        else: self.pepids = self._get_unique_pepids(par, cursor)
        MAX_UIS = par.max_uis
        common_filename = par.get_common_filename()
        self.allpeps = {}
        self.non_unique_count = 0
        self.total_count = 0
        self.start = time.time()
        progressm = progress.ProgressMeter(total=len(self.pepids), unit='peptides')
        self.coll_time = 0
        self.innermysql = 0
        for ii, pep in enumerate(self.pepids):
            p_id = pep['parent_id']
            q1 = pep['q1']
            non_unique = {}
            non_unique_q1 = {}
            non_unique_ppm = {}
            non_unique_clash = {}
            transitions = self._get_all_transitions(par, pep, cursor)
            nr_transitions = len( transitions )
            mystart = time.time()
            if nr_transitions == 0: continue #no transitions in this window
            if use_per_transition: collisions = self._get_all_collisions_per_transition(par, pep, transitions, cursor)
            elif calc_collisions: 
                collisions = list(self._get_all_collisions_calculate(par, pep, cursor))
            else: collisions = self._get_all_collisions(par, pep, cursor, transitions = transitions)
            #
            self.coll_time += time.time() - mystart
            if UIS_only:
                non_uis_list = get_non_UIS_from_transitions(transitions, 
                                            collisions, par, MAX_UIS)
                for i in range(1,min(MAX_UIS+1, nr_transitions+1)): 
                    restable = 'hroest.result_srmuis'
                    cursor.execute( "insert into %s (non_useable_UIS, total_UIS, \
                                   parent_key, uisorder, exp_key) values " % restable +\
                                   "(%s,%s,%s,%s,%s)" % (len(non_uis_list[i]) , 
                                       choose(nr_transitions, i), p_id , i, exp_key) )
                self.end = time.time()
                if not par.quiet: progressm.update(1)
                continue
            #here we loop through all possible combinations of transitions and
            #potential collisions and check whether we really have a collision
            q3_window_used = par.q3_window
            for t in transitions:
                if par.ppm: q3_window_used = par.q3_window * 10**(-6) * t[0]
                this_min = q3_window_used
                for c in collisions:
                    if abs( t[0] - c[0] ) <= this_min:
                        non_unique[ t[1] ] = t[0] - c[0]
            self.allpeps[ p_id ] = 1.0 - len( non_unique ) * 1.0  / nr_transitions
            self.non_unique_count += len( non_unique )
            self.total_count += len( transitions)
            self.end = time.time()
            if not par.quiet: progressm.update(1)
        self.total_time = self.end - self.start

    def _get_all_precursors(self, par, pep, cursor, 
         values="q1, modified_sequence, peptide_key, q1_charge", 
         bysequence=False):
        #we compare the parent ion against 4 different parent ions
        #thus we need to take the PEPTIDE key here
        vdict = { 'q1' : pep['q1'], 'ssrcalc' : pep['ssrcalc'],
                'peptide_key' : pep['peptide_key'], 'q1_window' : par.q1_window,
                'query_add' : par.query2_add, 'ssr_window' : par.ssrcalc_window,
                'pep' : par.peptide_table, 'values' : values, 
                'pepseq' : pep['mod_sequence']}
        if bysequence: selectby = "and %(pep)s.modified_sequence != '%(pepseq)s'" % vdict
        else: selectby = "and %(pep)s.peptide_key != %(peptide_key)d" % vdict
        vdict['selectby'] = selectby
        vdict['query_add'] += vdict['query_add'] + ' and isotope_nr <= %s ' % par.isotopes_up_to
        query2 = """
        select %(values)s
        from %(pep)s
        where ssrcalc > %(ssrcalc)s - %(ssr_window)s 
            and ssrcalc < %(ssrcalc)s + %(ssr_window)s
        and q1 > %(q1)s - %(q1_window)s and q1 < %(q1)s + %(q1_window)s
        %(selectby)s
        %(query_add)s
        """ % vdict
        if par.print_query: print query2
        cursor.execute( query2 )
        return cursor.fetchall()


    def _get_all_collisions_calculate(self, par, pep, cursor, 
      values="q1, modified_sequence, peptide_key, q1_charge, transition_group"):
        """
        Here we calculate all the fragments on the fly instead of getting them
        from the database.

        We do get the precursors which could interfere from the database and 
        then create a DDB.Peptide object for each precursor, and get the
        fragmentation patterns from those.
        """
        import Residues
        import DDB 
        R = Residues.Residues('mono')
        q3_low, q3_high = par.get_q3range_collisions()
        for c in self._get_all_precursors(par, pep, cursor, values=values):
            q1 = c[0]
            peptide_key = c[2]
            peptide = DDB.Peptide()
            peptide.set_sequence(c[1])
            peptide.charge = c[3]
            peptide.create_fragmentation_pattern(R)
            b_series = peptide.b_series
            y_series = peptide.y_series
            for ch in [1,2]:
                for pred in y_series:
                    q3 = ( pred + (ch -1)*R.mass_H)/ch
                    if q3 < q3_low or q3 > q3_high: continue
                    yield (q3, q1, 0, peptide_key)
                for pred in b_series:
                    q3 = ( pred + (ch -1)*R.mass_H)/ch
                    if q3 < q3_low or q3 > q3_high: continue
                    yield (q3, q1, 0, peptide_key)

    def _get_all_collisions_calculate_sub(self, precursors, par, R, q3_low, q3_high):
        for c in precursors:
            q1 = c[0]
            peptide_key = c[2]
            peptide = DDB.Peptide()
            peptide.set_sequence(c[1])
            peptide.charge = c[3]
            peptide.create_fragmentation_pattern( R, 
                aions      =  par.aions    ,
                aMinusNH3  =  par.aMinusNH3,
                bMinusH2O  =  par.bMinusH2O,
                bMinusNH3  =  par.bMinusNH3,
                bPlusH2O   =  par.bPlusH2O ,
                yMinusH2O  =  par.yMinusH2O,
                yMinusNH3  =  par.yMinusNH3,
                cions      =  par.cions    ,
                xions      =  par.xions    ,
                zions      =  par.zions )
            b_series = peptide.b_series
            y_series = peptide.y_series
            #TODO fix it
            peptide.allseries = b_series
            peptide.allseries.extend( y_series )
            for ch in [1,2]:
                for pred in peptide.allseries:
                    q3 = ( pred + (ch -1)*R.mass_H)/ch
                    if q3 < q3_low or q3 > q3_high: continue
                    yield (q3, q1, 0, peptide_key)

    def find_clashes(self, db, par, toptrans=False, pepids=None, 
                     use_per_transition=False, exp_key = None):
        #make sure we only get unique peptides
        MAX_UIS = par.max_uis
        cursor = db.cursor()
        if pepids is not None: self.pepids = pepids
        elif not toptrans: self.pepids = self._get_unique_pepids(par, cursor)
        else: self.pepids = self._get_unique_pepids_toptransitions(par, cursor)

        self.allpeps = {}
        self.allcollisions = []
        self.q3min_distr = []
        self.q3all_distr = []
        self.q1all_distr = []
        self.q1min_distr = []
        self.q3min_distr_ppm = []
        self.non_unique_count = 0
        self.count_pair_collisions = {}
        self.found3good = []
        self.total_count = 0
        self.colliding_peptides = {}
        self.UIS_redprob = [0 for i in range(MAX_UIS+1)]
        self.UIS_singleprob = [ {} for i in range(MAX_UIS+1)]
        self.count_analysed = 0
        common_filename = par.get_common_filename()
        self.start = time.time()
        progressm = progress.ProgressMeter(total=len(self.pepids), unit='peptides')

        self.q3min_hist = numpy.zeros(1000)
        self.q3all_hist = numpy.zeros(1000)
        self.q1all_hist = numpy.zeros(1000)
        self.q1min_hist = numpy.zeros(1000)

        self.coll_time = 0
        self.trans_time = 0

        for ii, pep in enumerate(self.pepids):

            one_runstart = time.time()
            p_id = pep['parent_id']
            q1 = pep['q1']
            non_unique = {}
            non_unique_q1 = {}
            non_unique_ppm = {}
            non_unique_clash = {}
            all_clashes = {}
            non_uis_list = [set() for i in range(MAX_UIS+1)]
            collisions_per_peptide = {}

            #END per loop variables
            mystart = time.time()
            if not toptrans: transitions = self._get_all_transitions(par, pep, cursor)
            else: transitions = self._get_all_transitions_toptransitions(par, pep, cursor)
            nr_transitions = len( transitions )
            self.trans_time += time.time() - mystart
            #
            mystart = time.time()
            if nr_transitions == 0: continue #no transitions in this window
            if use_per_transition: collisions = self._get_all_collisions_per_transition(par, pep, transitions, cursor)
            else: 
                #collisions = self._get_all_collisions(par, pep, cursor)
                collisions = self._get_all_collisions(par, pep, cursor, transitions = transitions)
            self.coll_time += time.time() - mystart
            #
            #collisions
            #q3, q1, srm_id, peptide_key
            #transitions
            #q3, srm_id
            #here we loop through all possible combinations of transitions and
            #potential collisions and check whether we really have a collision
            q3_window_used = par.q3_window
            for t in transitions:
                if par.ppm: q3_window_used = par.q3_window * 10**(-6) * t[0]
                this_min = q3_window_used
                q1min = par.q1_window
                for c in collisions:
                    if abs( t[0] - c[0] ) <= q3_window_used:
                        #gets all collisions
                        coll_key = "%s:%s" % (p_id, c[3])
                        if collisions_per_peptide.has_key(c[3]):
                            if not t[1] in collisions_per_peptide[c[3]]:
                                collisions_per_peptide[c[3]].append( t[1] )
                        else: collisions_per_peptide[c[3]] = [ t[1] ] 
                        #print coll_key
                        try: self.count_pair_collisions[ coll_key ] += 1
                        except KeyError: self.count_pair_collisions[ coll_key ] = 1
                        self.q1all_distr.append( q1 - c[1])
                        self.q3all_distr.append(t[0] - c[0] ) 
                        if abs( q1 - c[1] ) <= q1min:
                            non_unique_q1[ t[1] ] = q1 - c[1]
                            q1min = abs( q1 - c[1])
                    if abs( t[0] - c[0] ) <= this_min:
                        #gets the nearest collision
                        #print t[0] - c[0], c
                        this_min = abs( t[0] - c[0] )
                        non_unique[ t[1] ] = t[0] - c[0]
                        non_unique_ppm[ t[1] ] = (t[0] - c[0] ) * 10**6 / t[0]
                        all_clashes[ t[1] ] = c[2]

            #here we calculate the UIS for this peptide with the given RT-range
            for pepc in collisions_per_peptide.values():
                for i in range(1,MAX_UIS+1):
                    get_non_uis(pepc, non_uis_list[i], i)
            #for i in range(1,MAX_UIS+1): print len(non_uis_list[i]), choose(nr_transitions, i)
            for i in range(1,min(MAX_UIS+1, nr_transitions+1)): 
                #also in the UIS paper they compute the average of the probabilities
                self.UIS_redprob[i] += len(non_uis_list[i]) / choose(nr_transitions, i)
            #we could actually calculate the UIS combinations 
            #from the non-UIS combinations
            if False:
                srm_ids = [t[1] for t in transitions]
                get_uis(srm_ids, non_uis_list[2], 2)
                len(get_uis(srm_ids, non_uis_list[2], 2))
            #number of colliding peptides
            self.colliding_peptides[ p_id ] = len( collisions_per_peptide )
            self.allpeps[ p_id ] = 1.0 - len( non_unique ) * 1.0  / nr_transitions
            if len(transitions) - len(non_unique) > 3: self.found3good.append( p_id )
            if len(non_unique) > 0 and exp_key is not None:
                self.allcollisions.extend( all_clashes.items() )
                for v in non_unique.values(): self.q3min_distr.append( v )
                for v in non_unique_q1.values(): self.q1min_distr.append( v )
                for v in non_unique_ppm.values(): self.q3min_distr_ppm.append( v )
            self.non_unique_count += len( non_unique )
            self.total_count += len( transitions)
            self.count_analysed += 1
            self.end = time.time()
            if not par.quiet: progressm.update(1)
            #print "This run took %ss" % (time.time() - one_runstart )
            if ii > 0 and ii % 100 == 0 and exp_key is not None: 

                q1r = [-par.q1_window, par.q1_window] 
                q3r = [-par.q3_window, par.q3_window] 
                #print "pep ", ii
                #print len(self.q3all_distr), sum( self.q3all_hist)
                self.q1min_hist += numpy.histogram( self.q1min_distr, 1000, q1r)[0]
                self.q3min_hist += numpy.histogram( self.q3min_distr, 1000, q3r)[0]
                self.q1all_hist += numpy.histogram( self.q1all_distr, 1000, q1r)[0]
                self.q3all_hist += numpy.histogram( self.q3all_distr, 1000, q3r)[0]
                #print sum( self.q3all_hist)

                self.allpeps = {}
                self.allcollisions = []
                self.q3min_distr = []
                self.q3all_distr = []
                self.q1all_distr = []
                self.q1min_distr = []
                self.q3min_distr_ppm = []
                self.non_unique_count = 0
                self.count_pair_collisions = {}
                self.found3good = []
                self.colliding_peptides = {}
                self.UIS_redprob = [0 for i in range(MAX_UIS+1)]
                self.UIS_singleprob = [ {} for i in range(MAX_UIS+1)]

        self.total_time = self.end - self.start
        #the last one
        q1r = [-par.q1_window, par.q1_window] 
        q3r = [-par.q3_window, par.q3_window] 
        self.q1min_hist += numpy.histogram( self.q1min_distr, 1000, q1r)[0]
        self.q3min_hist += numpy.histogram( self.q3min_distr, 1000, q3r)[0]
        self.q1all_hist += numpy.histogram( self.q1all_distr, 1000, q1r)[0]
        self.q3all_hist += numpy.histogram( self.q3all_distr, 1000, q3r)[0]

    def find_clashes_toptrans_paola(self, db, par,
        use_per_transition=False, bgpar=None):
        """
        We find interferences the following way:
            Single for the top transition
            Pairs for the top 2 transitions
            Triplets for the top3 transitions
        We stop at par.max_uis
        """
        #make sure we only get unique peptides
        assert par.query1_add.find( ' order by Intensity DESC') != -1 
        assert not par.do_1vs 
        assert par.max_uis > 0
        cursor = db.cursor()
        background = False
        if bgpar == None: bgpar = par
        else: background = True
        self.pepids = self._get_unique_pepids_toptransitions(par, cursor)
        MAX_UIS = par.max_uis
        common_filename = par.get_common_filename()
        #self.peptide_specificity = []
        self.min_transitions = []
        start = time.time()
        progressm = progress.ProgressMeter(total=len(self.pepids), unit='peptides')
        for i, pep in enumerate(self.pepids):
            p_id = pep['parent_id']
            transitions = self._get_all_transitions_toptransitions(par, pep, cursor)
            nr_transitions = len( transitions )
            if nr_transitions == 0: continue #no transitions in this window
            if use_per_transition: collisions = self._get_all_collisions_per_transition(
                bgpar, pep, transitions, cursor)
            else: 
                if background:
                    #q =  "select peptide_key from %s where modified_sequence = '%s' " % (
                    #    bgpar.peptide_table, pep['mod_sequence'].split('/')[0]) 
                    #cursor.execute(q)
                    #oldpkey = pep['peptide_key']
                    #newpkey = cursor.fetchone()
                    #if newpkey is not None: pep['peptide_key'] = newpkey[0]
                    #collisions = self._get_all_collisions(bgpar, pep, cursor)
                    #pep['peptide_key'] = oldpkey
                    collisions = self._get_all_collisions(bgpar, pep, cursor, bysequence=True, transitions=transitions)
                else: 
                    collisions = self._get_all_collisions(bgpar, pep, cursor, transitions=transitions)
            #
            min_needed = self._getMinNeededTransitions(par, transitions, collisions)
            #
            self.min_transitions.append( [p_id, min_needed] )
            end = time.time()
            if not par.quiet: progressm.update(1)
        self.total_time = end - start

    def _getMinNeededTransitions(self, par, transitions, collisions):
        nr_transitions = len( transitions )
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
        return self._sub_getMinNeededTransitions(par, transitions, collisions_per_peptide)

    def _sub_getMinNeededTransitions(self, par, transitions, collisions_per_peptide):
        #take the top j transitions and see whether they, as a tuple, are
        #shared
        min_needed = -1
        for j in range(par.max_uis,0,-1): 
            #take top j transitions
            mytransitions = tuple(sorted([t[1] for t in transitions[:j]]))
            unuseable = False
            for k,v in collisions_per_peptide.iteritems():
                if tuple(sorted(v)[:j]) == mytransitions: unuseable=True
            if not unuseable: min_needed = j
        return min_needed

    def find_clashes_toptrans_3strike(self, db, par, pepids=None, 
                     use_per_transition=False):
        """
        The idea is to find UIS combinations that are globally UIS, locally
        clean and also there are no cases in which all the UIS coelute when
        they are in different peptides:
            * global UIS = whole RT
            * locally clean = no intereferences around the peptide
            * no coelution: find all peptides globally that share transitions
        This is a 3 strikes rule to find good UIS combinations.
        """
        #make sure we only get unique peptides
        VERY_LARGE_SSR_WINDOW = 9999999
        cursor = db.cursor()
        if pepids is not None: self.pepids = pepids
        elif not toptrans: self.pepids = self._get_unique_pepids(par, cursor)
        else: self.pepids = self._get_unique_pepids_toptransitions(par, cursor)
        self.count_analysed = 0
        common_filename = par.get_common_filename()
        start = time.time()
        progressm = progress.ProgressMeter(total=len(self.pepids), unit='peptides')
        myssrcalc = par.ssrcalc_window
        f = open('paola_uis.csv', 'w')
        ff = open('paola.log', 'w')
        for i, pep in enumerate(self.pepids):
            transitions_extend = self._get_all_transitions_toptransitions(par, pep,
                cursor, values = 'q3, m.id, m.sequence, m.Ion_description')
            transitions = [ (t[0], t[1]) for t in transitions_extend]
            nr_transitions = len( transitions )
            transitions_dic = [ [t[1], t[0]] for t in transitions]
            transitions_dic = dict(transitions_dic)
            if nr_transitions == 0: continue #no transitions in this window

            ###############################################################
            #strike 1: it has to be global UIS
            par.ssrcalc_window = VERY_LARGE_SSR_WINDOW
            collisions = self._get_all_collisions_per_transition(par, pep, transitions, cursor)
            pairs_1strike = get_UIS_from_transitions(transitions, collisions, par, 2)[2]
            par.ssrcalc_window = myssrcalc

            ###############################################################
            #strike 2: it has to be locally clean
            collisions = self._get_all_collisions_per_transition(par, pep, transitions, cursor)
            local_interferences = get_non_UIS_from_transitions(transitions, collisions, par, 1)
            pairs_2strike = []
            for pair in pairs_1strike:
                contaminated = False
                for dirty_t in local_interferences[1]:
                    if dirty_t[0] in pair: contaminated = True
                if not contaminated: pairs_2strike.append(pair)

            ###############################################################
            #strike 3: the two transitions shall not coelute elsewhere
            pairs_3strike = []
            for pair in pairs_2strike: 
                mass1 = transitions_dic[pair[0]]
                mass2 = transitions_dic[pair[1]]
                par.ssrcalc_window = VERY_LARGE_SSR_WINDOW
                globalc1 = self._get_collisions_per_transition(par, 
                               pep, mass1, cursor, values='ssrcalc')
                globalc2 = self._get_collisions_per_transition(par,
                               pep, mass2, cursor, values='ssrcalc')
                par.ssrcalc_window = myssrcalc
                contaminated = False
                for c1 in globalc1:
                    for c2 in globalc2:
                        if abs(c1[0]-c2[0]) < par.ssrcalc_window: contaminated = True
                if not contaminated: pairs_3strike.append(pair)

            ###############################################################
            #output
            ff.write( "%s,%s,%s,%s\n" % (transitions_extend[0][2], len(pairs_1strike), len(pairs_2strike), len(pairs_3strike)))
            for pair in pairs_3strike:
                p1 = [t for t in transitions_extend if pair[0] == t[1] ][0]
                p2 = [t for t in transitions_extend if pair[1] == t[1] ][0]
                f.write( "%s,%s,%s,%s,%s\n" % (p1[2], p1[0], p1[3], p2[0], p2[3]) )

            self.count_analysed += 1
            end = time.time()
            if not par.quiet: progressm.update(1)
        f.close()
        self.total_time = end - start

    def _get_unique_pepids(self, par, cursor, ignore_genomeoccurence=False):
        query = """
        select parent_id, q1, q1_charge, ssrcalc, peptide.id, modified_sequence
         from %s
         inner join
         ddb.peptide on peptide.id = %s.peptide_key
         inner join ddb.peptideOrganism on peptide.id = peptideOrganism.peptide_key 
         where genome_occurence = 1
         %s
        """ % (par.peptide_table, par.peptide_table, par.query_add )
        if ignore_genomeoccurence:
            query = """
            select parent_id, q1, q1_charge, ssrcalc, peptide_key, modified_sequence
             from %s
             where 4 = 4
             %s
            """ % (par.peptide_table, par.query_add )
        if par.print_query: print query
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
            }
            for r in res
        ]

    def _get_unique_pepids_toptransitions(self, par, cursor):
        #TODO we select all peptides that are in the peplink table
        #regardless whether their charge is actually correct
        query = """
        select parent_id, q1, q1_charge, ssrcalc, peptide.id, m.sequence
         from %s  srmPep
         inner join ddb.peptide on peptide.id = srmPep.peptide_key
         inner join ddb.peptideOrganism on peptide.id = peptideOrganism.peptide_key 
         inner join hroest.MRMPepLink_final l on l.peptide_key = peptide.id
         inner join hroest.MRMAtlas_qtrap_final_no_pyroGlu m on m.id = l.mrm_key
         where genome_occurence = 1
         and l.charge = q1_charge
         #make sure that the modifications are the same!
         and modified_sequence = left(m.sequence, length(m.sequence) -2)
         %s
         group by parent_id
        """ % (par.peptide_table, par.query_add )
        if par.print_query: print query
        cursor.execute( query )
        res = cursor.fetchall()
        return [
            {
                'parent_id' :  r[0],
                'q1' :         r[1],
                'q1_charge' :  r[2],
                'ssrcalc' :    r[3],
                'peptide_key' :r[4],
                'mod_sequence' :r[5].split('/')[0],
            }
            for r in res
        ]

    def _get_all_transitions_toptransitions(self, par, pep, cursor, values = 'q3, m.id'):
        q3_low, q3_high = par.get_q3range_transitions()
        query1 = """
        select %(values)s
        from %(pep)s p 
        inner join hroest.MRMPepLink_final l on p.peptide_key = l.peptide_key
        inner join hroest.MRMAtlas_qtrap_final_no_pyroGlu m on m.id = l.mrm_key
        where m.parent_charge = p.q1_charge
        and m.sequence = '%(mod_seq)s'
        and p.peptide_key = %(peptide_key)s
        and q3 > %(q3_low)s and q3 < %(q3_high)s         
        and m.parent_charge = %(q1_charge)s
        #make sure that the modifications are the same!
        and modified_sequence = left(m.sequence, length(m.sequence) -2)
        %(query_add)s
        """ % { 'peptide_key' : pep['peptide_key'], 'q3_low' : q3_low,
               'q3_high' : q3_high, 'query_add' : par.query1_add,
               'pep' : par.peptide_table, 'trans' : par.transition_table, 
               'values' : values, 'q1_charge' : pep['q1_charge'], 
               'mod_seq' : pep['mod_sequence'] + '/' + str(pep['q1_charge'])}
        if par.print_query: print query1
        cursor.execute( query1 )
        return cursor.fetchall()

    def _get_all_transitions(self, par, pep, cursor, values = "q3, srm_id"):
        q3_low, q3_high = par.get_q3range_transitions()
        query1 = """
        select %(values)s
        from %(pep)s
        inner join %(trans)s
          on %(pep)s.transition_group = %(trans)s.group_id
        where parent_id = %(parent_id)s
        and q3 > %(q3_low)s and q3 < %(q3_high)s         
        %(query_add)s
        """ % { 'parent_id' : pep['parent_id'], 'q3_low' : q3_low,
               'q3_high' : q3_high, 'query_add' : par.query1_add,
               'pep' : par.peptide_table, 'trans' : par.transition_table, 
               'values' : values}
        if par.print_query: print query1
        cursor.execute( query1 )
        return cursor.fetchall()

    def _get_all_collisions_per_transition(self, par, pep, transitions, cursor):
        assert False #not supported any more
        collisions = []
        for q3, id in transitions:
            collisions.extend( self._get_collisions_per_transition(par, pep, q3, cursor) )
        return  collisions

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
            cursor.execute( query2 + txt )
        return cursor.fetchall()

    def _get_collisions_per_transition(self, par, pep, q3, cursor, 
               values='q3, q1, srm_id, peptide_key'):
        q3_low, q3_high = par.get_q3_window_transitions(q3)
        #we compare the parent ion against 4 different parent ions
        #thus we need to take the PEPTIDE key here
        query2 = """
        select %(values)s
        from %(pep)s
        inner join %(trans)s
          on %(pep)s.transition_group = %(trans)s.group_id
        where ssrcalc > %(ssrcalc)s - %(ssr_window)s 
            and ssrcalc < %(ssrcalc)s + %(ssr_window)s
        and q1 > %(q1)s - %(q1_window)s and q1 < %(q1)s + %(q1_window)s
        and %(pep)s.peptide_key != %(peptide_key)d
        and q3 > %(q3_low)s and q3 < %(q3_high)s
        %(query_add)s
        """ % { 'q1' : pep['q1'], 'ssrcalc' : pep['ssrcalc'], 
                'peptide_key' : pep['peptide_key'],
               'q3_low':q3_low,'q3_high':q3_high, 'q1_window' : par.q1_window,
               'query_add' : par.query2_add, 'ssr_window' : par.ssrcalc_window,
               'pep' : par.peptide_table, 'trans' : par.transition_table, 
               'values': values}
        cursor.execute( query2 )
        return cursor.fetchall()

    def store_object(self, par):
        import pickle
        common_filename = par.get_common_filename()
        pickle.dump( self, open(common_filename + '.pkl' , 'w'))

    def store_in_file(self, par):
        import pickle
        common_filename = par.get_common_filename()
        try:
            pickle.dump( self.q1min_distr, open(common_filename + '_q1min_distr.pkl' , 'w'))
        except Exception: pass
        try: pickle.dump( self.q3min_distr, open(common_filename + '_q3min_distr.pkl' , 'w'))
        except Exception: pass
        try: pickle.dump( self.q3min_distr_ppm, open(common_filename + '_q3min_distr_ppm.pkl', 'w'))
        except Exception: pass
        try: pickle.dump( self.q3all_distr_ppm, open(common_filename + '_q3all_distr_ppm.pkl', 'w'))
        except Exception: pass
        try: pickle.dump( self.q1all_distr, open(common_filename + '_q1all_distr.pkl', 'w'))
        except Exception: pass
        try: pickle.dump( self.count_pair_collisions, open(common_filename + '_pair_collisions.pkl', 'w'))
        except Exception: pass
        pickle.dump( self.allpeps, open(common_filename + '_allpeps.pkl', 'w'))
        pickle.dump( [self.non_unique_count, self.total_count], open(common_filename + '_count.pkl', 'w'))

    def load_from_file(self, par, directory):
        import pickle
        #fname = directory + "yeast_True_False_20_250_100ppm_range300to2000"
        fname = directory + par.get_common_filename()
        try:
            self.q1min_distr = pickle.load( open(fname + "_q1min_distr.pkl"))
        except Exception: pass
        try:
            self.q3min_distr = pickle.load( open(fname + "_q3min_distr.pkl"))
        except Exception: pass
        try:
            self.q3min_distr_ppm = pickle.load( open(fname + "_q3min_distr_ppm.pkl"))
        except Exception: pass
        self.allpeps = pickle.load( open(fname + "_allpeps.pkl"))
        self.non_unique_count, self.total_count = pickle.load( open(fname + "_count.pkl"))
        self.pepids = [ (k,) for k in self.allpeps.keys()]

    def print_unique_histogram(self, par):
        mydist = [  self.allpeps[ p[0] ] for p in self.pepids]
        h, n = numpy.histogram( mydist , 10)
        n = [nn * 100.0 + 5 for nn in n]
        h = [hh * 100.0 / len(mydist) for hh in h]
        filename = par.get_common_filename()
        gnu = gnuplot.Gnuplot.draw_boxes_from_data( [h,n], filename + '.eps',
          'Unique transitions per peptide / %' , 'Occurence / %', keep_data = True,
                                             tmp_csv = filename + '.csv' )
        #gnu.add_to_body( "set yrange[0:80]" )
        #gnu.draw_boxes(keep_data = True)

    def print_cumm_unique(self, par):
        mydist = [  self.allpeps[ p[0] ] for p in self.pepids]
        h, n = numpy.histogram( mydist , 100)
        n = [nn * 100.0 for nn in n]
        h = [hh * 100.0 / len(mydist) for hh in h]
        cumh = get_cum_dist( h )
        h.append(100)
        cumh.append(100)
        filename = par.get_common_filename() + "_cum"
        plt_cmd = 'with lines lt -1 lw 2'
        gnu = gnuplot.Gnuplot.draw_from_data( [cumh,n], plt_cmd, filename + '.eps',
          'Unique transitions per peptide / %' , 'Cumulative Occurence / %', 
          keep_data=True, tmp_csv = filename + '.csv')
        gnu.add_to_body( "set yrange[0:100]" )
        gnu.draw(plt_cmd, keep_data=True)
        #title = '(2+,1+) transitions, Background of 4 combinations [Range 300 to 1500 Da]' )

    def print_cumm_unique_all(self, par, cursor):
        """ This is to plot the absolute distribition (transitions per peptide)
        not in percent
        """
        query = """
        select parent_id, count(*) as nr_transitions
        from %(peptable)s
        inner join %(tratable)s on transition_group = group_id 
        where 
        #q1 between 400 and 1200 
        q3 between 400 and 1200 
        and q3_charge = 1
        group by parent_id
        """ % {'peptable' : par.peptide_table ,
               'tratable' : par.transition_table,
                'query_add' : par.query1_add }
        cursor.execute( query)
        lines = cursor.fetchall()
        pep_nr_transitions = {}
        for l in lines: pep_nr_transitions[ l[0] ] = l[1]
        self.allpeps_abs = {}
        for p_id, pep in self.allpeps.items():
            self.allpeps_abs[ p_id ] = pep * pep_nr_transitions[p_id]
        #
        mydist = self.allpeps_abs.values()
        filename = par.get_common_filename() + "_cum_abs"
        #def print_cumm_unique(self, par, mydist, filename):
        h, n = numpy.histogram( mydist, 21, (0,21) )
        n = [nn for nn in n]
        h = [hh * 100.0 / len(mydist) for hh in h]
        cumh = collider.get_cum_dist( h )
        h.append(100)
        cumh.append(100)
        plt_cmd = 'with lines lt -1 lw 2'
        gnu = gnuplot.Gnuplot.draw_from_data( [cumh,n], plt_cmd, filename + '.eps',
          'Unique transitions per peptide' , 'Cumulative Occurence / %', 
          keep_data=True, tmp_csv = filename + '.csv')
        gnu.add_to_body( "set yrange[0:100]" )
        gnu.add_to_body( "set mxtics 5" )
        gnu.draw(plt_cmd, keep_data=True)
        #title = '(2+,1+) transitions, Background of 4 combinations [Range 300 to 1500 Da]' )

    def print_q3min(self, par):
        h, n = numpy.histogram( self.q3min_distr, 100)
        filename = par.get_common_filename() + '_q3distr' 
        n = [nn * 10**3 for nn in n]
        gnu  = gnuplot.Gnuplot.draw_boxes_from_data( [h,n], filename + '.eps',
          'Q3 difference / 10^{-3} Th', 'Number of transitions', keep_data = True )
        #gnu.add_to_body( "set xrange[%s:%s]"  % (-1, 1) )
        #gnu.add_to_body( "set yrange[0:100]" )
        #gnu.draw_boxes()

    def print_q3min_ppm(self, par):
        h, n = numpy.histogram( self.q3min_distr_ppm, 100, (-par.q3_window, par.q3_window) )
        filename = par.get_common_filename() + '_q3distr_ppm' 
        gnu  = gnuplot.Gnuplot.draw_boxes_from_data( [h,n], filename + '.eps',
          'Q3 difference / ppm', 'Number of transitions', keep_data = True )
        #gnu.add_to_body( "set xrange[%s:%s]"  % (-0.1, 0.1) )
        #gnu.add_to_body( "set xrange[%s:%s]"  % (-5, 5) )
        #gnu.draw_boxes()

    def print_q1all(self, par, bars = 50):
        window = par.q1_window
        filename = par.get_common_filename() + '_q1all_distr' 
        h, n = numpy.histogram( self.q1all_distr, bars, (-window, window) )
        shift = window * 1.0 / bars 
        n = [nn + shift for nn in n]
        gnu  = gnuplot.Gnuplot.draw_boxes_from_data( [h,n], filename + '.eps',
          'Q1 difference / Th', 'Number of transitions', keep_data = True )
        os.system( 'epstopdf %(a)s; rm %(a)s' % {'a' : filename + '.eps'})

    def print_q3all_ppm(self, par, bars = 50):
        window = par.q3_window
        if not par.ppm: window = par.q3_window  * 3000
        h, n = numpy.histogram( self.q3all_distr_ppm, bars, (-window, window) )
        shift = window * 1.0 / bars 
        n = [nn + shift for nn in n]
        #n = n[1500:4500]
        #h = h[1500:4500]
        filename = par.get_common_filename() + '_q3all_distr_ppm' 
        gnu  = gnuplot.Gnuplot.draw_boxes_from_data( [h,n], filename + '.eps',
          'Q3 difference / PPM', 'Number of transitions', keep_data = True )
        os.system( 'epstopdf %(a)s; rm %(a)s' % {'a' : filename + '.eps'})
        #os.system( 'convert %(a)s.pdf %(a)s.png' % {'a' : filename})

    def print_stats(self):
        print "Nonunique / Total transitions : %s / %s = %s" % (
            self.non_unique_count, self.total_count, 
            self.non_unique_count * 1.0 /self.total_count)
        below_1ppm = len( [q3 for q3 in self.q3min_distr_ppm if abs(q3) < 1] ) * 1.0/ len( self.q3min_distr_ppm)
        print('Percentage of collisions below 1 ppm: %02.2f~\%%' % (below_1ppm*100) )


def print_trans_collisions(par, db, p_id = 1, q3_low = 300, q3_high = 2000,
    query1_add = 'and q1_charge = 2 and q3_charge = 1'):
    res_str = ''
    cursor = db.cursor()
    non_unique = {}
    non_unique_q1 = {}
    non_unique_ppm = {}
    non_unique_clash = {}
    query1 = """
    select q3, srm_id, q1, ssrcalc, sequence, type, %(pep)s.peptide_key
    from %(pep)s
    inner join %(trans)s
      on %(pep)s.transition_group = %(trans)s.group_id
    inner join ddb.peptide on %(pep)s.peptide_key = peptide.id
    where parent_id = %(parent_id)s
    and q3 > %(q3_low)s and q3 < %(q3_high)s         %(query_add)s
    """ % { 'parent_id' : p_id, 'q3_low' : q3_low,
           'q3_high' : q3_high, 'query_add' : par.query1_add,
           'pep' : par.peptide_table, 'trans' : par.transition_table }
    nr_transitions = cursor.execute( query1 )
    transitions = cursor.fetchall()
    q1 = transitions[0][2]
    ssrcalc = transitions[0][3]
    sequence = transitions[0][4]
    peptide_key = transitions[0][6]
    prt = '\nAnalysing 2+ fragment ions of peptide nr %s\nWith sequence %s, q1 at %s and ssrcalc %s\n'
    res_str += prt % (p_id, sequence, q1, ssrcalc) 
    res_str += 'q3\tion      q3\tion q1   SSRcal dq3   sequence   q3charge\n'
    #
    query2 = """
    select q3, q1, ssrcalc, sequence, type, q3_charge, fragment_number
    from %(pep)s
    inner join %(trans)s
      on %(pep)s.transition_group = %(trans)s.group_id
    inner join ddb.peptide on %(pep)s.peptide_key = peptide.id
    where ssrcalc > %(ssrcalc)s - %(ssr_window)s 
        and ssrcalc < %(ssrcalc)s + %(ssr_window)s
    and q1 > %(q1)s - %(q1_window)s and q1 < %(q1)s + %(q1_window)s
    and %(pep)s.peptide_key != %(peptide_key)d
    and q3 > %(q3_low)s and q3 < %(q3_high)s
    %(query_add)s
    """ % { 'q1' : q1, 'ssrcalc' : ssrcalc, 'peptide_key' : peptide_key,
           'q3_low':q3_low,'q3_high':q3_high, 'q1_window' : par.q1_window,
           'query_add' : par.query2_add, 'ssr_window' : par.ssrcalc_window,
           'pep' : par.peptide_table, 'trans' : par.transition_table }
    tmp = cursor.execute( query2 )
    collisions = cursor.fetchall()
    #here we loop through all possible combinations of transitions and
    #potential collisions and check whether we really have a collision
    q3_window_used = par.q3_window
    isy = True; jj = len( transitions ) / 2 + 1
    jj = 0
    for ii,t in enumerate(transitions):
        ##before we didnt have y and b ions labelled
        if ii == len( transitions ) / 2: res_str += ( "=" * 75 + '\n'); isy = False; jj = 0
        #if isy: jj -= 1
        #else:  jj += 1
        jj += 1
        if par.ppm: q3_window_used = par.q3_window * 10**(-6) * t[0]
        this_min = q3_window_used + 1
        for c in collisions:
            if abs( t[0] - c[0] ) <= this_min:
                this_min = abs( t[0] - c[0] )
                non_unique[ t[1] ] = c
                non_unique_q1[ t[1] ] = q1 - c[1]
                non_unique_ppm[ t[1] ] = (t[0] - c[0] ) * 10**6 / t[0]
        if this_min < q3_window_used:
            c = non_unique[ t[1] ]
            res_str += ("""%.2f\t%s%s <==> %.2f\t%s%s %07.2f %.2f %.2f %s %s\n""" %
                  (t[0], t[5], jj, c[0], c[4], c[6], c[1], c[2], t[0] - c[0], c[3],
                   c[5] ) )
        else:
            res_str += ("""%.2f\t%s%s \n""" % (t[0], t[5], jj) )
    return res_str

class TransitionPeptide(object):
    def __init__(self):
        self.transitions = []
    def someprint(self):
        res_str = ''
        prt = '\nAnalysing peptide nr %s\nWith sequence %s, q1 at %s and ssrcalc %s\n'
        res_str += prt % (self.p_id, self.sequence, self.q1, self.ssrcalc) 
        res_str += 'q3\tion      q3\tion q1   SSRcal dq3   sequence   q3charge\n'
        return res_str
    def latex_describe(self):
        prt = '\nAnalysing peptide nr %s\nWith sequence %s, q1 at %s and ssrcalc %s\n'
        prt = 'Analysing peptide %s \\\\ With sequence %s, q1 at %s and ssrcalc %s'
        return prt % (self.p_id, self.sequence, self.q1, self.ssrcalc) 

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
        collisions_per_peptide = c_getnonuis.getnonuis(
            transitions, collisions, par.q3_window, par.ppm)
        for order in range(1,MAX_UIS+1):
            non_uis_list[order] = c_getnonuis.get_non_uis(
                collisions_per_peptide, order)

        if forceset: return [set(k.keys()) for k in non_uis_list]
        return non_uis_list

    except ImportError:
        #old way of doing it
        return get_non_UIS_from_transitions_old(transitions, collisions, par, MAX_UIS)


def _get_coll_per_peptide_sub(self, transitions, par, pep, cursor):

    try:
        #use range tree, really fast = 50
        #needs self.c_rangetree and self.parentid_lookup
        import c_getnonuis
        q3_low, q3_high = par.get_q3range_collisions()
        q1 = pep['q1']
        ssrcalc = pep['q1']
        q1_low = q1 - par.q1_window
        q1_high = q1 + par.q1_window
        ssrcalc_low = ssrcalc - par.ssrcalc_window
        ssrcalc_high = ssrcalc + par.ssrcalc_window
        #
        precursor_ids = tuple(self.c_rangetree.query_tree( 
            q1_low, -9999, q1_high,  9999 )  )
        precursors = tuple([self.parentid_lookup[myid[0]] for myid in precursor_ids
                            #dont select myself 
                           if parentid_lookup[myid[0]][2]  != pep['peptide_key']])
        return c_getnonuis.calculate_collisions_per_peptide( 
            transitions, precursors, q3_low, q3_high, par.q3_window, par.ppm)
    except AttributeError, ImportError:
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        #fastest = 100 
        import c_getnonuis
        precursors = self._get_all_precursors(par, pep, cursor)
        return c_getnonuis.calculate_collisions_per_peptide( 
            transitions, precursors, q3_low, q3_high, par.q3_window, par.ppm)

def get_coll_per_peptide(self, transitions, par, pep, cursor,
        do_not_calculate=False, forceNonCpp=False):
    if do_not_calculate:
        # slowest =  1000
        # only if we are forced to do that
        # e.g. if the precomputed transitions are somehow special
        # (from experimental data)
        #
        collisions = self._get_all_collisions(par, pep, cursor, transitions = transitions)
    else:
        try: 
            #try to use rangetree, really fast 50ms or less
            #if that doesnt work, try c++ calccollperpep (100ms)
            if forceNonCpp: import somedummymodulethatwillneverexist
            return _get_coll_per_peptide_sub(self, transitions, par, pep, cursor)
        except ImportError:
            #second-fastest = 522 
            #if we dont have any C++ code compiled
            collisions = list(self._get_all_collisions_calculate(par, pep, cursor))
    collisions_per_peptide = {}
    q3_window_used = par.q3_window
    for t in transitions:
        if par.ppm: q3_window_used = par.q3_window * 10**(-6) * t[0]
        for c in collisions:
            if abs( t[0] - c[0] ) <= q3_window_used:
                #gets all collisions
                if collisions_per_peptide.has_key(c[3]):
                    if not t[1] in collisions_per_peptide[c[3]]:
                        collisions_per_peptide[c[3]].append( t[1] )
                else: collisions_per_peptide[c[3]] = [ t[1] ] 
    return collisions_per_peptide 



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
            else: get_non_uis(pepc, non_uis_list[i], i)
    return non_uis_list

def get_UIS_from_transitions(transitions, collisions, par, MAX_UIS):
    """ Get all combinations that are UIS """
    uis_list = [ [] for i in range(MAX_UIS+1)]
    non_uis_list = get_non_UIS_from_transitions(transitions, collisions, par, MAX_UIS)
    srm_ids = [t[1] for t in transitions]
    for i in range(1,MAX_UIS+1):
        uis_list[i] = get_uis(srm_ids, non_uis_list[i], i)
    return uis_list

""
###########################################################################
#
# Setup Functions
#
def get_peptide_from_table(t, row):
    peptide = DDB.Peptide()
    peptide.set_sequence( t.row(row, 'sequence')  )
    peptide.genome_occurence = t.row( row, 'genome_occurence' )
    peptide.ssr_calc         = t.row( row, 'ssrcalc' )
    #peptide.gene_id          = t.row( row, 'gene_id' )
    peptide.mw               = t.row( row, 'molecular_weight' )
    peptide.id               = t.row( row, 'peptide_key' )
    peptide.pairs            = []
    peptide.non_unique = {}
    return peptide

def insert_peptide_in_db(self, db, peptide_table, transition_group):
    c = db.cursor()
    #insert peptide into db
    vals = "peptide_key, q1_charge, q1, ssrcalc, modified_sequence, isotope_nr, transition_group"
    q = "insert into %s (%s) VALUES (%s,%s,%s,%s,'%s', %s, %s)" % (
        peptide_table,
        vals, 
        self.id, self.charge, 
        self.charged_mass, self.ssr_calc, 
        self.get_modified_sequence(),
        0, #we only have the 0th isotope (0 C13 atoms)
        transition_group
    )
    c.execute(q)
    self.parent_id = db.insert_id()

def get_actual_mass(self):
    return self.peptide_mass / self.ion_charge

def insert_in_db(self, db, fragment_charge, transition_table):
    assert False #use fast now
    c = db.cursor()
    peptide = self.ass_peptide
    parent_id = peptide.parent_id
    vals = "type, parent_key, q3_charge, q3"
    ch = fragment_charge
    charged_y_series =  [ ( pred + (ch -1)*self.mass_H)/ch for pred in self.y_series ]
    charged_b_series =  [ ( pred + (ch -1)*self.mass_H)/ch for pred in self.b_series ]
    self.fragment_ids[ ch ] = [ [], [] ]
    for i, q3 in enumerate(charged_b_series):
        type = 'b'
        q = "insert into %s (%s) VALUES ('%s',%s,%s,%s)" % (
            transition_table,
            vals, 
            type, parent_id, fragment_charge, q3 )
        c.execute(q)
        self.fragment_ids[ch][0].append( db.insert_id()  )
    pepLen = len(peptide)
    for i, q3 in enumerate(charged_y_series):
        type = 'y'
        q = "insert into %s (%s) VALUES ('%s',%s,%s,%s)" % (
            transition_table,
            vals, 
            type, parent_id, fragment_charge, q3 )
        c.execute(q)
        self.fragment_ids[ch][1].append( db.insert_id()  )

def fast_insert_in_db(self, db, fragment_charge, transition_table, transition_group):
    c = db.cursor()
    vals = "type, fragment_number, group_id, q3_charge, q3 "
    q = "insert into %s (%s)" % (transition_table, vals)  + " VALUES (%s,%s,%s,%s,%s)" 
    ch = fragment_charge
    tr = len(self.y_series)
    charged_y =  [ ( pred + (ch -1)*self.mass_H)/ch for pred in self.y_series ]
    charged_b =  [ ( pred + (ch -1)*self.mass_H)/ch for pred in self.b_series ]
    many = [ ['y', i+1, transition_group, ch, q3] for i, q3 in enumerate(reversed(charged_y))] 
    manyb = [ ['b', i+1, transition_group, ch, q3] for i, q3 in enumerate(charged_b) ]
    many.extend( manyb )
    c.executemany( q, many)
    #
    #here we could recover the inserted ids but since they map to several 
    #peptides, this is not helpful
    return
    first_id = db.insert_id()
    self.fragment_ids[ ch ] = {
    'y' : [i for i in range(first_id, first_id + len(charged_y )) ] , 
    'b' : [i for i in range(first_id + len(charged_y ),  first_id + 2*len(charged_y )) ]
    }

def all_calculate_clashes_in_series_insert_db( S, S2, pairs_dic, 
                      MS2_bins, clash_distr, range_small, c, table, frag_range):
    for ch1 in frag_range:
        for ch2 in frag_range:
            calculate_clashes_in_series_insert_db( S, S2, ch1, ch2, pairs_dic, 
                                 MS2_bins, clash_distr, range_small, c, table )

class NonExecuteableDatabaseCursor:
    def __init__(self): pass
    def execute(self, t): print t

def calculate_clashes_in_series_insert_db(S, S2, charge1, charge2, pairs_dic, 
                        MS2_bins, clash_distr, myrange, c, table):
    """Compares all b/y ions of spectra S against those of peptide S2.
    If it finds fragments that are within the MS2_bins distance, it will 
    increase share (the return value) by one and flag the corresponding 
    fragment in BOTH associated peptides as non-unique.
    """
    #charge1 = ch1
    #charge2 = ch2
    vals = "coll_srm1, coll_srm2"
    #fragment_ids[ 0 ] = b, [1] = y
    share = 0
    pep1 = S.ass_peptide
    pep2 = S2.ass_peptide
    ch = charge1
    charged_y_series1=[ ( pred + (ch -1)*S.mass_H)/ch for pred in S.y_series ]
    charged_b_series1=[ ( pred + (ch -1)*S.mass_H)/ch for pred in S.b_series ]
    fragment_ids1 = S.fragment_ids[ ch ]
    ch = charge2
    charged_y_series2=[ ( pred + (ch -1)*S2.mass_H)/ch for pred in S2.y_series ]
    charged_b_series2=[ ( pred + (ch -1)*S2.mass_H)/ch for pred in S2.b_series ]
    fragment_ids2 = S2.fragment_ids[ ch ]
    #c = NonExecuteableDatabaseCursor() #for testing
    #every time we count share + 1 we cannot use one b or y
    #only exception: if it clashes with BOTH series, but how
    #often does that happen?
    for i, y in enumerate( charged_y_series1 ):
        if int(y) not in myrange: continue
        for kk, yy in enumerate( charged_y_series2 ):
            if abs( y - yy) < MS2_bins: 
                c.execute("insert into hroest.%s (%s) VALUES (%s,%s)" % ( 
                    table, vals, fragment_ids1[1][i], fragment_ids2[1][kk] )  )
        for kk, bb in enumerate( charged_b_series2 ):
            if abs( y - bb) < MS2_bins: 
                c.execute("insert into hroest.%s (%s) VALUES (%s,%s)" % ( 
                    table, vals, fragment_ids1[1][i], fragment_ids2[0][kk] )  )
    for i, b in enumerate( charged_b_series1 ):
        if int(b) not in myrange: continue
        for kk, yy in enumerate( charged_y_series2 ):
            if abs( b - yy) < MS2_bins: 
                c.execute("insert into hroest.%s (%s) VALUES (%s,%s)" % (
                    table, vals, fragment_ids1[0][i], fragment_ids2[1][kk] )  )
        for kk, bb in enumerate( charged_b_series2 ):
            if abs( b - bb) < MS2_bins: 
                c.execute("insert into hroest.%s (%s) VALUES (%s,%s)" % ( 
                    table, vals, fragment_ids1[0][i], fragment_ids2[0][kk] ) )
    if share > 0:
        pep1.shared_trans.append( pep2 )
        pep2.shared_trans.append( pep1 )
        #if mypeptide.spectrum in (S, S2): 
        #print "share", share, S.ass_peptide.sequence, S2.ass_peptide.sequence, sorted( mypeptide.non_unique)
    #if share > len( S.y_series):
    #    S.ass_peptide.pairs.append( S2 )
    #    S2.ass_peptide.pairs.append( S )
    #    if pairs_dic.has_key( S.peptide ): 
    #        pairs_dic[ S.peptide ].append( S2.peptide )
    #    else: pairs_dic[ S.peptide ] = [ S2.peptide ]
    ##clash_distr[ share ] += 1
    return share

def reset_pairs_unique(mass_bins):
    #reset pairs and non_unique entries to start level
    i = 0
    for b in mass_bins:
        i += len( b )
        for pep in b:
            pep.non_unique = {}
            pep.pairs = []
            pep.shared_trans = []

def calculate_clashes_in_series(S, S2, charge1, charge2, pairs_dic, 
                        MS2_bins, clash_distr, myrange):
    """Compares all b/y ions of spectra S against those of peptide S2.
    If it finds fragments that are within the MS2_bins distance, it will 
    increase share (the return value) by one and flag the corresponding 
    fragment in BOTH associated peptides as non-unique.
    """
    #charge1 = ch1
    #charge2 = ch2
    share = 0
    pep1 = S.ass_peptide
    pep2 = S2.ass_peptide
    ch = charge1
    charged_y_series1=[ ( pred + (ch -1)*S.mass_H)/ch for pred in S.y_series ]
    charged_b_series1=[ ( pred + (ch -1)*S.mass_H)/ch for pred in S.b_series ]
    ch = charge2
    charged_y_series2=[ ( pred + (ch -1)*S2.mass_H)/ch for pred in S2.y_series ]
    charged_b_series2=[ ( pred + (ch -1)*S2.mass_H)/ch for pred in S2.b_series ]
    #every time we count share + 1 we cannot use one b or y
    #only exception: if it clashes with BOTH series, but how
    #often does that happen?
    for i, y in enumerate( charged_y_series1 ):
        if int(y) not in myrange: continue
        for kk, yy in enumerate( charged_y_series2 ):
            if abs( y - yy) < MS2_bins: 
                share += 1
                pep1.non_unique[ 'y%s' % (i) ] = ''
                pep2.non_unique[ 'y%s' % (kk) ] = ''
        for kk, bb in enumerate( charged_b_series2 ):
            if abs( y - bb) < MS2_bins: 
                share += 1
                pep1.non_unique[ 'y%s' % (i) ] = ''
                pep2.non_unique[ 'b%s' % (kk) ] = ''
    for i, b in enumerate( charged_b_series1 ):
        if int(b) not in myrange: continue
        for kk, yy in enumerate( charged_y_series2 ):
            if abs( b - yy) < MS2_bins: 
                share += 1
                pep1.non_unique[ 'b%s' % (i) ] = ''
                pep2.non_unique[ 'y%s' % (kk) ] = ''
        for kk, bb in enumerate( charged_b_series2 ):
            if abs( b - bb) < MS2_bins: 
                share += 1
                pep1.non_unique[ 'b%s' % (i) ] = ''
                pep2.non_unique[ 'b%s' % (kk) ] = ''
    if share > 0:
        pep1.shared_trans.append( pep2 )
        pep2.shared_trans.append( pep1 )
        #if mypeptide.spectrum in (S, S2): 
        #print "share", share, S.ass_peptide.sequence, S2.ass_peptide.sequence, sorted( mypeptide.non_unique)
    #if share > len( S.y_series):
    #    S.ass_peptide.pairs.append( S2 )
    #    S2.ass_peptide.pairs.append( S )
    #    if pairs_dic.has_key( S.peptide ): 
    #        pairs_dic[ S.peptide ].append( S2.peptide )
    #    else: pairs_dic[ S.peptide ] = [ S2.peptide ]
    clash_distr[ share ] += 1
    return share



def _calculate_transitions_ch(peptides, charges, q3_low, q3_high):
    import Residues
    import DDB 
    R = Residues.Residues('mono')
    for p in peptides:
        q1 = p[0]
        peptide_key = p[2]
        peptide = DDB.Peptide()
        peptide.set_sequence(p[1])
        peptide.charge = 2 #dummy, wont need it later
        peptide.create_fragmentation_pattern(R)
        b_series = peptide.b_series
        y_series = peptide.y_series
        for ch in charges:
            for pred in y_series:
                q3 = ( pred + (ch -1)*R.mass_H)/ch
                if q3 < q3_low or q3 > q3_high: continue
                yield (q3, q1, 0, peptide_key)
            for pred in b_series:
                q3 = ( pred + (ch -1)*R.mass_H)/ch
                if q3 < q3_low or q3 > q3_high: continue
                yield (q3, q1, 0, peptide_key)

def calculate_transitions_ch(peptides, charges, q3_low, q3_high):
    try: 
        import c_getnonuis
        return c_getnonuis.calculate_transitions_ch(
            peptides, charges, q3_low, q3_high)
    except ImportError:
        return list(_calculate_transitions_ch(
            peptides, charges, q3_low, q3_high))


def get_nonuis_list(collisions_per_peptide, MAX_UIS):
    non_uis_list = [set() for i in range(MAX_UIS+1)]
    try:
        import c_getnonuis
        for order in range(1,MAX_UIS+1):
                non_uis_list[order] = c_getnonuis.get_non_uis(
                    collisions_per_peptide, order)
    except ImportError:
        for pepc in collisions_per_peptide.values():
            for i in range(1,MAX_UIS+1):
                get_non_uis(pepc, non_uis_list[i], i)
    return non_uis_list 

###UIS Code
def get_non_uis(pepc, non_uis, order):
    if len( pepc ) >= order: 
        non_uis.update( [tuple(sorted(p)) for p in combinations(pepc, order)] )

def get_non_uis_unsorted(pepc, non_uis, order):
    if len( pepc ) >= order: 
        non_uis.update( [tuple(p) for p in combinations(pepc, order)] )


def choose(i,r):
    assert i > 0
    assert r > 0
    assert i >= r
    if r == i: return 1
    return reduce( lambda x,y: x*y, range(1,i+1) ) * 1.0 / ( 
        reduce( lambda x,y: x*y, range(1,r+1) ) * 
        reduce( lambda x,y: x*y, range(1,i-r+1) ) ) 

def get_uis(srm_ids, non_uis, order):
    #prepare input, make sure its sorted
    non_uis = [ tuple(sorted(list(i))) for i in non_uis]
    srm_ids = sorted( srm_ids )
    #
    #the output of combinations is guaranteed to be sorted, if the input was sorted
    result = combinations(srm_ids, order)
    result = [r for r in result if not r in non_uis]
    return result

def permutations(iterable, r=None):
    #use itertools in 2.6
    #see
    #http://svn.python.org/view/python/trunk/Modules/itertoolsmodule.c?view=markup&pathrev=81889
    import itertools
    try:
        return itertools.permutations( iterable, r)
    except:
        return _permutations( iterable, r)


def _permutations(iterable, r=None):
    # permutations('ABCD', 2) --> AB AC AD BA BC BD CA CB CD DA DB DC
    # permutations(range(3)) --> 012 021 102 120 201 210
    #print iterable, r
    pool = tuple(iterable)
    n = len(pool)
    if r is None: r = n 
    else: r=r
    indices = range(n)
    cycles = range(n, n-r, -1)
    yield tuple(pool[i] for i in indices[:r])
    while n:
        for i in reversed(range(r)):
            cycles[i] -= 1
            if cycles[i] == 0:
                indices[i:] = indices[i+1:] + indices[i:i+1]
                cycles[i] = n - i
            else:
                j = cycles[i]
                indices[i], indices[-j] = indices[-j], indices[i]
                yield tuple(pool[i] for i in indices[:r])
                break
        else:
            return

def _combinations(N, M):
    """All index combinations of M elements drawn without replacement
     from a set of N elements.
    Order of elements does NOT matter."""
    index = range( M )
    while index[0] <= N-M:
        yield index[:]
        index[ M-1 ] += 1
        if index[ M-1 ] >= N:
            #now we hit the end, need to increment other positions than last
            #the last position may reach N-1, the second last only N-2 etc.
            j = M-1
            while j >= 0 and index[j] >= N-M+j: j -= 1
            #j contains the value of the index that needs to be incremented
            index[j] += 1
            k = j + 1
            while k < M: index[k] = index[k-1] + 1; k += 1; 

def combinations(iterable, r):
    """ All combinations of size r of the given iterable set. 
    Order of elements does NOT matter.
    """
    try:
        import itertools
        return itertools.combinations( iterable , r) 
    except:
        #we are before python 2.6
        return [ tuple([iterable[i] for i in indices]) for indices in _combinations( len(iterable) , r)] 



def _combinationsDiffLen(N):
    """All index combinations of M elements drawn without replacement
     from a set of N elements.
    Order of elements does NOT matter."""
    M = len(N)
    index = [0 for i in range( M )]
    while True:
        yield index[:]
        #print index[:] 
        #kk += 1
        #if index[2] == 73: break
        index[ M-1 ] += 1
        #if kk > 7: break
        if index[ M-1 ] >= N[ M-1 ]:
            #now we hit the end, need to increment other positions than last
            j = M-1
            while j >= 0 and index[j] >= N[j]-1: j -= 1
            #j contains the value of the index that needs to be incremented
            #when we are at the end of the interation, j will be -1
            if j <0: break
            index[j] += 1
            k = j + 1
            #set all other positions to zero again
            while k < M: index[k] = 0; k += 1; 

