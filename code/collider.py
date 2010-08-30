#!/usr/bin/python
import sys 
sys.path.append( '/home/hroest/msa/code/tppGhost' )
sys.path.append( '/home/hroest/lib/' )
import time
import numpy
import progress
import gnuplot

import DDB

#with this method, I can do 
##method 2
# 22.2 per second (1Da), thus do 500k in 6 hours ==> using per peptide methods
# 5 per second (25Da), thus do 500k in 24 hours ==> using per peptide methods
#
#10ppm, 1Da => 200 / s
#10ppm, 9Da => 20 / s
#10ppm, 25Da => 7 / s
#10ppm, 50Da => 4 / s
##method 1
# 0.5 per second, thus do 500k in 240 hours ==> using per transition methods
class SRM_parameters(object):
    def __init__(self): 
        self.do_1vs = True #check only one charge state?
        self.do_vs1 = False #background only one charge state?
        self.dontdo2p2f = True #do not look at 2+ parent / 2+ fragment ions
        self.ppm = True #measure q3 in ppm
        self.transition_table = 'hroest.srmTransitions_yeast'
        self.peptide_table = 'hroest.srmPeptides_yeast'
        self.q3_range = [400, 1200]
        self.ssrcalc_window = 2.0 / 2
        self.q1_window = 25.0 / 2.0
        self.q3_window = 10.0 / 2.0
        self.do_1_only = "and q1_charge = 2 and q3_charge = 1"
    def eval(self):
        self.query_add = ""
        self.query2_add = ""
        self.ppm_string = "Th"
        if self.do_1vs :
            self.query_add = "and q1_charge = 2"
            self.query1_add = self.do_1_only
        if self.do_vs1 : 
            self.query2_add = self.do_1_only
        elif self.dontdo2p2f:
            self.query2_add = "and not (q1_charge = 2 and q3_charge = 2)"
        if self.ppm: self.ppm_string = "PPM"
        self.experiment_type = """Experiment Type:
        check all four charge states [%s] vs all four charge states [%s] with
        thresholds of SSRCalc %s, Q1 %s (Th), Q3 %s (%s) and a range of %s to %s
        Da for the q3 transitions.  Ignore 2+ parent / 2+ fragment ions %s, 
        selecting from %s and %s.""" % ( 
            not self.do_1vs, not self.do_vs1,
          self.ssrcalc_window*2,  self.q1_window*2, self.q3_window*2, 
          self.ppm_string, self.q3_range[0], self.q3_range[1], self.dontdo2p2f, 
          self.peptide_table, self.transition_table)
    def get_q3_high(self, q1, q1_charge):
        q3_high = self.q3_range[1]
        if q3_high < 0: 
            if not self.query1_add == 'and q1_charge = 2 and q3_charge = 1':
                raise( 'relative q3_high not implemented for this combination')
            q3_high = q1 * q1_charge + q3_high
        return q3_high
    def get_q3_low(self):
        return self.q3_range[0]
    def get_common_filename(self):
        common_filename = '%s_%s_%s_%d_%d' % (self.peptide_tbl_identifier, 
                                              self.do_1vs, self.do_vs1, 
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
    def peptide_tbl(self): return self.peptide_table.split('.')[1]

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

    def find_clashes_small(self, db, par):
        #make sure we only get unique peptides
        cursor = db.cursor()
        self.pepids = self._get_unique_pepids(par, cursor)
        self.allpeps = {}
        self.non_unique_count = 0
        self.total_count = 0
        start = time.time()
        progressm = progress.ProgressMeter(total=len(self.pepids), unit='peptides')
        for i, pep in enumerate(self.pepids):
            p_id, q1, q1_charge, ssrcalc = pep
            non_unique = {}
            non_unique_q1 = {}
            non_unique_ppm = {}
            non_unique_clash = {}
            all_clashes = {}
            transitions = self._get_all_transitions(par, pep, cursor)
            nr_transitions = len( transitions )
            collisions = self._get_all_collisions(par, pep, cursor)
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
            end = time.time()
            progressm.update(1)
        self.total_time = end - start

    def find_clashes(self, db, par, toptrans=False):
        #make sure we only get unique peptides
        cursor = db.cursor()
        if not toptrans: self.pepids = self._get_unique_pepids(par, cursor)
        else: self.pepids = self._get_unique_pepids_toptransitions(par, cursor)
        self.allpeps = {}
        self.allcollisions = []
        self.q3min_distr = []
        self.q3all_distr_ppm = []
        self.q1all_distr = []
        self.q1min_distr = []
        self.q3min_distr_ppm = []
        self.non_unique_count = 0
        self.count_pair_collisions = {}
        self.found3good = []
        self.total_count = 0
        start = time.time()
        progressm = progress.ProgressMeter(total=len(self.pepids), unit='peptides')
        for i, pep in enumerate(self.pepids):
            p_id, q1, q1_charge, ssrcalc = pep
            non_unique = {}
            non_unique_q1 = {}
            non_unique_ppm = {}
            non_unique_clash = {}
            all_clashes = {}
            if not toptrans: transitions = self._get_all_transitions(par, pep, cursor)
            else: transitions = self._get_all_transitions_toptransitions(par, pep, cursor)
            nr_transitions = len( transitions )
            collisions = self._get_all_collisions(par, pep, cursor)
            #here we loop through all possible combinations of transitions and
            #potential collisions and check whether we really have a collision
            q3_window_used = par.q3_window
            for t in transitions:
                if par.ppm: q3_window_used = par.q3_window * 10**(-6) * t[0]
                this_min = q3_window_used
                for c in collisions:
                    self.count_pair_collisions[ p_id ] = 0
                    if abs( t[0] - c[0] ) <= q3_window_used:
                        #gets all collisions
                        coll_key = "%s:%s" % (p_id, c[3])
                        try: self.count_pair_collisions[ coll_key ] += 1
                        except KeyError: self.count_pair_collisions[ coll_key ] = 1
                        self.q1all_distr.append( q1 - c[1])
                        self.q3all_distr_ppm.append((t[0] - c[0] ) * 10**6 / t[0])
                    if abs( t[0] - c[0] ) <= this_min:
                        #gets the nearest collision
                        this_min = abs( t[0] - c[0] )
                        non_unique[ t[1] ] = t[0] - c[0]
                        non_unique_q1[ t[1] ] = q1 - c[1]
                        non_unique_ppm[ t[1] ] = (t[0] - c[0] ) * 10**6 / t[0]
                        all_clashes[ t[1] ] = c[2]
            self.allpeps[ p_id ] = 1.0 - len( non_unique ) * 1.0  / nr_transitions
            if len(transitions) - len(non_unique) > 3: self.found3good.append( p_id )
            if len(non_unique) > 0:
                self.allcollisions.extend( all_clashes.items() )
                for v in non_unique.values(): self.q3min_distr.append( v )
                for v in non_unique_q1.values(): self.q1min_distr.append( v )
                for v in non_unique_ppm.values(): self.q3min_distr_ppm.append( v )
            self.non_unique_count += len( non_unique )
            self.total_count += len( transitions)
            end = time.time()
            progressm.update(1)
        self.total_time = end - start

    def _get_unique_pepids(self, par, cursor):
        query = """
        select parent_id, q1, q1_charge, ssrcalc
         from %s
         inner join
         ddb.peptide on peptide.id = %s.peptide_key
         inner join ddb.peptideOrganism on peptide.id = peptideOrganism.peptide_key 
         where genome_occurence = 1
         %s
        """ % (par.peptide_table, par.peptide_table, par.query_add )
        cursor.execute( query )
        return cursor.fetchall()

    def _get_unique_pepids_toptransitions(self, par, cursor):
        query = """
            select parent_id, %(peptable)s.q1 as q1, q1_charge, ssrcalc 
            from %(db)s.%(peptable)s
            inner join
            ddb.peptide on peptide.id = %(peptable)s.peptide_key
            inner join ddb.peptideOrganism on peptide.id = peptideOrganism.peptide_key 
            inner join hroest.yeast_dp_data_light  on
            yeast_dp_data_light.modified_sequence = %(peptable)s.modified_sequence
            where genome_occurence = 1
            and q1_charge = prec_z  
            %(qadd)s
            group by parent_id, q1_charge
            """ % { 'db' : par.peptide_db, 'peptable' : par.peptide_tbl, 'qadd' : par.query_add  }
        #print query
        cursor.execute( query )
        return cursor.fetchall()


    def _get_all_transitions_toptransitions(self, par, pep, cursor):
        p_id, q1, q1_charge, ssrcalc = pep
        q3_high = par.get_q3_high(q1, q1_charge)
        q3_low = par.get_q3_low()
        query1 = """
        select distinct %(transtable)s.q3, srm_id, rank
        from %(pep)s
        inner join %(trans)s
          on parent_id = parent_key
          inner join hroest.yeast_dp_data_light  on
          yeast_dp_data_light.modified_sequence = %(peptable)s.modified_sequence
        where parent_id = %(parent_id)s
        and %(transtable)s.q3 > %(q3_low)s and %(transtable)s.q3 < %(q3_high)s         
        and rank <= 5
        and prec_z = q1_charge
        and frg_type = type
        and frg_z = q3_charge
        and frg_nr = fragment_number
        %(query_add)s
        """ % { 'parent_id' : p_id, 'q3_low' : q3_low,
               'q3_high' : q3_high, 'query_add' : par.query1_add,
               'pep' : par.peptide_table, 'trans' : par.transition_table, 
               'peptable' : par.peptide_tbl, 'transtable' : par.transition_tbl }
        cursor.execute( query1 )
        return cursor.fetchall()

    def _get_all_transitions(self, par, pep, cursor):
            p_id, q1, q1_charge, ssrcalc = pep
            q3_high = par.get_q3_high(q1, q1_charge)
            q3_low = par.get_q3_low()
            query1 = """
            select q3, srm_id
            from %(pep)s
            inner join %(trans)s
              on parent_id = parent_key
            where parent_id = %(parent_id)s
            and q3 > %(q3_low)s and q3 < %(q3_high)s         
            %(query_add)s
            """ % { 'parent_id' : p_id, 'q3_low' : q3_low,
                   'q3_high' : q3_high, 'query_add' : par.query1_add,
                   'pep' : par.peptide_table, 'trans' : par.transition_table }
            cursor.execute( query1 )
            return cursor.fetchall()

    def _get_all_collisions(self, par, pep, cursor):
            p_id, q1, q1_charge, ssrcalc = pep
            q3_high = par.get_q3_high(q1, q1_charge)
            q3_low = par.get_q3_low()
            query2 = """
            select q3, q1, srm_id, peptide_key
            from %(pep)s
            inner join %(trans)s
              on parent_id = parent_key
            where ssrcalc > %(ssrcalc)s - %(ssr_window)s 
                and ssrcalc < %(ssrcalc)s + %(ssr_window)s
            and q1 > %(q1)s - %(q1_window)s and q1 < %(q1)s + %(q1_window)s
            and parent_id != %(parent_id)d
            and q3 > %(q3_low)s and q3 < %(q3_high)s
            %(query_add)s
            """ % { 'q1' : q1, 'ssrcalc' : ssrcalc, 'parent_id' : p_id,
                   'q3_low':q3_low,'q3_high':q3_high, 'q1_window' : par.q1_window,
                   'query_add' : par.query2_add, 'ssr_window' : par.ssrcalc_window,
                   'pep' : par.peptide_table, 'trans' : par.transition_table }
            cursor.execute( query2 )
            return cursor.fetchall()

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
        pickle.dump( self.allpeps, open(common_filename + '_allpeps.pkl', 'w'))
        pickle.dump( [self.non_unique_count, self.total_count], open(common_filename + '_count.pkl', 'w'))

    def load_from_file(self, par, directory):
        import pickle
        #fname = directory + "yeast_True_False_20_250_100ppm_range300to2000"
        fname = directory + par.get_common_filename()
        try:
            self.q1min_distr = pickle.load( open(fname + "_q1min_distr.pkl"))
        except Exception: pass
        self.q3min_distr = pickle.load( open(fname + "_q3min_distr.pkl"))
        self.q3min_distr_ppm = pickle.load( open(fname + "_q3min_distr_ppm.pkl"))
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

    def print_q1all(self, par, bars = 50, window = par.q1_window):
        h, n = numpy.histogram( self.q1all_distr, bars, (-window, window) )
        shift = window * 1.0 / bars 
        n = [nn + shift for nn in n]
        filename = par.get_common_filename() + '_q1all_distr' 
        gnu  = gnuplot.Gnuplot.draw_boxes_from_data( [h,n], filename + '.eps',
          'Q1 difference / Th', 'Number of transitions', keep_data = True )

    def print_q3all_ppm(self, par, bars = 50, window = par.q3_window):
        if not par.ppm: window = par.q3_window  * 3000
        h, n = numpy.histogram( self.q3all_distr_ppm, bars, (-window, window) )
        shift = window * 1.0 / bars 
        n = [nn + shift for nn in n]
        filename = par.get_common_filename() + '_q3all_distr_ppm' 
        gnu  = gnuplot.Gnuplot.draw_boxes_from_data( [h,n], filename + '.eps',
          'Q3 difference / PPM', 'Number of transitions', keep_data = True )
        os.system( 'epstopdf %(a)s; rm %(a)s' % {'a' : filename + '.eps'})
        #os.system( 'convert %(a)s.pdf %(a)s.png' % {'a' : filename})

    def print_stats(self):
        print "Nonunique / Total transitions : %s / %s = %s" % (
            self.non_unique_count, self.total_count, 
            self.non_unique_count * 1.0 /self.total_count)
        mydist = [  self.allpeps[ p[0] ] for p in self.pepids]
        h, n = numpy.histogram( mydist , 100)
        print('Percentage of collisions below 1 ppm: %02.2f~\%%' % (sum( h[40:60] )* 100.0 / sum( h )))

def print_trans_collisions(par, p_id = 1, q3_low = 300, q3_high = 2000,
    query1_add = 'and q1_charge = 2 and q3_charge = 1'):
    cursor = db.cursor()
    non_unique = {}
    non_unique_q1 = {}
    non_unique_ppm = {}
    non_unique_clash = {}
    query1 = """
    select q3, srm_id, q1, ssrcalc, sequence, type
    from %(pep)s
    inner join %(trans)s
      on parent_id = parent_key
    inner join ddb.peptide on peptide_key = peptide.id
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
    prt = '\nAnalysing peptide nr %s\nWith sequence %s, q1 at %s and ssrcalc %s'
    print(prt % (p_id, sequence, q1, ssrcalc) )
    #
    query2 = """
    select q3, q1, ssrcalc, sequence, type, q3_charge
    from %(pep)s
    inner join %(trans)s
      on parent_id = parent_key
    inner join ddb.peptide on peptide_key = peptide.id
    where ssrcalc > %(ssrcalc)s - %(ssr_window)s 
        and ssrcalc < %(ssrcalc)s + %(ssr_window)s
    and q1 > %(q1)s - %(q1_window)s and q1 < %(q1)s + %(q1_window)s
    and parent_id != %(parent_id)d
    and q3 > %(q3_low)s and q3 < %(q3_high)s
    %(query_add)s
    """ % { 'q1' : q1, 'ssrcalc' : ssrcalc, 'parent_id' : p_id,
           'q3_low':q3_low,'q3_high':q3_high, 'q1_window' : par.q1_window,
           'query_add' : par.query2_add, 'ssr_window' : par.ssrcalc_window,
           'pep' : par.peptide_table, 'trans' : par.transition_table }
    tmp = cursor.execute( query2 )
    collisions = cursor.fetchall()
    #here we loop through all possible combinations of transitions and
    #potential collisions and check whether we really have a collision
    q3_window_used = par.q3_window
    isy = True; jj = len( transitions ) / 2 + 1
    for ii,t in enumerate(transitions):
        if ii == len( transitions ) / 2: print( "=" * 75); isy = False; jj = 0
        if isy: jj -= 1
        else:  jj += 1
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
            print("""%.2f\t%s%s <==> %.2f\t%s %07.2f %.2f %.2f %s %s""" % (t[0], t[5], jj, c[0], c[4], c[1], c[2], t[0] - c[0], c[3], c[5] ) )
        else:
            print("""%.2f\t%s%s """ % (t[0], t[5], jj) )

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

def insert_peptide_in_db(self, db, peptide_table):
    c = db.cursor()
    peptide = self.ass_peptide
    #insert peptide into db
    vals = "peptide_key, q1_charge, q1, ssrcalc, modified_sequence"
    q = "insert into %s (%s) VALUES (%s,%s,%s,%s,'%s')" % (
        peptide_table,
        vals, 
        peptide.id, peptide.charge, 
        get_actual_mass(self), peptide.ssr_calc, 
        peptide.get_modified_sequence()
    )
    c.execute(q)
    peptide.parent_id = db.insert_id()

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

def fast_insert_in_db(self, db, fragment_charge, transition_table):
    c = db.cursor()
    peptide = self.ass_peptide
    parent_id = peptide.parent_id
    vals = "type, fragment_number, parent_key, q3_charge, q3 "
    q = "insert into %s (%s)" % (transition_table, vals)  + " VALUES (%s,%s,%s,%s,%s)" 
    ch = fragment_charge
    tr = len(self.y_series)
    charged_y =  [ ( pred + (ch -1)*self.mass_H)/ch for pred in self.y_series ]
    charged_b =  [ ( pred + (ch -1)*self.mass_H)/ch for pred in self.b_series ]
    many = [ ['y', i+1, parent_id, ch, q3] for i, q3 in enumerate(reversed(charged_y))] 
    manyb = [ ['b', i+1, parent_id, ch, q3] for i, q3 in enumerate(charged_b) ]
    many.extend( manyb )
    c.executemany( q, many)
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


#make a list of all peptides

#old analysis

def read_fragment_ids(self, cursor, peptide, peptide_table, transition_table):
    #reads srmPeptide and srmTransitions in order to store the srm_ids of the 
    #fragments. Speed ~ 300 / s
    qq = """select 
    q3_charge, type, srm_id
    from %s 
    inner join %s on parent_id = parent_key
    where peptide_key = %s and q1_charge = %s
    order by q3_charge, type
    """ % (peptide_table, transition_table, peptide.id, peptide.charge)
    cursor.execute( qq ) 
    all_transitions = cursor.fetchall()
    self.fragment_ids = { 1 : [ [], [] ], 2 : [ [], [] ]}
    tmp = 0
    for r in all_transitions:
        if r[1] == 'b' : tmp =0
        else: tmp = 1
        self.fragment_ids[ r[0] ][tmp].append( r[2] )
