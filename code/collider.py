#!/usr/bin/python
import sys 
sys.path.append( '/home/hroest/msa/code/tppGhost' )
sys.path.append( '/home/hroest/lib/' )
import time
import numpy
import progress
import gnuplot

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
        if self.ppm: self.ppm_string = "PPM"
        self.experiment_type = """Experiment Type:
        check all four charge states [%s] vs all four charge states [%s] with
        thresholds of SSRCalc %s, Q1 %s (Th), Q3 %s (%s) and a range of %s to %s
        Da for the q3 transitions.  """ % ( not self.do_1vs, not self.do_vs1,
          self.ssrcalc_window*2,  self.q1_window*2, self.q3_window*2, 
          self.ppm_string, self.q3_range[0], self.q3_range[1])
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
        common_filename = 'yeast_%s_%s_%d_%d' % (self.do_1vs, self.do_vs1, 
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

def testcase():
    par = SRM_parameters()
    par.q1_window = 0.7 / 2
    par.q3_window = 1.0 / 2
    par.ppm = False
    par.transition_table = 'hroest.srmTransitions_test'
    par.peptide_table = 'hroest.srmPeptides_test'
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

    def find_clashes(self, db, par):
        #make sure we only get unique peptides
        cursor = db.cursor()
        self.pepids = self._get_unique_pepids(par, cursor)
        self.allpeps = {}
        self.allcollisions = []
        self.q3min_distr = []
        self.q1min_distr = []
        self.q3min_distr_ppm = []
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
                        this_min = abs( t[0] - c[0] )
                        non_unique[ t[1] ] = t[0] - c[0]
                        non_unique_q1[ t[1] ] = q1 - c[1]
                        non_unique_ppm[ t[1] ] = (t[0] - c[0] ) * 10**6 / t[0]
                        all_clashes[ t[1] ] = c[2]
            self.allpeps[ p_id ] = 1.0 - len( non_unique ) * 1.0  / nr_transitions
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
            and q3 > %(q3_low)s and q3 < %(q3_high)s         %(query_add)s
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
            select q3, q1, srm_id
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

    def store_in_file(self):
        import pickle
        common_filename = par.get_common_filename()
        try:
            pickle.dump( self.q1min_distr, open(common_filename + '_q1min_distr.pkl' , 'w'))
        except Exception: pass
        pickle.dump( self.q3min_distr, open(common_filename + '_q3min_distr.pkl' , 'w'))
        pickle.dump( self.q3min_distr_ppm, open(common_filename + '_q3min_distr_ppm.pkl', 'w'))
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

    def print_stats(self):
        print "Nonunique / Total transitions : %s / %s = %s" % (
            self.non_unique_count, self.total_count, 
            self.non_unique_count * 1.0 /self.total_count)
        mydist = [  self.allpeps[ p[0] ] for p in self.pepids]
        h, n = numpy.histogram( mydist , 100)
        print('Percentage of collisions below 1 ppm: %02.2f~\%%' % (sum( h[40:60] )* 100.0 / sum( h )))





