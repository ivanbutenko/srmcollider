class SRM_parameters(object):

    def __init__(self): 
        self.do_1vs = True #check only one charge state?
        self.do_vs1 = False #background only one charge state?
        self.dontdo2p2f = True #do not look at 2+ parent / 2+ fragment ions
        #self.considerIsotopes = False #do not consider the C13 isotopes
        self.isotopes_up_to = 0
        self.ppm = True #measure q3 in ppm
        self.transition_table = 'srmcollider.srmTransitions_yeast'
        self.peptide_table = 'srmcollider.srmPeptides_yeast'
        self.q3_range = [400, 1200]
        self.ssrcalc_window = 2.0 / 2
        self.q1_window = 25.0 / 2.0
        self.q3_window = 10.0 / 2.0
        self.max_uis = 5 #maximal order of UIS to calculate (no UIS => set to 0)
        self.do_1_only = "and q1_charge = 2 and q3_charge = 1"
        self.print_query = False
        # 
        self.parent_charges = [2,3] # parent ion charge states
        #self.select_floor = False
        self.quiet = False
        self.bions      =  True
        self.yions      =  True
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
        group.add_option("-c", "--config", dest="config_file", default='',
                          help="Configuration file (it overrides cmdline options!)" )
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
        group.add_option("-i", "--isotopes_up_to", dest="isotopes_up_to", 
            default=3, type='int', help="Consider isotopes of the precursors,"+\
            " (default up to 3amu) ")
        group.add_option("--q3_low", dest="q3_low", default=400, type="float",
                          help="Start of transition range to analyse (default 400)", metavar="Q3LOW")
        group.add_option("--q3_high", dest="q3_high", default=1400, type="float",
                          help="End of transition range to analyse (default 1400)", metavar="Q3HIGH")
        group.add_option("--max_uis", dest="max_uis", default=0, type='int',
                          help="maximal order of UIS to calculate " +
                          "(defaults to 0 == no UIS)" )
        group.add_option("--peptide_table", dest="peptide_table", default='srmcollider.srmPeptides_yeast',
                          help="MySQL table containing the peptides" )
        group.add_option("--transition_table", dest="transition_table", default='srmcollider.srmTransitions_yeast',
                          help="MySQL table containing the transitions" )
        group.add_option("--mysql_config", dest="mysql_config", default=default_mysql,
                          help="Location of mysql config file, defaults to %s" % default_mysql )
        group.add_option("-q", "--quiet", dest="quiet", default=False,
                          help="don't print status messages to stdout")
        parser.add_option_group(group)

    def parse_options(self, options):
        if options.config_file != '':
            self.read_parameter_file(options.config_file)
        else: self.__dict__.update( options.__dict__ )

        self.q3_range = [options.q3_low, options.q3_high]
        self.q1_window /= 2.0
        self.q3_window /= 2.0
        self.ssrcalc_window /= 2.0
        if self.ppm == 'True': self.ppm = True
        elif self.ppm == 'False': self.ppm = False
        elif self.ppm in [True, False]: pass
        else: 'wrong arg for ppm'; assert False

    def read_parameter_file(self, thefile):
        parameter = self
        execfile(thefile)

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
        #
        # is default since isotopes are not in the database any more
        self.query2_add += " and isotope_nr = 0 " 
        if self.ppm: self.ppm_string = "PPM"
        self.experiment_type = """Experiment Type:
        check all four charge states [%s] vs all four charge states [%s] with
        thresholds of SSRCalc %s, Q1 %s (Th), Q3 %s (%s) and a range of %s to %s
        Da for the q3 transitions.  Ignore 2+ parent / 2+ fragment ions %s, 
        selecting from %s and %s.
        Consider Isotopes up to: %s""" % ( 
            not self.do_1vs, not self.do_vs1,
          self.ssrcalc_window*2,  self.q1_window*2, self.q3_window*2, 
          self.ppm_string, self.q3_range[0], self.q3_range[1], self.dontdo2p2f, 
          self.peptide_table, self.transition_table, self.isotopes_up_to)
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

    def do_b_y_only(self):

        return (
        self.aions      ==  False and
        self.aMinusNH3  ==  False and
        self.bions      ==  True  and
        self.bMinusH2O  ==  False and
        self.bMinusNH3  ==  False and
        self.bPlusH2O   ==  False and
        self.cions      ==  False and
        self.xions      ==  False and
        self.yions      ==  True  and
        self.yMinusH2O  ==  False and
        self.yMinusNH3  ==  False and
        self.zions      ==  False 
        )

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

def testcase(testdatabase='srmcollider'):
    par = SRM_parameters()
    par.q1_window = 0.7 / 2
    par.q3_window = 1.0 / 2
    par.ppm = False
    par.transition_table = testdatabase + '.srmTransitions_test'
    par.peptide_table = testdatabase + '.srmPeptides_test'
    par.dontdo2p2f = False #do not look at 2+ parent / 2+ fragment ions
    #default 
    #par.considerIsotopes = False #do not consider the C13 isotopes
    par.isotopes_up_to = 0
    par.do_1vs = True #check only one charge state?
    par.do_vs1 = False #background only one charge state?
    par.q3_range = [400, 1200]
    par.ssrcalc_window = 2.0 / 2
    par.do_1_only = "and q1_charge = 2 and q3_charge = 1"
    #
    par.eval()
    return par

