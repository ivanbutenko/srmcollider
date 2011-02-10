#!/usr/bin/python
# -*- coding: utf-8  -*-
# vim:set fdm=marker:
import MySQLdb, sys; sys.path.append( '..')
import collider
                         
from optparse import OptionParser, OptionGroup
usage = "usage: %prog [options] experiment_key max_batch"
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "Prepare Rangetree Options",
                    "These are the options to prepare a rangetree run")

group.add_option("-f", "--file", dest="outfile", default='/tmp/tmp.sh',
                  help="Outfile (/tmp/tmp.sh is default)" , metavar="FILE")
parser.add_option_group(group)


###########################################################################
#Parse options
db = MySQLdb.connect(read_default_file="~/.my.cnf")
cursor = db.cursor()
par = collider.SRM_parameters()
par.parse_cmdl_args(parser)
options, args = parser.parse_args(sys.argv[1:])
filename = options.outfile
exp_key = args[0]
max_batch = float(args[1])
par.__dict__.update( options.__dict__ )
par.q3_range = [400, 1400]
par.eval()
mycollider = collider.SRMcollider()

cursor.execute( """
select q1
from %(peptide_table)s where q1 between %(lowq1)s and %(highq1)s
""" % {'peptide_table' : par.peptide_table, 
              'lowq1'  : par.q3_range[0],
              'highq1' : par.q3_range[1]
      } )

myrange=[]
for i,l in enumerate(cursor.fetchall()):
    if i % max_batch == 0: 
        myrange.append( int(l[0]))

myrange.append( int(l[0]))
myrange[0] = par.q3_range[0]
myrange[-1] = par.q3_range[1]

f = open(filename, 'w')
f.write("START=`uptime | cut -c-9`\n")
for i in range(1,len(myrange)):
    #print myrange[i-1], myrange[i]
    f.write("echo 'Starting CHUNK %s of %s at' `uptime | cut -c-9` '[START:' $START']'\n" 
            % (i, len(myrange)-1))
    f.write('python run_uis.py %s %s %s %s %s\n' % (exp_key, myrange[i-1], 
                                 myrange[i], par.ssrcalc_window, par.peptide_table))
f.close()

