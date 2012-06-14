#!/usr/bin/python
# -*- coding: utf-8  -*-

"""
 *
 * Program       : SRMCollider
 * Author        : Hannes Roest <roest@imsb.biol.ethz.ch>
 * Date          : 14.06.2012
 *
 *
 * Copyright (C) 2011 - 2012 Hannes Roest
 *
 * This library is free software; you can redistribute it and/or
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307, USA
 *
"""

import cStringIO
import cgi, sys

from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

import matplotlib
# chose a non-GUI backend
matplotlib.use( 'Agg' )
import pylab
import matplotlib.pyplot as pyplot

X,Y = 500, 275 #image width and height

# Pass data over filesystem
def get_data_from_file(datafile):
    #log = file(r"/tmp/data.txt" )
    log = file("/tmp/%s" % datafile)
    x = []
    y = []
    for (i, value) in enumerate(log):
        xx, yy = value.split('\t')
        x.append(float(xx))
        y.append(float(yy))
    return x,y

# Pass data over HTML, there is a limit 
def parse_data(data):
    x = []
    y = []
    for (i, value) in enumerate(data.split(",,")):
        xx, yy = value.split(',')
        x.append(float(xx))
        y.append(float(yy))
    return x,y

def graph(data, mz, rt, pepname):
    x,y = get_data_from_file(data)
    #x,y = parse_data(data)

    #write to file object
    pl1, = pylab.plot(x,y, 'bo', markersize=4, label="Interferences")
    pl2, = pylab.plot([mz],[rt], 'ro', label=pepname)
 
    pylab.xlabel("m/z", rotation="horizontal")
    pylab.ylabel("RT", rotation="horizontal")
    #pylab.xticks(rotation="vertical") #pylab.xticks(rotation=25)
    pylab.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
        ncol=2, mode="expand", borderaxespad=0.)

    # prevent a constant offset to be displayed
    # see http://stackoverflow.com/questions/3677368/matplotlib-format-axis-offset-values-to-whole-numbers-or-specific-number
    ax = pylab.gca()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%0.3f')) 

    #output to browser
    print "Content-type: image/png\n"
    # save the plot as a png and output directly to webserver
    pylab.savefig( sys.stdout, format='png' )

if __name__ == "__main__":
    form = cgi.FieldStorage()
    if "data" in form:
        data = form["data"].value
        mz = float(form["mz"].value)
        rt = float(form["rt"].value)
        pepname = form["pepname"].value
        graph(data, mz, rt, pepname)
    else:
        print "Content-type: text/html\n"
        print """<html><body>No input file given</body></html>"""

