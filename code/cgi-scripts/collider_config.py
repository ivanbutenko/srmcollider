#!/usr/bin/python
# -*- coding: utf-8  -*-

"""
 *
 * Program       : SRMCollider
 * Author        : Hannes Roest <roest@imsb.biol.ethz.ch>
 * Date          : 05.02.2011 
 *
 *
 * Copyright (C) 2011 - 2012 Hannes Roest
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

# This file contains some options that may need to be changed to suit your
# local settings. 
# Default installation uses a linux-user called "srmcollider" with the home at
# /home/srmcollider/ and a mysql conf file in the home directory.
# SSRCalc is expected to be installed under home in a ssrcalc/ folder.

# The website is assumed to live at /var/websites/srmcollider and the documents
# (e.g. the csv file) will be placed under documents and the script is expected
# to be placed at srmcollider/srmcollider.py. If you need a different setup,
# please change the lines at "Web setup"

# To configure the available genomes, please edit the point "Available genomes"
# and make sure that the genome_select HTML and the map_db_tables function are
# synchronized.

# [General Options]
SRMCOLLIDER_HOME = '/home/srmcollider/srmcollider/code' # direct to where collider.py lives

# [Database setup]
default_mysql = "/home/srmcollider/.srm.cnf" # custom .my.cnf
db_used = 'srmcollider' 
default_org_prefix = '.srmPeptides_' # common prefix of all database tables
# Database table containing SSRCalc values. Leave empty if not present.
default_ssrcalc = '' 

# [External programs]
# path to SSRCalc3.pl and SSRCalc3.par files
ssrcalc_path = '/home/srmcollider/ssrcalc/' 

# [Web setup]
myCSVFile_ = '/var/websites/srmcollider/documents/srmcollider_'
myUIS_CSVFile_ = '/var/websites/srmcollider/documents/uis_srmcollider_'
myUIS_CSVFile_rel_ = '/../documents/uis_srmcollider_'
myCSVFile_rel_ = '/../documents/srmcollider_'

# [Available genomes]
# If you want to add additional genomes, edit the genome_select HTML and the
# map_db_tables function
genomes_that_require_N15_data = ['yeastN15']

genome_select = """
    <option value="yeast">Yeast (tryptic)</option>
    <option value="human">Human (tryptic)</option>
"""

def map_db_tables(genome):
    # Determine which database should be used, needs to correspond to the
    # "genome_select" variable above.
    if genome == 'yeast':
        table_used =  'yeast_oxMetDeamid_miss1'
    elif genome == 'human':
        table_used =  'human_oxMetDeamid_miss1'
    else: 
        print "Genome not recognized";
        exit()
    return table_used


