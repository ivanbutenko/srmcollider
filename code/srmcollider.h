/*
 *
 * Program       : SRMCollider
 * Author        : Hannes Roest <roest@imsb.biol.ethz.ch>
 * Date          : 30.05.2012 
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
 */

#ifndef SRMCOLLIDER_H
#define SRMCOLLIDER_H

#include <iostream>
#include <vector>
#include <set>

/*
 * http://www.sisweb.com/referenc/source/exactmaa.htm
 * http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
*/
#define MASS_H 1.007825032
#define MASS_N 14.003074005
#define MASS_O 15.994914620
#define MASS_OH 17.002739651999999
#define MASS_NH3 17.026549101000001
#define MASS_H2O 18.010564683999998
#define MASS_CO 27.994914619999999

#define MASS_C   12.000000
#define MASS_C13 13.00335484
#define MASS_DIFFC13 (MASS_C13 - MASS_C)

//assume that there are never more than 32 transitions in an assay :-)
//please change when an error occurs
#define COMBINT uint32_t
#define COMBLIMIT 32

#define NOISOTOPEMODIFICATION 0
#define N15_ISOTOPEMODIFICATION 1

// Boost.Python headers
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
namespace python = boost::python;

// Function declarations of shared functions
//

// Function to calculate the b and y fragment series given a peptide sequence 
int _calculate_clashes(const char* sequence, double* b_series, double* y_series,
        double ch) ;

//functions to calculate a list of transitions for a given number of peptides
python::list _find_clashes_calculate_clashes_ch(python::tuple precursors,
        python::list charges, double q3_low, double q3_high );
python::list _find_clashes_calculate_clashes(python::tuple precursors,
        double q3_low, double q3_high ); //for charge states 1 and 2
python::list _calculate_clashes_wrapper(python::tuple &tlist, 
        double charge); //for a single peptide

//function to calculate the charged mass of a peptide
double calculate_charged_mass(python::tuple clist, int ch);


void _py_combinations(int M, int N, const python::list &mapping, python::dict &result) ;
void _combinations(int M, int N, std::vector<std::vector<int> >& result) ;
void _combinations_magic(int M, int N, COMBINT* mapping,
        python::dict &result) ;
python::dict get_non_uis(python::dict collisions_per_peptide, int order) ;
void get_non_uis_magic(std::vector<COMBINT>& newcollperpep, int max_tr, int
        order, std::set<COMBINT> & result);

#endif

