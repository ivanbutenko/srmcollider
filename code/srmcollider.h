#ifndef SRMCOLLIDER_H
#define SRMCOLLIDER_H

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

#define MASS_H 1.007825032
#define MASS_OH 17.002739651999999

//assume that there are never more than 32 transitions in an assay
#define COMBINT uint32_t
#define COMBLIMIT 32


/*
 *
 *
 * PYTHON interface
*/

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
namespace python = boost::python;

//definitions of shared functions
int _calculate_clashes(const char* sequence, double* b_series, double* y_series,
        double ch) ;
void _combinations(int M, int N, const python::list &mapping, python::dict &result) ;
python::dict get_non_uis(python::dict collisions_per_peptide, int order) ;


#endif

