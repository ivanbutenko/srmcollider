/*
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
 */

//include our own libraries
#ifndef SRMCOLLIDERLIB_H
#define SRMCOLLIDERLIB_H
#include "srmcollider.h"

struct SRMTransition
{
  double q1;
  double q3;
  long transition_id;
};

struct SRMPrecursor{
  char* sequence; 
  double q1;
  long transition_group;
  int q1_charge;
  int isotope_modification;
  int maximal_charge;
  double ssrcalc;
};

struct SRMParameters
{
  bool aions      ;
  bool aMinusNH3  ;
  bool bions      ;
  bool bMinusH2O  ;
  bool bMinusNH3  ;
  bool bPlusH2O   ;
  bool cions      ;
  bool xions      ;
  bool yions      ;
  bool yMinusH2O  ;
  bool yMinusNH3  ;
  bool zions      ;
  bool MMinusH2O  ;
  bool MMinusNH3  ;

  SRMParameters()
  {
      // default values
      aions      = false;
      aMinusNH3  = false;
      bions      = true;
      bMinusH2O  = false;
      bMinusNH3  = false;
      bPlusH2O   = false;
      cions      = false;
      xions      = false;
      yions      = true;
      yMinusH2O  = false;
      yMinusNH3  = false;
      zions      = false;
      MMinusH2O  = false;
      MMinusNH3  = false;
  }

};

/* 
 * Function to calculate the b and y fragment series given a peptide sequence 
 * Note that the y series is calculated "backwards" starting from left to right
 * instead of right to left. Whereas b_series[0] holds the b-1 ion, y_series[0]
 * holds the y-(n-1) ion and y_series[n-1] holds the y-1 ion.
 */

int _calculate_clashes_other_series_sub(const char* sequence, double* tmp, 
        double* series, double ch, 
        bool aions     , bool aMinusNH3 , bool bions     , bool bMinusH2O ,
        bool bMinusNH3 , bool bPlusH2O  , bool cions     , bool xions     ,
        bool yions     , bool yMinusH2O , bool yMinusNH3 , bool zions     ,
        bool MMinusH2O , bool MMinusNH3 ,
        int isotope_mod = NOISOTOPEMODIFICATION
        ) ;

int calculate_fragment_masses(const char* sequence, double* tmp, 
        double* series, double ch, SRMParameters& params, int isotope_mod);

int _calculate_clashes_other_series(const char* sequence, double* tmp, 
        double* series, double ch, python::object parameters) {

    bool aions      =  python::extract<bool>(parameters.attr("aions"));
    bool aMinusNH3  =  python::extract<bool>(parameters.attr("aMinusNH3"));
    bool bions      =  python::extract<bool>(parameters.attr("bions"));
    bool bMinusH2O  =  python::extract<bool>(parameters.attr("bMinusH2O"));
    bool bMinusNH3  =  python::extract<bool>(parameters.attr("bMinusNH3"));
    bool bPlusH2O   =  python::extract<bool>(parameters.attr("bPlusH2O"));
    bool cions      =  python::extract<bool>(parameters.attr("cions"));
    bool xions      =  python::extract<bool>(parameters.attr("xions"));
    bool yions      =  python::extract<bool>(parameters.attr("yions"));
    bool yMinusH2O  =  python::extract<bool>(parameters.attr("yMinusH2O"));
    bool yMinusNH3  =  python::extract<bool>(parameters.attr("yMinusNH3"));
    bool zions      =  python::extract<bool>(parameters.attr("zions"));
    bool MMinusH2O  =  python::extract<bool>(parameters.attr("MMinusH2O"));
    bool MMinusNH3  =  python::extract<bool>(parameters.attr("MMinusNH3"));

return _calculate_clashes_other_series_sub(sequence, tmp, 
        series, ch, aions, aMinusNH3, bions, bMinusH2O,
        bMinusNH3, bPlusH2O, cions, xions, yions, yMinusH2O,
        yMinusNH3, zions, MMinusH2O, MMinusNH3);

}

int _calculate_clashes_other_series_sub(const char* sequence, double* tmp, 
        double* series, double ch, 
        bool aions     , bool aMinusNH3 , bool bions     , bool bMinusH2O ,
        bool bMinusNH3 , bool bPlusH2O  , bool cions     , bool xions     ,
        bool yions     , bool yMinusH2O , bool yMinusNH3 , bool zions     ,
        bool MMinusH2O , bool MMinusNH3 ,
        int isotope_mod )
{

    SRMParameters params;
    params.aions      =  aions;
    params.aMinusNH3  =  aMinusNH3;
    params.bions      =  bions;
    params.bMinusH2O  =  bMinusH2O;
    params.bMinusNH3  =  bMinusNH3;
    params.bPlusH2O   =  bPlusH2O;
    params.cions      =  cions;
    params.xions      =  xions;
    params.yions      =  yions;
    params.yMinusH2O  =  yMinusH2O;
    params.yMinusNH3  =  yMinusNH3;
    params.zions      =  zions;
    params.MMinusH2O  =  MMinusH2O;
    params.MMinusNH3  =  MMinusNH3;
    return calculate_fragment_masses(sequence, tmp, 
        series, ch, params, isotope_mod);
}

int calculate_fragment_masses(const char* sequence, double* tmp, 
        double* series, double ch, SRMParameters& params, int isotope_mod) 
{

    int j, start, scounter;
    double acc_mass, res_mass;
    char c;
    bool inside;
    int frg_cnt = 0;

    acc_mass = 0.0;
    res_mass = 0.0;
    scounter = 0;

    inside = false;
    start = 0;
    j = 0; 
    if(isotope_mod == NOISOTOPEMODIFICATION) {
    //go through all characters in the sequence until 0 is hit
    while((c = sequence[j++])) {
        if(sequence[j] == '[') {
            start = j-1;
            inside = true;
        }
        else if(sequence[j-1] == ']') {
            //We found a modification
            switch(sequence[start]) {
                case 'M': 
                    if(!(sequence[start+2] == '1' && 
                         sequence[start+3] == '4' && 
                         sequence[start+4] == '7' )) {
                        PyErr_SetString(PyExc_ValueError, 
                            "Unknown modification for methionine");
                        boost::python::throw_error_already_set();
                        return -1;
                        }
                    res_mass = 147.03540462; break;
                case 'C': 
                    if(!(sequence[start+2] == '1' && 
                         sequence[start+3] == '6' && 
                         sequence[start+4] == '0' )) {
                        PyErr_SetString(PyExc_ValueError, 
                            "Unknown modification for cysteine");
                        boost::python::throw_error_already_set();
                        return -1;
                    }
                    res_mass = 160.030653721; break;
                case 'N': 
                    if(!(sequence[start+2] == '1' && 
                         sequence[start+3] == '1' && 
                         sequence[start+4] == '5' )) {
                        PyErr_SetString(PyExc_ValueError, 
                            "Unknown modification for asparagine");
                        boost::python::throw_error_already_set();
                        return -1;
                    }
                    res_mass = 115.026945583; break;

                default: 
                    PyErr_SetString(PyExc_ValueError, 
                        "Unknown modification ");
                    boost::python::throw_error_already_set();
                    return -1;
            }
            //'M[147]':  131.04049 + mass_O), # oxygen
            //'C[160]':  103.00919 + mass_CAM - mass_H ), # CAM replaces H
            //'N[115]':  114.04293 - mass_N - mass_H + mass_O

            acc_mass += res_mass;
            tmp[scounter] = acc_mass;
            scounter++;

            inside = false;
        }
        else if(inside) { }
        else {
            //We found a regular AA
            switch(c) {
                case 'A': res_mass = 71.03711; break;
                case 'C': res_mass = 103.00919; break;
                case 'D': res_mass = 115.02694; break;
                case 'E': res_mass = 129.04259; break;
                case 'F': res_mass = 147.06841; break;
                case 'G': res_mass = 57.02146; break;
                case 'H': res_mass = 137.05891; break;
                case 'I': res_mass = 113.08406; break;
                case 'K': res_mass = 128.09496; break;
                case 'L': res_mass = 113.08406; break;
                case 'M': res_mass = 131.04049; break;
                case 'N': res_mass = 114.04293; break;
                case 'P': res_mass = 97.05276; break;
                case 'Q': res_mass = 128.05858; break;
                case 'R': res_mass = 156.10111; break;
                case 'S': res_mass = 87.03203; break;
                case 'T': res_mass = 101.04768; break;
                case 'V': res_mass = 99.06841; break;
                case 'W': res_mass = 186.07931; break;
                case 'X': res_mass = 113.08406; break;
                case 'Y': res_mass = 163.06333; break;
                default: 
                    PyErr_SetString(PyExc_ValueError, 
                        "Unknown amino acid ");
                    boost::python::throw_error_already_set();
                    return -1;
            }

            acc_mass += res_mass;
            tmp[scounter] = acc_mass;
            scounter++;
        }
    }

    } else if(isotope_mod == N15_ISOTOPEMODIFICATION ) {
    // using values for N15 isotopic modification 
    while((c = sequence[j++])) {
        if(sequence[j] == '[') {
            start = j-1;
            inside = true;
        }
        else if(sequence[j-1] == ']') {
            //We found a modification
            switch(sequence[start]) {
                case 'M': 
                    if(!(sequence[start+2] == '1' && 
                         sequence[start+3] == '4' && 
                         sequence[start+4] == '7' )) {
                        PyErr_SetString(PyExc_ValueError, 
                            "Unknown modification for methionine");
                        boost::python::throw_error_already_set();
                        return -1;
                        }
                    res_mass = 148.0324344260; break;
                case 'C': 
                    if(!(sequence[start+2] == '1' && 
                         sequence[start+3] == '6' && 
                         sequence[start+4] == '0' )) {
                        PyErr_SetString(PyExc_ValueError, 
                            "Unknown modification for cysteine");
                        boost::python::throw_error_already_set();
                        return -1;
                    }
                    res_mass = 161.027683399; break;
                case 'N': 
                    if(!(sequence[start+2] == '1' && 
                         sequence[start+3] == '1' && 
                         sequence[start+4] == '5' )) {
                        PyErr_SetString(PyExc_ValueError, 
                            "Unknown modification for asparagine");
                        boost::python::throw_error_already_set();
                        return -1;
                    }
                    res_mass = 116.023977918; break;

                default: 
                    PyErr_SetString(PyExc_ValueError, 
                        "Unknown modification ");
                    boost::python::throw_error_already_set();
                    return -1;
            }
            //'M[147]':  131.04049 + mass_O), # oxygen
            //'C[160]':  103.00919 + mass_CAM - mass_H ), # CAM replaces H
            //'N[115]':  114.04293 - mass_N - mass_H + mass_O

            acc_mass += res_mass;
            tmp[scounter] = acc_mass;
            scounter++;

            inside = false;
        }
        else if(inside) { }
        else {
            //We found a regular AA
            switch(c) {
                case 'A': res_mass = 72.0341486780; break;
                case 'C': res_mass = 104.006219678; break;
                case 'D': res_mass = 116.023977918; break;
                case 'E': res_mass = 130.039627982; break;
                case 'F': res_mass = 148.065448806; break;
                case 'G': res_mass = 58.018498614; break;
                case 'H': res_mass = 140.050016538; break;
                case 'I': res_mass = 114.08109887; break;
                case 'K': res_mass = 130.0890328; break;
                case 'L': res_mass = 114.08109887; break;
                case 'M': res_mass = 132.037519806; break;
                case 'N': res_mass = 116.036997228; break;
                case 'P': res_mass = 98.049798742; break;
                case 'Q': res_mass = 130.052647292; break;
                case 'R': res_mass = 160.089250596; break;
                case 'S': res_mass = 88.029063298; break;
                case 'T': res_mass = 102.044713362; break;
                case 'V': res_mass = 100.065448806; break;
                case 'W': res_mass = 188.073382736; break;
                case 'X': res_mass = 114.08109887; break;
                case 'Y': res_mass = 164.060363426; break;

                default: 
                    PyErr_SetString(PyExc_ValueError, 
                        "Unknown amino acid ");
                    boost::python::throw_error_already_set();
                    return -1;
            }

            acc_mass += res_mass;
            tmp[scounter] = acc_mass;
            scounter++;
        }
    } 
    }


    // see also http://www.matrixscience.com/help/fragmentation_help.html
    // note that the b and y series only go up to y[n-1] and b[n-1] since
    // the last ion of the series would be the parent ion (y) or the parent
    // ion with a loss of water (b).
    //
    // One could add the last ion of the series (a,b,c = remove the -1)
    //
    if (params.yions) 
        // series[frg_cnt++] = acc_mass + 2*MASS_H + MASS_OH;
        for (int j=0; j<scounter-1; j++) series[frg_cnt++] = acc_mass - tmp[j] + 2*MASS_H + MASS_OH;
    if (params.bions) for (int j=0; j<scounter-1; j++) series[frg_cnt++] = tmp[j] + MASS_H ;

    if (params.aions) for (int j=0; j<scounter-1; j++) series[frg_cnt++] = tmp[j] + MASS_H - MASS_CO ;
    if (params.cions) for (int j=0; j<scounter-1; j++) series[frg_cnt++] = tmp[j] + MASS_H + MASS_NH3 ;
    if (params.xions) for (int j=0; j<scounter-1; j++) series[frg_cnt++] = acc_mass - tmp[j] + MASS_CO + MASS_OH;

    if (params.zions) 
        // series[frg_cnt++] = acc_mass + 2*MASS_H + MASS_OH - MASS_NH3;
        for (int j=0; j<scounter-1; j++) series[frg_cnt++] = acc_mass - tmp[j] + 2*MASS_H + MASS_OH - MASS_NH3; 

    if (params.aMinusNH3) for (int j=0; j<scounter-1; j++) series[frg_cnt++] = tmp[j] + MASS_H - MASS_CO - MASS_NH3; 
    if (params.bMinusH2O) for (int j=0; j<scounter-1; j++) series[frg_cnt++] = tmp[j] + MASS_H - MASS_H2O;
    if (params.bMinusNH3) for (int j=0; j<scounter-1; j++) series[frg_cnt++] = tmp[j] + MASS_H - MASS_NH3;
    if (params.bPlusH2O)  for (int j=0; j<scounter-1; j++) series[frg_cnt++] = tmp[j] + MASS_H + MASS_H2O;

    if (params.yMinusH2O) {
        // series[frg_cnt++] = acc_mass + MASS_H;
        for (int j=0; j<scounter-1; j++) 
            series[frg_cnt++] = acc_mass - tmp[j] + MASS_H; }
    if (params.yMinusNH3) {
        //series[frg_cnt++] = acc_mass + 2*MASS_H + MASS_OH - MASS_NH3;
        for (int j=0; j<scounter-1; j++) 
            series[frg_cnt++] = acc_mass - tmp[j] + 2*MASS_H + MASS_OH - MASS_NH3;
    }

    /*
    #these are all the more exotic ions
    self.b_series = [b + R.mass_H for b in fragment_series[:-1]] 
    self.y_series = [self.mass - y + 2*R.mass_H + R.mass_OH for y in fragment_series[:-1]]
    if aions: self.a_series = [f - R.mass_CO for f in self.b_series]
    if cions: self.c_series = [f + R.mass_NH3 for f in self.b_series]
    if xions: self.x_series = [f + R.mass_CO - 2*R.mass_H for f in self.y_series]
    if zions: self.z_series = [f - R.mass_NH3 for f in self.y_series]
    if aMinusNH3: self.a_minus_NH3 = [f - R.mass_CO - R.mass_NH3 for f in self.b_series]
    if bMinusH2O: self.b_minus_H2O = [b - R.mass_H2O for b in self.b_series]
    if bMinusNH3: self.b_minus_NH3 = [b - R.mass_NH3 for b in self.b_series]
    if bPlusH2O:  self.b_plus_H2O  = [b + R.mass_H2O for b in self.b_series]
    if yMinusH2O: self.y_minus_H2O = [y - R.mass_H2O for y in self.y_series]
    if yMinusNH3: self.y_minus_NH3 = [y - R.mass_NH3 for y in self.y_series]
    if MMinusH2O: self.waterloss = self.mass; self.allseries.append(self.waterloss)
    if MMinusNH3: self.nh3loss =   self.mass + R.mass_H2O - R.mass_NH3; self.allseries.append(self.nh3loss)
    */

    if (params.MMinusH2O) series[frg_cnt++] = acc_mass;
    if (params.MMinusNH3) series[frg_cnt++] = acc_mass + MASS_H2O - MASS_NH3;

    // calculate the charged mass of the fragments
    for (int j=0; j<frg_cnt; j++) series[j] = (series[j] + (ch-1)*MASS_H)/ch;
    return frg_cnt;
}

void calculate_transitions_with_charge(SRMPrecursor& p, std::vector<int>& charges, std::vector<SRMTransition>& result, 
  double* tmp, double* series, double & q3_low, double & q3_high, SRMParameters& param)
{

  int fragcount, k;
  double q3;

  for (std::vector<int>::iterator ch_it = charges.begin(); ch_it != charges.end(); ch_it++) {
    //fragcount = _calculate_clashes(p.sequence, b_series, y_series, (*ch_it) );
    fragcount = calculate_fragment_masses(p.sequence, tmp, series, (*ch_it), param, NOISOTOPEMODIFICATION);

    // go through all fragments of this precursor
    for (k=0; k<fragcount; k++) {
      q3 = series[k];
      if (q3 > q3_low && q3 < q3_high)
      {
        SRMTransition t = {p.q1, q3, p.transition_group};
        result.push_back(t);
      }
    }

  }

}

/* 
 * Function to calculate all transitions of a list of precursor peptides and
 * allows to select the charge states of these precursors.
 * Precursors are tuples of the form (q1, sequence, peptide_key).
 * It will return a list of tuples that are of the form 
 * (q3, q1, 0, peptide_key) 
 *
 * Input
 * precursors: (q1, sequence, peptide_key)
 *
 */
python::list _py_calculate_transitions_with_charge(python::tuple py_precursors,
        python::list py_charges, double q3_low, double q3_high ) {

    python::tuple clist;
    python::list py_result;

    std::vector<SRMTransition> result;
    std::vector<SRMPrecursor> precursors;
    std::vector<int> charges;

    SRMParameters param;

    long peptide_key;
    double q1;
    char* sequence;

    double* tmp = new double[1024];
    double* series = new double[1024];

    for (int kk=0; kk<python::extract<int>(py_charges.attr("__len__")()); kk++) {
      charges.push_back(python::extract< int >(py_charges[kk]));
    }
    for (int i=0; i<python::extract<int>(py_precursors.attr("__len__")()); i++) {
      SRMPrecursor p; 
      clist = python::extract< python::tuple >(py_precursors[i]);
      p.q1 = python::extract< double >(clist[0]);
      p.sequence = python::extract<char *>(clist[1]);
      p.transition_group = python::extract< long >(clist[2]);
      precursors.push_back(p);
    }

    for (int i=0; i< precursors.size(); i++) {
      calculate_transitions_with_charge(precursors[i], charges, result, tmp, series, q3_low, q3_high, param);
    }

    for (std::vector<SRMTransition>::iterator it = result.begin(); it != result.end(); it++)
    {
      py_result.append(python::make_tuple(it->q3, it->q1, 0, it->transition_id));
    }

    delete [] tmp;
    delete [] series;

    return py_result;
}        

/*
 * Calculate the charged mass for a sequence, supplied in a python tuple in the
 * first place.
*/
double calculate_charged_mass(python::tuple clist, int ch) {

    double charged_mass;
    double* tmp = new double[256];
    double* series = new double[256];

    char* sequence = python::extract<char *>(clist[1]);
    SRMParameters params;
    params.bions = true;
    params.yions = false;
    int fragcount = calculate_fragment_masses(sequence, tmp, series, ch, params, NOISOTOPEMODIFICATION);
 
    //In order to get the full mass, we need the "last" element of the b-series
    //(which is not technically part of the b series) and add water as well as
    //protons according to the charge to it. Then, the fragment has to be
    //divided by the charge.
    charged_mass = (tmp[fragcount] + MASS_H + MASS_H*ch + MASS_OH ) /ch  ;
 
    delete [] tmp;
    delete [] series;

    return charged_mass;
}        

/*
* This function calculates all combinations of M elements
* drawn without replacement from a set of N elements. Order of elements
* does NOT matter.
*
* M is the order
* N is the length of the input vector
* mapping is an array of integers of fixed length (see srmcollider.h)
* result is a python dictionary that holds all combinations, they are
*   implemented as integers with bitflags set or unset
*/
void _combinations_magic(int M, int N, COMBINT* mapping,
        python::dict &result) {
    // The basic idea is to create an index array of length M that contains
    // numbers between 0 and N. The indices denote the combination produced and
    // we always increment the rightmost index. If it goes above N, we try to
    // increase the one left to it until we find one that still can be
    // increased. We stop when the rightmost index hits N-M, we thus go from
    // (0,1,2,3...M-1) to (N-M,N-M+1,N-M+2...N-1)
    //
    // see http://svn.python.org/projects/python/trunk/Modules/itertoolsmodule.c
    // combinations_next fxn around line 2113

    int j, k;
    int* index = new int[M];
    COMBINT tmpres;

    python::tuple tmptuple;
    //initialize with numbers from 0 to M = range( M )
    for(int k=0;k<M;k++) index[k] = k;
    while (index[0] <= N-M) {

        //EVALUATE THE RESULT
        tmpres = 0;
        for(k=0;k<M;k++) 
            tmpres |= mapping[index[k]];
        result[tmpres] = 0;

        // We need to break if index[0] has the final value
        // The other option is to make the while condition (index[0] < N-M) 
        // and add an additional evaluation of the result to the end of the 
        // function (that would probably be faster).
        if(index[0] == N-M) break;

        index[ M-1 ] += 1;
        if (index[ M-1 ] >= N) {
            //#now we hit the end, need to increment other positions than last
            //#the last position may reach N-1, the second last only N-2 etc.
            j = M-1;
            while (j >= 0 and index[j] >= N-M+j) j -= 1;
            //#j contains the value of the index that needs to be incremented
            //TODO if(j<0) break;
            index[j] += 1;
            //now start from j+1 to the right and set all indices to their
            //lowest possible values
            k = j + 1;
            while (k < M) { index[k] = index[k-1] + 1; k += 1; } 
        }
    }

    delete [] index;
}

/*
* This function calculates all combinations of M elements
* drawn without replacement from a set of N elements. Order of elements
* does NOT matter.
*
* M is the order
* N is the length of the input vector
* mapping is a python list containing any object
* result is a python dictionary that holds all combinations, as tuples of size
*   M using elements that are found in the mapping pool
*/
void _combinations(int M, int N, const python::list &mapping, 
        python::dict &result) {
    // The basic idea is to create an index array of length M that contains
    // numbers between 0 and N. The indices denote the combination produced and
    // we always increment the rightmost index. If it goes above N, we try to
    // increase the one left to it until we find one that still can be
    // increased. We stop when the rightmost index hits N-M, we thus go from
    // (0,1,2,3...M-1) to (N-M,N-M+1,N-M+2...N-1)
    int j, k;
    int* index = new int[M];
    //int* mytuple = new int[M];
    python::tuple tmptuple;
    //initialize with numbers from 0 to M = range( M )
    for(int k=0;k<M;k++) index[k] = k;
    while (index[0] <= N-M) {

        //EVALUATE THE RESULT
        //only works for M up to order 5
        switch(M) {
            case 1: tmptuple = python::make_tuple( mapping[index[0]] ); break;
            case 2: tmptuple = python::make_tuple( mapping[index[0]], mapping[index[1]]); break;
            case 3: tmptuple = python::make_tuple( mapping[index[0]], mapping[index[1]], 
                            mapping[index[2]]); break;
            case 4: tmptuple = python::make_tuple( mapping[index[0]], mapping[index[1]], 
                            mapping[index[2]], mapping[index[3]]); break;
            case 5: tmptuple = python::make_tuple( mapping[index[0]], mapping[index[1]], 
                            mapping[index[2]], mapping[index[3]], mapping[index[4]]); break;
            case 6: tmptuple = python::make_tuple( mapping[index[0]],
                            mapping[index[1]], 
                            mapping[index[2]],
                            mapping[index[3]],
                            mapping[index[4]],
                            mapping[index[5]]
                    ); break;
            case 7: tmptuple = python::make_tuple( mapping[index[0]],
                            mapping[index[1]], 
                            mapping[index[2]],
                            mapping[index[3]],
                            mapping[index[4]],
                            mapping[index[5]],
                            mapping[index[6]]
                    ); break;
            case 8: tmptuple = python::make_tuple( mapping[index[0]],
                            mapping[index[1]], 
                            mapping[index[2]],
                            mapping[index[3]],
                            mapping[index[4]],
                            mapping[index[5]],
                            mapping[index[6]],
                            mapping[index[7]]
                    ); break;
            case 9: tmptuple = python::make_tuple( mapping[index[0]],
                            mapping[index[1]], 
                            mapping[index[2]],
                            mapping[index[3]],
                            mapping[index[4]],
                            mapping[index[5]],
                            mapping[index[6]],
                            mapping[index[7]],
                            mapping[index[8]]
                    ); break;
            case 10: tmptuple = python::make_tuple( mapping[index[0]],
                            mapping[index[1]], 
                            mapping[index[2]],
                            mapping[index[3]],
                            mapping[index[4]],
                            mapping[index[5]],
                            mapping[index[6]],
                            mapping[index[7]],
                            mapping[index[8]],
                            mapping[index[9]]
                    ); break;
            default:
                PyErr_SetString(PyExc_ValueError, 
                    "Order (M) larger than 5 is not implemented");
                boost::python::throw_error_already_set();
                return;
        }
        result[ tmptuple ] = 0; //append to result dict

        // We need to break if index[0] has the final value
        // The other option is to make the while condition (index[0] < N-M) 
        // and add an additional evaluation of the result to the end of the 
        // function (that would probably be faster).
        if(index[0] == N-M) break;

        index[ M-1 ] += 1;
        if (index[ M-1 ] >= N) {
            //#now we hit the end, need to increment other positions than last
            //#the last position may reach N-1, the second last only N-2 etc.
            j = M-1;
            while (j >= 0 and index[j] >= N-M+j) j -= 1;
            //#j contains the value of the index that needs to be incremented
            index[j] += 1;
            //TODO if(j<0) break;
            k = j + 1;
            while (k < M) { index[k] = index[k-1] + 1; k += 1; } 
        }
    }

    delete [] index;
}

/*
 * Function to calculate all non-UIS transitions from a dictionary
 * where for each key (colliding peptide key) the set of transitions
 * of the query peptide are stored that are interfered with is held
 * (collisions_per_peptide dictionary)
 * Note that this dictionary contains tuples that represent the combinations.
 *
 * It will return a list of all non UIS of the requested order.
 */
python::dict get_non_uis(python::dict collisions_per_peptide, int order) {

    python::list tmplist;
    python::list pepcollisions = python::extract< python::list >(
            collisions_per_peptide.values() );
    python::dict result;

    int tmp_length;
    int collision_length = python::extract<int>(pepcollisions.attr("__len__")());

    for (int i=0; i<collision_length; i++) {
        tmplist = python::extract< python::list >( pepcollisions[i] );
        tmp_length = python::extract<int>(tmplist.attr("__len__")());
        _combinations(order, tmp_length, tmplist, result);
    }

    return result;

    /* This above get_non_uis and _combinations functions 
     * correspond to the following python function (roughly)
     * 

        for pepc in collisions_per_peptide.values():
            for i in range(1,MAX_UIS+1):
                if len( pepc ) >= i: 
                    tmp = [ [pepc[j] for j in indices] for indices in 
                        _combinations( len(pepc) , i)] 
                    non_uis.update( [tuple(sorted(p)) for p in combinations(pepc, i)] )
        return non_uis_list

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

    */

}

/*
 * Function to calculate all non-UIS transitions from a dictionary
 * where for each key (colliding peptide key) the set of transitions
 * of the query peptide are stored that are interfered with is held
 * (collisions_per_peptide vector).
 * Note that this vector contains integers that represent the combinations.
 *
 * It will return a list of all non UIS of the requested order.
 */
python::dict get_non_uis_magic(std::vector<COMBINT>& newcollperpep, int max_tr, int
        order) {

    python::dict result;

    int onecounter;
    COMBINT tmparr;
    COMBINT mask;
    COMBINT* mapping = new COMBINT[max_tr];

    for (uint i=0; i<newcollperpep.size(); i++) {

        //count the number of binary ones in the bitarray
        mask = 1;
        onecounter = 0;
        tmparr = newcollperpep[i];

        //we know that there cannot be bits populated past max nr transitions
        for(int j=0; j<max_tr; j++) {
            //true if nonzero
            if(tmparr & mask) 
                mapping[onecounter++] = mask;
            mask <<=1;
        }

        /* The other way to do it
         * we shift the tmparr to the right
        for(int j=0; tmparr != 0; tmparr >>= 1) {
            onecounter += tmparr & mask;
        }
        */

        _combinations_magic(order, onecounter, mapping, result);
    }

    delete [] mapping;
    return result;
}




int __old_calculate_clashes(const char* sequence, double* b_series, double* y_series,
        double ch) {

    int j, start, scounter;
    double acc_mass, res_mass;
    char c;
    bool inside;

    acc_mass = 0.0;
    res_mass = 0.0;
    scounter = 0;

    inside = false;
    start = 0;
    j = 0; 
    //go through all characters in the sequence until 0 is hit
    while((c = sequence[j++])) {
        if(sequence[j] == '[') {
            start = j-1;
            inside = true;
        }
        else if(sequence[j-1] == ']') {
            //We found a modification
            switch(sequence[start]) {
                case 'M': 
                    if(!(sequence[start+2] == '1' && 
                         sequence[start+3] == '4' && 
                         sequence[start+4] == '7' )) {
                        PyErr_SetString(PyExc_ValueError, 
                            "Unknown modification for methionine");
                        boost::python::throw_error_already_set();
                        return -1;
                        }
                    res_mass = 147.03540462; break;
                case 'C': 
                    if(!(sequence[start+2] == '1' && 
                         sequence[start+3] == '6' && 
                         sequence[start+4] == '0' )) {
                        PyErr_SetString(PyExc_ValueError, 
                            "Unknown modification for cysteine");
                        boost::python::throw_error_already_set();
                        return -1;
                    }
                    res_mass = 160.030653721; break;
                case 'N': 
                    if(!(sequence[start+2] == '1' && 
                         sequence[start+3] == '1' && 
                         sequence[start+4] == '5' )) {
                        PyErr_SetString(PyExc_ValueError, 
                            "Unknown modification for cysteine");
                        boost::python::throw_error_already_set();
                        return -1;
                    }
                    res_mass = 115.026945583; break;

                default: 
                    PyErr_SetString(PyExc_ValueError, 
                        "Unknown modification ");
                    boost::python::throw_error_already_set();
                    return -1;
            }
            //'M[147]':  131.04049 + mass_O), # oxygen
            //'C[160]':  103.00919 + mass_CAM - mass_H ), # CAM replaces H
            //'N[115]':  114.04293 - mass_N - mass_H + mass_O

            acc_mass += res_mass;
            b_series[scounter] = acc_mass + MASS_H ;
            y_series[scounter] = - acc_mass + 2* MASS_H + MASS_OH ;
            scounter++;

            inside = false;
        }
        else if(inside) { }
        else {
            //We found a regular AA
            //TODO use hash map http://en.wikipedia.org/wiki/Hash_map_%28C%2B%2B%29
            switch(c) {
                case 'A': res_mass = 71.03711; break;
                case 'C': res_mass = 103.00919; break;
                case 'D': res_mass = 115.02694; break;
                case 'E': res_mass = 129.04259; break;
                case 'F': res_mass = 147.06841; break;
                case 'G': res_mass = 57.02146; break;
                case 'H': res_mass = 137.05891; break;
                case 'I': res_mass = 113.08406; break;
                case 'K': res_mass = 128.09496; break;
                case 'L': res_mass = 113.08406; break;
                case 'M': res_mass = 131.04049; break;
                case 'N': res_mass = 114.04293; break;
                case 'P': res_mass = 97.05276; break;
                case 'Q': res_mass = 128.05858; break;
                case 'R': res_mass = 156.10111; break;
                case 'S': res_mass = 87.03203; break;
                case 'T': res_mass = 101.04768; break;
                case 'V': res_mass = 99.06841; break;
                case 'W': res_mass = 186.07931; break;
                case 'X': res_mass = 113.08406; break;
                case 'Y': res_mass = 163.06333; break;
                default: 
                    PyErr_SetString(PyExc_ValueError, 
                        "Unknown amino acid ");
                    boost::python::throw_error_already_set();
                    return -1;
            }

            acc_mass += res_mass;
            b_series[scounter] = acc_mass + MASS_H ;
            y_series[scounter] = - acc_mass + 2* MASS_H + MASS_OH ;
            scounter++;
        }
    }

    // We computed the y series "backwards" using negative values and thus need
    // to add the total mass to all entries
    for (int j=0; j<scounter; j++) y_series[j] += acc_mass;

    for (int j=0; j<scounter-1; j++){

        y_series[j] = (y_series[j] + (ch-1)*MASS_H)/ch;
        b_series[j] = (b_series[j] + (ch-1)*MASS_H)/ch;

    }
    // Subtract one because we do not count the last mass as a fragment.
    // The complete peptide has not undergone fragmentation by definition and
    // is thus not a fragment. We still compute it but just dont consider it.
    return --scounter;
}

#endif
