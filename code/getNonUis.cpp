/*
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
 */

/*
 * The functions in this file allow to calculate the collisions_per_peptide
 * dictionary easily, either using precursors (e.g. peptide sequences) or
 * collisions (q1,q3 tuples) as input. The collisions_per_peptide dictionary
 * holds for each peptide in the background the exact transitions of the query
 * peptides that are shared.
 *
 * Furthermore, we provide interfaces for some of the shared library functions.
*/

#include <vector>
//include our own libraries
#include "srmcollider.h"
#include "srmcolliderLib.cpp"
#include "py_srmcolliderLib.h"
#include "combinatorics.h"

// Boost.Python headers
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
namespace python = boost::python;

#include "calculate_eUIS.cpp"

//using namespace std;

namespace SRMCollider
{
  namespace getNonUIS
  {
    // Function declarations

    // checks whether the current fragment has an allowed charge
    bool has_allowed_charge(int fragment_charge, int q1_charge, int maximal_charge)
    {

      // disallow doubly charged precursors and double charged fragments
      if(fragment_charge == 2 && q1_charge == 2)
      {return false;}
      // disallow precursors that exceed the maximal charge of the peptide
      if(q1_charge > maximal_charge )
      {return false;}
      // disallow doubly charged fragments where the maximal charge of the fragment is 1 or 2
      if(fragment_charge == 2 && maximal_charge < 3)
      {return false;}

      // TODO / implement: check each 2+ fragment whether it can hold the charge

      return true;
    }

    ///////////////////////////////////////////////////////////////////////////

    /*
     * Function to calculate the collisions_per_peptide out of a set of transitions
     * and precursor peptides that are colliding peptides.  It will return a
     * dictionary where for each key (colliding peptide key) the set of transitions
     * of the query peptide are stored that are interfered with is held.
     * Transitions are tuples of the form (q3, srm_id), precursors are tuples of
     * the form (q1, sequence, peptide_key).
     *
     * Calculate collisions per peptide
     *
        // go through all (potential) interfering precursors and store the
        // colliding SRM ids in a dictionary (they can be found at position 3 and 1
        // respectively).
     */
    python::dict _find_clashes_calculate_collperpeptide_other_ion_series(
            python::tuple transitions, python::list precursors, python::object par, 
            double q3_low, double q3_high, double q3window, bool ppm, bool forceChargeCheck) 
    {

        python::dict collisions_per_peptide, tmpdict;
        python::tuple tlist;
        python::list tmplist;

        std::vector<SRMCollider::Common::SRMPrecursor> c_precursors;
        std::vector<SRMCollider::Common::SRMTransition> c_transitions;

        int transitions_length = python::extract<int>(transitions.attr("__len__")());
        int precursor_length = python::extract<int>(precursors.attr("__len__")());
        int fragcount, i, j, k, ch, listmembers = 0;

        long t1, peptide_key;
        double t0, q3used = q3window;

        double* series = new double[1024];
        double* tmp_series = new double[1024];

        SRMCollider::Common::SRMParameters params;
        pyToC::initialize_param_obj(par, params);
        pyToC::initialize_precursors(precursors, c_precursors);
        pyToC::initialize_transitions(transitions, c_transitions);

        // go through all (potential) interfering precursors and store the
        // colliding SRM ids in a dictionary (they can be found at position 3 and 1
        // respectively).
        for (j=0; j<precursor_length; j++) {
            SRMCollider::Common::SRMPrecursor & precursor = c_precursors[j];

            for (ch=1; ch<=2; ch++) {
              try{
                fragcount = calculate_fragment_masses(precursor.sequence, tmp_series, series, ch,
                      params, precursor.isotope_modification);
              }
              catch (SRMCollider::Common::AANotFound& error)
              {
                PyErr_SetString(PyExc_ValueError, error.message.c_str());
                boost::python::throw_error_already_set();
                delete [] series;
                delete [] tmp_series;
                return collisions_per_peptide;
              }

                if(forceChargeCheck && !has_allowed_charge(ch, precursor.q1_charge, precursor.maximal_charge) )
                {continue;}

                for (i=0; i<transitions_length; i++) {
                    //ppm is 10^-6
                    t0 = c_transitions[i].q3;
                    if(ppm) {q3used = q3window / 1000000.0 * t0; } 

                        // go through all fragments of this precursor
                        for (k=0; k<fragcount; k++) {

                            if(fabs(t0-series[k]) < q3used ) {
                                // extract SRM_id from transition list and store it
                                // as dirty in a temporary dict
                                t1 = c_transitions[i].transition_id;
                                tmpdict[t1] = 0;
                                listmembers++; 
                            
                            }
                        }
                    } //loop over all transitions
                }

            //we keep one empty list around since we hope that most precursors dont 
            //use it. If it gets used, we store it in the result and create a new 
            //list to use from then on.
            if (listmembers>0) {
                peptide_key = c_precursors[j].transition_group;
                tmplist = tmpdict.keys();
                tmplist.sort();
                collisions_per_peptide[peptide_key] = tmplist;
                python::dict newlist;
                tmpdict = newlist;
            }
            listmembers = 0;

        } //end of loop over all precursors

        delete [] series;
        delete [] tmp_series;
        return collisions_per_peptide;
    }

    /*
     * Function to calculate the collision density out of a set of transitions
     * and precursor peptides that are colliding peptides.
     *
     * It will return an array where the number of interferences is recorded for
     * each transition.  Transitions are tuples of the form (q3, srm_id),
     * precursors are tuples of the form (q1, sequence, peptide_key).
     */
    python::list _find_clashes_calculate_colldensity(python::tuple transitions,
        python::list precursors, double q3_low, double q3_high, double q3window,
        bool ppm) 
    {

        python::dict collisions_per_peptide, tmpdict;
        python::tuple clist;
        python::tuple tlist;

        python::list tmplist;

        int transitions_length = python::extract<int>(transitions.attr("__len__")());
        int precursor_length = python::extract<int>(precursors.attr("__len__")());
        int fragcount, i, j, k, ch;

        SRMCollider::Common::SRMParameters param;

        double q3used = q3window;
        std::string sequence;

        double* tmp = new double[256];
        double* series = new double[256];

        int* cresult = new int[transitions_length];
        for (i=0; i<transitions_length; i++) cresult[i] = 0;

        double* tmptrans = new double[transitions_length];
        double* tmpq3used = new double[transitions_length];
        for (i=0; i<transitions_length; i++) {
            tlist = python::extract< python::tuple >(transitions[i]);
            //ppm is 10^-6
            tmptrans[i] = python::extract< double >(tlist[0]);
            if(ppm) {q3used = q3window / 1000000.0 * tmptrans[i]; } 
            tmpq3used[i] =q3used;
        }

        // go through all (potential) collisions
        // and store the colliding SRM ids in a dictionary (they can be found at
        // position 3 and 1 respectively)
        for (j=0; j<precursor_length; j++) {
            clist = python::extract< python::tuple >(precursors[j]);
            sequence = python::extract<std::string>(clist[1]);

            for (ch=1; ch<=2; ch++) {
                //fragcount = _calculate_clashes(sequence, b_series, y_series, ch);
              try
              {
                fragcount = calculate_fragment_masses(sequence, tmp, series, ch, param, NOISOTOPEMODIFICATION);
              }
              catch (SRMCollider::Common::AANotFound& error)
              {
                PyErr_SetString(PyExc_ValueError, error.message.c_str());
                boost::python::throw_error_already_set();
                delete [] series;
                delete [] tmp;
                python::list result;
                return result;
              }

                for (i=0; i<transitions_length; i++) {

                        // go through all fragments of this precursor
                        for (k=0; k<fragcount; k++) {
                            if(fabs(tmptrans[i]-series[k]) < tmpq3used[i]) cresult[i]++;
                            //if(fabs(tmptrans[i]-b_series[k]) < tmpq3used[i]) cresult[i]++; 
                        }
                    } //loop over all transitions
                }
        } //end of loop over all precursors

        delete [] tmp;
        delete [] series;

        // TODO replace with std::vector 
        // http://www.boost.org/doc/libs/1_41_0/libs/python/doc/v2/indexing.html
        // http://www.cplusplus.com/reference/stl/vector/reserve/
        //
        python::list result;
        for (i=0; i<transitions_length; i++) result.append( cresult[i] ) ;

        delete [] cresult ;
        delete [] tmptrans;
        delete [] tmpq3used;
        return result;
    }

    // we annotate the ion nr k that was produced by a call to calculate_fragment_masses, 
    // the annotated ion is of type "curr_ion"-l (e.g. y7 means curr_ion = "y" and l = 7). 
    void annotate_ion( int& l, int k, const std::string sequence, std::string& curr_ion, SRMCollider::Common::SRMParameters& params) 
    {
        int scounter, icounter;
        double* tmp = new double[256];
        double* series = new double[256];
        bool done = false;

        // get the number of amino acids in the sequence
        SRMCollider::Common::SRMParameters tmp_params;
        tmp_params.bions = false;
        scounter = calculate_fragment_masses(sequence, tmp, series, 1, tmp_params, NOISOTOPEMODIFICATION);
        scounter++;

        // default is ? (unknown)
        done = false;
        icounter = 0;
        curr_ion = "?";

        // if we use the same order as in calculate_fragment_masses, we can infer the annotation
        if (params.yions && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "y"; done = true; break;}; icounter++;}
        if (params.bions && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "b"; done = true; break;}; icounter++;}
        if (params.aions && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "a"; done = true; break;}; icounter++;}
        if (params.cions && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "c"; done = true; break;}; icounter++;}
        if (params.xions && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "x"; done = true; break;}; icounter++;}
        if (params.zions && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "z"; done = true; break;}; icounter++;}

        if (params.aMinusNH3 && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "aMinusNH3"; done = true; break;}; icounter++;}
        if (params.bMinusH2O && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "bMinusH2O"; done = true; break;}; icounter++;}
        if (params.bMinusNH3 && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "bMinusNH3"; done = true; break;}; icounter++;}
        if (params.bPlusH2O && !done)  for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "bPlusH2O"; done = true; break;}; icounter++;}

        if (params.yMinusH2O && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "yMinusH2O"; done = true; break;}; icounter++;}
        if (params.yMinusNH3 && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "yMinusNH3"; done = true; break;}; icounter++;}

        if (params.MMinusH2O && !done) {l=0;if(icounter==k){curr_ion = "MMinusH2O"; done = true;}; icounter++;}
        if (params.MMinusNH3 && !done) {l=0;if(icounter==k){curr_ion = "MMinusNH3"; done = true;}; icounter++;}


        l++; // ion series starts at 1, thus add one

        // need to reverse the backwards series x,y,z
        if (curr_ion == "y" || curr_ion == "x" || curr_ion == "z"
            || curr_ion == "yMinusNH3" || curr_ion == "yMinusH2O")
        {
            l = scounter -l;
        }

        delete [] tmp;
        delete [] series;
    }

    /*
     * Function to calculate the exact interfering transitions for each peptide.
     * It will return a Transitions are tuples of the form (q3, srm_id), precursors
     * are tuples of the form (q1, sequence, peptide_key).
     *
     * Used by the web-scripts to report the exact interfering transitions.
     */
    python::dict _find_clashes_forall_other_series(python::tuple transitions,
        python::list precursors, python::object par, double q3_low, double q3_high, double q3window,
        bool ppm, double q1_low, bool forceChargeCheck) 
    {

      python::dict result, tmpdict;
      python::tuple tlist;
      python::list tmplist;

      std::vector<SRMCollider::Common::SRMPrecursor> c_precursors;
      std::vector<SRMCollider::Common::SRMTransition> c_transitions;

      int transitions_length = python::extract<int>(transitions.attr("__len__")());
      int precursor_length = python::extract<int>(precursors.attr("__len__")());
      int fragcount, i, j, k, ch;
      int isotope_nr;

      long t1;
      double t0, q3used = q3window, q1_used;
      std::string sequence;
      int max_isotopes = python::extract<int>(par.attr("isotopes_up_to"));

      double* series = new double[10*256];
      double* tmp_series = new double[256];

      SRMCollider::Common::SRMParameters params;
      pyToC::initialize_param_obj(par, params);
      pyToC::initialize_precursors(precursors, c_precursors);
      pyToC::initialize_transitions(transitions, c_transitions);

      std::string curr_ion = "?";

      // go through all (potential) collisions
      // and store the colliding SRM ids in a dictionary (they can be found at
      // position 3 and 1 respectively)
      for (j=0; j<precursor_length; j++) {
          SRMCollider::Common::SRMPrecursor & precursor = c_precursors[j];
          sequence = c_precursors[j].sequence;

        for (ch=1; ch<=2; ch++) 
        {
          try
          {
            fragcount = calculate_fragment_masses(sequence, tmp_series, series, ch,
              params, NOISOTOPEMODIFICATION);
          }
          catch (SRMCollider::Common::AANotFound& error)
          {
            PyErr_SetString(PyExc_ValueError, error.message.c_str());
            boost::python::throw_error_already_set();
            delete [] series;
            delete [] tmp_series;
            return result;
          }

          if(forceChargeCheck && !has_allowed_charge(ch, precursor.q1_charge, precursor.maximal_charge) )
          {continue;}

          for (i=0; i<transitions_length; i++) {
            //ppm is 10^-6
            t0 = c_transitions[i].q3;
            if(ppm) {q3used = q3window / 1000000.0 * t0; } 

            // go through all fragments of this precursor
            for (k=0; k<fragcount; k++) {
              if(fabs(t0-series[k]) < q3used ) {

                t1 = c_transitions[i].transition_id;
                if( result.has_key(t1) ) 
                {
                  tmplist = python::extract<python::list>(result[t1]);
                }
                else
                {
                  python::list newlist;
                  tmplist = newlist;
                }

                int snumber = k; //ion number within series
                // We need to map back the ion number to a specific
                // ion series (this is for cosmetics only)
                annotate_ion(snumber, k, sequence, curr_ion, params);

                // Find the isotope with the least amount of C13
                // that is above the specified range (only for
                // cosmetics so that we can report whether the hit
                // was on a monoisotopic precursor or not)
                for(isotope_nr=0;isotope_nr<=max_isotopes;isotope_nr++)
                {
                  q1_used = precursor.q1 + (MASS_DIFFC13 * isotope_nr)/precursor.q1_charge;
                  if(q1_used > q1_low) {break;}
                }
                tmplist.append( python::make_tuple(series[k], q1_used, 0, precursor.transition_group,
                  curr_ion, snumber, (std::string)sequence, precursor.ssrcalc, isotope_nr, ch));
                result[t1] = tmplist;

              }
            }

          } //loop over all transitions
        } // loop over all charges

      } //end of loop over all precursors

      delete [] series;
      delete [] tmp_series;
      return result;
    }

    using namespace python;
    BOOST_PYTHON_MODULE(c_getnonuis)
    {

      /* 
       * Which functions are used where:
       *  calculate_collisions_per_peptide_other_ion_series - precursor.py  / collider.py
       *  calculate_transitions_ch - collider.py / precursor.py to calculate transitions (in Lib)
       *
       *  3strikes.py: 
       *   - get_non_uis  (in Lib)
       *   - calculate_eUIS
       *
       *  run_paola.py: 
       *   - calculate_transitions_ch  (in Lib)
       *   - uses c_integrated for minNeeded
       *
       *  run_swath.py: 
       *   - calculate_transitions_ch  (in Lib)
       *   - calculate_density 
       * 
       * 
       *  cgi-scripts/srmcollider.py
       *   - calculate_transitions_ch (in Lib) 
       *   - calculate_collisions_per_peptide_other_ion_series 
       *   - _find_clashes_forall_other_series
       *   - get_non_uis  (in Lib)
       *
       *
      */

        def("calculate_collisions_per_peptide_other_ion_series", _find_clashes_calculate_collperpeptide_other_ion_series, 
     "Function to calculate the collisions_per_peptide out of a set of transitions\n"
     "and precursor peptides that are colliding peptides.  It will return a\n"
     "dictionary where for each key (colliding peptide key) the set of transitions\n"
     "of the query peptide are stored that are interfered with is held.\n"
     "Transitions are tuples of the form (q3, srm_id), precursors are tuples of\n"
     "the form (q1, sequence, peptide_key).\n"
     "\n"
     "\n"
     " Signature\n"
     "dict calculate_collisions_per_peptide_other_ion_series(tuple transitions, tuple precursors, \n"
     " parameter object, double q3_low, double q3_high, double q3window, bool ppm, \n"
     " bool forceChargecheck)\n"
                ""); 

        def("get_non_uis", SRMCollider::Combinatorics::get_non_uis, 
     "Function to calculate all non-UIS transitions from a dictionary \n"
     "where for each key (colliding peptide key) the set of transitions\n"
     "of the query peptide are stored that are interfered with is held\n"
     "(collisions_per_peptide dictionary)\n"
     "It will return a list of all non UIS of the requested order.\n"
     "\n"
     "\n"
     " Signature\n"
     "list get_non_uis(dict collisions_per_peptide, int order)\n"
     );


        def("calculate_transitions_ch", SRMCollider::Common::_py_calculate_transitions_with_charge,
     "Function to calculate all transitions of a list of precursor peptides and \n"
     "allows to select the charge states of these precursors.\n"
     "Precursors are tuples of the form (q1, sequence, peptide_key).\n"
     "It will return a list of tuples that are of the form \n"
     "(q3, q1, 0, peptide_key) \n"
     "\n"
     "\n"
     " Signature\n"
     "list calculate_transitions_ch(tuple precursors, list charges, double q3_low, double q3_high ) \n"
     );

        def("calculate_charged_mass", SRMCollider::Common::_py_calculate_charged_mass, 
     "Calculate the charged mass for a sequence, supplied in a python tuple in the\n"
     "first place.\n"
     "\n"
     "\n"
     " Signature\n"
     "double calculate_charged_mass(python::tuple clist, int ch)"

                
                );
       def("calculate_density", _find_clashes_calculate_colldensity, "");

       /*
       * Used by the web-scripts to report the exact interfering transitions.
       */
       def("_find_clashes_forall_other_series", _find_clashes_forall_other_series, "");

       /*
       * Used by to calculate eUIs
       */
       def ("calculate_eUIS", py_calculate_eUIS);

    }

  }
}
