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
 * This file contains functions that allow the user to execute one calculation
 * of UIS or minimal number of transitions needed for a peptide with one call.
 * This means that instead of calling several functions after each other and
 * passing data structures between Python and Cpp, all can be done in Cpp which
 * is faster.
 *
 * Thus it combines the functionality of rangetree.cpp and getNonUis.cpp in one.
 *
*/

#ifndef SRMCOLLIDER_INTEGRATED_C
#define SRMCOLLIDER_INTEGRATED_C

//include our own libraries
#include "srmcollider.h"
#include "rangetree.cpp"
#include "combinatorics.h"
#include "srmcolliderLib.cpp"

// Boost.Python headers
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
namespace python = boost::python;

#include <vector>

using namespace SRMCollider::Common;
using namespace SRMCollider;
//using namespace std;

namespace SRMCollider 
{
  namespace IntegratedRun 
  {

    typedef SRMCollider::Common::SRMTransition Transition;

    void _pyToC_integratedrunTransitions(python::tuple& transitions, std::vector<Transition>& mytransitions)
    {
      /*
      * Transitions are tuples of the form (q3, transition_id)
      * convert to our struct.
      */
      int transitions_length = python::extract<int>(transitions.attr("__len__")());
      python::tuple tlist;
      for (int i=0; i<transitions_length; i++) {
          tlist = python::extract< python::tuple >(transitions[i]);
          double q3 = python::extract<double>(tlist[0]);
          long transition_id = python::extract<long>(tlist[1]);
          Transition entry = {q3, transition_id};
          mytransitions[i] = entry;
      }
    }

    inline void _charged_interference(const std::string& sequence, double* tmp_series, double* series, const int ch,
        const SRMParameters& params, const int isotope_modification, 
        const std::vector<Transition>& mytransitions, const bool ppm, const double q3window, COMBINT& currenttmp)
    {
      COMBINT one;
      double q3, q3used = q3window;
      Transition transition;
      int k;

      int fragcount = SRMCollider::Common::calculate_fragment_masses(
          sequence, tmp_series, series, ch, params, isotope_modification);
      //std::cout <<" with seq " << sequence << "/"<< ch << " i got " << fragcount << " vs " << mytransitions.size() << " wi " << q3window << " " << ppm << std::endl; 

      for (size_t i=0; i<mytransitions.size(); i++) {

          transition = mytransitions[i];
          q3 = transition.q3;
          //ppm is 10^-6
          if(ppm) {q3used = q3window / 1000000.0 * q3; } 

              // go through all fragments of this precursor
              for (k=0; k<fragcount; k++) {
                  if(fabs(q3-series[k]) < q3used ) {
                      //left bitshift == 2^i
                      one = 1;
                      currenttmp |= one << i;
                  }
              }
      } //loop over all transitions
    }

    inline void _charged_interference(const SRMPrecursor& precursor, const int ch, const SRMParameters& params, 
        const std::vector<Transition>& mytransitions, COMBINT& currenttmp)
    {
      COMBINT one;
      double q3, q3used = params.q3_window;
      Transition transition;
      size_t k;

      std::vector<double> result;
      result.reserve(1024);
      precursor.get_fragment_masses(result, ch, params);
      //std::cout <<" with seq " << sequence << "/"<< ch << " i got " << fragcount << " vs " << mytransitions.size() << " wi " << q3window << " " << ppm << std::endl; 

      for (size_t i=0; i<mytransitions.size(); i++) {
          transition = mytransitions[i];
          q3 = transition.q3;
          //ppm is 10^-6
          if(params.ppm) {q3used = params.q3_window / 1000000.0 * q3; } 

              // go through all fragments of this precursor
              //for (k=0; k<fragcount; k++) 
              for (k=0; k<result.size(); k++) 
              {
                  if(fabs(q3-result[k]) < q3used ) {
                      //left bitshift == 2^i
                      one = 1;
                      currenttmp |= one << i;
                  }
              }
      } //loop over all transitions

    }

  /* 
   * Calculate the minimally needed number of transitions needed to get a UIS
   * given transitions in their preferred order and interfering precursors.
   *
   * Uses integers with bitflags to store the combinations, thus the number of
   * transitions that can be considered is limited by COMBLIMIT.  Also note that
   * the leftmost bit needs to be set to zero all the time, otherwise an infinte
   * loop occurs, thus there is one bit we cannot use.
  */
  int min_needed(std::vector<Transition>& mytransitions, std::vector<SRMPrecursor>& precursors,
      int max_uis, double q3window, bool ppm, SRMParameters& params)
  {

      //use the defined COMBINT (default 32bit int) and some bitwise operation to do this :-)
      COMBINT one;
      COMBINT currenttmp = 0;
      std::vector<COMBINT> newcollperpep;

      int i, ch;
      size_t transitions_length = mytransitions.size();

      double* series = new double[1024];
      double* tmp_series = new double[1024];

      //Check whether we have more transitions than we have bits in our number
      if (transitions_length > COMBLIMIT) {
          // PyErr_SetString(PyExc_ValueError, 
          //     "Too many transitions, please adjust limit.");
          // boost::python::throw_error_already_set();
          return -1;
      }

      // Go through all (potential) collisions 
      //
      int maxoverlap = 0;
      for (size_t j=0; j< precursors.size(); j++) 
      {

        for (ch=1; ch<=2; ch++) 
        {
          // TODO do N15
          _charged_interference(precursors[j].sequence, tmp_series, series, ch,
            params, NOISOTOPEMODIFICATION, mytransitions, ppm, q3window, currenttmp);
        } //end loop over all charge states of this precursor
        if ( currenttmp ) 
        {
          //while we find the transitions from our relative order also in the
          //peptide we just looked at, increase i
          one = 1;
          i = 0;
          while( currenttmp & one << i ) i++;
          if( i > maxoverlap ) maxoverlap = i;
          currenttmp = 0;
        }
      }

      delete series;
      delete tmp_series;

      //we have counted from 0 to N-1, but now want the number of transitions
      ++maxoverlap;
      if(maxoverlap > (int)transitions_length) {return -1;}
      return maxoverlap;
  }

  // Python wrapper for min_needed
  int _py_min_needed(python::tuple transitions, python::tuple precursors,
      int max_uis, double q3window, bool ppm, python::object par)   
  {

      python::tuple clist;
      int i;
      std::string sequence;
      int transitions_length = python::extract<int>(transitions.attr("__len__")());
      int precursor_length = python::extract<int>(precursors.attr("__len__")());
      python::tuple tlist;
      
      //Check whether we have more transitions than we have bits in our number
      if (transitions_length > COMBLIMIT) {
          PyErr_SetString(PyExc_ValueError, 
              "Too many transitions, please adjust limit.");
          boost::python::throw_error_already_set();
          return -1;
      }

      SRMParameters params;
      pyToC::initialize_param_obj(par, params);

      std::vector<Transition> mytransitions(transitions_length);
      _pyToC_integratedrunTransitions(transitions, mytransitions);

      std::vector<SRMPrecursor> myprecursors(precursor_length);
      for (i=0; i<precursor_length; i++) {
          tlist = python::extract< python::tuple >(precursors[i]);
          sequence = python::extract<std::string>(tlist[1]);
          SRMPrecursor p;
          p.sequence = sequence;
          p.transition_group = p.q1 = p.maximal_charge = p.ssrcalc = -1;
          myprecursors[i] = p;
      }

      return min_needed(mytransitions, myprecursors, max_uis, q3window, ppm, params);
  }

  /*
   * Return the number of non-UIS for all orders up to max_uis
   * Given the transitions and then four numbers giving the coordinate window in
   * which the precursors should be found: 
   *   q1_low, ssrcalc_low, q1_high, ssrcalc_high, 
   * Also the peptide key of the current peptide has to be given (to exclude it)
   * as well as the maximal order of UIS to compute, the used q3 window and
   * whether the q3 window is given in ppm (parts per million), default unit is
   * Th which is m/z
   *
   * This function uses "bitwise operations" to speed up the calculation: each combination of
   * transitions is encoded in an integer, which allows easy storage and
   * computation by bitshits even though the code gets less easy to understand.
   * The default length of the integer is 32 bits thus we can at most handle
   * combinations for 31 transitions; if more are needed one needs to change
   * COMBINT to an integer format that usees more bits, e.g. unit64_t or even the
   * BigInt class for larger integers.
   * If there are more transitions provided than allowed, an error will be
   * thrown. 
  */
  void wrap_all_bitwise(std::vector<Transition>& mytransitions, double a, double b,
    double c, double d, long thistransitiongr, int max_uis, double q3window,
    bool ppm, int max_nr_isotopes, double isotope_correction, SRMParameters& params,
    SRMCollider::ExtendedRangetree::Rangetree_Q1_RT& rtree, std::vector<COMBINT>& newcollperpep)
  {
      //use the defined COMBINT (default 32bit int) and some bitwise operations to do this :-)
      COMBINT currenttmp = 0;
      int ch;

      SRMCollider::ExtendedRangetree::Precursor precursor;
      std::vector<SRMCollider::ExtendedRangetree::Key> OutputList;

      int iso;
      double q1_low = a; double q1_high = c;

      //Check whether we have more transitions than we have bits in our number
      int transitions_length = mytransitions.size();
      if (transitions_length > COMBLIMIT) {
          //PyErr_SetString(PyExc_ValueError, 
          //    "Too many transitions, please adjust limit.");
          //boost::python::throw_error_already_set();
      }

      // search for all matching peptides (including all potential isotopes)
      SRMCollider::ExtendedRangetree::Interval win(SRMCollider::ExtendedRangetree::Interval(
        SRMCollider::ExtendedRangetree::K::Point_2(a-isotope_correction,b),
        SRMCollider::ExtendedRangetree::K::Point_2(c,d)));
      rtree.my_rangetree->window_query(win, std::back_inserter(OutputList));
      std::vector<SRMCollider::ExtendedRangetree::Key>::iterator current=OutputList.begin();

      double* series = new double[1024];
      double* tmp_series = new double[1024];

      // Go through all (potential) collisions we just extracted from the rangetree
      //
      // This assumes that we do not have peptide_keys mapped to different
      // colliding transitions. In fact we do not have this situation even though
      // we have duplicate peptide_keys (from the isotopes). But they will
      // produce the same interefering transitions and thus the same entry in the
      // collisions per peptide table.

      while(current!=OutputList.end())
      {

        SRMPrecursor& p = current->second;
        p.q1 = current->first[0];

        bool proceed = false;
        // check whether there are any relevant isotopes, otherwise exclude this hit
        for (iso=0; iso<=max_nr_isotopes; iso++) {
            if (p.q1 + (MASS_DIFFC13 * iso)/p.q1_charge > q1_low && 
                p.q1 + (MASS_DIFFC13 * iso)/p.q1_charge < q1_high) 
            {
                proceed = true;
            }
        }

        if (!proceed || thistransitiongr == p.transition_group) {
          // go to next
          current++;
          continue;
        }

        for (ch=1; ch<=2; ch++) {

#if 0
          // 20 % speed penalty with the std::vector solution
          _charged_interference(p, ch, params, mytransitions, currenttmp);
#else
          _charged_interference(p.sequence, tmp_series, series, ch, params, 
            p.isotope_modification, mytransitions, ppm, q3window, currenttmp);
#endif
        }


        //Store current combination
        if ( currenttmp ) {
            newcollperpep.push_back(currenttmp);
            currenttmp = 0;
        }

        current++;
      }

      delete series;
      delete tmp_series;
  }

  // Python wrapper for wrap_all
  python::list _py_wrap_all_bitwise(python::tuple& py_transitions, double a, double b,
    double c, double d, long thistransitiongr, int max_uis, double q3window,
    bool ppm, int max_nr_isotopes, double isotope_correction, python::object& par,
    SRMCollider::ExtendedRangetree::Rangetree_Q1_RT& rtree)
  {

      //Check whether we have more transitions than we have bits in our number
      int transitions_length = python::extract<int>(py_transitions.attr("__len__")());
      if (transitions_length > COMBLIMIT) {
          PyErr_SetString(PyExc_ValueError, 
              "Too many transitions, please adjust limit.");
          boost::python::throw_error_already_set();
          python::list tlist;
          return tlist;
      }

      SRMParameters params;
      pyToC::initialize_param_obj(par, params);

      std::vector<Transition> mytransitions(transitions_length);
      _pyToC_integratedrunTransitions(py_transitions, mytransitions);

      std::vector<COMBINT> collisions_per_peptide; 
      wrap_all_bitwise(mytransitions, a, b, c, d,
          thistransitiongr, max_uis, q3window, ppm, max_nr_isotopes, isotope_correction, 
           params, rtree, collisions_per_peptide);

      std::vector<int> c_result;
      for(int i =1; i<= max_uis; i++) {
        std::set<COMBINT> combinations;
        SRMCollider::Combinatorics::get_non_uis_bitwise(collisions_per_peptide, transitions_length, i, combinations);
        c_result.push_back(combinations.size());
      }

      python::list result;
      for (size_t i = 0; i < c_result.size(); i++)
      {
        result.append(c_result[i]);
      }
      return result;
  }

  // Expose to Python
  using namespace python;
  BOOST_PYTHON_MODULE(c_integrated)
  {

      def("wrap_all_bitwise", _py_wrap_all_bitwise,
              
   "Return the number of non-UIS for all orders up to max_uis\n"
   "Given the transitions and then four numbers giving the coordinate window in\n"
   "which the precursors should be found: \n"
   "  q1_low, ssrcalc_low, q1_high, ssrcalc_high, \n"
   "Also the peptide key of the current peptide has to be given (to exclude it)\n"
   "as well as the maximal order of UIS to compute, the used q3 window and\n"
   "whether the q3 window is given in ppm (parts per million), default unit is\n"
   "Th which is m/z\n"
   "\n"
   "This function uses bitwise operations to speed up the calculation: each combination of\n"
   "transitions is encoded in an integer, which allows easy storage and\n"
   "computation by bitshits even though the code gets less easy to understand.\n"
   "The default length of the integer is 32 bits thus we can at most handle\n"
   "combinations for 31 transitions; if more are needed one needs to change\n"
   "COMBINT to an integer format that usees more bits, e.g. unit64_t or even the\n"
   "BigInt class for larger integers.\n"
   "If there are more transitions provided than allowed, an error will be\n"
   "thrown. \n"
   "\n"
   "\n"
   " Signature\n"
   "list wrap_all_bitwise(python::tuple transitions, double a, double b,\n"
          "double c, double d, long thistransitiongr, int max_uis, double q3window,\n"
          "bool ppm, int max_nr_isotopes, double isotope_correction)"
      );

      def("getMinNeededTransitions", _py_min_needed,
              
   "Calculate the minimally needed number of transitions needed to get a UIS\n"
   "given transitions in their preferred order and interfering precursors.\n"
   "\n"
   "Uses integers with bitflags to store the combinations, thus the number of\n"
   "transitions that can be considered is limited by COMBLIMIT.  Also note that\n"
   "the leftmost bit needs to be set to zero all the time, otherwise an infinte\n"
   "loop occurs, thus there is one bit we cannot use.\n"
   "\n"
   "\n"
   " Signature\n"
   "int min_needed(tuple transitions, tuple precursors,\n"
      "int max_uis, double q3window, bool ppm )\n"
              );
  }

  }
}
#endif
