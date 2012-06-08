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
#include "rangetree.h"
#include "combinatorics.h"
#include "srmcolliderLib.h"

#include <vector>

using namespace SRMCollider::Common;
//using namespace std;

namespace SRMCollider 
{
  namespace IntegratedRun 
  {

    typedef SRMCollider::Common::SRMTransition Transition;

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

  }
}
#endif
