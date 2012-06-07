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

#ifndef SRMCOLLIDER_INTEGRATED_H
#define SRMCOLLIDER_INTEGRATED_H

//include our own libraries
#include "srmcollider.h"
#include "rangetree.cpp"
#include "combinatorics.h"
#include "srmcolliderLib.cpp"

#include <iostream>
#include <vector>
#include <set>

#include <vector>

using namespace SRMCollider::Common;
using namespace SRMCollider;
//using namespace std;

namespace SRMCollider 
{
  namespace IntegratedRun 
  {

    typedef SRMCollider::Common::SRMTransition Transition;

    inline void _charged_interference(const std::string& sequence, double* tmp_series, double* series, const int ch,
        const SRMParameters& params, const int isotope_modification, 
        const std::vector<Transition>& mytransitions, const bool ppm, const double q3window, COMBINT& currenttmp);

    inline void _charged_interference(const SRMPrecursor& precursor, const int ch, const SRMParameters& params, 
        const std::vector<Transition>& mytransitions, COMBINT& currenttmp);

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
      int max_uis, double q3window, bool ppm, SRMParameters& params);

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
    SRMCollider::ExtendedRangetree::Rangetree_Q1_RT& rtree, std::vector<COMBINT>& newcollperpep);

  }
}
#endif
