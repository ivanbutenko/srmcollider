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

//include our own libraries
#ifndef SRMCOLLIDERLIB_H
#define SRMCOLLIDERLIB_H
#include "srmcollider.h"
#include <vector>
#include <string>

#include <boost/shared_ptr.hpp>

namespace SRMCollider 
{
  namespace Common 
  {
    struct AANotFound
    {

      AANotFound(std::string m)
      {
        message = m;
      }

      std::string message;

    };

    struct SRMTransition
    {
      double q3;
      long transition_id;
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

      bool ppm  ;
      double q3_window  ;

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

    struct SRMPrecursor
    {

      std::string sequence; 
      double q1;
      long transition_group;
      long parent_id;
      int q1_charge;
      int isotope_modification;
      double ssrcalc;

      SRMPrecursor()
      {
        isotope_modification = 0;
        q1_charge = 1;
      }

      SRMPrecursor(
        std::string sequence_, 
        double q1_,
        long transition_group_,
        long parent_id_,
        int q1_charge_,
        int isotope_modification_,
        double ssrcalc_
      )
      {
        sequence = sequence_;
        q1 = q1_;
        transition_group = transition_group_;
        parent_id =  parent_id_;
        q1_charge = q1_charge_;
        isotope_modification = isotope_modification_;
        ssrcalc = ssrcalc_;
      }

      int get_fragment_masses(double* tmp, double* series, double ch, const SRMParameters& params, int isotope_mod);

      /* 
       * Function to calculate the b and y fragment series given a peptide sequence 
       * Note that the y series is calculated "backwards" starting from left to right
       * instead of right to left. Whereas b_series[0] holds the b-1 ion, y_series[0]
       * holds the y-(n-1) ion and y_series[n-1] holds the y-1 ion.
       */
      // more OO, but signficantly slower (20% total) than the raw fragment
      void get_fragment_masses(std::vector<double>& result, double ch, const SRMParameters& params) const;

      double calculate_charged_mass();

      void set_maximal_charge(int mcharge);

      int get_maximal_charge();

      void calculate_transitions_with_charge(std::vector<int>& charges, std::vector<SRMTransition>& result, 
        double q3_low, double q3_high, const SRMParameters& param);

      //private:
      std::vector<boost::shared_ptr<SRMTransition> > transitions;

      int maximal_charge;
      private:

      int maximal_charge_;

    };

    /* 
     * Function to calculate the b and y fragment series given a peptide sequence 
     * Note that the y series is calculated "backwards" starting from left to right
     * instead of right to left. Whereas b_series[0] holds the b-1 ion, y_series[0]
     * holds the y-(n-1) ion and y_series[n-1] holds the y-1 ion.
     */
    int calculate_fragment_masses(const std::string& sequence, double* tmp, 
            double* series, double ch, const SRMParameters& params, int isotope_mod);

    void calculate_transitions_with_charge(SRMPrecursor& p, std::vector<int>& charges, std::vector<SRMTransition>& result, 
        double* tmp, double* series, double q3_low, double q3_high, const SRMParameters& param);
  }

}

#endif
