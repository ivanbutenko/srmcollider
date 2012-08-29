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
#ifndef SRMCOLLIDERLIB_C
#define SRMCOLLIDERLIB_C
#include "srmcolliderLib.h"

namespace SRMCollider 
{
  namespace Common 
  {

    int calculate_fragment_masses(const std::string& sequence, double* tmp, 
            double* series, double ch, const SRMParameters& params, int isotope_mod) 
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
                             sequence[start+4] == '7' )) 
                        {
                          throw AANotFound("Unknown modification for methionine");
                        }
                        res_mass = 147.03540462; break;
                    case 'C': 
                        if(!(sequence[start+2] == '1' && 
                             sequence[start+3] == '6' && 
                             sequence[start+4] == '0' )) 
                        {
                          throw AANotFound("Unknown modification for cystein");
                        }
                        res_mass = 160.030653721; break;
                    case 'N': 
                        if(!(sequence[start+2] == '1' && 
                             sequence[start+3] == '1' && 
                             sequence[start+4] == '5' )) 
                        {
                          throw AANotFound("Unknown modification for asparagine");
                        }
                        res_mass = 115.026945583; break;
                    case 'R': 
                        if(!(sequence[start+2] == '1' && 
                             sequence[start+3] == '6' && 
                             sequence[start+4] == '6' )) 
                        {
                          throw AANotFound("Unknown modification for arginine");
                        }
                        res_mass = 166.109379; break;
                    case 'K': 
                        if(!(sequence[start+2] == '1' && 
                             sequence[start+3] == '3' && 
                             sequence[start+4] == '6' )) 
                        {
                          throw AANotFound("Unknown modification for lysine");
                        }
                        res_mass = 136.109159; break;

                    default: 
                          throw AANotFound("Unknown modification");
                }
                //'M[147]':  131.04049 + mass_O), # oxygen
                //'C[160]':  103.00919 + mass_CAM - mass_H ), # CAM replaces H
                //'N[115]':  114.04293 - mass_N - mass_H + mass_O
                // 'K[136]' : ('heavy Lysine', 128.09496 + 8.014199),
                // 'R[166]' : ('heavy Arginine', 156.10111 + 10.008269),

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
                          throw AANotFound("Unknown amino acid");
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
                             sequence[start+4] == '7' )) 
                        {
                          throw AANotFound("Unknown modification for methionine");
                        }
                        res_mass = 148.0324344260; break;
                    case 'C': 
                        if(!(sequence[start+2] == '1' && 
                             sequence[start+3] == '6' && 
                             sequence[start+4] == '0' )) 
                        {
                          throw AANotFound("Unknown modification for cysteine");
                        }
                        res_mass = 161.027683399; break;
                    case 'N': 
                        if(!(sequence[start+2] == '1' && 
                             sequence[start+3] == '1' && 
                             sequence[start+4] == '5' )) 
                        {
                          throw AANotFound("Unknown modification for asparagine");
                        }
                        res_mass = 116.023977918; break;
                    // SILAC label and N15 is not very likely
                    // case 'R': 
                    // case 'K': 

                    default: 
                          throw AANotFound("Unknown modification");
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
                          throw AANotFound("Unknown amino acid");
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

    double SRMPrecursor::calculate_charged_mass() 
    {
        std::vector<double> result;
        SRMParameters params;
        params.MMinusH2O = true;
        params.yions = false;
        params.bions = false;
        get_fragment_masses(result, 1, params);
     
        // In order to get the full mass, we get the (M-H2O)+ ion and then add H2O
        // to it.
        return (result[0] + MASS_H + MASS_OH + MASS_H* this->q1_charge ) / this->q1_charge  ;
    }        

    void SRMPrecursor::set_maximal_charge(int mcharge) 
    {
      maximal_charge_ = mcharge;
    }

    int SRMPrecursor::get_maximal_charge() 
    {
      return \
        std::count(sequence.begin(), sequence.end(), 'R') + 
        std::count(sequence.begin(), sequence.end(), 'H') + 
        std::count(sequence.begin(), sequence.end(), 'K');
    }

    void SRMPrecursor::calculate_transitions_with_charge(std::vector<int>& charges, std::vector<SRMTransition>& result, 
      double q3_low, double q3_high, const SRMParameters& param)
  {

    size_t k;
    double q3;
    std::vector<double> series;

    for (std::vector<int>::iterator ch_it = charges.begin(); ch_it != charges.end(); ch_it++) {
      get_fragment_masses(series, (*ch_it), param);

      // go through all fragments of this precursor
      for (k=0; k<series.size(); k++) {
        q3 = series[k];
        if (q3 > q3_low && q3 < q3_high)
        {
          SRMTransition t = {q3, k};
          result.push_back(t);
        }
      }
    }

  }

    int SRMPrecursor::get_fragment_masses(double* tmp, double* series, double ch, const SRMParameters& params, int isotope_mod) 
    {
      return calculate_fragment_masses(sequence, tmp, series, ch, params, isotope_mod);
    }

    void SRMPrecursor::get_fragment_masses(std::vector<double>& result, double ch, const SRMParameters& params) const
    {

      double* series = new double[1024];
      double* tmp_series = new double[1024];
      int l = calculate_fragment_masses(sequence, tmp_series, series, ch, params, isotope_modification);
      result.assign(series, series+l);
      delete series;
      delete tmp_series;
    }

    void calculate_transitions_with_charge(SRMPrecursor& p, std::vector<int>& charges, std::vector<SRMTransition>& result, 
        double* tmp, double* series, double q3_low, double q3_high, const SRMParameters& param)
    {
      p.calculate_transitions_with_charge(charges, result, q3_low, q3_high, param);
    }

  }

}
#endif
