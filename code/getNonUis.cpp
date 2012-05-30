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

//include our own libraries
#include "srmcollider.h"
#include "srmcolliderLib.cpp"

//using namespace std;
using namespace SRMCollider::Common;

bool SortIntDoublePairSecond(const std::pair<int,double>& left, const std::pair<int,double>& right)
{
  return left.second < right.second;
}


// Function declarations

// functions to calculate collperpeptide dictionary
//// python::dict _find_clashes_calculate_collperpeptide(python::tuple transitions,
////         python::tuple precursors, double q3_low, double q3_high, double
////         q3window, bool ppm);
python::dict _find_clashes_calculate_collperpeptide_other_ion_series(
        python::tuple transitions, python::list precursors, python::object par, 
        double q3_low, double q3_high, double q3window, bool ppm, bool forceChargeCheck=false);

// calculate the number of collisions for each transition
python::list _find_clashes_calculate_colldensity(python::tuple transitions,
    python::list precursors, double q3_low, double q3_high, double q3window,
    bool ppm) ;

// Function to calculate the exact interfering transitions for each peptide for
// series other than b/y.
python::dict _find_clashes_forall_other_series(python::tuple transitions,
    python::list precursors, python::object par, double q3_low, double q3_high, double q3window,
    bool ppm, double q1_low, bool forceChargeCheck);
void _find_clashes_forall_other_series_sub( int& l, int ch, int k,
    const char* sequence, std::string& curr_ion, python::object par);

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

    std::vector<SRMPrecursor> c_precursors;
    std::vector<SRMTransition> c_transitions;

    int transitions_length = python::extract<int>(transitions.attr("__len__")());
    int precursor_length = python::extract<int>(precursors.attr("__len__")());
    int fragcount, i, j, k, ch, listmembers = 0;

    long t1, peptide_key;
    double t0, q3used = q3window;

    double* series = new double[1024];
    double* tmp_series = new double[1024];

    SRMParameters params;
    params.aions      =  python::extract<bool>(par.attr("aions"));
    params.aMinusNH3  =  python::extract<bool>(par.attr("aMinusNH3"));
    params.bions      =  python::extract<bool>(par.attr("bions"));
    params.bMinusH2O  =  python::extract<bool>(par.attr("bMinusH2O"));
    params.bMinusNH3  =  python::extract<bool>(par.attr("bMinusNH3"));
    params.bPlusH2O   =  python::extract<bool>(par.attr("bPlusH2O"));
    params.cions      =  python::extract<bool>(par.attr("cions"));
    params.xions      =  python::extract<bool>(par.attr("xions"));
    params.yions      =  python::extract<bool>(par.attr("yions"));
    params.yMinusH2O  =  python::extract<bool>(par.attr("yMinusH2O"));
    params.yMinusNH3  =  python::extract<bool>(par.attr("yMinusNH3"));
    params.zions      =  python::extract<bool>(par.attr("zions"));
    params.MMinusH2O  =  python::extract<bool>(par.attr("MMinusH2O"));
    params.MMinusNH3  =  python::extract<bool>(par.attr("MMinusNH3"));

    for (int i=0; i<python::extract<int>(precursors.attr("__len__")()); i++) {
      SRMPrecursor p;
      python::object precursor = python::extract< python::object >(precursors[i]);
      p.sequence = python::extract<char *>(precursor.attr("modified_sequence"));
      p.isotope_modification = python::extract<int>(precursor.attr("isotopically_modified"));
      p.q1_charge = python::extract<int>(precursor.attr("q1_charge"));
      p.maximal_charge = python::extract<int>(precursor.attr("to_peptide")().attr("get_maximal_charge")());
      p.transition_group = python::extract<long>(precursor.attr("transition_group"));
      c_precursors.push_back(p);
    }

    for (int i=0; i<python::extract<int>(transitions.attr("__len__")()); i++) {
      SRMTransition t;
      tlist = python::extract< python::tuple >(transitions[i]);
      t.q3 = python::extract< double >(tlist[0]);
      t.transition_id = python::extract<long>(tlist[1]);
      c_transitions.push_back(t);
    }

    // go through all (potential) interfering precursors and store the
    // colliding SRM ids in a dictionary (they can be found at position 3 and 1
    // respectively).
    for (j=0; j<precursor_length; j++) {
        SRMPrecursor & precursor = c_precursors[j];

        for (ch=1; ch<=2; ch++) {
            fragcount = calculate_fragment_masses(precursor.sequence, tmp_series, series, ch,
                  params, precursor.isotope_modification);

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
        //Sorting it costs something and could be done more efficiently since 
        //in fact, we only have to merge 2 presorted arrays. TODO Still the cost
        //is negligible.
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
    bool ppm) {

    python::dict collisions_per_peptide, tmpdict;
    python::tuple clist;
    python::tuple tlist;

    python::list tmplist;

    int transitions_length = python::extract<int>(transitions.attr("__len__")());
    int precursor_length = python::extract<int>(precursors.attr("__len__")());
    int fragcount, i, j, k, ch;

    SRMParameters param;

    double q3used = q3window;
    char* sequence;

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
        sequence = python::extract<char *>(clist[1]);

        for (ch=1; ch<=2; ch++) {
            //fragcount = _calculate_clashes(sequence, b_series, y_series, ch);
            fragcount = calculate_fragment_masses(sequence, tmp, series, ch, param, NOISOTOPEMODIFICATION);

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

/*
 * Function to calculate the exact interfering transitions for each peptide.
 * It will return a Transitions are tuples of the form (q3, srm_id), precursors
 * are tuples of the form (q1, sequence, peptide_key).
 *
 * Used by the web-scripts to report the exact interfering transitions.
 */
python::dict _find_clashes_forall_other_series(python::tuple transitions,
    python::list precursors, python::object par, double q3_low, double q3_high, double q3window,
    bool ppm, double q1_low, bool forceChargeCheck) {

    python::dict result, tmpdict;
    //python::object precursor;
    python::tuple tlist;
    python::list tmplist;

    std::vector<SRMPrecursor> c_precursors;
    std::vector<SRMTransition> c_transitions;

    int transitions_length = python::extract<int>(transitions.attr("__len__")());
    int precursor_length = python::extract<int>(precursors.attr("__len__")());
    int fragcount, i, j, k, ch;
    int isotope_nr;

    long t1;
    double t0, q3used = q3window, q1_used;
    char* sequence;
    int max_isotopes = python::extract<int>(par.attr("isotopes_up_to"));

    double* series = new double[10*256];
    double* tmp_series = new double[256];

    for (int i=0; i<python::extract<int>(precursors.attr("__len__")()); i++) {
      SRMPrecursor p;
      python::object precursor = python::extract< python::object >(precursors[i]);
      p.q1 = python::extract<double>(precursor.attr("q1"));
      p.sequence = python::extract<char *>(precursor.attr("modified_sequence"));
      p.transition_group = python::extract<long>(precursor.attr("transition_group"));

      p.ssrcalc = python::extract<double>(precursor.attr("ssrcalc"));
      p.q1_charge = python::extract<int>(precursor.attr("q1_charge"));
      p.maximal_charge = python::extract<int>(precursor.attr("to_peptide")().attr("get_maximal_charge")());
      c_precursors.push_back(p);
    }

    for (int i=0; i<python::extract<int>(transitions.attr("__len__")()); i++) {
      SRMTransition t;
      tlist = python::extract< python::tuple >(transitions[i]);
      t.q3 = python::extract< double >(tlist[0]);
      t.transition_id = python::extract<long>(tlist[1]);
      c_transitions.push_back(t);
    }

    std::string curr_ion = "?";

    // go through all (potential) collisions
    // and store the colliding SRM ids in a dictionary (they can be found at
    // position 3 and 1 respectively)
    for (j=0; j<precursor_length; j++) {
        SRMPrecursor & precursor = c_precursors[j];
        sequence = c_precursors[j].sequence;

        for (ch=1; ch<=2; ch++) {
            fragcount = _calculate_clashes_other_series(sequence, tmp_series,
                    series, ch, par);

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
                        if( result.has_key(t1) ) {
                            int snumber = k; //ion number within series
                            // We need to map back the ion number to a specific
                            // ion series (this is for cosmetics only)
                            _find_clashes_forall_other_series_sub(snumber, ch,
                                    k, sequence, curr_ion, par);

                            // Find the isotope with the least amount of C13
                            // that is above the specified range (only for
                            // cosmetics so that we can report whether the hit
                            // was on a monoisotopic precursor or not)
                            for(isotope_nr=0;isotope_nr<=max_isotopes;isotope_nr++)
                            {
                                q1_used = precursor.q1 + (MASS_DIFFC13 * isotope_nr)/precursor.q1_charge;
                                if(q1_used > q1_low) {break;}
                            }
                            tmplist = python::extract<python::list>(result[t1]);
                            tmplist.append( python::make_tuple(series[k], q1_used, 0, precursor.transition_group,
                              curr_ion, snumber, (std::string)sequence, precursor.ssrcalc, isotope_nr, ch));
                        }
                        else{
                            python::list newlist;
                            tmplist = newlist;
                            int snumber = k; //ion number within series
                            // We need to map back the ion number to a specific
                            // ion series (this is for cosmetics only)
                            _find_clashes_forall_other_series_sub(snumber, ch,
                                    k, sequence, curr_ion, par);

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
                }
            } //loop over all transitions
        }

    } //end of loop over all precursors

    delete [] series;
    delete [] tmp_series;
    return result;
}

void _find_clashes_forall_other_series_sub( int& l, int ch, int k,
    const char* sequence, std::string& curr_ion, python::object par) {

    bool aions      =  python::extract<bool>(par.attr("aions"));
    bool aMinusNH3  =  python::extract<bool>(par.attr("aMinusNH3"));
    bool bions      =  python::extract<bool>(par.attr("bions"));
    bool bMinusH2O  =  python::extract<bool>(par.attr("bMinusH2O"));
    bool bMinusNH3  =  python::extract<bool>(par.attr("bMinusNH3"));
    bool bPlusH2O   =  python::extract<bool>(par.attr("bPlusH2O"));
    bool cions      =  python::extract<bool>(par.attr("cions"));
    bool xions      =  python::extract<bool>(par.attr("xions"));
    bool yions      =  python::extract<bool>(par.attr("yions"));
    bool yMinusH2O  =  python::extract<bool>(par.attr("yMinusH2O"));
    bool yMinusNH3  =  python::extract<bool>(par.attr("yMinusNH3"));
    bool zions      =  python::extract<bool>(par.attr("zions"));
    bool MMinusH2O  =  python::extract<bool>(par.attr("MMinusH2O"));
    bool MMinusNH3  =  python::extract<bool>(par.attr("MMinusNH3"));

    int scounter, icounter;
    double* tmp = new double[256];
    double* series = new double[256];
    bool done = false;

    // get the number of amino acids in the sequence
    SRMParameters params;
    params.bions = false;
    scounter = calculate_fragment_masses(sequence, tmp, series, ch,  params, NOISOTOPEMODIFICATION);
    scounter++;
    done = false;
    icounter = 0;
    curr_ion = "?";

    if (yions && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "y"; done = true; break;}; icounter++;}
    if (bions && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "b"; done = true; break;}; icounter++;}
    if (aions && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "a"; done = true; break;}; icounter++;}
    if (cions && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "c"; done = true; break;}; icounter++;}
    if (xions && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "x"; done = true; break;}; icounter++;}
    if (zions && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "z"; done = true; break;}; icounter++;}

    if (aMinusNH3 && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "aMinusNH3"; done = true; break;}; icounter++;}
    if (bMinusH2O && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "bMinusH2O"; done = true; break;}; icounter++;}
    if (bMinusNH3 && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "bMinusNH3"; done = true; break;}; icounter++;}
    if (bPlusH2O && !done)  for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "bPlusH2O"; done = true; break;}; icounter++;}

    if (yMinusH2O && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "yMinusH2O"; done = true; break;}; icounter++;}
    if (yMinusNH3 && !done) for (l=0; l<scounter-1; l++) {if(icounter==k) {curr_ion = "yMinusNH3"; done = true; break;}; icounter++;}

    if (MMinusH2O && !done) {l=0;if(icounter==k){curr_ion = "MMinusH2O"; done = true;}; icounter++;}
    if (MMinusNH3 && !done) {l=0;if(icounter==k){curr_ion = "MMinusNH3"; done = true;}; icounter++;}


    l++; // ion series starts at 1, thus add one

    delete [] tmp;
    delete [] series;
}

/*
 * Function to calculate whether there exists any combination of values from M
 * arrays (one value from each array) such that the M values are within a
 * certain window. 
 *
 * The first argument is a list that contains the length of the lists to be
 * checked. The second argument a list of lists which contain the values to be
 * checked, the third argument the window size.
 *
 * The function returns a all the "forbidden" tuples, e.g. tuples of 
 * transitions that are interfering and can thus not be used for an eUIS.
 *
 */
void calculate_eUIS(std::vector<int>& N, std::vector<std::vector<double> >& c_ssrcalcvalues,
    double ssrwindow, std::vector<std::vector<int> >& all_nonuis) 
{

    int M = (int)N.size();

    int k, i;
    unsigned int m, n;
    int sumlen = 0;
    std::vector<int> index; index.resize(M);
    std::vector<int> sort_idx; sort_idx.resize(M);
    std::vector<bool> discarded_indices; discarded_indices.resize(M);
    std::vector<double> myssr; myssr.resize(M);
    std::vector< std::pair<int,double> > with_index;

    for(int k=0;k<M;k++) index[k] = 0;
    for(int k=0;k<M;k++) sumlen += N[k];
    for(int k=0;k<M;k++) discarded_indices[k] = false;

    double max_elem = 0;
    for(k = 0; k < M; k++) {
      for(i = 0; i < N[k]; i++) {
        if (c_ssrcalcvalues[k][i] > max_elem) {max_elem = c_ssrcalcvalues[k][i];}
      }
    }

    //# check whether there are any empty ssrcalcvalues
    int cnt =0;
    for(k=0; k < M; k++) {
      if(N[k] == 0)
      {
        discarded_indices[k] = true;
        cnt++;
      }
    }

    if(cnt==M) {return;}

    while(true) {

        for(k=0; k < M; k++) {
          if(!discarded_indices[k])
          {
            myssr[k] = c_ssrcalcvalues[k][ index[k] ];
          }
        }

        // find the pivot element
        double smin = max_elem;
        int piv_i = -1;
        for(k=0; k < M; k++) 
        {
          if(!discarded_indices[k])
            if(myssr[k] <= smin) {
              smin = myssr[k];
              piv_i = k;
            }
        }

        with_index.resize(0);
        //# we need to sort by we also need to have a map back to retrieve the original!
        // store them in a pair with the index, sort, retrieve values and index
        for(k=0; k < M; k++) {
          with_index.push_back(std::make_pair(k,myssr[k]));
        }
        std::stable_sort(with_index.begin(), with_index.end(), SortIntDoublePairSecond);
        for(k=0; k < M; k++) {
          sort_idx[k] = with_index[k].first;
          myssr[k] = with_index[k].second;
        }

        //# now find all N different combinations that are not UIS. Since they are
        //# sorted we only need to consider elements that have a higher index.
        for(k=0; k < M; k++) {
          if(discarded_indices[sort_idx[k]]) continue;

          std::vector<int> nonuis;
          nonuis.push_back(sort_idx[k]);

          for(i=k+1; i < M; i++) 
          {
              if(discarded_indices[sort_idx[i]]) continue;
              if(!(fabs(myssr[k] - myssr[i]) > ssrwindow))
              {
                nonuis.push_back(sort_idx[i]);
              }
          }

          // here we have to figure out whether we want to append it (e.g. if
          // it is not yet present).
          std::sort(nonuis.begin(), nonuis.end());
          bool is_present = false;
          bool this_not_present = true;
          for(m=0; m<all_nonuis.size(); m++)
          {
            this_not_present = false;
            for(n=0; n<all_nonuis[m].size() and n < nonuis.size(); n++)
            {
                if(nonuis[n] != all_nonuis[m][n])
                {
                  this_not_present = true; break;
                }
            }
            if(all_nonuis[m].size() != nonuis.size()) {this_not_present = true;}

            //cout << " compared " << m << endl;
            if(!this_not_present) {
              is_present = true; 
              break;
            }
          }

          if(!is_present)
          {
            all_nonuis.push_back(nonuis);

          }
        }

        //# Advance the pivot element
        index[piv_i] += 1;
        if(index[piv_i] >= N[piv_i])
        {
            discarded_indices[ piv_i ] = true;
            int dcount = 0;
            while(dcount < M && discarded_indices[dcount]) dcount++;
            if(dcount >= M)
                break;
        }

    }
    return;
}

// Python wrapper for calculate_eUIS
python::list py_calculate_eUIS(python::list myN, python::list py_ssrcalcvalues, double ssrwindow) 
{

    python::list result;
    std::vector<std::vector<int> > all_nonuis;

    int M = python::extract<int>(myN.attr("__len__")());
    int k, i;
    unsigned int m, n;

    // fill up N
    std::vector<int> N; N.resize(M);
    for(int k=0;k<M;k++) N[k] = python::extract<int>(myN[k]);

    // fil up c_ssrcalc values
    std::vector<std::vector<double> > c_ssrcalcvalues; c_ssrcalcvalues.resize(M);
    for(k = 0; k < M; k++) { c_ssrcalcvalues[k].resize(N[k]); }
    python::list tmplist;
    for(k = 0; k < M; k++) {
        tmplist = python::extract<python::list>(py_ssrcalcvalues[k]);
        for(i = 0; i < N[k]; i++) {
            c_ssrcalcvalues[k][i] = python::extract<double>(tmplist[i]); 
        }
    }

    calculate_eUIS(N, c_ssrcalcvalues, ssrwindow, all_nonuis);

    // convert to python datastructure
    for(m=0; m<all_nonuis.size(); m++)
    {
      python::list tmplist;
      for(n=0; n<all_nonuis[m].size(); n++)
      {
        tmplist.append(all_nonuis[m][n]);
      }
      result.append(tmplist);
    }

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

    def("get_non_uis", get_non_uis, 
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


    def("calculate_transitions_ch", _py_calculate_transitions_with_charge,
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

    def("calculate_charged_mass", SRMCollider::Common::calculate_charged_mass, 
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

