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

#ifndef SRMCOLLIDE_EUIS_H
#define SRMCOLLIDE_EUIS_H

#include <vector>

// Boost.Python headers
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
namespace python = boost::python;

bool SortIntDoublePairSecond(const std::pair<int,double>& left, const std::pair<int,double>& right)
{
  return left.second < right.second;
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
              if(!(std::fabs(myssr[k] - myssr[i]) > ssrwindow))
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

#endif
