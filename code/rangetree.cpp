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
 * This file has the abstraction layer to make calls to the CGAL rangetree and
 * retrieve precursor values from it. There are two functions: build tree and
 * query tree.
 *
 *
 * To install CGAL:
 *
 * sudo apt-get install libcgal-dev 
 * g++ -frounding-math -lCGAL -O3 rangetree.cpp 
 * http://www.cgal.org/Manual/latest/doc_html/cgal_manual/SearchStructures/Chapter_main.html
 * http://graphics.stanford.edu/courses/cs368-00-spring/TA/manuals/CGAL/ref-manual2/SearchStructures/Chapter_main.html
 * Q Public Licence http://www.gnu.org/licenses/license-list.html
 * http://www.cgal.org/Manual/latest/doc_html/cgal_manual/packages.html#part_XIII
*/

#ifndef SRMCOLLIDER_RANGETREE_C
#define SRMCOLLIDER_RANGETREE_C

//include our own libraries
#include "rangetree.h"

namespace SRMCollider
{
  namespace SimpleRangetree
  {
    void Rangetree_Q1_RT::new_rangetree  () 
    {
      my_rangetree = boost::shared_ptr<Range_tree_2_type>(new Range_tree_2_type); 
    }

    void Rangetree_Q1_RT::create_tree(python::tuple pepids) 
    {

        python::tuple tlist;
        std::vector<Key> InputList;
        int i, q1_charge;
        long parent_id;
        double ssrcalc, q1;

        int pepids_length = python::extract<int>(pepids.attr("__len__")());
        for (i=0; i<pepids_length; i++) {
            tlist = python::extract< python::tuple >(pepids[i]);

            parent_id = python::extract<long>(tlist[2]);
            q1_charge = python::extract<int>(tlist[3]);

            q1 = python::extract<double>(tlist[4]);
            ssrcalc = python::extract<double>(tlist[5]);

            struct Precursor entry = {parent_id, q1_charge};
            InputList.push_back(Key(K::Point_2(q1,ssrcalc), entry));
        }
        my_rangetree->make_tree(InputList.begin(),InputList.end());
    }

    python::list Rangetree_Q1_RT::query_tree(double a, double b, double c, double d, int max_nr_isotopes, double correction)   
    {
      std::vector<Key> OutputList;
      python::list result;

      double q1, q1_low = a, q1_high = c;
      int charge, iso;
      bool proceed;
      Interval win(Interval(K::Point_2(a-correction,b),K::Point_2(c,d)));
      my_rangetree->window_query(win, std::back_inserter(OutputList));
      std::vector<Key>::iterator current=OutputList.begin();
      while(current!=OutputList.end()){

          q1 = current->first[0];
          charge = current->second.q1_charge;
          proceed = false;
          // check whether there are any relevant isotopes, otherwise exclude this hit
          for (iso=0; iso<=max_nr_isotopes; iso++) {
              if (q1 + (MASS_DIFFC13 * iso)/charge > q1_low && 
                    q1 + (MASS_DIFFC13 * iso)/charge < q1_high) 
              {
                  proceed = true;
              }
          }

          if(proceed) {result.append(python::make_tuple( (*current).second.parent_id));}
          current++;
      }
      return result;

    }

    Range_tree_2_type *Range_tree_2 = new Range_tree_2_type;

    void create_tree(python::tuple pepids) 
    {

        python::tuple tlist;
        std::vector<Key> InputList;
        int i, q1_charge;
        long parent_id;
        double ssrcalc, q1;

        int pepids_length = python::extract<int>(pepids.attr("__len__")());
        for (i=0; i<pepids_length; i++) {
            tlist = python::extract< python::tuple >(pepids[i]);

            parent_id = python::extract<long>(tlist[2]);
            q1_charge = python::extract<int>(tlist[3]);

            q1 = python::extract<double>(tlist[4]);
            ssrcalc = python::extract<double>(tlist[5]);

            struct Precursor entry = {parent_id, q1_charge};
            InputList.push_back(Key(K::Point_2(q1,ssrcalc), entry));
        }
        Range_tree_2->make_tree(InputList.begin(),InputList.end());
    }

    python::list query_tree(double a, double b, double c, double d, int max_nr_isotopes, double correction)   
    {

      std::vector<Key> OutputList;
      python::list result;

      double q1, q1_low = a, q1_high = c;
      int charge, iso;
      bool proceed;
      Interval win(Interval(K::Point_2(a-correction,b),K::Point_2(c,d)));
      Range_tree_2->window_query(win, std::back_inserter(OutputList));
      std::vector<Key>::iterator current=OutputList.begin();
      while(current!=OutputList.end()){

          q1 = current->first[0];
          charge = current->second.q1_charge;
          proceed = false;
          // check whether there are any relevant isotopes, otherwise exclude this hit
          for (iso=0; iso<=max_nr_isotopes; iso++) {
              if (q1 + (MASS_DIFFC13 * iso)/charge > q1_low && 
                    q1 + (MASS_DIFFC13 * iso)/charge < q1_high) 
              {
                  proceed = true;
              }
          }

          if(proceed) {result.append(python::make_tuple( (*current).second.parent_id));}
          current++;
      }
      return result;

    }

  }

  namespace ExtendedRangetree
  {
    void Rangetree_Q1_RT::new_rangetree  () 
    {
      my_rangetree = boost::shared_ptr<Range_tree_2_type>(new Range_tree_2_type); 
    }

    void Rangetree_Q1_RT::create_tree(python::tuple pepids) 
    {

        python::tuple tlist;
        std::vector<Key> InputList;
        int i, q1_charge, isotope_modification;
        long parent_id, peptide_key;
        double ssrcalc, q1;
        std::string sequence;

        int pepids_length = python::extract<int>(pepids.attr("__len__")());
        for (i=0; i<pepids_length; i++) 
        {
            tlist = python::extract< python::tuple >(pepids[i]);

            sequence = python::extract<std::string>(tlist[0]);
            peptide_key = python::extract<long>(tlist[1]); // transition_group
            parent_id = python::extract<long>(tlist[2]);
            q1_charge = python::extract<int>(tlist[3]);

            q1 = python::extract<double>(tlist[4]);
            ssrcalc = python::extract<double>(tlist[5]);
            isotope_modification = python::extract<double>(tlist[8]);

            Precursor entry;
            entry.sequence = sequence;
            entry.transition_group = peptide_key; 
            entry.parent_id = parent_id;
            entry.q1_charge = q1_charge;
            //entry.q1 = q1;
            //entry.ssrcalc = ssrcalc; 
            entry.isotope_modification = isotope_modification; //= {sequence, peptide_key, parent_id, q1_charge, isotope_modification};
            InputList.push_back(Key(K::Point_2(q1,ssrcalc), entry));
        }
        my_rangetree->make_tree(InputList.begin(),InputList.end());
    }

    void Rangetree_Q1_RT::create_tree(std::vector<Precursor> precursors) 
    {
        std::vector<Key> InputList;
        for (size_t i=0; i<precursors.size(); i++) 
        {
            Precursor& entry = precursors[i];
            InputList.push_back(Key(K::Point_2(entry.q1,entry.ssrcalc), entry));
        }
        my_rangetree->make_tree(InputList.begin(),InputList.end());
    }

    python::list Rangetree_Q1_RT::query_tree(double a, double b, double c, double d, int max_nr_isotopes, double correction)
    {
      python::list result;
      std::vector< Precursor* > c_res = query_tree_c(a, b, c, d, max_nr_isotopes, correction);
      for (std::vector< Precursor* >::iterator it = c_res.begin(); it != c_res.end(); it++)
      {
        Precursor* p = *it;
        result.append(python::make_tuple(p->parent_id));
      }
      return result;
    }

    std::vector< Precursor* > Rangetree_Q1_RT::query_tree_c(double a, double b, double c, double d, int max_nr_isotopes, double correction)
    {
      std::vector<Key> OutputList;
      std::vector< Precursor* > result;

      int charge, iso;
      double q1, q1_low = a, q1_high = c;
      bool proceed;
      Interval win(Interval(K::Point_2(a-correction,b),K::Point_2(c,d)));
      my_rangetree->window_query(win, std::back_inserter(OutputList));
      std::vector<Key>::iterator current=OutputList.begin();
      while(current!=OutputList.end()){

          q1 = current->first[0];
          charge = current->second.q1_charge;
          proceed = false;
          // check whether there are any relevant isotopes, otherwise exclude this hit
          for (iso=0; iso<=max_nr_isotopes; iso++) {
              if (q1 + (MASS_DIFFC13 * iso)/charge > q1_low && 
                    q1 + (MASS_DIFFC13 * iso)/charge < q1_high) 
              {
                  proceed = true;
              }
          }

          if(proceed) {result.push_back( &current->second);}
          current++;
      }
      return result;
    }

  }
}

#endif
