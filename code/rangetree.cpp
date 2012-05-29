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
#include "srmcollider.h"
#include <boost/shared_ptr.hpp>


// Including headers from CGAL 
#include <CGAL/Cartesian.h>
#include <CGAL/Range_segment_tree_traits.h>
#include <CGAL/Range_tree_k.h>

namespace SRMCollider
{
  namespace SimpleRangetree
  {

struct Precursor{
  long parent_id;
  int q1_charge;
};


typedef CGAL::Cartesian<double> K;
typedef CGAL::Range_tree_map_traits_2<K, Precursor> Traits;
typedef CGAL::Range_tree_2<Traits> Range_tree_2_type;

typedef Traits::Key Key;                
typedef Traits::Interval Interval;    

struct Rangetree_Q1_RT 
{
  static boost::shared_ptr<Rangetree_Q1_RT> create () 
  { 
    return boost::shared_ptr<Rangetree_Q1_RT>(new Rangetree_Q1_RT); 
  }

  void new_rangetree  () 
  {
    my_rangetree = boost::shared_ptr<Range_tree_2_type>(new Range_tree_2_type); 
  }

  /* Create the rangetree that will be used throughout. This is essential. The
   * rangetree will stay in place while this module is loaded.
   * The tuples have the following structure:
   *   0
   *   1 
   *   2 - parent_id
   *   3 - q1_charge
   *   4 - q1
   *   5 - ssrcalc
  */
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
      my_rangetree->make_tree(InputList.begin(),InputList.end());
  }

  /* Query the rangetree. Format is (x1,y1,x2,y2), returns all entries that are
   * in the defined square defined by these four numbers.
   * Returns a list with keys that were stored in the tree.
   *
   * Requires the maximal number of isotopes to consider and a isotopic correction of the lower window.
   * The isotope correction of the lower window is an offset that is used to
   * include all the monoisotopic precursors that might have isotopes between x1
   * and x2. Since the isotopes are not actively stored, the monoisotopic
   * precursors have to be selected and then checked whether they have any
   * potential isotopes that could fall in the (x1,x2) window.
   * The isotope correction should be computed as:
   *    nr_isotopes_to_consider * mass_difference_of_C13 / minimal_parent_charge
  */
  python::list query_tree(double a, double b, double c, double d, int max_nr_isotopes, double correction)   
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

  boost::shared_ptr<Range_tree_2_type> my_rangetree;
};


Range_tree_2_type *Range_tree_2 = new Range_tree_2_type;

// Function declarations
void create_tree(python::tuple pepids) ;
python::list query_tree(double a, double b, double c, double d, int max_nr_isotopes, double correction);


/* Create the rangetree that will be used throughout. This is essential. The
 * rangetree will stay in place while this module is loaded.
 * The tuples have the following structure:
 *   0
 *   1 
 *   2 - parent_id
 *   3 - q1_charge
 *   4 - q1
 *   5 - ssrcalc
*/
void create_tree(python::tuple pepids) {

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

/* Query the rangetree. Format is (x1,y1,x2,y2), returns all entries that are
 * in the defined square defined by these four numbers.
 * Returns a list with keys that were stored in the tree.
 *
 * Requires the maximal number of isotopes to consider and a isotopic correction of the lower window.
 * The isotope correction of the lower window is an offset that is used to
 * include all the monoisotopic precursors that might have isotopes between x1
 * and x2. Since the isotopes are not actively stored, the monoisotopic
 * precursors have to be selected and then checked whether they have any
 * potential isotopes that could fall in the (x1,x2) window.
 * The isotope correction should be computed as:
 *    nr_isotopes_to_consider * mass_difference_of_C13 / minimal_parent_charge
*/
python::list query_tree(double a, double b, double c, double d, int max_nr_isotopes, double correction)   {

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

void pass_rtree(Rangetree_Q1_RT& t)
{
   std::cout << "got a tree " << std::endl;
   t.query_tree(2,1,0,0,0,1);
   std::cout << "got a tree " << std::endl;
}


  }
  namespace ExtendedRangetree
  {

    struct Precursor{
      char* sequence; 
      long peptide_key;
      long parent_id;
      int q1_charge;
      int isotope_modification;
    };

    struct Transition{
      double q3;
      long srm_id;
    };

    typedef CGAL::Cartesian<double> K;
    typedef CGAL::Range_tree_map_traits_2<K, Precursor> Traits;
    typedef CGAL::Range_tree_2<Traits> Range_tree_2_type;

    typedef Traits::Key Key;                
    typedef Traits::Interval Interval;    

    struct Rangetree_Q1_RT 
    {
      static boost::shared_ptr<Rangetree_Q1_RT> create () 
      { 
        return boost::shared_ptr<Rangetree_Q1_RT>(new Rangetree_Q1_RT); 
      }

      void new_rangetree  () 
      {
        my_rangetree = boost::shared_ptr<Range_tree_2_type>(new Range_tree_2_type); 
      }

      // Create the rangetree that will be used throughout. This is essential
      void create_tree(python::tuple pepids) {

          python::tuple tlist;
          std::vector<Key> InputList;
          int i, q1_charge, isotope_modification;
          long parent_id, peptide_key;
          double ssrcalc, q1;
          char* sequence;

          int pepids_length = python::extract<int>(pepids.attr("__len__")());
          for (i=0; i<pepids_length; i++) 
          {
              tlist = python::extract< python::tuple >(pepids[i]);

              sequence = python::extract<char *>(tlist[0]);
              peptide_key = python::extract<long>(tlist[1]);
              parent_id = python::extract<long>(tlist[2]);
              q1_charge = python::extract<int>(tlist[3]);

              q1 = python::extract<double>(tlist[4]);
              ssrcalc = python::extract<double>(tlist[5]);
              isotope_modification = python::extract<double>(tlist[8]);

              struct Precursor entry = {sequence, peptide_key, parent_id, q1_charge, isotope_modification};
              InputList.push_back(Key(K::Point_2(q1,ssrcalc), entry));
          }
          my_rangetree->make_tree(InputList.begin(),InputList.end());
      }

      /* Query the rangetree. Format is (x1,y1,x2,y2), returns all entries that are
       * in the defined square defined by these four numbers.
       * Returns a list with keys that were stored in the tree.
       *
       * Requires the maximal number of isotopes to consider and a isotopic correction of the lower window.
       * The isotope correction of the lower window is an offset that is used to
       * include all the monoisotopic precursors that might have isotopes between x1
       * and x2. Since the isotopes are not actively stored, the monoisotopic
       * precursors have to be selected and then checked whether they have any
       * potential isotopes that could fall in the (x1,x2) window.
       * The isotope correction should be computed as:
       *    nr_isotopes_to_consider * mass_difference_of_C13 / minimal_parent_charge
      */
      python::list query_tree(double a, double b, double c, double d, int max_nr_isotopes, double correction)   
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

      boost::shared_ptr<Range_tree_2_type> my_rangetree;
    };

  }
}


// Expose to Python
using namespace python;
BOOST_PYTHON_MODULE(c_rangetree)
{

    def("create_tree", SRMCollider::SimpleRangetree::create_tree, 
 "Query the rangetree. Format is (x1,y1,x2,y2), returns all entries that are \n"
 "in the defined square defined by these four numbers. \n"
 "Returns a list with keys that were stored in the tree. \n"
 "\n"
 "\n"
 " Signature\n"
 "list query_tree(double a, double b, double c, double d)\n"
            "");
    def("query_tree", SRMCollider::SimpleRangetree::query_tree,
            
 "Create the rangetree that will be used throughout. This is essential. The \n"
 "rangetree will stay in place while this module is loaded. \n"
 "\n"
 "\n"
 " Signature\n"
 "void create_tree(tuple pepids) \n"
            "");

class_<SRMCollider::SimpleRangetree::Rangetree_Q1_RT,
  boost::shared_ptr<SRMCollider::SimpleRangetree::Rangetree_Q1_RT> >("Rangetree_Q1_RT",init<>())
        .def("create",&SRMCollider::SimpleRangetree::Rangetree_Q1_RT::create )
        .staticmethod("create")
        .def("new_rangetree",&SRMCollider::SimpleRangetree::Rangetree_Q1_RT::new_rangetree)
        .def("create_tree",&SRMCollider::SimpleRangetree::Rangetree_Q1_RT::create_tree)
        .def("query_tree",&SRMCollider::SimpleRangetree::Rangetree_Q1_RT::query_tree)
    ;

class_<SRMCollider::ExtendedRangetree::Rangetree_Q1_RT,
  boost::shared_ptr<SRMCollider::ExtendedRangetree::Rangetree_Q1_RT> >("ExtendedRangetree_Q1_RT",init<>())
        .def("create",&SRMCollider::ExtendedRangetree::Rangetree_Q1_RT::create )
        .staticmethod("create")
        .def("new_rangetree",&SRMCollider::ExtendedRangetree::Rangetree_Q1_RT::new_rangetree)
        .def("create_tree",&SRMCollider::ExtendedRangetree::Rangetree_Q1_RT::create_tree)
        .def("query_tree",&SRMCollider::ExtendedRangetree::Rangetree_Q1_RT::query_tree)
    ;

  def("pass_rtree", SRMCollider::SimpleRangetree::pass_rtree);

}
