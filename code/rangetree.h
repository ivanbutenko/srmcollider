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

#ifndef SRMCOLLIDER_RANGETREE_H
#define SRMCOLLIDER_RANGETREE_H

// Boost.Python headers
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
namespace python = boost::python;

#include <boost/shared_ptr.hpp>

// Including headers from CGAL 
#include <CGAL/Cartesian.h>
#include <CGAL/Range_segment_tree_traits.h>
#include <CGAL/Range_tree_k.h>

#include <vector>

//include our own libraries
#include "srmcollider.h"
#include "srmcolliderLib.h"

namespace SRMCollider
{
  namespace SimpleRangetree
  {

    struct Precursor
    {
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

      void new_rangetree ();

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
      void create_tree(python::tuple pepids);

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
      python::list query_tree(double a, double b, double c, double d, int max_nr_isotopes, double correction);

      boost::shared_ptr<Range_tree_2_type> my_rangetree;
    };

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
    void create_tree(python::tuple pepids);

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
    python::list query_tree(double a, double b, double c, double d, int max_nr_isotopes, double correction);
  }

  namespace ExtendedRangetree
  {

    typedef SRMCollider::Common::SRMPrecursor Precursor;

    struct Transition
    {
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

      void new_rangetree ();

      // Create the rangetree that will be used throughout. This is essential
      void create_tree(python::tuple pepids);

      void create_tree(std::vector<Precursor> precursors);

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
      python::list query_tree(double a, double b, double c, double d, int max_nr_isotopes, double correction);
      std::vector< Precursor* > query_tree_c(double a, double b, double c, double d, int max_nr_isotopes, double correction);

      boost::shared_ptr<Range_tree_2_type> my_rangetree;
    };

  }
}
#endif
