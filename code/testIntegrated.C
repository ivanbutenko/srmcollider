/*
 *
 * Program       : SRMCollider
 * Author        : Hannes Roest <roest@imsb.biol.ethz.ch>
 * Date          : 30.05.2012 
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MyTest
 
// Including headers from CGAL 
#include <CGAL/Cartesian.h>
#include <CGAL/Range_segment_tree_traits.h>
#include <CGAL/Range_tree_k.h>

#include <boost/test/unit_test.hpp>
#include <iostream>

//include our own libraries
#include "srmcollider.h"
#include "srmcolliderLib.h"
#include "integratedrun.h"
#include "rangetree.h"

#define EPS_05 boost::test_tools::fraction_tolerance(1.e-5) 
 
void test_helper_get_transitions(std::string sequence, std::vector<SRMCollider::IntegratedRun::Transition>& transitions)
{
  // get the transition
  std::vector<int> charges; charges.push_back(1);
  //char* sequence = (char*)"YYLLDYR";
  double* series = new double[1024];
  double* tmp_series = new double[1024];
  SRMPrecursor p;
  p.sequence = sequence;
  p.transition_group = 1;
  p.isotope_modification = 0;
  std::vector<SRMTransition> srm_transitions;
  SRMParameters tmp_param;
  calculate_transitions_with_charge(p, charges, srm_transitions, series, tmp_series, 400, 1500, tmp_param);
  for (size_t i = 0; i < srm_transitions.size(); i++) {
    SRMCollider::IntegratedRun::Transition t;
    t.q3 = srm_transitions[i].q3;
    t.transition_id = srm_transitions[i].transition_id;
    transitions.push_back(t);
  }
  delete tmp_series;
  delete series;
}

// this tests wrap_all_bitwise;
BOOST_AUTO_TEST_CASE( wrap_all_bitwise_3peptides_TEST)
{
  /*
  """
  The target is YYLLDYR with these transitions and numbers

    (842.4412197, 0), y6+
    (679.3778897, 1), y5+
    (566.2938297, 2), y4+
    (453.2097697, 3), y3+
    (440.2185450, 4), b3+
    (553.3026050, 5), b4+
    (668.3295450, 6), b5+
    (831.3928750, 7)  b6+ 

  The peptides GGLIVELGDK b5+ ion interferes with the targets b3+ ion which leads to 665: [4]

  The peptides NGTDGGLQVAIDAMR b9+ ion (842.4008) interferes with the targets y6+ ion
  and also the y11++ ion (565.8035) interferes with the targets y4+ ion which leads to 618: [0, 2].

              sequence = python::extract<char *>(tlist[0]);
              peptide_key = python::extract<long>(tlist[1]);
              parent_id = python::extract<long>(tlist[2]);
              q1_charge = python::extract<int>(tlist[3]);

              q1 = python::extract<double>(tlist[4]);
              ssrcalc = python::extract<double>(tlist[5]);
              isotope_modification = python::extract<double>(tlist[8]);
  """
  */

  double q3_window = 1.0 / 2.0;
  double ppm = false;
  double isotopes_up_to = 3;
  double isotope_correction = 1.0;
  int peptide_key = 3;
  int max_uis = 5;
  SRMCollider::Common::SRMParameters param;
  param.ppm = ppm;
  param.q3_window = q3_window;

  Py_Initialize();
  boost::python::tuple t1 = boost::python::make_tuple("GGLIVELGDK", 1, 1, 2, 500.0, 25.0, -1,-1, 0);
  boost::python::tuple t2 = boost::python::make_tuple("NGTDGGLQVAIDAMR", 2, 2, 2, 500.0, 25.0, -1,-1, 0);
  boost::python::tuple t3 = boost::python::make_tuple("YYLLDYR", 3, 3, 2, 500.0, 25.0, -1,-1, 0);
  boost::python::tuple alltuples = boost::python::make_tuple(t1,t2,t3);
  SRMCollider::ExtendedRangetree::Rangetree_Q1_RT tree;
  tree.new_rangetree();
  tree.create_tree(alltuples);

  // get the transition
  std::vector<SRMCollider::IntegratedRun::Transition> transitions;
  std::string sequence = "YYLLDYR";
  test_helper_get_transitions(sequence, transitions);
  BOOST_CHECK_EQUAL(transitions.size(), 8);

  // get the collisions per peptide
  std::vector<COMBINT> collisions_per_peptide; 
  SRMCollider::IntegratedRun::wrap_all_bitwise(transitions, 499, 24.0, 501, 26.0, 
      peptide_key, max_uis, q3_window, ppm, isotopes_up_to, isotope_correction, param, tree,
      collisions_per_peptide);

  // the interference with transition 4 leads to 2^^4 = 16
  // and the interference with transition 0 and 2 leads to 2^^0 + 2^^2 = 5
  static const int arr1[] = {5, 16};
  std::vector<int> cmp_collisions_per_peptide (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  std::sort(collisions_per_peptide.begin(), collisions_per_peptide.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(collisions_per_peptide.begin(), collisions_per_peptide.end(),
      cmp_collisions_per_peptide.begin(), cmp_collisions_per_peptide.end());

}

// this tests get_non_uis_bitwise(std::vector<COMBINT>, int, int, std::vector<COMBINT>);
BOOST_AUTO_TEST_CASE( get_non_uis_bitwise_3peptides_TEST)
{
  /*
  """
  The target is YYLLDYR with these transitions and numbers

    (842.4412197, 0), y6+
    (679.3778897, 1), y5+
    (566.2938297, 2), y4+
    (453.2097697, 3), y3+
    (440.2185450, 4), b3+
    (553.3026050, 5), b4+
    (668.3295450, 6), b5+
    (831.3928750, 7)  b6+ 

  The peptides GGLIVELGDK b5+ ion interferes with the targets b3+ ion which leads to 665: [4]

  The peptides NGTDGGLQVAIDAMR b9+ ion (842.4008) interferes with the targets y6+ ion
  and also the y11++ ion (565.8035) interferes with the targets y4+ ion which leads to 618: [0, 2].

              sequence = python::extract<char *>(tlist[0]);
              peptide_key = python::extract<long>(tlist[1]);
              parent_id = python::extract<long>(tlist[2]);
              q1_charge = python::extract<int>(tlist[3]);

              q1 = python::extract<double>(tlist[4]);
              ssrcalc = python::extract<double>(tlist[5]);
              isotope_modification = python::extract<double>(tlist[8]);
  """
  */

  double q3_window = 1.0 / 2.0;
  double ppm = false;
  double isotopes_up_to = 3;
  double isotope_correction = 1.0;
  int peptide_key = 3;
  int max_uis = 5;
  SRMCollider::Common::SRMParameters param;
  param.ppm = ppm;
  param.q3_window = q3_window;

  Py_Initialize();
  boost::python::tuple t1 = boost::python::make_tuple("GGLIVELGDK", 1, 1, 2, 500.0, 25.0, -1,-1, 0);
  boost::python::tuple t2 = boost::python::make_tuple("NGTDGGLQVAIDAMR", 2, 2, 2, 500.0, 25.0, -1,-1, 0);
  boost::python::tuple t3 = boost::python::make_tuple("YYLLDYR", 3, 3, 2, 500.0, 25.0, -1,-1, 0);
  boost::python::tuple alltuples = boost::python::make_tuple(t1,t2,t3);
  SRMCollider::ExtendedRangetree::Rangetree_Q1_RT tree;
  tree.new_rangetree();
  tree.create_tree(alltuples);

  // get the transition
  std::vector<SRMCollider::IntegratedRun::Transition> transitions;
  std::string sequence = "YYLLDYR";
  test_helper_get_transitions(sequence, transitions);
  BOOST_CHECK_EQUAL(transitions.size(), 8);

  // get the collisions per peptide
  std::vector<COMBINT> collisions_per_peptide; 
  SRMCollider::IntegratedRun::wrap_all_bitwise(transitions, 499, 24.0, 501, 26.0, 
      peptide_key, max_uis, q3_window, ppm, isotopes_up_to, isotope_correction, param, tree,
      collisions_per_peptide);

  // reduce them to a specific order
  std::vector<int> result;
  for(int i =1; i<= max_uis; i++) 
  {
    std::set<COMBINT> combinations;
    SRMCollider::Combinatorics::get_non_uis_bitwise(collisions_per_peptide, transitions.size(), i, combinations);
    result.push_back(combinations.size());
  }

  static const int arr2[] = {3,1,0,0,0};
  std::vector<int> cmp_result (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

  BOOST_CHECK_EQUAL_COLLECTIONS(result.begin(), result.end(), cmp_result.begin(), cmp_result.end());
}

BOOST_AUTO_TEST_CASE( get_min_needed_3peptides_TEST)
{

  double q3_window = 1.0 / 2.0;
  double ppm = false;
  //double isotopes_up_to = 3;
  //double isotope_correction = 1.0;
  //int peptide_key = 3;
  int max_uis = 5;
  SRMCollider::Common::SRMParameters param;
  param.ppm = ppm;
  param.q3_window = q3_window;

  Py_Initialize();
  SRMCollider::Common::SRMPrecursor prec;
  std::vector<SRMCollider::Common::SRMPrecursor> precursors;
  prec.sequence = "GGLIVELGDK";
  prec.transition_group = prec.q1 = prec.maximal_charge = prec.ssrcalc = -1;
  precursors.push_back(prec);
  prec.sequence = "NGTDGGLQVAIDAMR";
  prec.transition_group = prec.q1 = prec.maximal_charge = prec.ssrcalc = -1;
  precursors.push_back(prec);

  std::vector<SRMCollider::IntegratedRun::Transition> transitions;
  std::string sequence = "YYLLDYR";
  test_helper_get_transitions(sequence, transitions);
  BOOST_CHECK_EQUAL(transitions.size(), 8);

  // get the collisions per peptide
  std::vector<COMBINT> collisions_per_peptide; 
  int min_tr = SRMCollider::IntegratedRun::min_needed(transitions, precursors, max_uis, q3_window, ppm, param);
  BOOST_CHECK_EQUAL(min_tr, 2);

  // Now lets give the transitions in a different order...we should now need 3 transitions to measure this peptide
  std::vector<SRMCollider::IntegratedRun::Transition> new_transitions;
  new_transitions.push_back(transitions[0]);
  new_transitions.push_back(transitions[2]);
  new_transitions.push_back(transitions[4]);
  new_transitions.push_back(transitions[1]);
  new_transitions.push_back(transitions[3]);
  new_transitions.push_back(transitions[5]);
  // get the collisions per peptide
  collisions_per_peptide.clear(); 
  min_tr = SRMCollider::IntegratedRun::min_needed(new_transitions, precursors, max_uis, q3_window, ppm, param);
  BOOST_CHECK_EQUAL(min_tr, 3);

}

BOOST_AUTO_TEST_CASE( wrap_all_bitwise_TEST)
{
}

