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
#include <boost/test/unit_test.hpp>
#define EPS_05 boost::test_tools::fraction_tolerance(1.e-5) 

#include "rangetree.h"

// this tests the results returned by the tree
void test_tree_c(SRMCollider::ExtendedRangetree::Rangetree_Q1_RT tree)
{

  int nr_isotopes = 0;
  std::vector< SRMCollider::ExtendedRangetree::Precursor* > result;

  result = tree.query_tree_c(490, 20, 510, 30, nr_isotopes, 10);
  BOOST_CHECK_EQUAL(result.size(), 3);
  // also check that the pointers are not null
  BOOST_CHECK(result[0] != NULL);
  BOOST_CHECK(result[1] != NULL);
  BOOST_CHECK(result[2] != NULL);

  result = tree.query_tree_c(500.9, 20, 510, 30, nr_isotopes, 10);
  BOOST_CHECK_EQUAL(result.size(), 2);

  result = tree.query_tree_c(501.5, 20, 510, 30, nr_isotopes, 10);
  BOOST_CHECK_EQUAL(result.size(), 1);

  // now with isotopes, all should fall in the window again
  nr_isotopes = 3;
  result = tree.query_tree_c(500.9, 20, 510, 30, nr_isotopes, 10);
  BOOST_CHECK_EQUAL(result.size(), 3);

  result = tree.query_tree_c(501.5, 20, 510, 30, nr_isotopes, 10);
  BOOST_CHECK_EQUAL(result.size(), 3);

  result = tree.query_tree_c(501.6, 20, 510, 30, nr_isotopes, 10);
  BOOST_CHECK_EQUAL(result.size(), 2);

  // with 2 isotopes, the GGLIVELGDK will produce 500.0, 500.5, 501.0 
  nr_isotopes = 2;
  result = tree.query_tree_c(501.8, 20, 510, 30, nr_isotopes, 10);
  BOOST_CHECK_EQUAL(result.size(), 2);

  result = tree.query_tree_c(500.9, 20, 510, 30, nr_isotopes, 10);
  BOOST_CHECK_EQUAL(result.size(), 3);
}

// this tests ExtendedRangetree using C code only
BOOST_AUTO_TEST_CASE( ExtendedRangetree_TEST_C )
{
  int q1_charge = 2;

  SRMCollider::ExtendedRangetree::Precursor p1("GGLIVELGDK",      500.0, 1, 1, q1_charge, 0, 25.0);
  SRMCollider::ExtendedRangetree::Precursor p2("NGTDGGLQVAIDAMR", 501.0, 2, 2, q1_charge, 0, 24.0);
  SRMCollider::ExtendedRangetree::Precursor p3("YYLLDYR",         502.0, 3, 3, q1_charge, 0, 23.0);
  std::vector< SRMCollider::ExtendedRangetree::Precursor > precursors; 
  precursors.push_back(p1);
  precursors.push_back(p2);
  precursors.push_back(p3);
  SRMCollider::ExtendedRangetree::Rangetree_Q1_RT tree;
  tree.new_rangetree();
  tree.create_tree(precursors);
  test_tree_c(tree);
}

// this tests the results returned by the tree, using the python wrapper for
// querying the tree
template <class RTreeType>
void test_tree(RTreeType tree)
{
  int nr_isotopes = 0;
  int length = 5;
  python::list result;

  result = tree.query_tree(490, 20, 510, 30, nr_isotopes, 10);
  length = python::extract<int>(result.attr("__len__")());
  BOOST_CHECK_EQUAL(length, 3);

  result = tree.query_tree(500.9, 20, 510, 30, nr_isotopes, 10);
  length = python::extract<int>(result.attr("__len__")());
  BOOST_CHECK_EQUAL(length, 2);

  result = tree.query_tree(501.5, 20, 510, 30, nr_isotopes, 10);
  length = python::extract<int>(result.attr("__len__")());
  BOOST_CHECK_EQUAL(length, 1);

  // now with isotopes, all should fall in the window again
  nr_isotopes = 3;
  result = tree.query_tree(500.9, 20, 510, 30, nr_isotopes, 10);
  length = python::extract<int>(result.attr("__len__")());
  BOOST_CHECK_EQUAL(length, 3);

  result = tree.query_tree(501.5, 20, 510, 30, nr_isotopes, 10);
  length = python::extract<int>(result.attr("__len__")());
  BOOST_CHECK_EQUAL(length, 3);

  result = tree.query_tree(501.6, 20, 510, 30, nr_isotopes, 10);
  length = python::extract<int>(result.attr("__len__")());
  BOOST_CHECK_EQUAL(length, 2);

  // with 2 isotopes, the GGLIVELGDK will produce 500.0, 500.5, 501.0 
  nr_isotopes = 2;
  result = tree.query_tree(501.8, 20, 510, 30, nr_isotopes, 10);
  length = python::extract<int>(result.attr("__len__")());
  BOOST_CHECK_EQUAL(length, 2);

  result = tree.query_tree(500.9, 20, 510, 30, nr_isotopes, 10);
  length = python::extract<int>(result.attr("__len__")());
  BOOST_CHECK_EQUAL(length, 3);
}

// this tests ExtendedRangetree;
BOOST_AUTO_TEST_CASE( ExtendedRangetree_TEST )
{
  int q1_charge = 2;
  Py_Initialize();
  boost::python::tuple t1 = boost::python::make_tuple("GGLIVELGDK", 1, 1, q1_charge,      500.0, 25.0, -1,-1, 0);
  boost::python::tuple t2 = boost::python::make_tuple("NGTDGGLQVAIDAMR", 2, 2, q1_charge, 501.0, 24.0, -1,-1, 0);
  boost::python::tuple t3 = boost::python::make_tuple("YYLLDYR", 3, 3, q1_charge,         502.0, 23.0, -1,-1, 0);
  boost::python::tuple alltuples = boost::python::make_tuple(t1,t2,t3);
  SRMCollider::ExtendedRangetree::Rangetree_Q1_RT tree;
  tree.new_rangetree();
  tree.create_tree(alltuples);
  test_tree(tree);
  test_tree_c(tree);
}

BOOST_AUTO_TEST_CASE( SimpleRangetree_TEST )
{
  int q1_charge = 2;

  Py_Initialize();
  boost::python::tuple t1 = boost::python::make_tuple("GGLIVELGDK", 1, 1, q1_charge,      500.0, 25.0, -1,-1, 0);
  boost::python::tuple t2 = boost::python::make_tuple("NGTDGGLQVAIDAMR", 2, 2, q1_charge, 501.0, 24.0, -1,-1, 0);
  boost::python::tuple t3 = boost::python::make_tuple("YYLLDYR", 3, 3, q1_charge,         502.0, 23.0, -1,-1, 0);
  boost::python::tuple alltuples = boost::python::make_tuple(t1,t2,t3);
  SRMCollider::SimpleRangetree::Rangetree_Q1_RT tree;
  tree.new_rangetree();
  tree.create_tree(alltuples);
  test_tree(tree);
  //test_tree_c(tree); // doesnt exist yet for simple range tree
}

