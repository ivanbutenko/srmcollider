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
#include <iostream>

//include our own libraries
#include "srmcollider.h"
#include "combinatorics.h"

#define EPS_05 boost::test_tools::fraction_tolerance(1.e-5) 

using namespace SRMCollider::Combinatorics;
 
BOOST_AUTO_TEST_CASE( _py_combinations_TEST)
{
  Py_Initialize();
  python::dict result;
  python::list mapping;
  for (int i = 0; i < 5; i++)
  {
    mapping.append(i);
  }
  _py_combinations(2, 5, mapping, result);

  python::list result_list = python::extract<python::list>(result.attr("keys")());
  int result_l = python::extract<int>(result_list.attr("__len__")());

  BOOST_CHECK_EQUAL(result_l, 10);
}

BOOST_AUTO_TEST_CASE( _combinations_TEST_1)
{
  Py_Initialize();
  std::vector<std::vector<int> > result;
  std::vector<std::vector<int> > cmp_result;
  int select_nr = 2;
  int collection_size = 5;
  _combinations(select_nr, collection_size, result);

  for (int i = 0; i < collection_size; i++)
  {
    for (int j = i+1; j < collection_size; j++)
    {
      std::vector<int> tmp;
      tmp.push_back(i);
      tmp.push_back(j);
      cmp_result.push_back(tmp);
    }
  }

  for (int i = 0; i < collection_size; i++)
  {
    BOOST_CHECK_EQUAL_COLLECTIONS(result[i].begin(), result[i].end(), cmp_result[i].begin(), cmp_result[i].end());
  }
}

BOOST_AUTO_TEST_CASE( _combinations_TEST_2)
{
  Py_Initialize();
  std::vector<std::vector<int> > result;
  std::vector<std::vector<int> > cmp_result;
  int select_nr = 3;
  int collection_size = 15;
  _combinations(select_nr, collection_size, result);

  for (int i = 0; i < collection_size; i++)
  {
    for (int j = i+1; j < collection_size; j++)
    {
      for (int k = j+1; k < collection_size; k++)
      {
        std::vector<int> tmp;
        tmp.push_back(i);
        tmp.push_back(j);
        tmp.push_back(k);
        cmp_result.push_back(tmp);
      }
    }
  }

  for (int i = 0; i < collection_size; i++)
  {
    BOOST_CHECK_EQUAL_COLLECTIONS(result[i].begin(), result[i].end(), cmp_result[i].begin(), cmp_result[i].end());
  }
}

BOOST_AUTO_TEST_CASE( _combinations_bitwise_TEST)
{
}
