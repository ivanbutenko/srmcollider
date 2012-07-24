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
#include "srmcolliderLib.h"

void calculate_eUIS(std::vector<int>& N, std::vector<std::vector<double> >& c_ssrcalcvalues,
    double ssrwindow, std::vector<std::vector<int> >& all_nonuis);

#include "calculate_eUIS.cpp"

#define EPS_05 boost::test_tools::fraction_tolerance(1.e-5) 
 
BOOST_AUTO_TEST_CASE( eUIS_test)
{

  /*
   *
   *
        self.ssrcalcvalues_four_example  = [  
            [ 1,   5,    8,   8.5], 
            [ 1.5,         8.2], 
            [ 2.1, 5.1,    8.1], 
            [1.0]  ]



  * The expected result
  *
  * [0, 3]
  * [3]
  * [1]
  * [2]
  * [0]
  * [0, 2]
  * [0, 1]
  * [0, 1, 2]
  * [1, 2]
  *
  */
  std::vector<std::vector<int> > all_nonuis;

  static const int arr1[] = {4, 2, 3, 1};
  std::vector<int> N (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );

  std::vector<std::vector<double> > c_ssrcalcvalues;
  std::vector<double> tmp;
  tmp.clear(); 
  tmp.push_back(1); tmp.push_back(5); tmp.push_back(8); tmp.push_back(8.5);
  c_ssrcalcvalues.push_back(tmp);
  tmp.clear(); 
  tmp.push_back(1.5); tmp.push_back(8.2);
  c_ssrcalcvalues.push_back(tmp);
  tmp.clear(); 
  tmp.push_back(2.1); tmp.push_back(5.1); tmp.push_back(8.1);
  c_ssrcalcvalues.push_back(tmp);
  tmp.clear(); 
  tmp.push_back(1); 
  c_ssrcalcvalues.push_back(tmp);

  double strike3_ssrcalcwindow = 0.3;
  calculate_eUIS(N, c_ssrcalcvalues, strike3_ssrcalcwindow, all_nonuis);

  BOOST_CHECK_EQUAL(all_nonuis.size(), 9);

  BOOST_CHECK_EQUAL(all_nonuis[0].size(), 2);
  BOOST_CHECK_EQUAL(all_nonuis[1].size(), 1);
  BOOST_CHECK_EQUAL(all_nonuis[2].size(), 1);
  BOOST_CHECK_EQUAL(all_nonuis[3].size(), 1);
  BOOST_CHECK_EQUAL(all_nonuis[4].size(), 1);
  BOOST_CHECK_EQUAL(all_nonuis[5].size(), 2);
  BOOST_CHECK_EQUAL(all_nonuis[6].size(), 2);
  BOOST_CHECK_EQUAL(all_nonuis[7].size(), 3);
  BOOST_CHECK_EQUAL(all_nonuis[8].size(), 2);

  BOOST_CHECK_EQUAL(all_nonuis[0][0], 0);
  BOOST_CHECK_EQUAL(all_nonuis[0][1], 3);

  BOOST_CHECK_EQUAL(all_nonuis[5][0], 0);
  BOOST_CHECK_EQUAL(all_nonuis[5][1], 2);

  BOOST_CHECK_EQUAL(all_nonuis[6][0], 0);
  BOOST_CHECK_EQUAL(all_nonuis[6][1], 1);

  BOOST_CHECK_EQUAL(all_nonuis[7][0], 0);
  BOOST_CHECK_EQUAL(all_nonuis[7][1], 1);
  BOOST_CHECK_EQUAL(all_nonuis[7][2], 2);

}

