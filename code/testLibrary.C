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
#include "py_srmcolliderLib.h"

#define EPS_05 boost::test_tools::fraction_tolerance(1.e-5) 

using namespace SRMCollider::Common;
using namespace std;

BOOST_AUTO_TEST_CASE( _calculate_fragment_masses_PEPTIDE_N14 )
{

  const char* sequence = "PEPTIDE";
  double* series = new double[1024];
  double* tmp_series = new double[1024];
  double ch = 2;
  SRMCollider::Common::SRMParameters params;
  int fragcount = calculate_fragment_masses(sequence, tmp_series, series, ch, params, NOISOTOPEMODIFICATION);

  BOOST_CHECK_EQUAL(fragcount, 12);
  BOOST_CHECK ( boost::test_tools::check_is_close( 352.161417374, series[0], EPS_05 ) ) ;
  BOOST_CHECK ( boost::test_tools::check_is_close( 287.640122374, series[1], EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close( 239.113742374, series[2], EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close( 188.589902374, series[3], EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close( 132.047872374, series[4], EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close( 74.5344023740, series[5], EPS_05 ) );

  BOOST_CHECK ( boost::test_tools::check_is_close( 49.53420503, series[6], EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close( 114.0555000, series[7], EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close( 162.5818800, series[8], EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close( 213.1057200, series[9], EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close( 269.6477500, series[10],EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close( 327.1612200, series[11],EPS_05 ) );
                                                                
}

BOOST_AUTO_TEST_CASE( _calculate_fragment_masses_PEPTIDE_N15 )
{

  const char* sequence = "PEPTIDE";
  double* series = new double[1024];
  double* tmp_series = new double[1024];
  double ch = 2;
  SRMCollider::Common::SRMParameters params;
  int fragcount = calculate_fragment_masses(sequence, tmp_series, series, ch, params, N15_ISOTOPEMODIFICATION);

  BOOST_CHECK_EQUAL(fragcount, 12);
  BOOST_CHECK ( boost::test_tools::check_is_close(355.153 , series[0], EPS_05 ) ) ;
  BOOST_CHECK ( boost::test_tools::check_is_close(290.133 , series[1], EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close(241.108 , series[2], EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close(190.085 , series[3], EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close(133.045 , series[4], EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close(75.0329 , series[5], EPS_05 ) );

  BOOST_CHECK ( boost::test_tools::check_is_close(50.0327 , series[6], EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close(115.053 , series[7], EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close(164.077 , series[8], EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close(215.1   , series[9], EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close(272.14  , series[10],EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close(330.152 , series[11],EPS_05 ) );

}

BOOST_AUTO_TEST_CASE( calculate_transitions_with_charge_three_peptide_test )
{
    /*
    The target is YYLLDYR with these transitions and numbers

      (842.4412197, 0), y6+
      (679.3778897, 1), y5+
      (566.2938297, 2), y4+
      (453.2097697, 3), y3+
      (440.2185450, 4), b3+
      (553.3026050, 5), b4+
      (668.3295450, 6), b5+
      (831.3928750, 7)  b6+ 

    */

  std::string sequence = "YYLLDYR";
  double* b_series = new double[256];
  double* y_series = new double[256];

  SRMPrecursor p;
  p.sequence = sequence;
  p.transition_group = 1;
  std::vector<SRMTransition> result;
  std::vector<int> charges;
  charges.push_back(1);
  double q3_low = 400;
  double q3_high = 1500;

  SRMParameters param;
  param.yions = true;
  param.bions = true;
  calculate_transitions_with_charge(p, charges, result, b_series, y_series, q3_low, q3_high, param);

  // y series
  BOOST_CHECK ( boost::test_tools::check_is_close( 842.4412197 , result[0].q3, EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close( 679.3778897 , result[1].q3, EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close( 566.2938297 , result[2].q3, EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close( 453.2097697 , result[3].q3, EPS_05 ) );

  // b series
  BOOST_CHECK ( boost::test_tools::check_is_close( 440.2185450 , result[4].q3, EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close( 553.3026050 , result[5].q3, EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close( 668.3295450 , result[6].q3, EPS_05 ) );
  BOOST_CHECK ( boost::test_tools::check_is_close( 831.3928750 , result[7].q3, EPS_05 ) );
                                                                
}

BOOST_AUTO_TEST_CASE( _py_calculate_charged_mass_TEST )
{

  Py_Initialize();
  double mass;
  python::tuple t = python::make_tuple(0, "YYLLDYR", 0);

  mass = SRMCollider::Common::_py_calculate_charged_mass(t, 1);
  BOOST_CHECK ( boost::test_tools::check_is_close( 1005.50460 , mass, EPS_05 ) );

  mass = SRMCollider::Common::_py_calculate_charged_mass(t, 2);
  BOOST_CHECK ( boost::test_tools::check_is_close( 503.25623 , mass, EPS_05 ) );

  mass = SRMCollider::Common::_py_calculate_charged_mass(t, 3);
  BOOST_CHECK ( boost::test_tools::check_is_close( 335.84011 , mass, EPS_05 ) );
                                                                
}

