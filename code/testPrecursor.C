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

#include <iostream>

//include our own libraries
#include "srmcollider.h"
#include "srmcolliderLib.h"

using namespace SRMCollider::Common;
using namespace std;

BOOST_AUTO_TEST_CASE( _calculate_fragment_masses_PEPTIDE_N14 )
{

  std::vector<double> series;
  SRMPrecursor p;
  p.sequence = "PEPTIDE";
  p.isotope_modification = NOISOTOPEMODIFICATION;
  p.q1_charge = 2;
  double ch = 2;
  SRMCollider::Common::SRMParameters params;
  params.bions = true;
  params.yions = true;
  p.get_fragment_masses(series, ch, params);

  BOOST_CHECK_EQUAL(series.size(), 12);
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

  std::vector<double> series;
  SRMPrecursor p;
  p.sequence = "PEPTIDE";
  p.isotope_modification = N15_ISOTOPEMODIFICATION;
  p.q1_charge = 2;
  double ch = 2;
  SRMCollider::Common::SRMParameters params;
  params.bions = true;
  params.yions = true;
  p.get_fragment_masses(series, ch, params);

  BOOST_CHECK_EQUAL(series.size(), 12);
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

BOOST_AUTO_TEST_CASE( calculate_charged_mass_TEST )
{

  double mass;
  SRMPrecursor p;
  p.sequence = "YYLLDYR";
  p.isotope_modification = NOISOTOPEMODIFICATION;

  p.q1_charge = 1;
  mass = p.calculate_charged_mass();
  BOOST_CHECK ( boost::test_tools::check_is_close( 1005.50460 , mass, EPS_05 ) );

  p.q1_charge = 2;
  mass = p.calculate_charged_mass();
  BOOST_CHECK ( boost::test_tools::check_is_close( 503.25623 , mass, EPS_05 ) );

  p.q1_charge = 3;
  mass = p.calculate_charged_mass();
  BOOST_CHECK ( boost::test_tools::check_is_close( 335.84011 , mass, EPS_05 ) );
                                                                
}

BOOST_AUTO_TEST_CASE( get_maximal_charge_TEST )
{

  SRMPrecursor p;
  p.sequence = "YYLLDYR";
  BOOST_CHECK_EQUAL ( p.get_maximal_charge(), 1);

  p.sequence = "YHHYLLDYR";
  BOOST_CHECK_EQUAL ( p.get_maximal_charge(), 3);

  p.sequence = "YHHYLKKLDYR";
  BOOST_CHECK_EQUAL ( p.get_maximal_charge(), 5);

                                                                
}

