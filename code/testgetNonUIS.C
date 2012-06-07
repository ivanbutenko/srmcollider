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

#include "getNonUis.cpp"

#define EPS_05 boost::test_tools::fraction_tolerance(1.e-5) 

using namespace SRMCollider::Common;
 
BOOST_AUTO_TEST_CASE(has_allowed_charge_TEST)
{
  // TODO
}

BOOST_AUTO_TEST_CASE(_find_clashes_calculate_colldensity_TEST)
{
  // TODO
}
BOOST_AUTO_TEST_CASE(_find_clashes_calculate_collperpeptide_other_ion_series_TEST)
{
  // TODO
}

BOOST_AUTO_TEST_CASE(_find_clashes_forall_other_series_TEST)
{
  // TODO
}


