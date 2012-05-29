
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MyTest
 
#include <boost/test/unit_test.hpp>
#include <iostream>

//include our own libraries
#include "srmcollider.h"
#include "combinatorics.h"

#define EPS_05 boost::test_tools::fraction_tolerance(1.e-5) 
 
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

BOOST_AUTO_TEST_CASE( _combinations_magic_TEST)
{
}
