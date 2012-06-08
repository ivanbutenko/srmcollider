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

//include our own libraries
#include "integratedrun.cpp"
#include "srmcolliderLib.cpp"
#include "rangetree.cpp"
#include "py_srmcolliderLib.h"

// Boost.Python headers
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
namespace python = boost::python;

namespace SRMCollider 
{
  namespace IntegratedRun 
  {

    void _pyToC_integratedrunTransitions(python::tuple& transitions, std::vector<Transition>& mytransitions)
    {
      /*
      * Transitions are tuples of the form (q3, transition_id)
      * convert to our struct.
      */
      int transitions_length = python::extract<int>(transitions.attr("__len__")());
      python::tuple tlist;
      for (int i=0; i<transitions_length; i++) {
          tlist = python::extract< python::tuple >(transitions[i]);
          double q3 = python::extract<double>(tlist[0]);
          long transition_id = python::extract<long>(tlist[1]);
          Transition entry = {q3, transition_id};
          mytransitions[i] = entry;
      }
    }

  // Python wrapper for min_needed
  int _py_min_needed(python::tuple transitions, python::tuple precursors,
      int max_uis, double q3window, bool ppm, python::object par)   
  {

      python::tuple clist;
      int i;
      std::string sequence;
      int transitions_length = python::extract<int>(transitions.attr("__len__")());
      int precursor_length = python::extract<int>(precursors.attr("__len__")());
      python::tuple tlist;
      
      //Check whether we have more transitions than we have bits in our number
      if (transitions_length > COMBLIMIT) {
          PyErr_SetString(PyExc_ValueError, 
              "Too many transitions, please adjust limit.");
          boost::python::throw_error_already_set();
          return -1;
      }

      SRMParameters params;
      pyToC::initialize_param_obj(par, params);

      std::vector<Transition> mytransitions(transitions_length);
      _pyToC_integratedrunTransitions(transitions, mytransitions);

      std::vector<SRMPrecursor> myprecursors(precursor_length);
      for (i=0; i<precursor_length; i++) {
          tlist = python::extract< python::tuple >(precursors[i]);
          sequence = python::extract<std::string>(tlist[1]);
          SRMPrecursor p;
          p.sequence = sequence;
          p.transition_group = p.q1 = p.maximal_charge = p.ssrcalc = -1;
          myprecursors[i] = p;
      }

      int result = -1;
      try
      {
        result = min_needed(mytransitions, myprecursors, max_uis, q3window, ppm, params);
      }
      catch (SRMCollider::Common::AANotFound& error)
      {
        PyErr_SetString(PyExc_ValueError, error.message.c_str());
        boost::python::throw_error_already_set();
        return result;
      }
      return result;
  }


  // Python wrapper for wrap_all
  python::list _py_wrap_all_bitwise(python::tuple& py_transitions, double a, double b,
    double c, double d, long thistransitiongr, int max_uis, double q3window,
    bool ppm, int max_nr_isotopes, double isotope_correction, python::object& par,
    SRMCollider::ExtendedRangetree::Rangetree_Q1_RT& rtree)
  {

      //Check whether we have more transitions than we have bits in our number
      int transitions_length = python::extract<int>(py_transitions.attr("__len__")());
      if (transitions_length > COMBLIMIT) {
          PyErr_SetString(PyExc_ValueError, 
              "Too many transitions, please adjust limit.");
          boost::python::throw_error_already_set();
          python::list tlist;
          return tlist;
      }

      SRMParameters params;
      pyToC::initialize_param_obj(par, params);

      std::vector<Transition> mytransitions(transitions_length);
      _pyToC_integratedrunTransitions(py_transitions, mytransitions);

      std::vector<COMBINT> collisions_per_peptide; 

      try
      {
        wrap_all_bitwise(mytransitions, a, b, c, d,
            thistransitiongr, max_uis, q3window, ppm, max_nr_isotopes, isotope_correction, 
             params, rtree, collisions_per_peptide);
      }
      catch (SRMCollider::Common::AANotFound& error)
      {
        PyErr_SetString(PyExc_ValueError, error.message.c_str());
        boost::python::throw_error_already_set();
        python::list result;
        return result;
      }

      std::vector<int> c_result;
      for(int i =1; i<= max_uis; i++) {
        std::set<COMBINT> combinations;
        SRMCollider::Combinatorics::get_non_uis_bitwise(collisions_per_peptide, transitions_length, i, combinations);
        c_result.push_back(combinations.size());
      }

      python::list result;
      for (size_t i = 0; i < c_result.size(); i++)
      {
        result.append(c_result[i]);
      }
      return result;
  }


  // Expose to Python
  using namespace python;
  BOOST_PYTHON_MODULE(c_integrated)
  {

      def("wrap_all_bitwise", _py_wrap_all_bitwise,
              
   "Return the number of non-UIS for all orders up to max_uis\n"
   "Given the transitions and then four numbers giving the coordinate window in\n"
   "which the precursors should be found: \n"
   "  q1_low, ssrcalc_low, q1_high, ssrcalc_high, \n"
   "Also the peptide key of the current peptide has to be given (to exclude it)\n"
   "as well as the maximal order of UIS to compute, the used q3 window and\n"
   "whether the q3 window is given in ppm (parts per million), default unit is\n"
   "Th which is m/z\n"
   "\n"
   "This function uses bitwise operations to speed up the calculation: each combination of\n"
   "transitions is encoded in an integer, which allows easy storage and\n"
   "computation by bitshits even though the code gets less easy to understand.\n"
   "The default length of the integer is 32 bits thus we can at most handle\n"
   "combinations for 31 transitions; if more are needed one needs to change\n"
   "COMBINT to an integer format that usees more bits, e.g. unit64_t or even the\n"
   "BigInt class for larger integers.\n"
   "If there are more transitions provided than allowed, an error will be\n"
   "thrown. \n"
   "\n"
   "\n"
   " Signature\n"
   "list wrap_all_bitwise(python::tuple transitions, double a, double b,\n"
          "double c, double d, long thistransitiongr, int max_uis, double q3window,\n"
          "bool ppm, int max_nr_isotopes, double isotope_correction)"
      );

      def("getMinNeededTransitions", _py_min_needed,
              
   "Calculate the minimally needed number of transitions needed to get a UIS\n"
   "given transitions in their preferred order and interfering precursors.\n"
   "\n"
   "Uses integers with bitflags to store the combinations, thus the number of\n"
   "transitions that can be considered is limited by COMBLIMIT.  Also note that\n"
   "the leftmost bit needs to be set to zero all the time, otherwise an infinte\n"
   "loop occurs, thus there is one bit we cannot use.\n"
   "\n"
   "\n"
   " Signature\n"
   "int min_needed(tuple transitions, tuple precursors,\n"
      "int max_uis, double q3window, bool ppm )\n"
              );
  }

  }
}
