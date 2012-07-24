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

#ifndef PY_SRMCOLLIDERLIB_H
#define PY_SRMCOLLIDERLIB_H

// Boost.Python headers
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
namespace python = boost::python;

namespace SRMCollider 
{
  namespace Common
  {
    /* 
     * Function to calculate all transitions of a list of precursor peptides and
     * allows to select the charge states of these precursors.
     * Precursors are tuples of the form (q1, sequence, peptide_key).
     * It will return a list of tuples that are of the form 
     * (q3, q1, 0, peptide_key) 
     *
     * Input
     * precursors: (q1, sequence, peptide_key)
     *
     */
    python::list _py_calculate_transitions_with_charge(python::tuple py_precursors,
            python::list py_charges, double q3_low, double q3_high ) 
    {

        python::tuple clist;
        python::list py_result;

        //std::vector<SRMTransition> result;
        std::vector<SRMPrecursor> precursors;
        std::vector<int> charges;

        SRMParameters param;

        double* tmp = new double[1024];
        double* series = new double[1024];

        for (int kk=0; kk<python::extract<int>(py_charges.attr("__len__")()); kk++) {
          charges.push_back(python::extract< int >(py_charges[kk]));
        }
        for (int i=0; i<python::extract<int>(py_precursors.attr("__len__")()); i++) {
          SRMPrecursor p; 
          clist = python::extract< python::tuple >(py_precursors[i]);
          p.q1 = python::extract< double >(clist[0]);
          p.sequence = python::extract<std::string>(clist[1]);
          p.transition_group = python::extract< long >(clist[2]);
          precursors.push_back(p);
        }

        for (size_t i=0; i< precursors.size(); i++) {
          std::vector<SRMTransition> result;
          precursors[i].calculate_transitions_with_charge(charges, result, q3_low, q3_high, param);
          for (std::vector<SRMTransition>::iterator it = result.begin(); it != result.end(); it++)
          {
            py_result.append(python::make_tuple(it->q3, precursors[i].q1, 1, precursors[i].transition_group));
          }
        }

        /*
         * we cannot do this because of the way we test this function in Python...
        for (std::vector<SRMTransition>::iterator it = result.begin(); it != result.end(); it++)
        {
          py_result.append(python::make_tuple(it->q3, -1, 0, it->transition_id));
        }
        */

        delete [] tmp;
        delete [] series;

        return py_result;
    }        

    double _py_calculate_charged_mass(python::tuple clist, int ch) 
    {
        std::string sequence = python::extract<std::string>(clist[1]);
        SRMPrecursor p;
        p.sequence = sequence;
        p.isotope_modification = NOISOTOPEMODIFICATION;
        p.q1_charge = ch;
        return p.calculate_charged_mass();
    }        
  }

  namespace pyToC
  {

    void initialize_param_obj(python::object& par, Common::SRMParameters& params)
    {
        params.aions      =  python::extract<bool>(par.attr("aions"));
        params.aMinusNH3  =  python::extract<bool>(par.attr("aMinusNH3"));
        params.bions      =  python::extract<bool>(par.attr("bions"));
        params.bMinusH2O  =  python::extract<bool>(par.attr("bMinusH2O"));
        params.bMinusNH3  =  python::extract<bool>(par.attr("bMinusNH3"));
        params.bPlusH2O   =  python::extract<bool>(par.attr("bPlusH2O"));
        params.cions      =  python::extract<bool>(par.attr("cions"));
        params.xions      =  python::extract<bool>(par.attr("xions"));
        params.yions      =  python::extract<bool>(par.attr("yions"));
        params.yMinusH2O  =  python::extract<bool>(par.attr("yMinusH2O"));
        params.yMinusNH3  =  python::extract<bool>(par.attr("yMinusNH3"));
        params.zions      =  python::extract<bool>(par.attr("zions"));
        params.MMinusH2O  =  python::extract<bool>(par.attr("MMinusH2O"));
        params.MMinusNH3  =  python::extract<bool>(par.attr("MMinusNH3"));

        params.ppm  =  python::extract<bool>(par.attr("ppm"));
        params.q3_window  =  python::extract<double>(par.attr("q3_window"));
    }

    void initialize_precursors(python::list& precursors, std::vector<Common::SRMPrecursor>& c_precursors)
    {
      for (int i=0; i<python::extract<int>(precursors.attr("__len__")()); i++) {
        Common::SRMPrecursor p;
        python::object precursor = python::extract< python::object >(precursors[i]);
        python::object q1 = python::extract< python::object >(precursor.attr("q1"));
        python::object ssrcalc = python::extract< python::object >(precursor.attr("ssrcalc"));
        // check for None
        if(q1.ptr() != python::object().ptr() )
        {
          p.q1 = python::extract<double>(precursor.attr("q1"));
        }
        if(ssrcalc.ptr() != python::object().ptr() )
        {
          p.ssrcalc = python::extract<double>(precursor.attr("ssrcalc"));
        }
        p.sequence = python::extract<std::string>(precursor.attr("modified_sequence"));
        p.isotope_modification = python::extract<int>(precursor.attr("isotopically_modified"));
        p.q1_charge = python::extract<int>(precursor.attr("q1_charge"));
        p.maximal_charge = python::extract<int>(precursor.attr("to_peptide")().attr("get_maximal_charge")());
        p.transition_group = python::extract<long>(precursor.attr("transition_group"));
        c_precursors.push_back(p);
      }
    }

    void initialize_transitions(python::tuple& transitions, std::vector<Common::SRMTransition>& c_transitions)
    {
      for (int i=0; i<python::extract<int>(transitions.attr("__len__")()); i++) {
        Common::SRMTransition t;
        python::tuple tlist = python::extract< python::tuple >(transitions[i]);
        t.q3 = python::extract< double >(tlist[0]);
        t.transition_id = python::extract<long>(tlist[1]);
        c_transitions.push_back(t);
      }
    }

  }
}


#endif
