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

// Including headers from CGAL 
#include <CGAL/Cartesian.h>
#include <CGAL/Range_segment_tree_traits.h>
#include <CGAL/Range_tree_k.h>

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>

//include our own libraries
#include "srmcollider.h"
#include "srmcolliderLib.h"
#include "combinatorics.h"
#include "integratedrun.h"
#include "rangetree.h"

#include "getNonUis.cpp"
#include "py_integratedrun.cpp"

// this should be equivalent to running
// python run_integrated.py 1999 700 715 --q1_window=1.0 --q3_window=1.0 --ssrcalc_window=10 -p hroest.srmpeptides_yeast --max_uis 5
//
int main(int argc, const char ** argv)
{

  Py_Initialize();

  // now time to insert the current working directory into the python path so module search can take advantage
  // this must happen after python has been initialised
  boost::filesystem::path workingDir = boost::filesystem3::complete("./").normalize();
  // TODO use filesystem or filesystem3
  PyObject* sysPath = PySys_GetObject("path");
  PyList_Insert( sysPath, 0, PyString_FromString(workingDir.string().c_str()));

  python::object ignored;

  std::cout << "Will try to import all the modules" << std::endl;
  python::object mMainModule = python::import("__main__");
  python::object precursor = boost::python::import("precursor");
  python::object collider = boost::python::import("collider");
  python::object MySQLdb = boost::python::import("MySQLdb");
  python::object par = collider.attr("SRM_parameters")();
  std::cout << "Imported all the modules successfully, trying to get a db cusor..." << std::endl;

  double min_q1 = 700; 
  double max_q1 = 715; 
  double q1_window = 1.0/2.0;
  double q3_window = 1.0/2.0;
  double ssrcalc_window = 10/2.0;
  int max_uis = 5;
  int max_nr_isotopes = 3;
  bool ppm = false;

  ignored = par.attr("set_single_peptide_table")("hroest.srmPeptides_yeast");
  ignored = par.attr("set_q3_range")(400, 1400);
  ignored = par.attr("set_default_vars")();
  ignored = par.attr("set_mysql_config")("~/.my.cnf");
  python::object db = par.attr("get_db")();
  python::object cursor = db.attr("cursor")();
  std::cout << "Got a db cursor, trying to get precursors from the db..." << std::endl;

  python::object myprecursors = precursor.attr("Precursors")();
  myprecursors.attr("getFromDB")(par, cursor, min_q1 - q1_window, max_q1 + q1_window);
  std::cout << "Got all precursors from the db, trying to build a rangetree..." << std::endl;

  python::list precursors_to_evaluate = python::extract<python::list>(myprecursors.attr("getPrecursorsToEvaluate")(min_q1, max_q1));
  double isotope_correction = python::extract<double>(par.attr("calculate_isotope_correction")());
  //python::object r_tree = myprecursors.attr("build_extended_rangetree")();

  python::tuple alltuples = python::extract<python::tuple>(myprecursors.attr("get_alltuples_extended_rangetree")());
  SRMCollider::ExtendedRangetree::Rangetree_Q1_RT rangetree;
  rangetree.new_rangetree();
  rangetree.create_tree(alltuples);
  std::cout << "Built a rangetree, starting calculations..." << std::endl;

  SRMParameters params;
  SRMCollider::pyToC::initialize_param_obj(par, params);
  params.q3_window = q3_window;
  params.ppm = false;

  for (int i = 0; i < precursors_to_evaluate.attr("__len__")(); i++)
  {
    python::object precursor = precursors_to_evaluate[i];
    if (i % 1000 == 0) std::cout << "doing " << i << std::endl;

    double ssrcalc= python::extract<double>(precursor.attr("ssrcalc"));
    double q1= python::extract<double>(precursor.attr("q1"));
    python::object o_transition_group = precursor.attr("get_tr")();
    long transition_group = python::extract<long>(o_transition_group);
    //python::object o_transition_group = precursor.attr("transition_group ");
    //long transition_group = python::extract<long>(precursor.attr("transition_group "));
    double ssrcalc_low = ssrcalc - ssrcalc_window + 0.001;
    double ssrcalc_high = ssrcalc + ssrcalc_window - 0.001;
    python::tuple py_transitions = python::extract<python::tuple>(precursor.attr("calculate_transitions_from_param")(par));

    //python::object transitions = precursor.attr("calculate_transitions_from_param")(par);

#if 0
    SRMCollider::IntegratedRun::_py_wrap_all_bitwise(py_transitions,
        q1 - q1_window, ssrcalc_low, q1 + q1_window, ssrcalc_high,
        transition_group, max_uis, q3_window, ppm, max_nr_isotopes,
        isotope_correction, par, rangetree);
#else

    double transitions_length = python::extract<double>(py_transitions.attr("__len__")());
    std::vector<SRMCollider::Common::SRMTransition> mytransitions(transitions_length);
    SRMCollider::IntegratedRun::_pyToC_integratedrunTransitions(py_transitions, mytransitions);

    std::vector<COMBINT> collisions_per_peptide; 
    SRMCollider::IntegratedRun::wrap_all_bitwise(mytransitions,
      q1 - q1_window, ssrcalc_low, q1 + q1_window, ssrcalc_high,
      transition_group, max_uis, q3_window, ppm, max_nr_isotopes, isotope_correction,
         params, rangetree, collisions_per_peptide);

      std::vector<int> c_result;
      for(int i =1; i<= max_uis; i++) {
        std::set<COMBINT> combinations;
        SRMCollider::Combinatorics::get_non_uis_bitwise(collisions_per_peptide, transitions_length, i, combinations);
        c_result.push_back(combinations.size());
      }

#endif



            /*
            precursor.q1 - par.q1_window, ssrcalc_low, precursor.q1 + par.q1_window,  ssrcalc_high,
            precursor.transition_group, min(par.max_uis,len(transitions)), par.q3_window, par.ppm,
            par.isotopes_up_to, isotope_correction, par, r_tree)


  // Python wrapper for wrap_all
  python::list _py_wrap_all_bitwise(python::tuple py_transitions, double a, double b,
    double c, double d, long thistransitiongr, int max_uis, double q3window,
    bool ppm, int max_nr_isotopes, double isotope_correction, python::object par,
    SRMCollider::ExtendedRangetree::Rangetree_Q1_RT& rtree)
  {

            */


  }


/*
# python::exec( "for order in range(1,6):
#   sum_all = sum([p[0]*1.0/p[1] for p in prepare if p[3] == order]) 
#   nr_peptides = len([p for p in prepare if p[3] == order])
#   if not par.quiet and not nr_peptides ==0: print 'Order %s, Average non useable UIS %s' % (order, sum_all *1.0/ nr_peptides)")
*/




/*
 *
 *
 *
 *
  void wrap_all_bitwise(std::vector<Transition> mytransitions, double a, double b,
    double c, double d, long thistransitiongr, int max_uis, double q3window,
    bool ppm, int max_nr_isotopes, double isotope_correction, SRMParameters params,
    SRMCollider::ExtendedRangetree::Rangetree_Q1_RT& rtree, std::vector<COMBINT>& newcollperpep)
 *
 *
  python::list _py_wrap_all_bitwise(python::tuple py_transitions, double a, double b,
    double c, double d, long thistransitiongr, int max_uis, double q3window,
    bool ppm, int max_nr_isotopes, double isotope_correction, python::object par,
    SRMCollider::ExtendedRangetree::Rangetree_Q1_RT& rtree)
 *
progressm = progress.ProgressMeter(total=len(precursors_to_evaluate), unit='peptides')
prepare  = []
for precursor in precursors_to_evaluate:
    transitions = precursor.calculate_transitions_from_param(par)
    #correct rounding errors, s.t. we get the same results as before!
    ssrcalc_low = precursor.ssrcalc - par.ssrcalc_window + 0.001
    ssrcalc_high = precursor.ssrcalc + par.ssrcalc_window - 0.001
    try:
        result = c_integrated.wrap_all_bitwise(transitions,
            precursor.q1 - par.q1_window, ssrcalc_low, precursor.q1 + par.q1_window,  ssrcalc_high,
            precursor.transition_group, min(par.max_uis,len(transitions)), par.q3_window, par.ppm,
            par.isotopes_up_to, isotope_correction, par, r_tree)
    except ValueError: 
      print "Too many transitions for", precursor.modification
      continue
    for order in range(1,min(par.max_uis+1, len(transitions)+1)): 
        prepare.append( (result[order-1], collider.choose(len(transitions), 
            order), precursor.parent_id , order, exp_key)  )
    #//break;
    progressm.update(1)
*/



/*
db = par.get_db()
cursor = db.cursor()

mycollider = collider.SRMcollider()
par = collider.SRM_parameters()

*/



  /*


myprecursors = Precursors()
myprecursors.getFromDB(par, db.cursor(), min_q1 - par.q1_window, max_q1 + par.q1_window)
precursors_to_evaluate = myprecursors.getPrecursorsToEvaluate(min_q1, max_q1)
isotope_correction = par.calculate_isotope_correction()
r_tree = myprecursors.build_extended_rangetree ()

        bp::object ignored = bp::exec("hello = file('hello.txt', 'w')\n"
              "hello.write('Hello world!')\n"
              "hello.close()", mMainNamespace);
  */

}


// #local arguments
// exp_key = sys.argv[1]
// min_q1 = float(sys.argv[2])
// max_q1 = float(sys.argv[3])
// db = par.get_db()
// cursor = db.cursor()
// 
// if par.max_uis ==0: 
//     print "Please change --max_uis option, 0 does not make sense here"
//     sys.exit()

/*
# Get the precursors
###########################################################################
myprecursors = Precursors()
myprecursors.getFromDB(par, db.cursor(), min_q1 - par.q1_window, max_q1 + par.q1_window)
precursors_to_evaluate = myprecursors.getPrecursorsToEvaluate(min_q1, max_q1)
isotope_correction = par.calculate_isotope_correction()
r_tree = myprecursors.build_extended_rangetree ()

progressm = progress.ProgressMeter(total=len(precursors_to_evaluate), unit='peptides')
prepare  = []
for precursor in precursors_to_evaluate:
    transitions = precursor.calculate_transitions_from_param(par)
    #correct rounding errors, s.t. we get the same results as before!
    ssrcalc_low = precursor.ssrcalc - par.ssrcalc_window + 0.001
    ssrcalc_high = precursor.ssrcalc + par.ssrcalc_window - 0.001
    try:
        result = c_integrated.wrap_all_bitwise(transitions,
            precursor.q1 - par.q1_window, ssrcalc_low, precursor.q1 + par.q1_window,  ssrcalc_high,
            precursor.transition_group, min(par.max_uis,len(transitions)), par.q3_window, par.ppm,
            par.isotopes_up_to, isotope_correction, par, r_tree)
    except ValueError: 
      print "Too many transitions for", precursor.modification
      continue

    for order in range(1,min(par.max_uis+1, len(transitions)+1)): 
        prepare.append( (result[order-1], collider.choose(len(transitions), 
            order), precursor.parent_id , order, exp_key)  )
    #//break;
    progressm.update(1)

for order in range(1,6):
    sum_all = sum([p[0]*1.0/p[1] for p in prepare if p[3] == order]) 
    nr_peptides = len([p for p in prepare if p[3] == order])
    if not par.quiet and not nr_peptides ==0: print "Order %s, Average non useable UIS %s" % (order, sum_all *1.0/ nr_peptides)
    #cursor.execute("insert into hroest.result_completegraph_aggr (sum_nonUIS, nr_peptides, uisorder, experiment) VALUES (%s,%s,%s,'%s')" % (sum_all, nr_peptides, order, exp_key))

*/
