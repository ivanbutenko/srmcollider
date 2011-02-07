
/*
 * The functions in this 
 *
 * Signature
 * void combinations(int,int,string,list<string>,list<list<int> >)
 *
 *
 * sudo apt-get install libcgal-dev 
 * g++ -frounding-math -lCGAL -O3 rangetree.cpp 
 * http://www.cgal.org/Manual/latest/doc_html/cgal_manual/SearchStructures/Chapter_main.html
 * http://graphics.stanford.edu/courses/cs368-00-spring/TA/manuals/CGAL/ref-manual2/SearchStructures/Chapter_main.html
 * Q Public Licence http://www.gnu.org/licenses/license-list.html
 * http://www.cgal.org/Manual/latest/doc_html/cgal_manual/packages.html#part_XIII
*/

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

/*
 *
 *
 * PYTHON interface
*/

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

/*
CGAL
*/
#include <CGAL/Cartesian.h>
#include <CGAL/Range_segment_tree_traits.h>
#include <CGAL/Range_tree_k.h>

typedef CGAL::Cartesian<double> K;
typedef CGAL::Range_tree_map_traits_2<K, long> Traits;
typedef CGAL::Range_tree_2<Traits> Range_tree_2_type;

typedef Traits::Key Key;                
typedef Traits::Interval Interval;    
Range_tree_2_type *Range_tree_2 = new Range_tree_2_type;

namespace python = boost::python;

void create_tree(python::tuple pepids) {

    python::tuple tlist;
    std::vector<Key> InputList;
    int i;
    long parent_id;
    double ssrcalc, q1;

    int pepids_length = python::extract<int>(pepids.attr("__len__")());
    for (i=0; i<pepids_length; i++) {
        tlist = python::extract< python::tuple >(pepids[i]);

        parent_id = python::extract<long>(tlist[2]);
        q1 = python::extract<double>(tlist[4]);
        ssrcalc = python::extract<double>(tlist[5]);

        InputList.push_back(Key(K::Point_2(q1,ssrcalc), parent_id));
    }
    Range_tree_2->make_tree(InputList.begin(),InputList.end());
}


python::list query_tree(double a, double b, double c, double d)   {

  std::vector<Key> OutputList;
  python::list result;

  Interval win(Interval(K::Point_2(a,b),K::Point_2(c,d)));
  Range_tree_2->window_query(win, std::back_inserter(OutputList));
  std::vector<Key>::iterator current=OutputList.begin();
  while(current!=OutputList.end()){
      result.append(python::make_tuple(
                  (*current).second));
                  /*
                  (*current).first.x(), 
                  (*current).second.peptide_key, 
                  (*current).second.q1_charge));
                  */
      current++;
    }
  return result;

}

void initcreate_tree() {;}
void initquery_tree() {;}

using namespace python;
BOOST_PYTHON_MODULE(c_rangetree)
{

    def("create_tree", create_tree, "");
    def("query_tree", query_tree, "");

}
