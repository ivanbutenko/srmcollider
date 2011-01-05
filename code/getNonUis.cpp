
/*
 * The functions in this 
 *
 * Signature
 *void combinations(int,int,string,list<string>,list<list<int> >)
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

namespace python = boost::python;

python::dict _getnonuis_wrapper(python::tuple transitions,
        python::tuple collisions, int max_uis, double q3window, 
        bool ppm) {

    python::dict collisions_per_peptide;
    python::tuple clist;
    python::tuple tlist;
    python::list tmplist;

    int collision_length = python::extract<int>(collisions.attr("__len__")());
    int transitions_length = python::extract<int>(transitions.attr("__len__")());
    int tmplen;
    long c3, t1, tmplong;
    bool tmpbool;
    double t0, q3used = q3window;


    for (int i=0; i<transitions_length; i++) {
        tlist = python::extract< python::tuple >(transitions[i]);
        for (int j=0; j<collision_length; j++) {
            clist = python::extract< python::tuple >(collisions[j]);

            //ppm is 10^-6
            t0 = python::extract< double >(tlist[0]);
            if(ppm) {q3used = q3window / 1000000.0 * t0; } 
            if(fabs(t0-python::extract< double >(clist[0]) ) <  q3used) {

                c3 = python::extract<long>(clist[3]);
                t1 = python::extract<long>(tlist[1]);

                if(collisions_per_peptide.has_key(c3)) {
                    //append to the list in the dictionary
                    tmplist = python::extract<python::list>( 
                            collisions_per_peptide[c3] );
                    tmplen = python::extract<int>(tmplist.attr("__len__")());
                    tmpbool = false;
                    for (int k=0; k<tmplen; k++) {
                        tmplong = python::extract<long>(tmplist[k]);
                        if(tmplong == t1) {tmpbool = true;}
                    }
                    if(not tmpbool) { tmplist.append(t1); }
                }
                else {
                    //create new list in the dictionary
                    python::list newlist;
                    newlist.append(t1);
                    collisions_per_peptide[c3] = newlist;
                }
            }
        }
    }

    /*
     * In Python, this corresponds to the following function
     *

        collisions_per_peptide = {}
        q3_window_used = par.q3_window
        for t in transitions:
            if par.ppm: q3_window_used = par.q3_window * 10**(-6) * t[0]
            this_min = q3_window_used
            for c in collisions:
                if abs( t[0] - c[0] ) <= q3_window_used:
                    #gets all collisions
                    if collisions_per_peptide.has_key(c[3]):
                        if not t[1] in collisions_per_peptide[c[3]]:
                            collisions_per_peptide[c[3]].append( t[1] )
                    else: collisions_per_peptide[c[3]] = [ t[1] ] 


    */

    return collisions_per_peptide;
}

void initgetnonuis() {;}

using namespace python;
BOOST_PYTHON_MODULE(c_getnonuis)
{
    def("getnonuis", _getnonuis_wrapper, 

 "void tesat()\n"

           );
}

