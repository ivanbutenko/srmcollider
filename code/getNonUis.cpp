
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

boost::python::dict _getnonuis_wrapper(boost::python::list transitions,
        boost::python::list collisions, int max_uis, double q3window, 
        bool ppm) {

    boost::python::dict collisions_per_peptide;
    boost::python::list clist;
    boost::python::list tlist;
    boost::python::list tmplist;

    int collision_length = boost::python::extract<int>(collisions.attr("__len__")());
    int transitions_length = boost::python::extract<int>(transitions.attr("__len__")());
    int c3, t1, tmplen, tmpint, t0;
    bool tmpbool;
    double q3used = q3window;


    for (int i=0; i<transitions_length; i++) {
        tlist = boost::python::extract< boost::python::list >(transitions[i]);
        for (int j=0; j<collision_length; j++) {
            clist = boost::python::extract< boost::python::list >(collisions[j]);

            //ppm is 10^-6
            t0 = boost::python::extract< double >(tlist[0]);
            if(ppm) {q3used = q3window / 1000000.0 * t0; } 
            if(fabs(t0-boost::python::extract< double >(clist[0]) ) <  q3used) {

                c3 = boost::python::extract<int>(clist[3]);
                t1 = boost::python::extract<int>(tlist[1]);

                if(collisions_per_peptide.has_key(c3)) {
                    //append to the list in the dictionary
                    tmplist = boost::python::extract<boost::python::list>( 
                            collisions_per_peptide[c3] );
                    tmplen = boost::python::extract<int>(tmplist.attr("__len__")());
                    tmpbool = false;
                    for (int k=0; k<tmplen; k++) {
                        tmpint = boost::python::extract<int>(tmplist[k]);
                        if(tmpint == t1) {tmpbool = true;}
                    }
                    if(not tmpbool) { tmplist.append(t1); }
                }
                else {
                    //create new list in the dictionary
                    boost::python::list newlist;
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

using namespace boost::python;
BOOST_PYTHON_MODULE(c_getnonuis)
{
    def("getnonuis", _getnonuis_wrapper, 

 "void tesat()\n"

           );
}

