/*
 *
 * Program       : SRMCollider
 * Author        : Hannes Roest <roest@imsb.biol.ethz.ch>
 * Date          : 05.02.2011 
 *
 *
 * Copyright (C) 2011 Hannes Roest
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

/*
 * The functions in this file allow to calculate the collisions_per_peptide
 * dictionary easily, either using precursors (e.g. peptide sequences) or
 * collisions (q1,q3 tuples) as input. The collisions_per_peptide dictionary
 * holds for each peptide in the background the exact transitions of the query
 * peptides that are shared.
 *
 * Furthermore, we provide interfaces for some of the shared library functions.
*/

//include our own libraries
#include "srmcollider.h"
#include "srmcolliderLib.cpp"

// Function declarations

//functions to calculate collperpeptide dictionary
python::dict _getnonuis_wrapper(python::tuple transitions, python::tuple
        collisions, double q3window, bool ppm);
python::dict _find_clashes_calculate_collperpeptide(python::tuple transitions,
        python::tuple precursors, double q3_low, double q3_high, double
        q3window, bool ppm);

//calculate closest collision in q3 space for each transition
python::dict _find_clashes_core_non_unique(python::tuple transitions,
        python::tuple collisions, double q3window, bool ppm);

/*
 * Function to calculate the collisions_per_peptide out of a set of 
 * transitions and collisions. It will return a dictionary 
 * where for each key (colliding peptide key) the set of transitions
 * of the query peptide are stored that are interfered with is held.
 * Transitions are tuples of the form (q3, srm_id), collisions are tuples of the
 * form (q3, q1, srm_id, peptide_key).
 *
 * Input
 * transitions: (q3, srm_id)
 * collisions: (q3, q1, srm_id, peptide_key)
 *
 */
python::dict _getnonuis_wrapper(python::tuple transitions,
        python::tuple collisions, double q3window, bool ppm) {

    python::dict collisions_per_peptide;
    python::tuple clist;
    python::tuple tlist;
    python::list tmplist;

    int collision_length = python::extract<int>(collisions.attr("__len__")());
    int transitions_length = python::extract<int>(transitions.attr("__len__")());
    int tmplen;
    long c3, t1, tmplong;
    bool already_in_list;
    double t0, q3used = q3window;

    for (int i=0; i<transitions_length; i++) {
        tlist = python::extract< python::tuple >(transitions[i]);

        //ppm is 10^-6
        t0 = python::extract< double >(tlist[0]);
        if(ppm) {q3used = q3window / 1000000.0 * t0; } 

        // go through all (potential) collisions
        // and store the colliding SRM ids in a dictionary (they can be found at
        // position 3 and 1 respectively)
        for (int j=0; j<collision_length; j++) {
            clist = python::extract< python::tuple >(collisions[j]);

            if(fabs(t0-python::extract< double >(clist[0]) ) <  q3used) {

                c3 = python::extract<long>(clist[3]);
                t1 = python::extract<long>(tlist[1]);

                if(collisions_per_peptide.has_key(c3)) {
                    //append to the list in the dictionary
                    ///unless its already in the list
                    tmplist = python::extract<python::list>( 
                            collisions_per_peptide[c3] );
                    tmplen = python::extract<int>(tmplist.attr("__len__")());
                    already_in_list = false;
                    for (int k=0; k<tmplen; k++) {
                        tmplong = python::extract<long>(tmplist[k]);
                        if(tmplong == t1) {already_in_list = true;}
                    }
                    if(not already_in_list) { tmplist.append(t1); }
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

/*
 * Function to calculate the collision density out of a set of transitions
 * and precursor peptides that are colliding peptides.
 *
 * It will return a
 * dictionary where for each key (colliding peptide key) the set of transitions
 * of the query peptide are stored that are interfered with is held.
 * Transitions are tuples of the form (q3, srm_id), precursors are tuples of
 * the form (q1, sequence, peptide_key).
 */
python::list _find_clashes_calculate_colldensity(python::tuple transitions,
    python::list precursors, double q3_low, double q3_high, double q3window,
    bool ppm) {

    python::dict collisions_per_peptide, tmpdict;
    python::tuple clist;
    python::tuple tlist;

    python::list tmplist;

    int transitions_length = python::extract<int>(transitions.attr("__len__")());
    int precursor_length = python::extract<int>(precursors.attr("__len__")());
    int fragcount, i, j, k, ch;

    double q3used = q3window;
    char* sequence;

    double* b_series = new double[256];
    double* y_series = new double[256];

    int* cresult = new int[transitions_length];
    for (i=0; i<transitions_length; i++) cresult[i] = 0;

    double* tmptrans = new double[transitions_length];
    double* tmpq3used = new double[transitions_length];
    for (i=0; i<transitions_length; i++) {
        tlist = python::extract< python::tuple >(transitions[i]);
        //ppm is 10^-6
        tmptrans[i] = python::extract< double >(tlist[0]);
        if(ppm) {q3used = q3window / 1000000.0 * tmptrans[i]; } 
        tmpq3used[i] =q3used;
    }

    // go through all (potential) collisions
    // and store the colliding SRM ids in a dictionary (they can be found at
    // position 3 and 1 respectively)
    for (j=0; j<precursor_length; j++) {
        clist = python::extract< python::tuple >(precursors[j]);
        sequence = python::extract<char *>(clist[1]);

        for (ch=1; ch<=2; ch++) {
            fragcount = _calculate_clashes(sequence, b_series, y_series, ch);

            for (i=0; i<transitions_length; i++) {

                    // go through all fragments of this precursor
                    for (k=0; k<fragcount; k++) {
                        if(fabs(tmptrans[i]-y_series[k]) < tmpq3used[i]) cresult[i]++;
                        if(fabs(tmptrans[i]-b_series[k]) < tmpq3used[i]) cresult[i]++; 
                    }
                } //loop over all transitions
            }
    } //end of loop over all precursors

    delete [] b_series;
    delete [] y_series;

    // TODO replace with std::vector 
    // http://www.boost.org/doc/libs/1_41_0/libs/python/doc/v2/indexing.html
    // http://www.cplusplus.com/reference/stl/vector/reserve/
    //
    python::list result;
    for (i=0; i<transitions_length; i++) result.append( cresult[i] ) ;

    delete [] cresult ;
    delete [] tmptrans;
    delete [] tmpq3used;
    return result;
}



/*
 * Function to calculate the collisions_per_peptide out of a set of transitions
 * and precursor peptides that are colliding peptides.  It will return a
 * dictionary where for each key (colliding peptide key) the set of transitions
 * of the query peptide are stored that are interfered with is held.
 * Transitions are tuples of the form (q3, srm_id), precursors are tuples of
 * the form (q1, sequence, peptide_key).
 */
python::dict _find_clashes_calculate_collperpeptide(python::tuple transitions,
    python::tuple precursors, double q3_low, double q3_high, double q3window,
    bool ppm) {

    python::dict collisions_per_peptide, tmpdict;
    python::tuple clist;
    python::tuple tlist;
    python::list tmplist;

    int transitions_length = python::extract<int>(transitions.attr("__len__")());
    int precursor_length = python::extract<int>(precursors.attr("__len__")());
    int fragcount, i, j, k, ch, listmembers = 0;

    long t1, peptide_key;
    double t0, q3used = q3window;
    char* sequence;

    double* b_series = new double[256];
    double* y_series = new double[256];

    // go through all (potential) collisions
    // and store the colliding SRM ids in a dictionary (they can be found at
    // position 3 and 1 respectively)
    for (j=0; j<precursor_length; j++) {
        clist = python::extract< python::tuple >(precursors[j]);
        sequence = python::extract<char *>(clist[1]);

        for (ch=1; ch<=2; ch++) {
            fragcount = _calculate_clashes(sequence, b_series, y_series, ch);

            for (i=0; i<transitions_length; i++) {
                tlist = python::extract< python::tuple >(transitions[i]);
                //ppm is 10^-6
                t0 = python::extract< double >(tlist[0]);
                if(ppm) {q3used = q3window / 1000000.0 * t0; } 

                    // go through all fragments of this precursor
                    for (k=0; k<fragcount; k++) {

                        if(fabs(t0-y_series[k]) < q3used || 
                           fabs(t0-b_series[k]) < q3used) {
                            // extract SRM_id from transition list and store it
                            // as dirty in a temporary dict
                            t1 = python::extract<long>(tlist[1]);
                            tmpdict[t1] = 0;
                            listmembers++; 
                        
                        }
                    }
                } //loop over all transitions
            }

        //we keep one empty list around since we hope that most precursors dont 
        //use it. If it gets used, we store it in the result and create a new 
        //list to use from then on.
        //Sorting it costs something and could be done more efficiently since 
        //in fact, we only have to merge 2 presorted arrays. TODO Still the cost
        //is negligible.
        if (listmembers>0) {
            peptide_key = python::extract< long >(clist[2]);
            tmplist = tmpdict.keys();
            tmplist.sort();
            collisions_per_peptide[peptide_key] = tmplist;
            python::dict newlist;
            tmpdict = newlist;
        }
        listmembers = 0;

    } //end of loop over all precursors

    delete [] b_series;
    delete [] y_series;
    return collisions_per_peptide;
}

/* 
 * Function to calculate the closest collision in q3-space for each transition.
 * Transitions are tuples of the form (q3, srm_id), collisions are tuples of
 * the form (q3, q1, srm_id, peptide_key).
 * It will return a dictionary that contains for each srm_id the distance to
 * the closest hit.
 *
 * Input
 * transitions: (q3, srm_id)
 * collisions: (q3, q1, srm_id, peptide_key)
 *
 */
python::dict _find_clashes_core_non_unique(python::tuple transitions,
        python::tuple collisions, double q3window, bool ppm) {

    python::dict non_unique;
    python::tuple clist;
    python::tuple tlist;

    int collision_length = python::extract<int>(collisions.attr("__len__")());
    int transitions_length = python::extract<int>(transitions.attr("__len__")());
    long t1 ;
    double t0, c0, q3used = q3window, this_min;

    for (int i=0; i<transitions_length; i++) {
        tlist = python::extract< python::tuple >(transitions[i]);

        //ppm is 10^-6
        t0 = python::extract< double >(tlist[0]);
        if(ppm) {q3used = q3window / 1000000.0 * t0; } 
        this_min = q3used;

        // go through all (potential) collisions
        // and store closest collision in Q3 as a dictionary entry
        for (int j=0; j<collision_length; j++) {
            clist = python::extract< python::tuple >(collisions[j]);
            c0 = python::extract< double >(clist[0]);

            if(fabs(t0-c0) < this_min) {

                t1 = python::extract<long>(tlist[1]);
                this_min = fabs(t0-c0);
                non_unique[t1] = this_min;
            }
        }
    }

    /*
     * In Python, this corresponds to the following function
     *
            for t in transitions:
                if par.ppm: q3_window_used = par.q3_window * 10**(-6) * t[0]
                this_min = q3_window_used
                for c in collisions:
                    if abs( t[0] - c[0] ) <= this_min:
                        non_unique[ t[1] ] = t[0] - c[0]
            */

    return non_unique;
}


















/*
 * Function to calculate the exact interfering transitions for each peptide.
 * .  It will return a
 * Transitions are tuples of the form (q3, srm_id), precursors are tuples of
 * the form (q1, sequence, peptide_key).
 */
python::dict _find_clashes_forall(python::tuple transitions,
    python::tuple precursors, double q3_low, double q3_high, double q3window,
    bool ppm) {

    python::dict result, tmpdict;
    python::tuple clist;
    python::tuple tlist;
    python::list tmplist;

    int transitions_length = python::extract<int>(transitions.attr("__len__")());
    int precursor_length = python::extract<int>(precursors.attr("__len__")());
    int fragcount, i, j, k, ch, listmembers = 0;
    int isotope_nr;

    long t1, peptide_key;
    double t0, q1, ssrcalc, q3used = q3window;
    char* sequence;

    double* b_series = new double[256];
    double* y_series = new double[256];

    // go through all (potential) collisions
    // and store the colliding SRM ids in a dictionary (they can be found at
    // position 3 and 1 respectively)
    for (j=0; j<precursor_length; j++) {
        clist = python::extract< python::tuple >(precursors[j]);
        q1 = python::extract<double>(clist[0]);
        sequence = python::extract<char *>(clist[1]);
        peptide_key = python::extract<long>(clist[2]);

        ssrcalc = python::extract<double>(clist[3]);
        isotope_nr = python::extract<int>(clist[4]);

        for (ch=1; ch<=2; ch++) {
            fragcount = _calculate_clashes(sequence, b_series, y_series, ch);

            for (i=0; i<transitions_length; i++) {
                tlist = python::extract< python::tuple >(transitions[i]);
                //ppm is 10^-6
                t0 = python::extract< double >(tlist[0]);
                if(ppm) {q3used = q3window / 1000000.0 * t0; } 

                    // go through all fragments of this precursor
                    for (k=0; k<fragcount; k++) {

                        if(fabs(t0-y_series[k]) < q3used || 
                           fabs(t0-b_series[k]) < q3used) {
                            // extract SRM_id from transition list and store it
                            // as dirty in a temporary dict
                            t1 = python::extract<long>(tlist[1]);
                            if( result.has_key(t1) ) {
                                tmplist = python::extract<python::list>(result[t1]);
                                if(fabs(t0-y_series[k]) < q3used)
                                    // for the y series, the ion number is backwards
                                    tmplist.append( python::make_tuple(y_series[k],
                                    q1, 0, peptide_key, "y", fragcount - k, clist[1], ssrcalc, isotope_nr, ch));
                                else if(fabs(t0-b_series[k]) < q3used)
                                    tmplist.append( python::make_tuple(b_series[k],
                                    q1, 0, peptide_key, "b", k+1, clist[1], ssrcalc, isotope_nr, ch));
                            }
                            else{
                                python::list newlist;
                                tmplist = newlist;
                                if(fabs(t0-y_series[k]) < q3used)
                                    tmplist.append( python::make_tuple(y_series[k],
                                    q1, 0, peptide_key, "y", fragcount - k, clist[1], ssrcalc, isotope_nr, ch));
                                else if(fabs(t0-b_series[k]) < q3used)
                                    tmplist.append( python::make_tuple(b_series[k], 
                                    q1, 0, peptide_key, "b", k+1, clist[1], ssrcalc, isotope_nr, ch));
                                result[t1] = tmplist;
                            }
                        }
                    }
                } //loop over all transitions
            }

        /*
         * In Python, this corresponds to the following function
         *
        for t in transitions:
            #if par.ppm: q3_window_used = par.q3_window * 10**(-6) * t[0]
            this_min = q3_window_used
            for c in collisions:
                if abs( t[0] - c[0] ) <= this_min:
                    try: 
                        non_unique[ t[1] ].append( c )
                    except:
                        non_unique[ t[1] ] = [c]
        */

    } //end of loop over all precursors

    delete [] b_series;
    delete [] y_series;
    return result;
}


/*
 * Function to calculate whether there exists any combination of values from M
 * arrays (one value from each array) such that the M values are within a
 * certain window. 
 *
 * The first argument is a list that contains te length of the lists to be
 * checked. The second argument a list of lists which contain the values to be
 * checked, the third argument the window size. The function returns true if an
 * M-tuple of such values exists, otherwise false.
 */
bool thirdstrike(python::list myN, python::list py_ssrcalcvalues, double
        ssrwindow ) {

    python::list result;
    int M = python::extract<int>(myN.attr("__len__")());

    int j, k, i;
    int* index = new int[M];
    int* N = new int[M];
    double* myssr = new double[M];
    double avg;
    bool contaminationfree;
    bool contaminated = false;

    //check whether allocation was successfull
    if (! (index && N && myssr)) {
        PyErr_SetString(PyExc_ValueError, 
            "Memory allocation failed. Sorry.");
        python::throw_error_already_set();
        return false; }

    for(int k=0;k<M;k++) index[k] = 0;
    for(int k=0;k<M;k++) N[k] = python::extract<int>(myN[k]);
    int sumlen = 0;
    for(int k=0;k<M;k++) sumlen += N[k];

    /*
     * Memory layout
     *
     * First we allocate memory for an array that contains pointers to arrays.
     * We need M arrays, so this array has M entries. Then we allocate enough
     * consecutive memory to hold all our data, e.g. sumlen data. After this we
     * can fill our first array with the points to the correct place in the
     * second memory allocation space.
     * 
     * If our memory looks like this, 
     *
     * [ . . . . . . x . . . . . . y . . . . . . . . . . . z . . . ]
     *
     * and x,y,z are the starting points of an array, then c_ssrcalcvalues[0]
     * will point to the start, c_ssrcalcvalues[1] will point to x, 2 to y etc.
     * We thus have a matrix with different row lengths.
     * 
    */

    //allocate and check whether allocation was successfull
    double **c_ssrcalcvalues = (double **) malloc(M * sizeof(double *));
    if (!c_ssrcalcvalues) {
        PyErr_SetString(PyExc_ValueError, 
            "Memory allocation failed. Sorry.");
        python::throw_error_already_set();
        return false; }

    //allocate and check whether allocation was successfull
    c_ssrcalcvalues[0] = (double *) malloc(sumlen * sizeof(double));
    if (!c_ssrcalcvalues[0]) {
        PyErr_SetString(PyExc_ValueError, 
            "Memory allocation failed. Sorry.");
        python::throw_error_already_set();
        return false; }

    //calculate start position for each array in allocated memory and then fill 
    for(i = 1; i < M; i++)  
        c_ssrcalcvalues[i] = c_ssrcalcvalues[i-1] + N[i-1];
    python::list tmplist;
    for(k = 0; k < M; k++) {
        tmplist = python::extract<python::list>(py_ssrcalcvalues[k]);
        for(i = 0; i < N[k]; i++) {
            c_ssrcalcvalues[k][i] = python::extract<double>(tmplist[i]); }
    }

    while(true) {
        //EVALUATE THE RESULT
        //
        avg = 0;
        for(k=0; k < M; k++) {
            //get ssrcalc value from precursor index[k] of transition k
            myssr[k] = c_ssrcalcvalues[k][ index[k] ];
            avg += myssr[k];
        }
        avg /= M;
        contaminationfree = false;
        for(k=0; k < M; k++) {
            //if one of them deviates more than ssrwindow from the avg, 
            //there is no contamination, we are done with this index combination
            if( fabs(myssr[k] - avg) > ssrwindow) {
                contaminationfree = true; break;}
        }
        if(not contaminationfree) {contaminated = true; break;}

        //CALCULATE NEW INDEX
        // go through all combinations of M-tuples from the different arrays
        index[ M-1 ] += 1;
        if (index[ M-1 ] >= N[ M-1 ]) {
            //#now we hit the end, need to increment other positions than last
            j = M-1;
            while (j >= 0 and index[j] >= N[j]-1) {j -= 1;}
            //#j contains the value of the index that needs to be incremented
            //#when we are at the end of the interation, j will be -1
            if (j <0) break;
            index[j] += 1;
            k = j + 1;
            //#set all other positions to zero again
            while (k < M) {index[k] = 0; k += 1;  }
        }

    }

    free((void *)c_ssrcalcvalues);
    delete [] index;
    delete [] N;
    delete [] myssr;

    return contaminated;
}



using namespace python;
BOOST_PYTHON_MODULE(c_getnonuis)
{

    def("getnonuis", _getnonuis_wrapper, 
 "Function to calculate the collisions_per_peptide out of a set of \n"
 "transitions and collisions. It will return a dictionary \n"
 "where for each key (colliding peptide key) the set of transitions\n"
 "of the query peptide are stored that are interfered with is held.\n"
 "Transitions are tuples of the form (q3, srm_id), collisions are tuples of the\n"
 "form (q3, q1, srm_id, peptide_key)\n"
 "\n"
 "\n"
 " Signature\n"
 "dict getnonuis(tuple transitions, tuple collisions, double q3window, bool ppm)\n"
 );

    def("calculate_collisions_per_peptide", _find_clashes_calculate_collperpeptide, 
 "Function to calculate the collisions_per_peptide out of a set of transitions\n"
 "and precursor peptides that are colliding peptides.  It will return a\n"
 "dictionary where for each key (colliding peptide key) the set of transitions\n"
 "of the query peptide are stored that are interfered with is held.\n"
 "Transitions are tuples of the form (q3, srm_id), precursors are tuples of\n"
 "the form (q1, sequence, peptide_key).\n"
 "\n"
 "\n"
 " Signature\n"
 "dict calculate_collisions_per_peptide(tuple transitions, tuple precursors, \n"
 "       double q3_low, double q3_high, double q3window, bool ppm)\n"
 );


    def("get_non_uis", get_non_uis, 
 "Function to calculate all non-UIS transitions from a dictionary \n"
 "where for each key (colliding peptide key) the set of transitions\n"
 "of the query peptide are stored that are interfered with is held\n"
 "(collisions_per_peptide dictionary)\n"
 "It will return a list of all non UIS of the requested order.\n"
 "\n"
 "\n"
 " Signature\n"
 "list get_non_uis(dict collisions_per_peptide, int order)\n"
 );


    def("core_non_unique", _find_clashes_core_non_unique, 
 "Function to calculate the closest collision in q3-space for each transition.\n"
 "Transitions are tuples of the form (q3, srm_id), collisions are tuples of\n"
 "the form (q3, q1, srm_id, peptide_key).\n"
 "It will return a dictionary that contains for each srm_id the distance to\n"
 "the closest hit.\n"
 "\n"
 "\n"
 " Signature\n"
 "dict core_non_unique(tuple transitions, tuple collisions, double q3window, bool ppm)\n"
 );

    def("calculate_transitions", _find_clashes_calculate_clashes, 
 "Function to calculate all transitions of a list of precursor peptides.\n"
 "Precursors are tuples of the form (q1, sequence, peptide_key).\n"
 "It will return a list of tuples that are of the form \n"
 "(q3, q1, 0, peptide_key) \n"
 "\n"
 "\n"
 " Signature\n"
 "list calculate_transitions(tuple precursors, double q3_low, double q3_high ) \n"
 );
    def("calculate_transitions_ch", _find_clashes_calculate_clashes_ch,
 "Function to calculate all transitions of a list of precursor peptides and \n"
 "allows to select the charge states of these precursors.\n"
 "Precursors are tuples of the form (q1, sequence, peptide_key).\n"
 "It will return a list of tuples that are of the form \n"
 "(q3, q1, 0, peptide_key) \n"
 "\n"
 "\n"
 " Signature\n"
 "list calculate_transitions_ch(tuple precursors, list charges, double q3_low, double q3_high ) \n"
 );

    def("calculate_transitions_inner", _calculate_clashes_wrapper, 
 "Function to calculate all transitions of a precursor peptide ion of a"
 "defined charge state.  The input is a tuple of the form (q1, sequence,"
 "peptide_key) and the charge state.  It will return a list containing the b"
 "and y fragments."
 "\n"
 "\n"
 " Signature\n"
 "list calculate_transitions_inner(tuple precursor, double charge) \n"
            "");

    def("calculate_charged_mass", calculate_charged_mass, 
 "Calculate the charged mass for a sequence, supplied in a python tuple in the\n"
 "first place.\n"
 "\n"
 "\n"
 " Signature\n"
 "double calculate_charged_mass(python::tuple clist, int ch)"

            
            );
 def("calculate_density", _find_clashes_calculate_colldensity, "");
 def("_find_clashes_forall", _find_clashes_forall, "");

    def ("thirdstrike", thirdstrike);
}

