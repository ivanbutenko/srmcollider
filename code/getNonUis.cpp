
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
#include <time.h>

using namespace std;

/*
 *
 *
 * PYTHON interface
*/

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

#define MASS_H 1.007825032
#define MASS_OH 17.002739651999999

namespace python = boost::python;

python::dict _getnonuis_wrapper(python::tuple transitions, python::tuple
        collisions, double q3window, bool ppm);
python::list _find_clashes_calculate_clashes(python::tuple precursors,
        double q3_low, double q3_high );
double calculate_charged_mass(python::tuple clist, int ch);
python::dict _find_clashes_calculate_collperpeptide(python::tuple transitions,
        python::tuple precursors, double q3_low, double q3_high, double q3window, bool ppm);
python::list _calculate_clashes_wrapper(python::tuple &tlist, double charge);
int _calculate_clashes(const char* sequence, double* b_series, double* y_series, double ch);
python::dict _find_clashes_core_non_unique(python::tuple transitions,
        python::tuple collisions, double q3window, bool ppm);
void _combinations(int M, int N, const python::list &mapping, python::dict &result) ;
python::dict get_non_uis(python::dict collisions_per_peptide, int order) ;



/* Input
 * transitions: (q3, srm_id)
 * collisions: (q3, q1, srm_id, peptide_key)
 *
 * Function to calculate the collisions_per_peptide out of a set of 
 * transitions and collisions. It will return a dictionary 
 * where for each key (colliding peptide key) the set of transitions
 * of the query peptide are stored that are interfered with is held.
 * Transitions are tuples of the form (q3, srm_id), collisions are tuples of the
 * form (q3, q1, srm_id, peptide_key).
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


/* Input
 * precursors: (q1, sequence, peptide_key)
 *
 * Function to calculate all transitions of a list of precursor peptides and
 * allows to select the charge states of these precursors.
 * Precursors are tuples of the form (q1, sequence, peptide_key).
 * It will return a list of tuples that are of the form 
 * (q3, q1, 0, peptide_key) 
 */
python::list _find_clashes_calculate_clashes_ch(python::tuple precursors,
        python::list charges, double q3_low, double q3_high ) {

    python::tuple clist;
    python::list result;

    int precursor_length = python::extract<int>(precursors.attr("__len__")());
    int charges_length = python::extract<int>(charges.attr("__len__")());
    long peptide_key;
    int ch, fragcount, k;
    double q3, q1;
    char* sequence;

    double* b_series = new double[256];
    double* y_series = new double[256];

    for (int i=0; i<precursor_length; i++) {
        clist = python::extract< python::tuple >(precursors[i]);
        q1 = python::extract< double >(clist[0]);
        sequence = python::extract<char *>(clist[1]);
        peptide_key = python::extract< long >(clist[2]);

        for (int kk=0; kk<charges_length; kk++) {
            ch = python::extract< int >(charges[kk]);
            fragcount = _calculate_clashes(sequence, b_series, y_series, ch);
            // go through all fragments of this precursor
            for (k=0; k<fragcount; k++) {
                q3 = y_series[k];
                if (q3 > q3_low && q3 < q3_high)
                    result.append(python::make_tuple(q3, q1, 0, peptide_key) );
            }
            for (k=0; k<fragcount; k++) {
                q3 = b_series[k];
                if (q3 > q3_low && q3 < q3_high)
                    result.append(python::make_tuple(q3, q1, 0, peptide_key) );
            }
        }
    }

    delete [] b_series;
    delete [] y_series;

    return result;
}        



double calculate_charged_mass(python::tuple clist, int ch) {

    double charged_mass;
    double* b_series = new double[256];
    double* y_series = new double[256];

    char* sequence = python::extract<char *>(clist[1]);
    int fragcount = _calculate_clashes(sequence, b_series, y_series, ch);

    //In order to get the full mass, we need the "last" element of the b-series
    //(which is not technically part of the b series) and add water as well as
    //protons according to the charge to it. Then, the fragment has to be
    //divided by the charge.
    charged_mass = (b_series[fragcount] + MASS_H*ch + MASS_OH ) /ch  ;

    delete [] b_series;
    delete [] y_series;

    return charged_mass;
}        



/* Input
 * precursors: (q1, sequence, peptide_key)
 *
 * Function to calculate all transitions of a list of precursor peptides.
 * Precursors are tuples of the form (q1, sequence, peptide_key).
 * It will return a list of tuples that are of the form 
 * (q3, q1, 0, peptide_key) 
 */
python::list _find_clashes_calculate_clashes(python::tuple precursors,
        double q3_low, double q3_high ) {

    python::tuple clist;
    python::list result;

    int precursor_length = python::extract<int>(precursors.attr("__len__")());
    long peptide_key;
    int ch, fragcount, k;
    double q3, q1;
    char* sequence;

    double* b_series = new double[256];
    double* y_series = new double[256];

    for (int i=0; i<precursor_length; i++) {
        clist = python::extract< python::tuple >(precursors[i]);
        q1 = python::extract< double >(clist[0]);
        sequence = python::extract<char *>(clist[1]);
        peptide_key = python::extract< long >(clist[2]);

        for (ch=1; ch<=2; ch++) {
            fragcount = _calculate_clashes(sequence, b_series, y_series, ch);
            // go through all fragments of this precursor
            for (k=0; k<fragcount; k++) {
                q3 = y_series[k];
                if (q3 > q3_low && q3 < q3_high)
                    result.append(python::make_tuple(q3, q1, 0, peptide_key) );
            }
            for (k=0; k<fragcount; k++) {
                q3 = b_series[k];
                if (q3 > q3_low && q3 < q3_high)
                    result.append(python::make_tuple(q3, q1, 0, peptide_key) );
            }
        }
    }

    /*
     * Python code
     *

        for c in self.__get_all_precursors(par, pep, cursor, values=values):
            q1 = c[0]
            peptide_key = c[2]
            peptide = DDB.Peptide()
            peptide.set_sequence(c[1])
            peptide.charge = c[3]
            peptide.create_fragmentation_pattern(R)
            b_series = peptide.b_series
            y_series = peptide.y_series
            for ch in [1,2]:
                for pred in y_series:
                    q3 = ( pred + (ch -1)*R.mass_H)/ch
                    if q3 < q3_low or q3 > q3_high: continue
                    yield (q3, q1, 0, peptide_key)
                for pred in b_series:
                    q3 = ( pred + (ch -1)*R.mass_H)/ch
                    if q3 < q3_low or q3 > q3_high: continue
                    yield (q3, q1, 0, peptide_key)
    */
    delete [] b_series;
    delete [] y_series;

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

    python::dict collisions_per_peptide;
    python::tuple clist;
    python::tuple tlist;
    python::list tmplist;

    int transitions_length = python::extract<int>(transitions.attr("__len__")());
    int precursor_length = python::extract<int>(precursors.attr("__len__")());
    int fragcount, i, j, k, ch, listmembers = 0;
    long t1, peptide_key;
    bool append_to_list;
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

        for (i=0; i<transitions_length; i++) {
            append_to_list = false;
            tlist = python::extract< python::tuple >(transitions[i]);
            //ppm is 10^-6
            t0 = python::extract< double >(tlist[0]);
            if(ppm) {q3used = q3window / 1000000.0 * t0; } 

            for (ch=1; ch<=2; ch++) {
                fragcount = _calculate_clashes(sequence, b_series, y_series, ch);
                // go through all fragments of this precursor
                for (k=0; k<fragcount; k++) {
                    if(fabs(t0-y_series[k])<q3used || 
                       fabs(t0-b_series[k])<q3used) {append_to_list = true; }
                }
            }
            //append if necessary
            if(append_to_list) {
                t1 = python::extract<long>(tlist[1]);
                tmplist.append(t1);
                listmembers++; 
            }
        }

        if (listmembers>0) {
            peptide_key = python::extract< long >(clist[2]);
            collisions_per_peptide[peptide_key] = tmplist;
            //create new list
            python::list newlist;
            tmplist = newlist;
        }
        listmembers = 0;

    }

    delete [] b_series;
    delete [] y_series;

    return collisions_per_peptide;
}



/* Input
 * precursors: (q1, sequence, peptide_key)
 *
 * Function to calculate all transitions of a precursor peptide ion of a
 * defined charge state.  The input is a tuple of the form (q1, sequence,
 * peptide_key) and the charge state.  It will return a list containing the b
 * and y fragments.
 */

python::list _calculate_clashes_wrapper(python::tuple &tlist, double charge) {

    python::list result;
    double* b_series = new double[256];
    double* y_series = new double[256];
    int k;

    char* sequence = python::extract<char *>(tlist[1]);
    int fragcount = _calculate_clashes(sequence, b_series, y_series, charge);

    for (k=0; k<fragcount; k++) result.append(y_series[k]);
    for (k=0; k<fragcount; k++) result.append(b_series[k]);

    return result;

}

/*
 * Input is a tuple of (q1, sequence, peptide_key)
*/
int _calculate_clashes(const char* sequence, double* b_series, double* y_series,
        double ch) {

    int j, start, scounter;
    double acc_mass, res_mass;
    char c;
    bool inside;

    acc_mass = 0.0;
    res_mass = 0.0;
    scounter = 0;

    inside = false;
    start = 0;
    j = 0; 
    while(c = sequence[j++]) {
        if(sequence[j] == '[') {
            start = j-1;
            inside = true;
        }
        else if(sequence[j-1] == ']') {
            //for(k=start; k < j; k++) cout << sequence[k];
            //We found a modification
            switch(sequence[start]) {
                case 'M': 
                    if(!(sequence[start+2] == '1' && 
                         sequence[start+3] == '4' && 
                         sequence[start+4] == '7' )) {
                        PyErr_SetString(PyExc_ValueError, 
                            "Unknown modification for methionine");
                        boost::python::throw_error_already_set();
                        return -1;
                        }
                    res_mass = 147.03540462; break;
                case 'C': 
                    if(!(sequence[start+2] == '1' && 
                         sequence[start+3] == '6' && 
                         sequence[start+4] == '0' )) {
                        PyErr_SetString(PyExc_ValueError, 
                            "Unknown modification for cysteine");
                        boost::python::throw_error_already_set();
                        return -1;
                    }
                    res_mass = 160.030653721; break;
                default: 
                    PyErr_SetString(PyExc_ValueError, 
                        "Unknown modification ");
                    boost::python::throw_error_already_set();
                    return -1;
            }
            //'M[147]':  131.04049 + mass_O), # oxygen
            //'C[160]':  103.00919 + mass_CAM - mass_H ), # CAM replaces H

            acc_mass += res_mass;
            b_series[scounter] = acc_mass + MASS_H ;
            y_series[scounter] = - acc_mass + 2* MASS_H + MASS_OH ;
            scounter++;

            inside = false;
        }
        else if(inside) { }
        else {
            //We found a regular AA
            //TODO use hash map http://en.wikipedia.org/wiki/Hash_map_%28C%2B%2B%29
            switch(c) {
                case 'A': res_mass = 71.03711; break;
                case 'C': res_mass = 103.00919; break;
                case 'D': res_mass = 115.02694; break;
                case 'E': res_mass = 129.04259; break;
                case 'F': res_mass = 147.06841; break;
                case 'G': res_mass = 57.02146; break;
                case 'H': res_mass = 137.05891; break;
                case 'I': res_mass = 113.08406; break;
                case 'K': res_mass = 128.09496; break;
                case 'L': res_mass = 113.08406; break;
                case 'M': res_mass = 131.04049; break;
                case 'N': res_mass = 114.04293; break;
                case 'P': res_mass = 97.05276; break;
                case 'Q': res_mass = 128.05858; break;
                case 'R': res_mass = 156.10111; break;
                case 'S': res_mass = 87.03203; break;
                case 'T': res_mass = 101.04768; break;
                case 'V': res_mass = 99.06841; break;
                case 'W': res_mass = 186.07931; break;
                case 'X': res_mass = 113.08406; break;
                case 'Y': res_mass = 163.06333; break;
                default: 
                    PyErr_SetString(PyExc_ValueError, 
                        "Unknown amino acid ");
                    boost::python::throw_error_already_set();
                    return -1;
            }

            acc_mass += res_mass;
            b_series[scounter] = acc_mass + MASS_H ;
            y_series[scounter] = - acc_mass + 2* MASS_H + MASS_OH ;
            scounter++;
        }
    }

    for (int j=0; j<scounter; j++) y_series[j] += acc_mass;

    for (int j=0; j<scounter-1; j++){

        y_series[j] = (y_series[j] + (ch-1)*MASS_H)/ch;
        b_series[j] = (b_series[j] + (ch-1)*MASS_H)/ch;

    }
    //subtract one because we do not count the last mass as a fragment
    //the complete peptide can not undergo fragmentation and is thus not
    //a fragment. We still compute it but just dont consider it
    return --scounter;
}


/* Input
 * transitions: (q3, srm_id)
 * collisions: (q3, q1, srm_id, peptide_key)
 *
 * Function to calculate the closest collision in q3-space for each transition.
 * Transitions are tuples of the form (q3, srm_id), collisions are tuples of
 * the form (q3, q1, srm_id, peptide_key).
 * It will return a dictionary that contains for each srm_id the distance to
 * the closest hit.
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


//M is the order
//N is the length of the input vector

void _combinations(int M, int N, const python::list &mapping, 
        python::dict &result) {
    // The basic idea is to create an index array of length M that contains
    // numbers between 0 and N. The indices denote the combination produced and
    // we always increment the rightmost index. If it goes above N, we try to
    // increase the one left to it until we find one that still can be
    // increased. We stop when the rightmost index hits N-M, we thus go from
    // (0,1,2,3...M-1) to (N-M,N-M+1,N-M+2...N-1)
    int j, k;
    int* index = new int[M];
    //int* mytuple = new int[M];
    python::tuple tmptuple;
    //initialize with numbers from 0 to M = range( M )
    for(int k=0;k<M;k++) index[k] = k;
    while (index[0] <= N-M) {

        //EVALUATE THE RESULT
        //only works for M up to order 5
        switch(M) {
            case 1: tmptuple = python::make_tuple( mapping[index[0]] ); break;
            case 2: tmptuple = python::make_tuple( mapping[index[0]], mapping[index[1]]); break;
            case 3: tmptuple = python::make_tuple( mapping[index[0]], mapping[index[1]], 
                            mapping[index[2]]); break;
            case 4: tmptuple = python::make_tuple( mapping[index[0]], mapping[index[1]], 
                            mapping[index[2]], mapping[index[3]]); break;
            case 5: tmptuple = python::make_tuple( mapping[index[0]], mapping[index[1]], 
                            mapping[index[2]], mapping[index[3]], mapping[index[4]]); break;
            default:
                PyErr_SetString(PyExc_ValueError, 
                    "Order (M) larger than 5 is not implemented");
                boost::python::throw_error_already_set();
                return;
        }
        result[ tmptuple ] = 0; //append to result dict

        // We need to break if index[0] has the final value
        // The other option is to make the while condition (index[0] < N-M) 
        // and add an additional evaluation of the result to the end of the 
        // function (that would probably be faster).
        if(index[0] == N-M) break;

        index[ M-1 ] += 1;
        if (index[ M-1 ] >= N) {
            //#now we hit the end, need to increment other positions than last
            //#the last position may reach N-1, the second last only N-2 etc.
            j = M-1;
            while (j >= 0 and index[j] >= N-M+j) j -= 1;
            //#j contains the value of the index that needs to be incremented
            index[j] += 1;
            k = j + 1;
            while (k < M) { index[k] = index[k-1] + 1; k += 1; } 
        }
    }

    delete [] index;
}



/*
 * Function to calculate all non-UIS transitions from a dictionary
 * where for each key (colliding peptide key) the set of transitions
 * of the query peptide are stored that are interfered with is held
 * (collisions_per_peptide dictionary)
 * It will return a list of all non UIS of the requested order.
 */
python::dict get_non_uis(python::dict collisions_per_peptide, int order) {

    python::list tmplist;
    python::list pepcollisions = python::extract< python::list >(
            collisions_per_peptide.values() );
    python::dict result;

    int tmp_length;
    int collision_length = python::extract<int>(pepcollisions.attr("__len__")());

    for (int i=0; i<collision_length; i++) {
        tmplist = python::extract< python::list >( pepcollisions[i] );
        tmp_length = python::extract<int>(tmplist.attr("__len__")());
        _combinations(order, tmp_length, tmplist, result);
    }

    return result;
}

/* This above get_non_uis and _combinations functions 
 * correspond to the following python function (roughly)
 * 

    for pepc in collisions_per_peptide.values():
        for i in range(1,MAX_UIS+1):
            if len( pepc ) >= i: 
                tmp = [ [pepc[j] for j in indices] for indices in 
                    _combinations( len(pepc) , i)] 
                non_uis.update( [tuple(sorted(p)) for p in combinations(pepc, i)] )
    return non_uis_list

    def _combinations(N, M):
        """All index combinations of M elements drawn without replacement
         from a set of N elements.
        Order of elements does NOT matter."""
        index = range( M )
        while index[0] <= N-M:
            yield index[:]
            index[ M-1 ] += 1
            if index[ M-1 ] >= N:
                #now we hit the end, need to increment other positions than last
                #the last position may reach N-1, the second last only N-2 etc.
                j = M-1
                while j >= 0 and index[j] >= N-M+j: j -= 1
                #j contains the value of the index that needs to be incremented
                index[j] += 1
                k = j + 1
                while k < M: index[k] = index[k-1] + 1; k += 1; 

*/

int main() {
/* Input
 * transitions: (q3, srm_id)
 * precursors: (q1, sequence, peptide_key)
 */
    cout << "main" << endl;
    python::tuple transitions = python::make_tuple( 
            python::make_tuple(500.2, 1),
            python::make_tuple(600.2, 1)
            );
    python::tuple precursors = python::make_tuple( 
            python::make_tuple(500.2, "PEPTIDE", 101),
            python::make_tuple(600.2, "DEPTIDE", 102)
            );
    python::dict res =  _find_clashes_calculate_collperpeptide(transitions,
        precursors, 400, 1400, 2.0, false);

}

void initgetnonuis() {;}
void initcore_non_unique() {;}
void initget_non_uis() {;}
void initcalculate_transitions_inner() {;}

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
 "Function to return the mass of a peptide."
 "\n"
 "\n"
 " Signature\n"
 "double calculate_charged_mass(python::tuple clist, int ch)"

            
            );
}


