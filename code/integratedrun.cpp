
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

#define MASS_H 1.007825032
#define MASS_OH 17.002739651999999

//assume that there are never more than 32 transitions in an assay
#define COMBINT uint32_t


/*
CGAL
*/
#include <CGAL/Cartesian.h>
#include <CGAL/Range_segment_tree_traits.h>
#include <CGAL/Range_tree_k.h>

struct Precursor{
  char* sequence; 
  long peptide_key;
  long parent_id;
  int q1_charge;
};

struct Transition{
  double q3;
  long srm_id;
};




typedef CGAL::Cartesian<double> K;
typedef CGAL::Range_tree_map_traits_2<K, Precursor> Traits;
typedef CGAL::Range_tree_2<Traits> Range_tree_2_type;

typedef Traits::Key Key;                
typedef Traits::Interval Interval;    
Range_tree_2_type *Range_tree_2 = new Range_tree_2_type;


namespace python = boost::python;
int _calculate_clashes(const char* sequence, double* b_series, double* y_series,
        double ch) ;
void _combinations(int M, int N, const python::list &mapping, python::dict &result) ;
python::dict get_non_uis(python::dict collisions_per_peptide, int order) ;

void create_tree(python::tuple pepids) {

    python::tuple tlist;
    std::vector<Key> InputList;
    int i, q1_charge;
    long parent_id, peptide_key;
    double ssrcalc, q1;
    char* sequence;

    int pepids_length = python::extract<int>(pepids.attr("__len__")());
    for (i=0; i<pepids_length; i++) {
        tlist = python::extract< python::tuple >(pepids[i]);

        sequence = python::extract<char *>(tlist[0]);
        peptide_key = python::extract<long>(tlist[1]);
        parent_id = python::extract<long>(tlist[2]);
        q1_charge = python::extract<int>(tlist[3]);

        q1 = python::extract<double>(tlist[4]);
        ssrcalc = python::extract<double>(tlist[5]);

        struct Precursor entry = {sequence, peptide_key, parent_id, q1_charge};
        InputList.push_back(Key(K::Point_2(q1,ssrcalc), entry));

        /*
        parent_id = python::extract<long>(tlist[2]);
        q1 = python::extract<double>(tlist[4]);
        ssrcalc = python::extract<double>(tlist[5]);

        InputList.push_back(Key(K::Point_2(q1,ssrcalc), parent_id));
        */
    }
    Range_tree_2->make_tree(InputList.begin(),InputList.end());
}



//M is the order
//N is the length of the input vector

void _combinations_magic(int M, int N, int* mapping,
        python::dict &result) {
    // The basic idea is to create an index array of length M that contains
    // numbers between 0 and N. The indices denote the combination produced and
    // we always increment the rightmost index. If it goes above N, we try to
    // increase the one left to it until we find one that still can be
    // increased. We stop when the rightmost index hits N-M, we thus go from
    // (0,1,2,3...M-1) to (N-M,N-M+1,N-M+2...N-1)

    int j, k;
    int* index = new int[M];
    COMBINT tmpres;

    python::tuple tmptuple;
    //initialize with numbers from 0 to M = range( M )
    for(int k=0;k<M;k++) index[k] = k;
    while (index[0] <= N-M) {

        //EVALUATE THE RESULT
        tmpres = 0;
        for(k=0;k<M;k++) 
            tmpres |= mapping[index[k]];
        result[tmpres] = 0;

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
python::dict get_non_uis_magic(vector<COMBINT>& newcollperpep, int max_tr, int order) {

    python::dict result;

    int onecounter;
    COMBINT tmparr;
    COMBINT mask;
    int* mapping = new int[max_tr];

    for (uint i=0; i<newcollperpep.size(); i++) {

        //count the number of binary ones in the bitarray
        mask = 1;
        onecounter = 0;
        tmparr = newcollperpep[i];

        //we know that there cannot be bits populated past max nr transitions
        for(int j=0; j<max_tr; j++) {
            //true if nonzero
            if(tmparr & mask) 
                mapping[onecounter++] = mask;
            mask <<=1;
        }

        /* The other way to do it
         * we shift the tmparr to the right
        for(int j=0; tmparr != 0; tmparr >>= 1) {
            onecounter += tmparr & mask;
        }
        */

        _combinations_magic(order, onecounter, mapping, result);
    }

    delete [] mapping;
    return result;
}


int min_needed(python::tuple transitions, python::tuple precursors,
    int max_uis, double q3window, bool ppm )   {

    //use the defined COMBINT (default 32bit int) and some magic to do this :-)
    COMBINT one;
    COMBINT currenttmp = 0;
    std::vector<COMBINT> newcollperpep;

    python::tuple clist;
    int fragcount, i, k, ch, j;
    long srm_id;
    double q3, q3used = q3window;
    char* sequence;

    Transition transition;
    std::vector<Key> OutputList;

    double* b_series = new double[256];
    double* y_series = new double[256];

    int precursor_length = python::extract<int>(precursors.attr("__len__")());
    int transitions_length = python::extract<int>(transitions.attr("__len__")());

    /*
    * Transitions are tuples of the form (q3, srm_id)
    * convert to our struct.
    */
    python::tuple tlist;
    vector<Transition> mytransitions(transitions_length);
    for (i=0; i<transitions_length; i++) {
        tlist = python::extract< python::tuple >(transitions[i]);
        q3 = python::extract<double>(tlist[0]);
        srm_id = python::extract<long>(tlist[1]);
        struct Transition entry = {q3, srm_id};
        mytransitions[i] = entry;
    }

    // Go through all (potential) collisions 
    //
    int maxoverlap = 0;
    for (j=0; j<precursor_length; j++) {
        clist = python::extract< python::tuple >(precursors[j]);
        sequence = python::extract<char *>(clist[1]);
        for (ch=1; ch<=2; ch++) {
            fragcount = _calculate_clashes(sequence, b_series, y_series, ch);
            for(int i = 0; i != transitions_length; i++) {
                //ppm is 10^-6
                transition = mytransitions[i];
                q3 = transition.q3;
                if(ppm) {q3used = q3window / 1000000.0 * q3; } 
                    // go through all fragments of this precursor
                    for (k=0; k<fragcount; k++) {
                        if(fabs(q3-y_series[k]) < q3used || 
                           fabs(q3-b_series[k]) < q3used) {
                            //left bitshift == 2^i
                            one = 1;
                            currenttmp |= one << i;
                        }
                    }
                } //loop over all transitions
            } //end loop over all charge states of this precursor
        if ( currenttmp ) {
            //while we find the transitions from our relative order also in the peptide
            //we just looked at, increase i
            one = 1;
            i = 0;
            while( currenttmp & one << i ) i++;
            if( i > maxoverlap ) maxoverlap = i;
            currenttmp = 0;
        }
    }
    //we have counted from 0 to N-1  but now want the number of transitions
    return ++maxoverlap;
}


python::list wrap_all_magic(python::tuple transitions, double a, double b,
        double c, double d, long thispeptide_key, int max_uis, double q3window,
        bool ppm )   {

    //use the defined COMBINT (default 32bit int) and some magic to do this :-)
    COMBINT one;
    COMBINT currenttmp = 0;
    std::vector<COMBINT> newcollperpep;


    int fragcount, i, k, ch;
    long srm_id;
    double q3, q3used = q3window;
    char* sequence;

    Transition transition;
    Precursor precursor;
    std::vector<Key> OutputList;

    double* b_series = new double[256];
    double* y_series = new double[256];

    /*
    * Transitions are tuples of the form (q3, srm_id)
    * convert to our struct.
    */
    python::tuple tlist;
    int transitions_length = python::extract<int>(transitions.attr("__len__")());
    vector<Transition> mytransitions(transitions_length);
    for (i=0; i<transitions_length; i++) {
        tlist = python::extract< python::tuple >(transitions[i]);
        q3 = python::extract<double>(tlist[0]);
        srm_id = python::extract<long>(tlist[1]);
        struct Transition entry = {q3, srm_id};
        mytransitions[i] = entry;
    }

    Interval win(Interval(K::Point_2(a,b),K::Point_2(c,d)));
    Range_tree_2->window_query(win, std::back_inserter(OutputList));
    std::vector<Key>::iterator current=OutputList.begin();

    // Go through all (potential) collisions we just extracted from the rangetree
    //
    // This assumes that we do not have peptide_keys mapped to different
    // colliding transitions. In fact we do not have this situation even though
    // we have duplicate peptide_keys (from the isotopes). But they will
    // produce the same interefering transitions and thus the same entry in the
    // collisions per peptide table.
    while(current!=OutputList.end()){
        if (thispeptide_key != (*current).second.peptide_key) {

            precursor =  (*current).second;
            sequence = precursor.sequence;
            for (ch=1; ch<=2; ch++) {
                fragcount = _calculate_clashes(sequence, b_series, y_series, ch);

                //for(std::vector<int>::size_type i = 0; i != transitions_length; i++) {
                for(int i = 0; i != transitions_length; i++) {
                    //ppm is 10^-6
                    transition = mytransitions[i];
                    q3 = transition.q3;
                    if(ppm) {q3used = q3window / 1000000.0 * q3; } 
                        // go through all fragments of this precursor
                        for (k=0; k<fragcount; k++) {
                            if(fabs(q3-y_series[k]) < q3used || 
                               fabs(q3-b_series[k]) < q3used) {

                                //left bitshift == 2^i
                                one = 1;
                                currenttmp |= one << i;
                            }
                        }
                    } //loop over all transitions
                }

            //Store current combination
            if ( currenttmp ) {
                newcollperpep.push_back(currenttmp);
                currenttmp = 0;
            }

        }//end of if
      current++;
    }

    //this takes about 50% or more of the time if we have many collisions_per_pep (10k)
    //and below 10% if we have few collision_per_pep (0.1k)
    python::list result;
    for(i =1; i<= max_uis; i++) {
        result.append( 
                get_non_uis_magic(newcollperpep, transitions_length, i).attr("__len__")() );
    }

    delete [] b_series;
    delete [] y_series;

    return result;
}

python::list wrap_all(python::tuple transitions, 
        double a, double b, double c, double d,
        long thispeptide_key, 
        int max_uis, double q3window,
        bool ppm
        )   {

    std::vector<Key> OutputList;

    python::tuple tlist;

    python::dict collisions_per_peptide, tmpdict;
    python::list tmplist;
    python::list precursors;
    int fragcount, i, k, ch, listmembers = 0;
    long srm_id, peptide_key;
    double q3, q3used = q3window;
    char* sequence;

    Transition transition;
    Precursor precursor;

    double* b_series = new double[256];
    double* y_series = new double[256];

    /*
    * Transitions are tuples of the form (q3, srm_id)
    * convert to our struct.
    */
    int transitions_length = python::extract<int>(transitions.attr("__len__")());
    vector<Transition> mytransitions(transitions_length);
    for (i=0; i<transitions_length; i++) {
        tlist = python::extract< python::tuple >(transitions[i]);
        q3 = python::extract<double>(tlist[0]);
        srm_id = python::extract<long>(tlist[1]);
        struct Transition entry = {q3, srm_id};
        mytransitions[i] = entry;
    }

    Interval win(Interval(K::Point_2(a,b),K::Point_2(c,d)));
    Range_tree_2->window_query(win, std::back_inserter(OutputList));
    std::vector<Key>::iterator current=OutputList.begin();

    // go through all (potential) collisions we just extracted from the rangetree
    // and store the colliding SRM ids in a dictionary (they can be found at
    // position 3 and 1 respectively)
    //
    // This assumes that we do not have peptide_keys mapped to different
    // colliding transitions. In fact we do not have this situation even though
    // we have duplicate peptide_keys (from the isotopes). But they will
    // produce the same interefering transitions and thus the same entry in the
    // collisions per peptide table.
    while(current!=OutputList.end()){
        if (thispeptide_key != (*current).second.peptide_key) {

            precursor =  (*current).second;
            sequence = precursor.sequence;
            for (ch=1; ch<=2; ch++) {
                fragcount = _calculate_clashes(sequence, b_series, y_series, ch);

                for(int i = 0; i != transitions_length; i++) {
                    //ppm is 10^-6
                    transition = mytransitions[i];
                    q3 = transition.q3;
                    if(ppm) {q3used = q3window / 1000000.0 * q3; } 
                        // go through all fragments of this precursor
                        for (k=0; k<fragcount; k++) {
                            if(fabs(q3-y_series[k]) < q3used || 
                               fabs(q3-b_series[k]) < q3used) {

                                srm_id = transition.srm_id;
                                tmpdict[srm_id] = 0; //dummy dict
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
            if (listmembers>0 ) {
                peptide_key = precursor.peptide_key;
                tmplist = tmpdict.keys();
                tmplist.sort();
                collisions_per_peptide[peptide_key] = tmplist;
                python::dict newlist;
                tmpdict = newlist;
                listmembers = 0;
            }

        }//end of if
      current++;
    }

    //this takes about 50% or more of the time if we have many collisions_per_pep (10k)
    //and below 10% if we have few collision_per_pep (0.1k)
    python::list result;
    for(i =1; i<= max_uis; i++) {
        result.append( get_non_uis(collisions_per_peptide, i).attr("__len__")() );
    }

    delete [] b_series;
    delete [] y_series;

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
    while((c = sequence[j++])) {
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

using namespace python;
BOOST_PYTHON_MODULE(c_integrated)
{

    def("create_tree", create_tree, "");
    def("wrap_all", wrap_all, "");
    def("wrap_all_magic", wrap_all_magic, "");
    def("getMinNeededTransitions", min_needed, "");
}


