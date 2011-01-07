
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

namespace python = boost::python;

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


python::list _find_clashes_calculate_clashes(python::tuple precursors,
        double q3_low, double q3_high ) {

    python::tuple tlist;
    python::list result;

    int precursor_length = python::extract<int>(precursors.attr("__len__")());
    long peptide_key;
    int j, k, start;
    double q1, acc_mass, res_mass;
    char* sequence;
    char c;
    bool inside;

    for (int i=0; i<precursor_length; i++) {
        tlist = python::extract< python::tuple >(precursors[i]);

        acc_mass = 0.0;
        res_mass = 0.0;
        python::list fragment_series;
        python::list b_series;
        python::list y_series;

        q1 = python::extract< double >(tlist[0]);
        peptide_key = python::extract< long >(tlist[2]);
        sequence = python::extract<char *>(tlist[1]);

        /* Python code
         *
            for q in re.finditer( '([A-Z]\[\d*\]|[A-Z])', seq):
                element = q.group(0)
                res_mass = R.residues[element][1]
                self.mass += res_mass
                fragment_series.append( self.mass )

            self.b_series = [b + R.mass_H for b in fragment_series[:-1]] 
            self.y_series = [self.mass - y + 2*R.mass_H + R.mass_OH for y 
                in fragment_series[:-1]]
        */
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
                            return result;
                            }
                        res_mass = 147.03540462; break;
                    case 'C': 
                        if(!(sequence[start+2] == '1' && 
                             sequence[start+3] == '6' && 
                             sequence[start+4] == '0' )) {
                            PyErr_SetString(PyExc_ValueError, 
                                "Unknown modification for cysteine");
                            boost::python::throw_error_already_set();
                            return result;
                        }
                        res_mass = 160.030653721; break;
                    default: 
                        PyErr_SetString(PyExc_ValueError, 
                            "Unknown modification ");
                        boost::python::throw_error_already_set();
                        return result;
                }
                //'M[147]':  131.04049 + mass_O), # oxygen
                //'C[160]':  103.00919 + mass_CAM - mass_H ), # CAM replaces H

                acc_mass += res_mass;
                fragment_series.append( acc_mass );
                b_series.append( acc_mass + MASS_H );
                y_series.append( - acc_mass + 2* MASS_H + MASS_OH );
                inside = false;
            }
            else if(inside) { }
            else {
                //We found a regular AA
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
                        return result;
                }

                fragment_series.append( res_mass );
                acc_mass += res_mass;
                b_series.append( acc_mass + MASS_H );
                y_series.append( - acc_mass + 2* MASS_H + MASS_OH );
            }
        }

        int y_length = python::extract<int>(y_series.attr("__len__")());
        for (int j=0; j<y_length; j++) y_series[j] += acc_mass;

        for(int ch=1; ch<=2; ch++) {
            for (int j=0; j<y_length-1; j++){
                /*
                cout << python::extract<double>(y_series[j]) << " ";
                cout << python::extract<double>(b_series[j]) << endl;
                */

                double q3 = (python::extract<double>(y_series[j]) + (ch-1)*MASS_H)/ch;
                if (q3 > q3_low && q3 < q3_high)
                    result.append(python::make_tuple(q3, q1, 0, peptide_key) );
            }

            for (int j=0; j<y_length-1; j++){
                double q3 = (python::extract<double>(b_series[j]) + (ch-1)*MASS_H)/ch;
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

    return result;
}        


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

            if(fabs(t0-c0) <  this_min) {

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
    int j, k, tmpint;
    int* index = new int[M];
    //int* mytuple = new int[M];
    python::tuple tmptuple;
    //initialize with numbers from 0 to M = range( M )
    for(int k=0;k<M;k++) index[k] = k;
    while (index[0] <= N-M) {

        //EVALUATE THE RESULT
        python::list tmplist;
        for(int i=0; i<M; ++i) tmplist.append( mapping[index[i]] ); 
        tmplist.sort();
        //this doesnt work for some reason
        //tmptuple = python::make_tuple( tmplist );
        switch(M) {
            case 1: tmptuple = python::make_tuple( tmplist[0] ); break;
            case 2: tmptuple = python::make_tuple( tmplist[0], tmplist[1] ); break;
            case 3: tmptuple = python::make_tuple( tmplist[0], tmplist[1],
                            tmplist[2]); break;
            case 4: tmptuple = python::make_tuple( tmplist[0], tmplist[1], 
                            tmplist[2], tmplist[3]); break;
            case 5: tmptuple = python::make_tuple( tmplist[0], tmplist[1], 
                            tmplist[2], tmplist[3], tmplist[4]); break;
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


python::dict get_non_uis(python::dict collisions_per_peptide, int order) {

    python::list tmplist;
    python::list pepcollisions = python::extract< python::list >(
            collisions_per_peptide.values() );
    python::dict result;

    int tmp_length;
    int collision_length = python::extract<int>(pepcollisions.attr("__len__")());

    for (int i=0; i<collision_length; i++) {
        tmplist = python::extract< python::list >( pepcollisions[i] );
        int tmp_length = python::extract<int>(tmplist.attr("__len__")());
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

void initgetnonuis() {;}
void initcore_non_unique() {;}
void initget_non_uis() {;}

using namespace python;
BOOST_PYTHON_MODULE(c_getnonuis)
{

    def("getnonuis", _getnonuis_wrapper, 
 "getnonuis(tuple transitions, tuple collisions, double q3window, bool ppm)\n"
           );

    def("get_non_uis", get_non_uis, 
 "void get_non_uis(dict collisions_per_peptide, int order)\n"
           );

    def("core_non_unique", _find_clashes_core_non_unique, 
 "core_non_unique(tuple transitions, tuple collisions, double q3window, bool ppm)\n"
           );

    def("calculate_transitions", _find_clashes_calculate_clashes, 
 ""
           );

}

