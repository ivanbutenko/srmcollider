//include our own libraries
#include "srmcollider.h"

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

//M is the order
//N is the length of the input vector

/*
* This function calculates all combinations of M elements
* drawn without replacement from a set of N elements. Order of elements
* does NOT matter.
*
* M is the order
* N is the length of the input vector
* mapping is an array of 
*/
void _combinations_magic(int M, int N, COMBINT* mapping,
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

