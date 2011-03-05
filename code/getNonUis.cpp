/*
 * The functions in this 
 *
 * Signature
 *void combinations(int,int,string,list<string>,list<list<int> >)
*/

//include our own libraries
#include "srmcollider.h"
#include "srmcolliderLib.cpp"


python::dict _getnonuis_wrapper(python::tuple transitions, python::tuple
        collisions, double q3window, bool ppm);
python::list _find_clashes_calculate_clashes(python::tuple precursors,
        double q3_low, double q3_high );
double calculate_charged_mass(python::tuple clist, int ch);
python::dict _find_clashes_calculate_collperpeptide(python::tuple transitions,
        python::tuple precursors, double q3_low, double q3_high, double q3window, bool ppm);
python::list _calculate_clashes_wrapper(python::tuple &tlist, double charge);
python::dict _find_clashes_core_non_unique(python::tuple transitions,
        python::tuple collisions, double q3window, bool ppm);

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

    python::dict collisions_per_peptide, tmpdict;
    python::tuple clist;
    python::tuple tlist;
    python::list tmplist;

    int transitions_length = python::extract<int>(transitions.attr("__len__")());
    int precursor_length = python::extract<int>(precursors.attr("__len__")());
    int fragcount, i, j, k, ch, listmembers = 0, tmplen;

    long t1, peptide_key, tmplong;
    bool already_in_list;
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
                            t1 = python::extract<long>(tlist[1]);
                            tmpdict[t1] = 0; //dummy dict
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


