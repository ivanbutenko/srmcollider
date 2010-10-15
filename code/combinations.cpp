/*
 * The functions in this file calculate all combinations of M elements drawn
 * without replacement from a set of N elements.  Order of elements does NOT
 * matter.
 *
 * The function Display is used to write to a ofstream which was defined
 * before. Instead of writing the indices, a string mapping can be provided.
 *
 * The function _combinations is used to calculate the indices and call Display
 *
 * The function _combinations_wrapper can be called from Python. It opens the
 * file and converts the python string-list mapping into a C++ vector<string>
 *
 * Example: M = 2, N = 4, mapping = ['Spam', 'Eggs', 'Bacon', 'Sausage' ]
 * Result:
 * Spam, Eggs
 * Spam, Bacon
 * Spam, Sausage
 * Eggs, Bacon
 * Eggs, Sausage
 * Bacon, Sausage
*/

#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

void Display(int vi[], int size, ofstream &myfile, const vector<string>&
        string_mapping) {
    for(int i=0; i<size-1; ++i)
        myfile << string_mapping[vi[i]] << ",";
    myfile << string_mapping[vi[size-1]];
    myfile << endl;
}

void _combinations(int M, int N, ofstream &myfile, const vector<string>&
        string_mapping) {
    // The basic idea is to create an index array of length M that contains
    // numbers between 0 and N. The indices denote the combination produced and
    // we always increment the rightmost index. If it goes above N, we try to
    // increase the one left to it until we find one that still can be
    // increased. We stop when the rightmost index hits N-M, we thus go from
    // (0,1,2,3...M-1) to (N-M,N-M+1,N-M+2...N-1)
    int j, k;
    int* index = new int[M];
    //initialize with numbers from 0 to M = range( M )
    for(int k=0;k<M;k++) index[k] = k;
    while (index[0] < N-M) {
        Display( index, M, myfile, string_mapping);
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
    Display( index, M, myfile, string_mapping);
    delete [] index;
}



/* int main()
{
    cout << "starting \n";

    clock_t start, finish;
    start = clock();

    vector<string> mystrings(4);
    mystrings[0] =  "Spam";
    mystrings[1] = "Eggs";
    mystrings[2] =  "Bacon";
    mystrings[3] =  "Sausage";

    ofstream myfile;
    myfile.open ("test.out");
    _combinations( 2, 4, myfile, mystrings);
    myfile.close();

    cout << "done \n";
    finish = clock();
    cout <<  finish - start << endl;
    cout << CLOCKS_PER_SEC  << endl;
    cout << ( (finish - start) * 1.0 /CLOCKS_PER_SEC ) << endl;

return 0;
} */






/*
 *
 *
 * PYTHON interface
*/

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

void _combinations_wrapper(int M, int N, const char* filename,
        boost::python::list mapping) {
    ofstream myfile;
    myfile.open (filename);

    int len = boost::python::extract<int>(mapping.attr("__len__")());
    //Do Error checking, the mapping needs to be at least as long as N 
    if (len < N) {
        PyErr_SetString(PyExc_ValueError, 
            "The string mapping must be at least of length N");
        boost::python::throw_error_already_set();
        return;
    }
    //Also M should not be bigger than N (not much point) and not smaller than 0
    if(M>N) {
        PyErr_SetString(PyExc_ValueError, 
            "M cannot be larger than N");
        boost::python::throw_error_already_set();
        return;
    }
    if(N < 0 || M < 0) {
        PyErr_SetString(PyExc_ValueError, 
            "M and N need to be larger than 0");
        boost::python::throw_error_already_set();
        return;
    }

    vector<string> mystrings(len);
    for (int i=0; i<len; i++) {
        mystrings[i] = boost::python::extract<char const *>(mapping[i]);
    }

    _combinations( M, N, myfile, mystrings);
    myfile.close();
}

void initcombinations() {;}

using namespace boost::python;
BOOST_PYTHON_MODULE(c_combinations)
{
    def("combinations", _combinations_wrapper);
}

