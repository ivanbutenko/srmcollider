/*
 *
 * Program       : SRMCollider
 * Author        : Hannes Roest <roest@imsb.biol.ethz.ch>
 * Date          : 05.02.2011 
 *
 *
 * Copyright (C) 2011 - 2012 Hannes Roest
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
 * The functions in this file calculate all combinations of M elements
 * drawn without replacement from a set of N elements. Order of elements
 * does NOT matter.
 *
 * The function Display is used to write to a ofstream which was defined
 * before. Instead of writing the indices, a string mapping can be
 * provided.
 *
 * The function _combinations is used to calculate the indices and call
 * Display
 *
 * The function _combinations_wrapper can be called from Python. It
 * opens the file and converts the Python arguments:
 * # M: number of elements to be picked
 * # N: pool of elements
 * # filename: output filename
 * # mapping: (list of strings) needs to be supplied to map the output
 *   indices to meaningful strings.
 * # exclude: (list of lists of integers) a list with 'exceptions'. Each
 *   list consists of integers, has length M and defines a combination
 *   to be skipped. The indices start at 0 and end at M-1 and need to be
 *   ordered.
 *
 * Example: M = 2, N = 4, mapping = ['Spam', 'Eggs', 'Bacon', 'Sausage' ] 
 * Call: combinations(2, 4,'outfile', ['Spam', 'Eggs', 'Bacon', 'Sausage'], [])
 * Result:
 * Spam, Eggs
 * Spam, Bacon
 * Spam, Sausage
 * Eggs, Bacon
 * Eggs, Sausage
 * Bacon, Sausage
 *
 * Example: M = 2, N = 4, mapping = ['Spam', 'Eggs', 'Bacon', 'Sausage' ] 
 *          exclude = [ [0,1], [1,2] ]
 * Call: combinations(2, 4,'outfile', ['Spam', 'Eggs', 'Bacon', 'Sausage'],  
 *                    [ [0,1], [1,2] ] )
 * Result:
 * Spam, Bacon
 * Spam, Sausage
 * Eggs, Sausage
 * Bacon, Sausage
 *
 *
 * Signature
 *void combinations(int,int,string,list<string>,list<list<int> >)
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
        string_mapping, const vector<vector<int> >& exclude) {
    // The basic idea is to create an index array of length M that contains
    // numbers between 0 and N. The indices denote the combination produced and
    // we always increment the rightmost index. If it goes above N, we try to
    // increase the one left to it until we find one that still can be
    // increased. We stop when the rightmost index hits N-M, we thus go from
    // (0,1,2,3...M-1) to (N-M,N-M+1,N-M+2...N-1)
    int j, k;
    bool found = false;
    int* index = new int[M];
    //initialize with numbers from 0 to M = range( M )
    for(int k=0;k<M;k++) index[k] = k;
    while (index[0] < N-M) {
        /* Do we need to skip this one, e.g. is it in the exclude?  */
        for(std::vector<int>::size_type i = 0; i != exclude.size(); i++) {
            found = true;
            for(std::vector<int>::size_type j = 0; j != exclude[i].size(); j++) {
                /* std::cout << someVector[i]; ... */
                if( exclude[i][j] != index[j] ) found = false;
            }
            if(found) break;
        }
        if(!found) Display( index, M, myfile, string_mapping);

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
    //TODO also check whether the last one is here!
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

    vector<vector<int> > exclude(0);

    ofstream myfile;
    myfile.open ("test.out");
    _combinations( 2, 4, myfile, mystrings, exclude);
    myfile.close();

    cout << "done \n";
    finish = clock();
    cout <<  finish - start << endl;
    cout << CLOCKS_PER_SEC  << endl;
    cout << ( (finish - start) * 1.0 /CLOCKS_PER_SEC ) << endl;

return 0;
}  */






/*
 *
 *
 * PYTHON interface
*/

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

void _combinations_wrapper(int M, int N, const char* filename,
        boost::python::list mapping,
        boost::python::list exclude) {
    /*
     * M and and N are the combinations parameter and will produce M choose N results
     * filename is the name of the output file
     * mapping is a list of strings of size N which maps the indices to output
     * exclude is a list of lists (a list of indices to exclude)
     * 
    */
    ofstream myfile;
    boost::python::list inner_list;
    myfile.open (filename);

    int mapping_length = boost::python::extract<int>(mapping.attr("__len__")());
    //Do Error checking, the mapping needs to be at least as long as N 
    if (mapping_length < N) {
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

    /* Convert the mapping */
    vector<string> mystrings(mapping_length);
    for (int i=0; i<mapping_length; i++) {
        mystrings[i] = boost::python::extract<char const *>(mapping[i]);
    }

    /* Convert the exclude list */
    int exclude_length = boost::python::extract<int>(exclude.attr("__len__")());
    vector< vector<int> > exclude_indices(exclude_length, vector<int>(M,0));    
    for (int i=0; i<exclude_length; i++) {
        inner_list = boost::python::extract< boost::python::list >(exclude[i]);
        int inner_length = boost::python::extract<int>(inner_list.attr("__len__")());

        /* exclude list needs to be the lenght of the indices we produce */
        if(inner_length != M) {
            PyErr_SetString(PyExc_ValueError, 
                "In c_combinations.combinations module... \n"
                "The exclude list needs to consist of a list of indices you \n"
                "want to exclude. It is thus a list of lists where each inner\n"
                "element has length M.");
            boost::python::throw_error_already_set();
            return;
        }

        for (int j=0; j<inner_length; j++) {
            exclude_indices[i][j] = boost::python::extract<int>(inner_list[j]);
        }
    }

    /*

    import c_combinations, string, profile
    #c_combinations.combinations(7, 30, "c.out", [l for l in string.letters[:30]], [ range(7),range(7) ])
    c_combinations.combinations(7, 30, "c.out", [l for l in string.letters[:30]], range(7))
    c_combinations.combinations(7, 30, "c.out", [l for l in string.letters[:30]], [ [ [1,2,4] ] ])
    c_combinations.combinations(  2,  4, 'outfile', ['Spam', 'Eggs', 'Bacon', 'Sausage' ],  [  [0,1], [1,2] ] )
    print c_combinations.combinations.__doc__

    profile.run('c_combinations.combinations(7, 30, "c.out", [l for l in string.letters[:30]], [ range(7) for i in range(100) ])')

    c_combinations.combinations(2, 5, "c.out", [l for l in string.letters[:30]],
        [ [0,1],  [1,2], [2,3 ], [1,3] ] )

    */

    _combinations( M, N, myfile, mystrings, exclude_indices);
    myfile.close();
}

void initcombinations() {;}

using namespace boost::python;
BOOST_PYTHON_MODULE(c_combinations)
{
    def("combinations", _combinations_wrapper, 

 " The functions in this file calculate all combinations of M elements\n"
 " drawn without replacement from a set of N elements. Order of elements\n"
 " does NOT matter.\n"
 "\n"
 " The function combinations can be called from Python. It\n"
 " opens the file and converts the Python arguments:\n"
 " # M: number of elements to be picked\n"
 " # N: pool of elements\n"
 " # filename: output filename\n"
 " # mapping: (list of strings) needs to be supplied to map the\n"
 "   output to meaningful strings. \n"
 " # exclude: (list of lists of integers) a list with 'exceptions'. Each\n"
 "   list consists of integers, has length M and defines a combination\n"
 "   to be skipped. The indices start at 0 and end at M-1 and need to be\n"
 "   ordered.\n"
 "\n"
 " Example: M = 2, N = 4, mapping = ['Spam', 'Eggs', 'Bacon', 'Sausage' ] \n"
 " Call: combinations(2, 4,'outfile', ['Spam', 'Eggs', 'Bacon', 'Sausage'], [])\n"
 " Result:\n"
 " Spam, Eggs\n"
 " Spam, Bacon\n"
 " Spam, Sausage\n"
 " Eggs, Bacon\n"
 " Eggs, Sausage\n"
 " Bacon, Sausage\n"
 "\n"
 " Example: M = 2, N = 4, mapping = ['Spam', 'Eggs', 'Bacon', 'Sausage' ] \n"
 "          exclude = [ [0,1], [1,2] ]\n"
 " Call: combinations(2, 4,'outfile', ['Spam', 'Eggs', 'Bacon', 'Sausage'],  \n"
 "                    [ [0,1], [1,2] ] )\n"
 " Result:\n"
 " Spam, Bacon\n"
 " Spam, Sausage\n"
 " Eggs, Sausage\n"
 " Bacon, Sausage\n"
 "\n"
 "\n"
 " Signature\n"
 "void combinations(int,int,string,list<string>,list<list<int> >)\n"

           );
}

