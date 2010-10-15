#include <Python.h>
#include <ctime>
#include <iostream>
#include <fstream>

using namespace std;

void Display(int vi[], int size, ofstream &myfile)
{
    for(size_t i=0; i<size; ++i)
        myfile<<vi[i]<<",";
    myfile<<endl;
}



void _combinations(int M, int N, ofstream &myfile) {
    int j, k;
    int* index = new int[M];
    //initialize with number from 0 to M = range( M )
    for(int k=0;k<M;k++) index[k] =  k ;
    while (index[0] < N-M) {
        Display( index , M, myfile);
        //cout << index[0] << N-M << endl;
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
    Display( index , M, myfile);
    delete [] index;
}

void _combinations_wrapper(int M, int N, const char* filename) {
    ofstream myfile;
    myfile.open (filename);
    _combinations( M, N, myfile);
    myfile.close();
}


/*
int main()
{
    cout << "starting \n";

    clock_t start, finish;
    start = clock();

    _combinations_wrapper(3,6, "test.out");

    cout << "done \n";
    finish = clock();
    cout <<  finish - start << endl;
    cout << CLOCKS_PER_SEC  << endl;
    cout << ( (finish - start) * 1.0 /CLOCKS_PER_SEC ) << endl;

return 0;
}

*/

void initcombinations() {;}



#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(h_combinations)
{
    def("combinations", _combinations_wrapper);
}

