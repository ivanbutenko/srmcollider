#include <iostream>
#include <ctime>
#include <vector>
#include <algorithm>

using namespace std;

void Display(const vector<int>& vi)
{
    for(size_t i=0; i<vi.size(); ++i)
        cout<<vi[i]<<",";
    cout<<endl;
}



void _combinations(int M, int N) {
    vector<int> index;
    //initialize with number from 0 to M = range( M )
    int j, k;
    for(int k=0;k<M;k++) index.push_back( k );
    while (index[0] < N-M) {
        Display( index );
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
    Display( index );
}


int main()
{
    cout << "starting \n";

    clock_t start, finish;
    start = clock();

    _combinations(3,6);
    cout << "done \n";
    finish = clock();
    cout <<  finish - start << endl;
    cout << CLOCKS_PER_SEC  << endl;
    cout << ( (finish - start) * 1.0 /CLOCKS_PER_SEC ) << endl;
return 0;
}
