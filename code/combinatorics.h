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

//include our own libraries
#ifndef COMBINATORICS_H
#define COMBINATORICS_H
#include "srmcollider.h"

//using namespace std;

namespace SRMCollider 
{

  namespace Combinatorics 
  {
    /*
    * This function calculates all combinations of M elements
    * drawn without replacement from a set of N elements. Order of elements
    * does NOT matter.
    *
    * M is the order
    * N is the length of the input vector
    * mapping is an array of integers of fixed length (see srmcollider.h)
    * result is a python dictionary that holds all combinations, they are
    *   implemented as integers with bitflags set or unset
    */
    void _combinations_bitwise(int M, int N, COMBINT* mapping,
            std::set<COMBINT>& result) 
    {
        // The basic idea is to create an index array of length M that contains
        // numbers between 0 and N. The indices denote the combination produced and
        // we always increment the rightmost index. If it goes above N, we try to
        // increase the one left to it until we find one that still can be
        // increased. We stop when the rightmost index hits N-M, we thus go from
        // (0,1,2,3...M-1) to (N-M,N-M+1,N-M+2...N-1)
        //
        // see http://svn.python.org/projects/python/trunk/Modules/itertoolsmodule.c
        // combinations_next fxn around line 2113

        int j, k;
        int* index = new int[M];
        COMBINT tmpres;

        // Speed measurements
        //python::dict t; // 0m28.892s
        //std::map<COMBINT, int> t; // 0m37.133s
        //std::set t; // 0m19.013s -- thus the std::set.insert is the slowest part here
        //std::vector<COMBINT> t; // 0m2.831s
        // do not add = 0m1.227s

        //initialize with numbers from 0 to M = range( M )
        for(int k=0;k<M;k++) index[k] = k;
        while (index[0] <= N-M) {

            //EVALUATE THE RESULT
            tmpres = 0;
            for(k=0;k<M;k++) 
                tmpres |= mapping[index[k]];
            result.insert(tmpres);

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
                //TODO if(j<0) break;
                index[j] += 1;
                //now start from j+1 to the right and set all indices to their
                //lowest possible values
                k = j + 1;
                while (k < M) { index[k] = index[k-1] + 1; k += 1; } 
            }
        }

        delete [] index;
    }

    /*
    * This function calculates all combinations of M elements
    * drawn without replacement from a set of N elements. Order of elements
    * does NOT matter.
    *
    * M is the order
    * N is the length of the input vector
    * mapping is an array of integers of fixed length (see srmcollider.h)
    * result is a python dictionary that holds all combinations, they are
    *   implemented as integers with bitflags set or unset
    */
    void _combinations(int M, int N, std::vector<std::vector<int> > &result) 
    {
        // The basic idea is to create an index array of length M that contains
        // numbers between 0 and N. The indices denote the combination produced and
        // we always increment the rightmost index. If it goes above N, we try to
        // increase the one left to it until we find one that still can be
        // increased. We stop when the rightmost index hits N-M, we thus go from
        // (0,1,2,3...M-1) to (N-M,N-M+1,N-M+2...N-1)
        //
        // see http://svn.python.org/projects/python/trunk/Modules/itertoolsmodule.c
        // combinations_next fxn around line 2113

        int j, k;
        int* index = new int[M];
        std::vector<int> tmpres;

        //initialize with numbers from 0 to M = range( M )
        for(int k=0;k<M;k++) index[k] = k;
        while (index[0] <= N-M) {

            //EVALUATE THE RESULT
            tmpres.clear();
            for(k=0;k<M;k++) 
            {
              tmpres.push_back(index[k]);
            }
            result.push_back(tmpres);

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
                //TODO if(j<0) break;
                index[j] += 1;
                //now start from j+1 to the right and set all indices to their
                //lowest possible values
                k = j + 1;
                while (k < M) { index[k] = index[k-1] + 1; k += 1; } 
            }
        }

        delete [] index;

    }

    /*
    * This function calculates all combinations of M elements
    * drawn without replacement from a set of N elements. Order of elements
    * does NOT matter.
    *
    * the result will be stored as tuple-keys in a python dict
    */
    void _py_combinations(int M, int N, const python::list &mapping, 
            python::dict &result) 
    {

      std::vector<std::vector<int> > result_vec;
      _combinations(M, N, result_vec);

      // convert to python tuple
      python::tuple tmptuple;
      for (size_t i = 0; i < result_vec.size(); i++)
      {

        std::vector<int> index = result_vec[i];

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
            case 6: tmptuple = python::make_tuple( mapping[index[0]],
                            mapping[index[1]], 
                            mapping[index[2]],
                            mapping[index[3]],
                            mapping[index[4]],
                            mapping[index[5]]
                    ); break;
            case 7: tmptuple = python::make_tuple( mapping[index[0]],
                            mapping[index[1]], 
                            mapping[index[2]],
                            mapping[index[3]],
                            mapping[index[4]],
                            mapping[index[5]],
                            mapping[index[6]]
                    ); break;
            case 8: tmptuple = python::make_tuple( mapping[index[0]],
                            mapping[index[1]], 
                            mapping[index[2]],
                            mapping[index[3]],
                            mapping[index[4]],
                            mapping[index[5]],
                            mapping[index[6]],
                            mapping[index[7]]
                    ); break;
            case 9: tmptuple = python::make_tuple( mapping[index[0]],
                            mapping[index[1]], 
                            mapping[index[2]],
                            mapping[index[3]],
                            mapping[index[4]],
                            mapping[index[5]],
                            mapping[index[6]],
                            mapping[index[7]],
                            mapping[index[8]]
                    ); break;
            case 10: tmptuple = python::make_tuple( mapping[index[0]],
                            mapping[index[1]], 
                            mapping[index[2]],
                            mapping[index[3]],
                            mapping[index[4]],
                            mapping[index[5]],
                            mapping[index[6]],
                            mapping[index[7]],
                            mapping[index[8]],
                            mapping[index[9]]
                    ); break;
            default:
                PyErr_SetString(PyExc_ValueError, 
                    "Order (M) larger than 5 is not implemented");
                boost::python::throw_error_already_set();
                return;
        }
        result[ tmptuple ] = 0; //append to result dict

      }
    }

    /*
     * Function to calculate all non-UIS transitions from a dictionary
     * where for each key (colliding peptide key) the set of transitions
     * of the query peptide are stored that are interfered with is held
     * (collisions_per_peptide dictionary)
     * Note that this dictionary contains tuples that represent the combinations.
     *
     * It will return a list of all non UIS of the requested order.
     */
    python::dict get_non_uis(python::dict collisions_per_peptide, int order) 
    {

        python::list tmplist;
        python::list pepcollisions = python::extract< python::list >(
                collisions_per_peptide.values() );
        python::dict result;

        int tmp_length;
        int collision_length = python::extract<int>(pepcollisions.attr("__len__")());

        for (int i=0; i<collision_length; i++) {
            tmplist = python::extract< python::list >( pepcollisions[i] );
            tmp_length = python::extract<int>(tmplist.attr("__len__")());
            _py_combinations(order, tmp_length, tmplist, result);
        }

        return result;

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

    }

    /*
     * Function to calculate all non-UIS transitions from a dictionary
     * where for each key (colliding peptide key) the set of transitions
     * of the query peptide are stored that are interfered with is held
     * (collisions_per_peptide vector).
     * Note that this vector contains integers that represent the combinations.
     *
     * It will return a list of all non UIS of the requested order.
     */
    void get_non_uis_bitwise(std::vector<COMBINT>& newcollperpep, int max_tr, int
            order, std::set<COMBINT> & result) 
    {
        int onecounter;
        COMBINT tmparr;
        COMBINT mask;
        COMBINT* mapping = new COMBINT[max_tr];

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

            _combinations_bitwise(order, onecounter, mapping, result);
        }

        delete [] mapping;
    }
  }
}

#endif
