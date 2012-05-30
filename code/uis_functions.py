"""
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
"""

# Get a list of all non-UIS combinbations (using c++ if possible)
def get_nonuis_list(collisions_per_peptide, MAX_UIS):
    non_uis_list = [set() for i in range(MAX_UIS+1)]
    try:
        import c_getnonuis
        for order in range(1,MAX_UIS+1):
                non_uis_list[order] = c_getnonuis.get_non_uis(
                    collisions_per_peptide, order)
    except ImportError:
        for pepc in collisions_per_peptide.values():
            for i in range(1,MAX_UIS+1):
                get_non_uis(pepc, non_uis_list[i], i)
    return non_uis_list 

def get_non_uis(pepc, non_uis, order):
    if len( pepc ) >= order: 
        non_uis.update( [tuple(sorted(p)) for p in combinations(pepc, order)] )

def choose(i,r):
    assert i > 0
    assert r > 0
    assert i >= r
    if r == i: return 1
    return reduce( lambda x,y: x*y, range(1,i+1) ) * 1.0 / ( 
        reduce( lambda x,y: x*y, range(1,r+1) ) * 
        reduce( lambda x,y: x*y, range(1,i-r+1) ) ) 

# get the UIS (unique ion signatures) from a set of non-UIS
def get_uis(srm_ids, non_uis, order):
    #prepare input, make sure its sorted
    non_uis = [ tuple(sorted(list(i))) for i in non_uis]
    srm_ids = sorted( srm_ids )
    #
    #the output of combinations is guaranteed to be sorted, if the input was sorted
    result = combinations(srm_ids, order)
    result = [r for r in result if not r in non_uis]
    return result

def _combinations(N, M):
    """All index combinations of M elements drawn without replacement
     from a set of N elements.
    Order of elements does NOT matter."""
    # TODO test
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

def combinations(iterable, r, force_non_native=False):
    """ All combinations of size r of the given iterable set. 
    Order of elements does NOT matter.
    """
    try:
        import itertools
        if force_non_native: import itertools_dummy_module
        return itertools.combinations( iterable , r) 
    except:
        #we are before python 2.6
        return [ tuple([iterable[i] for i in indices]) for indices in _combinations( len(iterable) , r)] 

def _combinationsDiffLen(N):
    """All index combinations of M elements drawn without replacement
     from a set of N elements.
    Order of elements does NOT matter."""
    M = len(N)
    index = [0 for i in range( M )]
    while True:
        yield index[:]
        #print index[:] 
        #kk += 1
        #if index[2] == 73: break
        index[ M-1 ] += 1
        #if kk > 7: break
        if index[ M-1 ] >= N[ M-1 ]:
            #now we hit the end, need to increment other positions than last
            j = M-1
            while j >= 0 and index[j] >= N[j]-1: j -= 1
            #j contains the value of the index that needs to be incremented
            #when we are at the end of the interation, j will be -1
            if j <0: break
            index[j] += 1
            k = j + 1
            #set all other positions to zero again
            while k < M: index[k] = 0; k += 1; 

