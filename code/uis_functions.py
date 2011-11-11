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

def permutations(iterable, r=None):
    #use itertools in 2.6
    #see
    #http://svn.python.org/view/python/trunk/Modules/itertoolsmodule.c?view=markup&pathrev=81889
    import itertools
    try:
        return itertools.permutations( iterable, r)
    except:
        return _permutations( iterable, r)

def _permutations(iterable, r=None):
    # permutations('ABCD', 2) --> AB AC AD BA BC BD CA CB CD DA DB DC
    # permutations(range(3)) --> 012 021 102 120 201 210
    #print iterable, r
    pool = tuple(iterable)
    n = len(pool)
    if r is None: r = n 
    else: r=r
    indices = range(n)
    cycles = range(n, n-r, -1)
    yield tuple(pool[i] for i in indices[:r])
    while n:
        for i in reversed(range(r)):
            cycles[i] -= 1
            if cycles[i] == 0:
                indices[i:] = indices[i+1:] + indices[i:i+1]
                cycles[i] = n - i
            else:
                j = cycles[i]
                indices[i], indices[-j] = indices[-j], indices[i]
                yield tuple(pool[i] for i in indices[:r])
                break
        else:
            return

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

def combinations(iterable, r):
    """ All combinations of size r of the given iterable set. 
    Order of elements does NOT matter.
    """
    try:
        import itertools
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

