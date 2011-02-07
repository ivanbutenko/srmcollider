import unittest
import time
import test_shared

import sys
sys.path.append( '..')
sys.path.append( '/home/hroest/projects/' )
import collider
import silver
import DDB 

try:
    import c_getnonuis
except ImportError:
    print "=" * 75, """
Module c_getnonuis is not available. Please compile it if you want to use it.
""", "=" * 75


def runcpp(self, pep, transitions, precursors):
    #first we run the C++ version
    st = time.time()
    for kk in range(10):
        tmp = c_getnonuis.calculate_collisions_per_peptide( 
            transitions, tuple(precursors), self.q3_low, self.q3_high, self.par.q3_window, self.par.ppm)
    ctime = (time.time() - st)/10.0
    return tmp, ctime

def runpy(self, pep, transitions, precursors):
    #now we run the pure Python version
    st = time.time()
    collisions = list(collider.SRMcollider._get_all_collisions_calculate_sub(
            collider.SRMcollider(), precursors, self.R, self.q3_low, self.q3_high))
    collisions_per_peptide = {}
    q3_window_used = self.par.q3_window
    for t in transitions:
        if self.par.ppm: q3_window_used = self.par.q3_window * 10**(-6) * t[0]
        for c in collisions:
            if abs( t[0] - c[0] ) <= q3_window_used:
                #gets all collisions
                if collisions_per_peptide.has_key(c[3]):
                    if not t[1] in collisions_per_peptide[c[3]]:
                        collisions_per_peptide[c[3]].append( t[1] )
                else: 
                    collisions_per_peptide[c[3]] = [ t[1] ] ; 
    oldtime = time.time() - st
    return collisions_per_peptide, oldtime

class Test_speed(unittest.TestCase):

    def setUp(self):
        try:
            import c_getnonuis
            self.transitions = test_shared.transitions_def1
            self.collisions = test_shared.collisions_def1
            self.pep1 = test_shared.peptide1
            self.pep2 = test_shared.peptide2

            self.transitions_12_between300_1500 = test_shared.transitions_12_between300_1500
            self.pep1_yseries = test_shared.pep1_yseries
            self.pep1_bseries = test_shared.pep1_bseries

            tuples = []
            tuples.append(self.pep1)
            tuples.append(self.pep2)
            self.tuples = tuple( tuples)

            self.R = silver.Residues.Residues('mono')

            class Minimal: pass
            self.par = Minimal()
            self.par.q3_window = 4.0
            self.par.ppm = False
            self.q3_high = 1500
            self.q3_low = 300
            self.MAX_UIS = 5

        except ImportError:
            pass

    def test_collperpeptide1(self):
        print '\nTesting collisions_per_peptide'
        print "old vs new time"

        pep = test_shared.runpep1
        transitions = test_shared.runtransitions1
        precursors = test_shared.runprecursors1
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])

        tmp, ctime = runcpp(self, pep, transitions, precursors)
        collisions_per_peptide, oldtime = runpy(self, pep, transitions, precursors)
        self.assertEqual(collisions_per_peptide, tmp )

        print ctime, oldtime
        print "Speedup:", oldtime / ctime

    def test_collperpeptide2(self):
        print '\nTesting collisions_per_peptide'
        print "old vs new time"

        pep = test_shared.runpep2
        transitions = test_shared.runtransitions2
        precursors = test_shared.runprecursors2
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])

        tmp, ctime = runcpp(self, pep, transitions, precursors)
        collisions_per_peptide, oldtime = runpy(self, pep, transitions, precursors)
        self.assertEqual(collisions_per_peptide, tmp )

        print ctime, oldtime
        print "Speedup:", oldtime / ctime
        
    def test_get_non_uis1(self):
        print '\nTesting non_uis_list'
        print "old vs new time"

        pep = test_shared.runpep1
        transitions = test_shared.runtransitions1
        precursors = test_shared.runprecursors1
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        collisions_per_peptide, ctime = runcpp(self, pep, transitions, precursors)

        MAX_UIS = self.MAX_UIS

        non_uis_list = [set() for i in range(MAX_UIS+1)]
        #here we calculate the UIS for this peptide with the given RT-range
        st = time.time()
        for pepc in collisions_per_peptide.values():
            for i in range(1,MAX_UIS+1):
                collider.get_non_uis(pepc, non_uis_list[i], i)
        oldtime = time.time() - st

        st = time.time()
        for kk in range(10):
            non_uis_list_new = [{} for i in range(MAX_UIS+1)]
            for order in range(1,MAX_UIS+1):
                non_uis_list_new[order] = c_getnonuis.get_non_uis(
                    collisions_per_peptide, order)

        non_uis_list_new = [set( res.keys() ) for res in non_uis_list_new]
        ctime = (time.time() - st)/10

        self.assertEqual(non_uis_list_new, non_uis_list)
        print ctime, oldtime
        print "Speedup:", oldtime / ctime

    def test_get_non_uis2(self):
        print '\nTesting non_uis_list'
        print "old vs new time"

        pep = test_shared.runpep2
        transitions = test_shared.runtransitions2
        precursors = test_shared.runprecursors2
        transitions = tuple([ (t[0], i) for i,t in enumerate(transitions)])
        collisions_per_peptide, ctime = runcpp(self, pep, transitions, precursors)

        MAX_UIS = self.MAX_UIS

        non_uis_list = [set() for i in range(MAX_UIS+1)]
        #here we calculate the UIS for this peptide with the given RT-range
        st = time.time()
        for pepc in collisions_per_peptide.values():
            for i in range(1,MAX_UIS+1):
                collider.get_non_uis(pepc, non_uis_list[i], i)
        oldtime = time.time() - st

        st = time.time()
        for kk in range(10):
            non_uis_list_new = [{} for i in range(MAX_UIS+1)]
            for order in range(1,MAX_UIS+1):
                non_uis_list_new[order] = c_getnonuis.get_non_uis(
                    collisions_per_peptide, order)

        non_uis_list_new = [set( res.keys() ) for res in non_uis_list_new]
        ctime = (time.time() - st)/10

        self.assertEqual(non_uis_list_new, non_uis_list)
        print ctime, oldtime
        print "Speedup:", oldtime / ctime

    def test_calculatetrans(self):
        pep = {'q1_charge': 2, 'q1': 417.21177749999998, 'sequence': 'HSISFDK', 'peptide_key': 11036419L, 'parent_id': 214997L, 'ssrcalc': 15.279999999999999}

        q3_high = self.q3_high
        q3_low = self.q3_low
        par = self.par 

        p_id = pep['parent_id']
        q1 = pep['q1']
        transitions = c_getnonuis.calculate_transitions(
                ((q1, pep['sequence'], p_id),), q3_low, q3_high)

        tnew = [tt[0] for tt in transitions]
        tnew.sort()
        mycollider = collider.SRMcollider()
        transitions = mycollider._get_all_transitions(par, pep, cursor)
        told = [tt[0] for tt in transitions]
        told.sort()

        #old and new should be similar
        for o,n in zip(told, tnew):
            assert abs(o -n) < 1e-5


#f = open('/tmp/tst.py', 'w')
#f.write('collpepresult = {')
#for k,v in collisions_per_peptide.iteritems():
#    f.write('%s : %s,\n' % (k,repr(v)))
#
#f.write('}')
#f.close()
#

if __name__ == '__main__':
    unittest.main()

