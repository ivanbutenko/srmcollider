import random

class Bins(object):
    def __init__(self, nr, min, max):
        self.bin_nr = nr
        self.max = max
        self.min = min
        self.range = max - min
        evenspace = (max - min) * 1.0 / (self.bin_nr) 
        self.bin_sizes = [ evenspace for i in range(self.bin_nr) ]
    def change_bin_size_locally(self, bin, amount, go_left = True):
        self.bin_sizes[bin] += amount
        if bin == 0:
            self.bin_sizes[bin +1] -= amount 
        elif bin == len(self) -1:
            self.bin_sizes[bin -1] -= amount 
        else: 
            if go_left:
                self.bin_sizes[bin +1] -= amount 
            else:
                self.bin_sizes[bin -1] -= amount
        assert sum(self) - self.max + self.min < 10**(-6)
    def change_bin_size(self, bin, amount):
        tot = sum([s for i,s in enumerate(self) if i != bin] )
        for i,s in enumerate(self):
            if i == bin: continue
            self.bin_sizes[i] -= s/tot * amount
        self.bin_sizes[bin] += amount
        assert sum(self) - self.max + self.min < 10**(-6)
    def __len__(self):
        return self.bin_nr
    def __getitem__(self, x):
        return self.bin_sizes[x]
    def get_boundaries(self):
        res = [self.min]
        for s in self:
            res.append( res[-1] + s)
        return res
    #def find_bin(self, n):
    #    for i, boundary in enumerate(self.distr):
    #        #print i, boundary
    #        if boundary > n: return i-1
    #    return self.bin_nr-1
    #def bin_sizes(self):
    #    res = [self.distr[i+1] - self.distr[i] for i in range(self.bin_nr-1)]
    #    res.append( self.max - self.distr[-1])
    #    return res
    #
    #def change_bin_size(self, bin, amount):
    #    #self.distr[ bin ] -= amount /2.0
    #    old_size = self.distr[bin+1]  - self.distr[bin]
    #    change_others = amount * 1.0 / (self.bin_nr -1)
    #    #for i in range(bin+2, self.bin_nr):
    #    #for i in range(self.bin_nr-1, bin+1, -1):
    #    #    self.distr[i] += amount/2.0 - change_others  * (i-bin-1)
    #    for i in range(1, bin+1):
    #        self.distr[i] -= change_others * (i)
    #    #now change the actual bin we want to resize
    #    self.distr[ bin +1 ] = self.distr[bin] + amount + old_size

class Window_finder(object):
    def __init__(self, values, windows, min, max):
        self.bins = Bins(windows, min, max)
        self.values = values
    def run(self, iterations=1000,overlap=2,max_change_amount=100,eval_every=100,
            temperature=1):
        e = 2.7182818284590451
        ord_values = self.values
        ord_values.sort()
        b = self.bins
        f = self.evaluate_f( b, ord_values, overlap)
        self.f =f
        #factor_t = 1.0 - temperature
        if not self.__dict__.has_key('min'):
            self.min = f
            self.best = b[:]
        for i in range(iterations):
            old = b[:]
            old_f = f
            if i % 3 == 0:
                #select a random bin and change it a random amount, distribute
                #the remainder equally over all other bins
                b.change_bin_size( int(random.random() * len(b) ), (random.random()-0.5) * max_change_amount)
            elif i % 3 == 1: 
                #randomly give to neighbor bin
                #we need this to flatten out hills and values
                direction = random.random()
                if direction < 0.5: direction = True
                else: direction = False
                bin = int( random.random() * len(b)   )
                #f = self.evaluate_f( b, ord_values, overlap)
                #mean = sum( self.elements_perbin) / len( self.elements_perbin)
                #diff = [e - mean for e in self.elements_perbin]
                #print bin
                #print diff
                #print b.bin_sizes
                b.change_bin_size_locally( bin, 
                     random.random() * max_change_amount, go_left=direction)
                #print b.bin_sizes
                #f = self.evaluate_f( b, ord_values, overlap)
                #mean = sum( self.elements_perbin) / len( self.elements_perbin)
                #diff = [e - mean for e in self.elements_perbin]
                #print diff
            else: 
                #find the border between valley and hill and try to exchange
                f = self.evaluate_f( b, ord_values, overlap)
                mean = sum( self.elements_perbin) / len( self.elements_perbin)
                diff = [e - mean for e in self.elements_perbin]
                #print "=" * 75
                k = 0
                for k in range(int(random.random() * len(b) ), len(b)-1):
                    #print k, diff[k] * diff[k+1] 
                    if diff[k] * diff[k+1] < 0: break
                bin = k
                amount = random.random() * max_change_amount
                #print "change ", bin, amount
                if diff[k] > 0: amount = -1.0 * amount
                #print "change ", bin, amount
                #print f
                #print b.bin_sizes
                #mean = sum( self.elements_perbin) / len( self.elements_perbin)
                #print [e - mean for e in self.elements_perbin]
                b.change_bin_size_locally( bin, amount, go_left=True)
                #f = self.evaluate_f( b, ord_values, overlap)
                #print b.bin_sizes
                #print f
                #print [e - mean for e in self.elements_perbin]
            f = self.evaluate_f( b, ord_values, overlap)
            self.f = f
            if f < self.min:
                self.best = b[:]
                self.min = f
            #print old_f, f
            if f > old_f: 
                #revert
                #if random.random() <  old_f *1.0 / f - factor_t:
                #if random.random() > 2.7182818284590451**(- old_f *0.3 / f):
                if random.random() < e**( (old_f *1.0 - f) / temperature):
                    #print "dont change ", f, old_f, 2.7182818284590451**( (old_f *1.0 - f) / temperature)
                    pass
                else:
                    b.bin_sizes = old
                    f = old_f
                    self.f = f
            if i % eval_every == 0:
                f = self.evaluate_f( b, ord_values, overlap)
                mean = sum( self.elements_perbin) / len( self.elements_perbin)
                print "Iteration %s, function is %s, best is %s" % (i, f, self.min), \
                [e - mean for e in self.elements_perbin]

    def set_start(self, start):
        self.bins.bin_sizes = start

    def meta_run(self, start):
        self.bins.bin_sizes = start
        self.results = []
        self.amount = 0.1
        for i in range(100): 
            self.run(max_change_amount=0.1, overlap = 2, iterations=1000)
            self.results.append( self.bins.bin_sizes )

    def evaluate_f(self, b, ord_values, overlap=0):
        """The objective function calculates how evenly distributed the bins are.
        It calculates the squared standard deviation times N of the number of
        elements in each bin. It will take the overlap into account."""
        j = 0 #current index on b
        self.elements_perbin  = [0 for i in range(len(b))]
        last_boundary = b.min
        next_boundary = last_boundary + b[j] - overlap / 2.0
        current_overlap = 0
        elemenents_last_overlap = 0
        in_overlap = False
        prev_i = 0
        for i,v in enumerate(ord_values):
            if v > next_boundary:
                if in_overlap:
                    #we were in overlap, and leave now
                    in_overlap = False
                    #print "leave ov ", next_boundary, j, prev_i, v
                    j += 1
                    if j > len(b): print "end"
                    next_boundary = next_boundary + b[j] - overlap 
                    #add the overlap entries to BOTH sides
                    self.elements_perbin[j] += i - prev_i
                    self.elements_perbin[j-1] += i - prev_i
                    prev_i = i
                else: 
                    in_overlap = True
                    #we enter the overlap
                    #print "enter ov ", next_boundary, j
                    last_boundary = next_boundary 
                    next_boundary = next_boundary + overlap 
                    self.elements_perbin[j] += i - prev_i
                    prev_i = i
        assert(in_overlap == True)  #the last overlap
        j += 1
        self.elements_perbin[j-1] += i - prev_i
        prev_i = i
        #we want to know how equally they are distributed
        #return the squared std deviation times N
        mean = sum( self.elements_perbin) / len( self.elements_perbin)
        return sum([ (e-mean) * (e-mean) for e in self.elements_perbin])

