import MySQLdb
db = MySQLdb.connect(read_default_file="~/.my.cnf")
c = db.cursor()

#from 28 we use 2 and 4 (Q1 min / all)
#from 30 we can use all 
exp_key = 30
c.execute("select bin,value from hroest.result_diffhist where exp_key = %s and type = %s" % 
          (exp_key, 1))


#from 36, use 2 and 4
#from 37, use 1 and 3
#from 38, use 1 and 3
nr = c.execute("select bin , sum(value) from hroest.result_diffhist where exp_key = 38 and type = 3 and bin > -10 and bin < 10 group by floor(bin*4) ")
nr = c.execute("select bin , sum(value) from hroest.result_diffhist where exp_key = 38 and type = 1 group by floor(bin*2.0) ")
res = c.fetchall()
#bar( [r[0] for r in res] , [r[1] for r in res])
sum( [r[1] for r in res])
n = [r[1] for r in res]
bins = [r[0] for r in res]
bins.append( -bins[0])
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
fig = plt.figure()
ax = fig.add_subplot(111)
# get the corners of the rectangles for the histogram
left = np.array(bins[:-1])
right = np.array(bins[1:])
bottom = np.zeros(len(left))
top = bottom + n
# we need a (numrects x numsides x 2) numpy array for the path helper
# function to build a compound path
XY = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T
# get the Path object
barpath = path.Path.make_compound_path_from_polys(XY)
# make a patch out of it
patch = patches.PathPatch(barpath, facecolor='black', edgecolor='black', alpha=0.8)
ax.add_patch(patch)
# update the view limits
ax.set_xlim(left[0], right[-1])
ax.set_ylim(bottom.min(), top.max())
#title('Q3 distribution (best) %s bins'  % nr)
xlabel( 'Difference / Th')
plt.show()

plt.savefig( 'q1mindistr_exp38_50.pdf', format='pdf')







nr = c.execute("select total_UIS, uisorder from hroest.result_srmuis where exp_key = 32")

nr = c.execute(
    """
    select total_UIS-non_useable_UIS  
    from hroest.result_srmuis  r
    inner join hroest.srmPeptides_yeast sp on sp.parent_id = r.parent_key
    inner join ddb.peptide p on p.id = sp.peptide_key
    inner join ddb.peptideOrganism o on p.id = o.peptide_key
    where exp_key = 32 and uisorder = 4
    #and genome_occurence = 1
    and length(p.sequence) >= 8
    """
)
uisall = c.fetchall()

clf()
#hist(uisall, 1000)
tmp = hist(uisall, 100, [0,8000], color='black')
#tmp = hist(uisall, 100, [0,8000])
#hist(uisall, 1000, [0,20000])
mean(uisall)
median(uisall)


savefig('UIS3distr_go1.pdf' , format='pdf')

savefig('UIS4distr_32_length7.pdf' , format='pdf')


clf()
uis5 = [u[0] for u in uisall if u[1] == 5]
hist(uis5, 100, [0,20000])

hist(uis5, 1000, [15000,16000])

for i in range(1,5):

    myuis = [u[0] for u in uisall if u[1] == i]
    myhist = numpy.histogram( myuis, 101, [0,100])
    len(myuis)
    #plt.plot(  log(myhist[1][:-1]), myhist[0] )
    plt.plot(  myhist[1][:-1], myhist[0] )

plt.xlim(0, 10000)





#from 26 we can use 2 and 4 (Q3 min / all)
exp_key = 26
c.execute("select value from hroest.result_diffdistr where exp_key = %s and type = %s" % 
          (exp_key, 4))
res = c.fetchall()

clf()
hist(res, 400)

clf()
hist(res, 500, [-5,5])


title('Q3 distribution (all)')
xlabel('Difference / Th')
savefig( 'q3alldistr_smallq1_all400.pdf', format='pdf')



#from 27 we can use 1 and 3 (Q1 min / all)
exp_key = 26



c.execute("select value from hroest.result_diffdistr where exp_key = 25 and type = 2")
res = c.fetchall()
clf()
hist(res, 800)

hist(res, 800, [-2,2])





c.execute("select value from hroest.result_diffdistr where exp_key = 25 and type = 3")
res = c.fetchall()
clf()
hist(res, 800)

clf()
hist(res, 150)
title('Q1 distribution (all)')
xlabel('Difference / Th')
savefig( 'q1alldistr.pdf', format='pdf')



c.execute("select value from hroest.result_diffdistr where exp_key = 25 and type = 4")
res = c.fetchall()

clf()
hist(res, 200)

title('Q3 distribution (all)')
xlabel('Difference / Th')

savefig( 'q3alldistr200.pdf', format='pdf')
































n = 5
N = 13


for n in range(6):
    #all fragments (it could be that the last fragment
    #has all the C13s)
    s = 0
    #for m in range(4,N): s += (1-1.0*n/N)**m
    p = 1.0 * n/N
    for m in range(4,N): s += (1-p)**(m-1)*p*m
    print n, s / (N-1-3)


n = 2
p = 1.0 * n/N
for m in range(4,N): 
    binoms = [(1-p)**(m-i)*(1.0*p)**i*choose(m,i) for i in range(m)]
    for b in binoms[:5]: b
    sum(  binoms )
    print "888"


for m in range(1,N-n+1): m





def choose(i,r):
    if r == i: return 1
    if r == 0: return 1
    assert i > 0
    assert r > 0
    assert i >= r
    return reduce( lambda x,y: x*y, range(1,i+1) ) * 1.0 / ( reduce( lambda x,y: x*y, range(1,r+1) ) * reduce( lambda x,y: x*y, range(1,i-r+1) ) ) 



