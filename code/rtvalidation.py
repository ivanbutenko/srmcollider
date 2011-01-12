import MySQLdb
#from numpy import *
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
db = MySQLdb.connect(read_default_file="~/.my.cnf")
cursor = db.cursor()
cursor.execute( """
select sequence, Tr, file_id from ddb.peptide  inner join hroest.orinner_scores on
peptide.sequence = orinner_scores.pepseq where experiment_key =
3394 and pg_rank = 1
order by sequence, Tr;
               """)
rtpeps_byfile = cursor.fetchall()


cursor.execute( """ select peptide.sequence, ssrcalc from ddb.peptide  
inner join compep.ssrcalc_prediction on ssrcalc_prediction.sequence = peptide.sequence 
where experiment_key = 3394 """)
allpeps = cursor.fetchall()
pepdic = dict( allpeps )


#p = allpeps[0]
for p in allpeps:
    tr = [ r[1] for r in rtpeps_byfile if r[0] == p[0] ]
    if len(tr) == 0: continue
    #n, bins, patches = plt.hist(tr, 20, normed=1, facecolor='green', alpha=0.75)
    n, bins, patches = plt.hist(tr, 10, alpha=0.75)
    plt.xlabel('Retention time')
    plt.ylabel('Occurence')
    plt.title(r'RT distribution for %s' % p[0])
    plt.grid(True)
    #savefig( 'fig_%s.pdf' % p, format='pdf' )
    #plt.clf()


#use the RT peptides to get a regression curve for each MS run
regvalues = [ 0 for i in range(28) ]
for i in range(1,28):
    rtpeps = [ [r[0], r[1], r[1] / 60.0] for r in rtpeps_byfile if r[2] == i]
    corr = np.array([ [r[1] / 60.0, pepdic[r[0]]] for r in rtpeps])
    x = corr[:,0]
    y = corr[:,1]
    A = np.vstack([x, np.ones(len(x))]).T
    #we fit y = m*x + c
    #Rsquare is 1-SSE/SST 
    #where SSE is the squared residuals
    #and SST is the squared deviation from the mean
    [m, c], SSE, d, dummy = np.linalg.lstsq( A, y )
    SST = sum( ( y-np.mean( y )) * ( y-np.mean( y )) )
    Rsquared = 1- SSE/SST
    regvalues[i] = [m,c]
    #print Rsquared



#boxplot of the m and intersect values
plt.clf()
regvalues[0] = [0,0]
plt.boxplot( np.array(regvalues)[:,0][1:] )
plt.savefig( 'boxplot_m.pdf', format='pdf' )
plt.clf()
plt.boxplot( np.array(regvalues)[:,1][1:] )
plt.savefig( 'boxplot_intersect.pdf', format='pdf' )








#select 1% FDR
#which is at 0.08
#but now only get the non-decoy peptides
cursor.execute( """
select pepseq, ssrcalc, Tr_min , run_id
from hroest.orinner_peakgroups       
inner join compep.ssrcalc_prediction on ssrcalc_prediction.sequence = pepseq      
where m_score < 0.08 and pg_rank = 1 
and decoy = 1
               """)
measured_peps = cursor.fetchall()


#without correction
x = [pep[1] for pep in measured_peps]
y = [pep[2] for pep in measured_peps]
A = np.vstack([x, np.ones(len(x))]).T
#we fit y = m*x + c
#Rsquare is 1-SSE/SST 
#where SSE is the squared residuals
#and SST is the squared deviation from the mean
[m, c], SSE, d, dummy = np.linalg.lstsq( A, y )
SST = sum( ( y-np.mean( y )) * ( y-np.mean( y )) )
Rsquared = 1- SSE/SST
plt.clf()
plt.plot(x, y, 'o', label='Original data') 
plt.plot(x, m*x + c, 'r', label='Fitted line')
plt.xlabel( 'SSRCalc')
plt.ylabel( 'RT')
plt.title( 'Correlation without correction: $R^2 = %s$' % Rsquared[0] )
plt.legend()
#plt.show()
plt.savefig( 'corr_uncorrected.pdf', format='pdf' )
print "Created corr_uncorrected.pdf"



#now start the correction
#first find all that belong together
by_runnr = [[] for i in range(28)]
for mpep in measured_peps:
    filename = mpep[3][30:]
    run_nr = int(filename.split('_')[0])
    by_runnr[run_nr].append( mpep)


x = []
y = []
for i in range(1,28):
    tr = np.array( [t[2] for t in by_runnr[i]] )
    ssr = np.array( [t[1] for t in by_runnr[i]] )
    m, c = regvalues[i]
    x.extend( tr * m + c)
    y.extend( ssr )


x = np.array( x )
y = np.array( y )
A = np.vstack([x, np.ones(len(x))]).T
#we fit y = m*x + c
#Rsquare is 1-SSE/SST 
#where SSE is the squared residuals
#and SST is the squared deviation from the mean
[m, c], SSE, d, dummy = np.linalg.lstsq( A, y )
SST = sum( ( y-np.mean( y )) * ( y-np.mean( y )) )
Rsquared = 1- SSE/SST
plt.clf()
plt.plot(x, y, 'o', label='Original data') 
plt.plot(x, m*x + c, 'r', label='Fitted line')
plt.xlabel( 'SSRCalc')
plt.ylabel( 'RT')
plt.title( 'Correlation with correction: $R^2 = %s$' % Rsquared[0] )
plt.legend()
#plt.show()
plt.savefig( 'corr_corrected.pdf', format='pdf' )
print "Created corr_corrected.pdf"





#get a probability distribution from SSRCalc to RT
count  = 0
rdiffs = []
sdiffs = []
z = []
rt_window = [-5,5]
for i in range(len(measured_peps)-1):
    for j in range(i+1, len(measured_peps) ):
        real_diff = measured_peps[i][2] - measured_peps[j][2]
        #if abs(real_diff) > diffrange: continue
        #if real_diff < rt_window[0] or real_diff > rt_window[1]:continue
        ssrcalc_diff = measured_peps[i][1] - measured_peps[j][1]
        rdiffs.append( real_diff )
        sdiffs.append( ssrcalc_diff )
        z.append( (measured_peps[i][2] + measured_peps[j][2])/2.0 )


#answer 2 questions
#how big does the SSRCalc window have to be, in order to find
#at least x % of all true interactions 
len( [yy for yy in y if yy > rt_window[0] and yy < rt_window[1] ])
len( [yy for yy in y if yy > rt_window[0] and yy < rt_window[1] ]) * 1.0 / len( y)

def count_rates(win_size, ssrcalc_size, rdiffs, sdiffs):
    TP = 0
    FP = 0
    TN = 0
    FN = 0
    for i in range(len(rdiffs)):
        if abs(rdiffs[i]) < win_size:
            if abs(sdiffs[i]) < ssrcalc_size: TP += 1
            else: FN += 1
        else:
            if abs(sdiffs[i]) < ssrcalc_size: FP += 1
            else: TN += 1
    return TP, FP, TN, FN

def get_ROC_area(fprs, tprs):
    sum = 0
    for x in range(len(fprs)-1):
        sum += (fprs[x+1] - fprs[x]) * (tprs[x+1] + tprs[x]) / 2.0
    return sum


#calculate the ROC curves for different real window sizes
#The Question:
#   How large does the SSRCalc window need to be get at least x
#   peptides inside the real window and get at most y FP
plt.clf()
win_sizes = [ 0.25 / 2.0, 0.5 / 2.0, 1 / 2.0, 2 / 2.0,  3 / 2.0 ,  5 / 2.0, 
             10.0 / 2.0, 
            # 20.0 / 2,
            # 50.0 / 2,
            ]

allfprs = []
alltprs = []
allssrcalc = []
alllabel = []
for win_size in win_sizes:
    ssrcalcs  = []
    fprs  = []
    tprs = []
    for ssrcalc_size in [i*0.1 for i in range(101)]:
        TP, FP, TN, FN = count_rates(win_size, ssrcalc_size, rdiffs, sdiffs)
        FPR = FP*1.0 / (FP + TN)
        TPR = TP*1.0 / (TP + FN)
        fprs.append( FPR )
        ssrcalcs.append( ssrcalc_size  )
        tprs.append( TPR )
    #for the correct ROC plot we need to have it end at (1,1)
    fprs.append( 1.0 )
    tprs.append( 1.0 )
    allfprs.append( fprs)
    alltprs.append( tprs[:])
    allssrcalc.append( ssrcalcs )
    ROC_area = get_ROC_area(fprs, tprs)
    mylabel = 'RT +/- %0.1f min (AUC = %0.4f)' % (win_size, ROC_area )
    alllabel.append( mylabel )
    print ROC_area, win_size

#make the first plot
for i in range(len(allfprs)):
    fprs = allfprs[i]
    tprs = alltprs[i]
    ssrcals = allssrcalc[i]
    mylabel = alllabel[i]
    #plot ROC curves
    x = np.array( fprs )
    y = np.array( tprs )
    plt.plot(x, y, label = mylabel )

plt.plot(x, x, label='Random (AUC = 0.5)')
plt.xlabel( 'False Positive Rate')
plt.ylabel( 'True Positive Rate')
plt.title( 'SSRCalc predicting two peptides to be in a RT window')
plt.axis([0, 1, 0, 1])
font = font_manager.FontProperties(size=12)
plt.legend(loc='best', prop=font)
plt.grid(True)
plt.savefig( 'SSRCalc_ROC.pdf', format='pdf' )
plt.clf()
print "Created SSRCalc_ROC.pdf"

#make the second plot
for i in range(len(allfprs)):
    fprs = allfprs[i]
    tprs = alltprs[i]
    ssrcals = allssrcalc[i]
    mylabel = alllabel[i]
    tprs.pop()
    x = np.array( ssrcalcs )
    y = np.array( tprs )
    plt.plot(x, y, label = mylabel )


plt.xlabel( 'SSRCalc value')
plt.ylabel( 'True Positive Rate')
plt.title( 'SSRCalc predicting two peptides to be in a RT window')
font = font_manager.FontProperties(size=12)
plt.legend(loc='best', prop=font)
plt.grid(True)
plt.savefig( 'SSRCalc_TPR.pdf', format='pdf' )
plt.clf()
print "Created SSRCalc_TPR.pdf"



#calculate the correlation between SSRcalc and real RT
x = np.array( rdiffs  )
y = np.array( sdiffs )
A = np.vstack([x, np.ones(len(x))]).T
#we fit y = m*x + c
#Rsquare is 1-SSE/SST 
#where SSE is the squared residuals
#and SST is the squared deviation from the mean
[m, c], SSE, d, dummy = np.linalg.lstsq( A, y )
SST = sum( ( y-np.mean( y )) * ( y-np.mean( y )) )

plt.clf()
plt.plot(x, y, 'o', label='Original data') 
plt.plot( [min(x), max(x) ], [min(x)*m + c, max(x)*m + c ], 'r', label='Fitted line')
plt.xlabel( '$\Delta$ RT / min')
plt.ylabel( r'$\Delta$ SSRCalc')
plt.title( 'Difference-Plot: $R^2 = %0.4f$' % Rsquared[0] )
#plt.legend(loc='best')
#plt.show()
plt.savefig( 'difference_correlation.pdf', format='pdf' )
plt.savefig( 'difference_correlation.png', format='png' )
print "Created difference_correlation.pdf"

#now we normalize the SSRcalc values so that they correspond to minutes
#calculate a histogram of the difference between real and predicted value
plt.clf()
norm_y = (y-c)/m
norm_diff = (norm_y - x) #/ z
#hist( norm_diff , 50, [-15,15]) 
n, bins, patches = plt.hist( norm_diff , 50, normed=True )
plt.title('Difference between RT and SSRCalc')
# add a 'best fit' line
mu = np.mean( norm_diff )
sigma = np.std( norm_diff )
line_y = mlab.normpdf( bins, mu, sigma)
l = plt.plot(bins, line_y, 'r--', linewidth=1)
plt.xlabel('$\Delta$ RT - norm($\Delta$ SSRCalc)')
plt.axis([-15, 15, 0, round(max(n)+0.01, 2) ])
plt.grid(True)
#plt.show()
plt.savefig( 'difference_histogram.pdf', format='pdf' )
print "Created difference_histogram.pdf"


#compare all against all
rsvalues = []
for i in range(1,27):
    rtpeps1 = [ [r[0] , r[1] / 60.0] for r in rtpeps_byfile if r[2] == i]
    rtpeps1 = dict( rtpeps1)
    for j in range(i+1,28):
        corr = np.array([ [rtpeps1[ r[0] ],  r[1] / 60.0] for r in rtpeps_byfile if r[2] == j])
        #
        x = corr[:,0]
        y = corr[:,1]
        A = np.vstack([x, np.ones(len(x))]).T
        #we fit y = m*x + c
        #Rsquare is 1-SSE/SST 
        #where SSE is the squared residuals
        #and SST is the squared deviation from the mean
        [m, c], SSE, d, dummy = np.linalg.lstsq( A, y )
        SST = sum( ( y-np.mean( y )) * ( y-np.mean( y )) )
        Rsquared = 1- SSE/SST
        rsvalues.append( Rsquared )

# the histogram of the data
plt.clf()
n, bins, patches = plt.hist(rsvalues, 50, facecolor='green', alpha=0.75)
plt.xlabel(r'$R^2$')
plt.ylabel('Occurence')
plt.title(r'$R^2$ values of RT peptides in 27 MS runs')
plt.grid(True)
plt.savefig( 'Rsqaured.pdf', format='pdf' )
print "Created Rsquared.pdf"





"""

###########################################################################
mu, sigma = 100, 15
x = mu + sigma*np.random.randn(10000)

# the histogram of the data
n, bins, patches = plt.hist(x, 50, normed=1, facecolor='green', alpha=0.75)

# add a 'best fit' line
y = mlab.normpdf( bins, mu, sigma)
l = plt.plot(bins, y, 'r--', linewidth=1)

plt.xlabel('Smarts')
plt.ylabel('Probability')
plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
plt.axis([40, 160, 0, 0.03])
plt.grid(True)

#plt.show()
plt.savefig( 'test.pdf', format='pdf' )
"""

