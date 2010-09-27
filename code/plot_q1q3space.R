library(RMySQL)
dbh <- dbConnect( MySQL(),default.file='/etc/httpd/conf/.my.cnf')
library(plotrix, lib.loc = '/projects/compep/R/')
library(hexbin)
result_dir <- "/home/hroest/srm_clashes/results/yeast/"

df <- dbGetQuery(dbh, "
select q1, q3 
from hroest.srmPeptides_new 
inner join hroest.srmTransitions_new on parent_id = parent_key
where 
q1 between 400 and 1500 
and q3 between 400 and 1500 
")

#plot the q1-q3 space using hexbins -- first the small one
pdf(file=paste(result_dir,"q1_q3_hexbins_120_small.pdf", sep="") )
hxbinobj = hexbin(x=df[,1], y=df[,2], xbins=110 )
plot( hxbinobj, xlab="", ylab="" )
dev.off()

##pdf(file=paste(result_dir,"q1_q3_hexbins_1200_small.pdf", sep="") )

hxbinobj = hexbin(x=df[,1], y=df[,2], xbins=1100 )

#now the big one with 1200x1200 bins
png(file=paste(result_dir,"q1_q3_hexbins_1200_small_bigfile.png", sep=""),
    width = 3500, height = 3500, units = "px", res = 500 )
plot( hxbinobj , xlab="", ylab="")
dev.off()

pdf(file=paste(result_dir,"q1_q3_hexbins_1200_small.pdf", sep="") )
plot( hxbinobj , xlab="", ylab="")
dev.off()




###########################################################################
###########################################################################
# plot the Q1 / RT space
#nice colourfull 2D histogram
df <- dbGetQuery(dbh, "
select ssrcalc, q1
from hroest.srmPeptides_yeast_allCAM 
where
q1 between 400 and 1500 
and ssrcalc between 0 and 80
")


#first calculate...
library(gplots)
ybins = round( (1500 - 400) / 15 ) 
q1rthist = hist2d(df, same.scale=FALSE, nbins=c(80,  ybins), 
    col=gray((4:0)/4), 
    xlab='SSRCalc', ylab ='Q1 (m/z)', 
    main='RT-Q1 Histogram, 15 Da bins ', labels=TRUE)


#then draw...
pdf(file=paste(result_dir,"q1_rt_hist_15Da.pdf", sep="") )
image( q1rthist$x,q1rthist$y, q1rthist$counts , col=rainbow(50) ,
    xlab='SSRCalc', ylab ='Q1 (m/z)', 
    main='RT-Q1 Histogram, 15 Da bins ' )
dev.off()

pdf(file=paste(result_dir,"q1_rt_hist_15Da_gray.pdf", sep="") )
image( q1rthist$x,q1rthist$y, q1rthist$counts , 
    xlab='SSRCalc', ylab ='Q1 (m/z)', 
    main='RT-Q1 Histogram, 15 Da bins ',
col=gray((50:0)/50) )
dev.off()

#with legend but countours
h2d = q1rthist
pdf(file=paste(result_dir,"q1_rt_hist_15Da_cont.pdf", sep="") )
filled.contour( h2d$x, h2d$y, h2d$counts, nlevels=50,
    xlab='SSRCalc', ylab ='Q1 (m/z)', 
    main='RT-Q1 Histogram, 15 Da bins ', 
    col=rainbow(50) 
    )
dev.off()




