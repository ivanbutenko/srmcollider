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
ssr <- dbGetQuery(dbh, "
select ssrcalc
from hroest.srmPeptides_yeast_allCAM 
")

ssrhist <-hist(ssr$ssrcalc)

###########################################################################
###########################################################################
# plot the Q1 / RT space
#nice colourfull 2D histogram
df <- dbGetQuery(dbh, "
select ssrcalc, q1
from hroest.srmPeptides_yeast_allCAM 
where
q1 between 400 and 1200 
and ssrcalc between 0 and 80
")






#first calculate...
library(gplots)
ybins = round( (1200 - 400) / 25 ) 
xbins = round( 80 / 0.5 )
pdf(file=paste(result_dir,"test.pdf", sep="") )
q1rthist = hist2d(df, same.scale=FALSE, nbins=c(xbins,  ybins), 
    #col=gray((4:0)/4), 
    col=rainbow(50),
    xlab='SSRCalc', ylab ='Q1 (m/z)', 
    main='RT-Q1 Histogram, 25 Da bins ', labels=TRUE)
dev.off()


#then draw...
pdf(file=paste(result_dir,"q1_rt_hist_25Da.pdf", sep="") )
image( q1rthist$x,q1rthist$y, q1rthist$counts , col=rainbow(50) ,
    xlab='SSRCalc', ylab ='Q1 (m/z)', 
    main='RT-Q1 Histogram, 25 Da bins ' )
dev.off()


count_hist = hist(q1rthist$counts)
#then draw...
pdf(file=paste(result_dir,"q1_rt_counts_25Da.pdf", sep="") )
hist(q1rthist$counts, breaks=50)
dev.off()


pdf(file=paste(result_dir,"rt_hist.pdf", sep="") )
ssrhist <-hist(ssr$ssrcalc)
dev.off()

pdf(file=paste(result_dir,"test.pdf", sep="") )
plot(hh$breaks, tt)

dev.off()




pdf(file=paste(result_dir,"q1_rt_hist_25Da_gray.pdf", sep="") )
image( q1rthist$x,q1rthist$y, q1rthist$counts , 
    xlab='SSRCalc', ylab ='Q1 (m/z)', 
    main='RT-Q1 Histogram, 25 Da bins ',
col=gray((50:0)/50) )
dev.off()

#with legend but countours
h2d = q1rthist
pdf(file=paste(result_dir,"q1_rt_hist_25Da_cont.pdf", sep="") )
filled.contour( h2d$x, h2d$y, h2d$counts, nlevels=50,
    xlab='SSRCalc', ylab ='Q1 (m/z)', 
    main='RT-Q1 Histogram, 25 Da bins ', 
    col=rainbow(50) 
    )
dev.off()







df.data1 <- dbGetQuery(dbh, "
select nr_interfering from hroest.result_tmp_bins_2ssrcalc_25Da_yeast 
")
df.data2 <- dbGetQuery(dbh, "
select nr_interfering from hroest.result_tmp_bins_2ssrcalc_25Da_yeast_pa
")
df.data3 <- dbGetQuery(dbh, "
select nr_interfering from hroest.result_tmp_bins_0_5_ssrcalc_25Da_yeast 
")
df.data4 <- dbGetQuery(dbh, "
select nr_interfering from hroest.result_tmp_bins_0_5_ssrcalc_25Da_yeast_pa
")

pdf(file=paste(result_dir,"boxplot.pdf", sep="") )
boxplot(df.data1$nr_interfering, df.data2$nr_interfering,  df.data3$nr_interfering, df.data4$nr_interfering,
names=c("2 (yeast)", "2 (yeast PA)", "0.5 (yeast)", "0.5 (yeast PA)"),
 main="Interfering peptides per precursor")
dev.off()

pdf(file=paste(result_dir,"boxplot_05.pdf", sep="") )
boxplot( df.data3$nr_interfering, df.data4$nr_interfering,
names=c("0.5 (yeast)", "0.5 (yeast PA)"),
 main="Interfering peptides per precursor")
dev.off()

pdf(file=paste(result_dir,"histInter_2ssrcalc_yeast.pdf", sep="") )
hist(df.data1$nr_interfering, breaks = 20, xlab ="Number of Interfering Peptides",
 main="Interfering peptides per precursor (2 SSRCalc, yeast)")
dev.off()

pdf(file=paste(result_dir,"histInter_2ssrcalc_yeast_pa.pdf", sep="") )
hist(df.data2$nr_interfering, breaks = 20, xlab ="Number of Interfering Peptides",
 main="Interfering peptides per precursor (2 SSRCalc, yeast PA)")
dev.off()


pdf(file=paste(result_dir,"histInter_05ssrcalc_yeast.pdf", sep="") )
hist(df.data3$nr_interfering, breaks = 20, xlab ="Number of Interfering Peptides",
 main="Interfering peptides per precursor (0.5 SSRCalc, yeast)")
dev.off()

pdf(file=paste(result_dir,"histInter_05ssrcalc_yeast_pa.pdf", sep="") )
hist(df.data4$nr_interfering, breaks = 20, xlab ="Number of Interfering Peptides",
 main="Interfering peptides per precursor (0.5 SSRCalc, yeast PA)")
dev.off()
