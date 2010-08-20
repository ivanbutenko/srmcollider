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

