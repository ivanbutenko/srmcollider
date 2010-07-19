library(RMySQL)
dbh <- dbConnect( MySQL(),default.file='/etc/httpd/conf/.my.cnf')
library(plotrix, lib.loc = '/projects/compep/R/')
library(hexbin)

df <- dbGetQuery(dbh, "
select * from hroest.srmClashesHuman
where experiment_key = 3061 
and q1 between 400 and 1500 
and q3 between 400 and 1500 
")

#plot the q1-q3 space using hexbins
pdf(file="./results/q1_q3_hexbins_120_human_small.pdf" )
hxbinobj = hexbin(x=df[,7], y=df[,8], xbins=110 )
plot( hxbinobj, xlab="", ylab="" )
dev.off()

pdf(file="./results/q1_q3_hexbins_1200_human_small.pdf" )

hxbinobj = hexbin(x=df[,7], y=df[,8], xbins=1100 )

png(file="./results/q1_q3_hexbins_1200_human_small_bigfile.png", 
    width = 5000, height = 5000, units = "px", res = 500 )
plot( hxbinobj , xlab="", ylab="")
dev.off()
pdf(file="./results/q1_q3_hexbins_1200_human_small.pdf" )
plot( hxbinobj , xlab="", ylab="")
dev.off()

