library(RMySQL)
dbh <- dbConnect( MySQL(),default.file='/etc/httpd/conf/.my.cnf')
library(plotrix, lib.loc = '/projects/compep/R/')
library(hexbin)

df <- dbGetQuery(dbh, "select * from hroest.srmClashes
where experiment_key = 3061 
and q1 between 300 and 1500 
and q3 between 300 and 1500 
#limit 100000
" )

#plot the q1-q3 space using hexbins
pdf(file="./results/q1_q3_hexbins_120_human.pdf" )
hxbinobj = hexbin(x=df[,7], y=df[,8], xbins=120 )
plot( hxbinobj, xlab="", ylab="" )
dev.off()
pdf(file="./results/q1_q3_hexbins_1200_human.pdf" )
hxbinobj = hexbin(x=df[,7], y=df[,8], xbins=1200 )
plot( hxbinobj , xlab="", ylab="")
dev.off()

