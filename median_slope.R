###################################
# median calculation of the slope
###################################
uf=commandArgs()

mission  <- uf[2]

file=paste("./",mission,"_SLOPE.txt",sep="")

SLOPE=read.table(file,sep=" ")

JULD=SLOPE$V2
SLOPE_value=SLOPE$V3
QC=SLOPE$V4

ind_good=which(QC=="2" | QC=="1" )
#ind_good=which( QC=="1" )

SLOPE_OUT=median(SLOPE_value[ind_good])
SLOPE_SD=sd(SLOPE_value[ind_good])
print(length(SLOPE_value[ind_good]))

SLOPE_table=data.frame(mission,SLOPE_OUT,SLOPE_SD)
write.table(file="SLOPE_stats.txt",SLOPE_table,col.names=F,row.names=F,append=TRUE)


