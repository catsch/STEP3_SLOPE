########################################################################
## SLOPE estimation for CHLA based on radiometry profiles at 490
## with Xiaogang Xing routines  
##
## Inputs : Dark and quenching corrected CHLA profiles and Ed490 profiles
## 
## Outputs: Estimation of the slope of CHLA
## 
## Version : 26 January 2021
########################################################################

library(ncdf4)
library(stringr)

source("./Xing11.R")
source("./RunningFilter.R")

uf=commandArgs()

mission  <- uf[2]

#### Creating the list of files for which we need to recompute  
liste_to_do=read.table("liste_all_B",header=FALSE, as.is=TRUE)

# List of the file to process
LIST_nc=liste_to_do$V1

for (IDnc in LIST_nc) {

### Initialising the Ed490 FLAG

	FLAG_NO_Ed490=FALSE

#### Getting the name of the core file
	file_in_C=str_replace(IDnc,"/B","/")

# if B and C are not in the same mode 
	if (!file.exists(file_in_C)) file_in_C=str_replace(file_in_C,"profiles/R","profiles/D")
	if (!file.exists(file_in_C)) file_in_C=str_replace(file_in_C,"profiles/D","profiles/R")

#########################
# Open the nc File
#########################

# Open the C file
	filenc_C=nc_open(file_in_C,readunlim=FALSE,write=FALSE)

# Open the B file
	filenc_B=nc_open(IDnc,readunlim=FALSE,write=FALSE)

############################################################
## work on index 
#############################################################

#### Get the list of parameters in the profile

	STATION_PARAMETERS=ncvar_get(filenc_B,"STATION_PARAMETERS")

# Stations parameters has a fixed length 64 characters 

	CHLA_STRING=str_pad("CHLA",64,"right")

	Ed490_STRING=str_pad("DOWN_IRRADIANCE490",64,"right")

# Find the profile containing CHLA and PAR  

	index_chla=which(STATION_PARAMETERS == CHLA_STRING, arr.ind=TRUE)

	index_Ed490=which(STATION_PARAMETERS == Ed490_STRING, arr.ind=TRUE)	

	if (length(index_chla)==0 | length(index_Ed490)==0) {

		next # jump to the next profile if no chla or no Ed490 in the profile

	} else {

		iprof_chla=index_chla[,2]

		iprof_Ed490=index_Ed490[,2]

	} 

###############################################################################
######## Jump to next profile or dont take PAR into account if PRES_QC is too bad
###############################################################################

	PROFILE_PRES_CTD_QC=strsplit(ncvar_get(filenc_C,"PROFILE_PRES_QC"),split="")

	if (PROFILE_PRES_CTD_QC[[1]][iprof_chla]=="E" | PROFILE_PRES_CTD_QC[[1]][iprof_chla]=="F" | PROFILE_PRES_CTD_QC[[1]][iprof_chla]=="D") next


###########################################################
# PRESSURE
###########################################################

	PRES_CTD=ncvar_get(filenc_B,"PRES")

	JULD=unique(ncvar_get(filenc_B,"JULD"))

	LATITUDE=unique(ncvar_get(filenc_B,"LATITUDE"))

	LONGITUDE=unique(ncvar_get(filenc_B,"LONGITUDE"))

	CYCLE_NUMBER=unique(ncvar_get(filenc_B,"CYCLE_NUMBER"))

###########################################################
# Enough light measurements
###########################################################

	Ed490=ncvar_get(filenc_B,"DOWN_IRRADIANCE490")
		
	if (length(which(!is.na(Ed490[,iprof_Ed490]))) < 2 ) FLAG_NO_Ed490=TRUE

	PROFILE_ED_QC=str_squish(ncvar_get(filenc_B,"PROFILE_DOWN_IRRADIANCE490_QC"))

############################################################
# Get Ed490 and CHLA on the same levels 
############################################################

#### Working on PAR
	if (!FLAG_NO_Ed490) {
 
		Ed490_CHLA=approx(PRES_CTD[,iprof_Ed490],Ed490[,iprof_Ed490],PRES_CTD[,iprof_chla], rule=1)$y

#### Filtering a bit the PAR 
		MED_Ed490_CHLA=RunningFilter(2,Ed490_CHLA,na.fill=T, ends.fill=F, Method="Median")

	} 

################################################################
#### Getting CHLA 
################################################################ 

###  Getting FLUO_CHLA

	FLUO=ncvar_get(filenc_B,"CHLA_ADJUSTED")

###  Initialising the FLUO_CHLA derived variables  

	MED_FLUO=FLUO # unspiked FLUO

###  PROFILE_CHLA

	PROFILE_CHLA_QC=str_squish(ncvar_get(filenc_B,"PROFILE_CHLA_QC"))

############################################################################################
# Filtering the signal To remove the spikes for BBP and CHLA (Briggs et al.) 
#	1. median filter 5 
#	2. Mean filter 7 
############################################################################################

### 1st median filter 5

	MED_FLUO[,iprof_chla]=RunningFilter(2,FLUO[,iprof_chla],na.fill=T, ends.fill=T, Method="Median")

### 2nd Mean filter 7 for BBP and CHLA 

	MED_FLUO[,iprof_chla]=RunningFilter(3,MED_FLUO[,iprof_chla],na.fill=T, ends.fill=T, Method="Mean")

#### Slope estimation 

	XING=Xing11(PRES_CTD[,iprof_chla], MED_FLUO[,iprof_chla], 100*MED_Ed490_CHLA)

	SLOPE=XING$F490

	FLAG=XING$FlagF490

	if ( PROFILE_CHLA_QC == "F" | PROFILE_ED_QC == "F") FLAG=4

###########################################
# Write the estimation in the output file 
###########################################

	slope=data.frame(CYCLE_NUMBER,JULD,SLOPE,FLAG,LATITUDE,LONGITUDE)

	write.table(file=paste(mission,"_SLOPE.txt",sep=""),slope,col.names=F,row.names=F,append=TRUE)
	
	nc_close(filenc_B)



}# end loop on B files  
