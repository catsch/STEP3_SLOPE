
Xing11 <- function(Depth, Fluo, Ed490,  FluoRange=c(0.1,1000), DepthRange=c(0,1000), Ed490Range=c(0.1,10000), F490Lim=5, Method="Old") {

KW2 <- 0.0166

X2 <- 0.08253

e2 <- 0.62588

if(Method=="New") {

X2 <-  0.0886

e2 <- 0.8281

}

ReEd490 <- NA*Depth

ReKd490 <- NA*Depth

FlagF490 <- NA

F490 <- NA

CEd490 <- log(Ed490)

CEd490[is.infinite(CEd490)] <- NA

len <- length(Depth)

FlagRegression <- c(1:len*0+1)

C2n <- c(1:len*NA)

W2 <- Depth*KW2

A2n <- CEd490 + W2

ResDC <- list(NA,NA,NA,NA)

if(!all(is.na(CEd490)) & !all(is.na(Fluo))) {

Fluo[is.na(Fluo)] <- 0

Fluo[Fluo<0] <- 0



C2n[1] <- Depth[1]*Fluo[1]^e2

for(im in 2:len) {

C2n[im] <- C2n[im-1]+(Depth[im]-Depth[im-1])*Fluo[im]^e2

}

# Intergration of Kw



FlagRegression[is.na(CEd490)] <- 0

FlagRegression[is.na(C2n)] <- 0

FlagRegression[!is.na(Fluo) & Fluo<FluoRange[1]] <- 0

FlagRegression[!is.na(Fluo) & Fluo>FluoRange[2]] <- 0

FlagRegression[Depth<DepthRange[1]] <- 0

FlagRegression[Depth>DepthRange[2]] <- 0

FlagRegression[!is.na(CEd490) & Ed490<Ed490Range[1]] <- 0

FlagRegression[!is.na(CEd490) & Ed490>Ed490Range[2]] <- 0

ResDC <- CombinedX(A2n,C2n,FlagRegression)

F490 <- round((-1*ResDC[[1]]/X2)^(1/e2),3)

FlagF490 <- ResDC[[3]]

}

if(!is.na(F490) & F490>F490Lim) { 

F490 <- NA

FlagF490 <- 4

}

if(!is.na(F490)) ReEd490 <- exp(ResDC[[1]]*C2n + ResDC[[2]] - W2)

if(!is.na(F490)) ReKd490 <- KW2 + X2*(F490*Fluo)^e2

return(list(F490=F490, FlagF490=FlagF490, ReEd490=ReEd490, ReKd490=ReKd490, Cn=C2n, An=A2n, Slope=ResDC[[1]], Intercept=ResDC[[2]], FlagReg=FlagRegression))

}

########################################

CombinedX <- function(An,Cn,FlagAn) {

# CombinedX is a function to detect cloud part in the Ed(490) profile and find the Slope between An and Cn through linear regression analysis

# An is the lnEd(z) + Intergration of Kw

# Cn is the intergration of Fluo^e

# FlagAn is the flag of Ed(z), All points with FlagAn=0 are removed when applying regression analysis

Flag1 <- FlagAn 

QCFLAG <- 1

SLOPE <- NA

NP <- NA

DIFFL <- c(1:length(An)*NA)

try(NP <- coef(lm(An[FlagAn==1] ~ Cn[FlagAn==1])))

if(is.na(NP)) {

SLOPE <- NA 

QCFLAG <- 4

} else {

SP1 <- 1

SP2 <- 1

while(1>0) {

P <- NP

# DIFFL is the difference between A and regressed A derived from a regressed linear function of C

DIFFL[FlagAn==1] <- An[FlagAn==1]-P[1]-P[2]*Cn[FlagAn==1] 

# If the An is much lower than the regressed line (-0.05), the point is treated as cloud-contaminated and flag it as 0.

FlagAn[!is.na(DIFFL) & DIFFL<(-0.05)] <- 0

try(NP <- coef(lm(An[FlagAn==1] ~ Cn[FlagAn==1])))

if (is.na(NP)) break

#SP1 and SP2 is two parameters to examine the difference between the last regressed results (slope and intercept) and the former one

SP1 <- abs((NP[1]-P[1])/P[1])

SP2 <- abs((NP[2]-P[2])/P[2])

if (is.na(SP1) || is.na(SP2) ) break

# If both two relative differences are lower than 5%, the function is considered to finish the cloud detection and the linear regression is already stable

if(SP1<=0.05 && SP2<=0.05) break

}

}

WP <- NA

try(WP <- coef(lm(An[FlagAn==1] ~ Cn[FlagAn==1])))

#If the regressed slope is positive (which should be negative) or the good points (FlagA==1) is less than 10, the regression analysis is treated as failed, then the SLOPE is NA.

if(is.na(WP) || is.na(WP[2]) || WP[2]>=0 || length(An[FlagAn==1])<10) {

SLOPE <- NA

QCFLAG <- 4

} else {

sqr <- 0

sqr <- summary(lm(An[FlagAn==1] ~ Cn[FlagAn==1]))$r.squared

#If the r2 is lower than 0.98, the regression analysis is treated as failed 

if(sqr<0.98) {

SLOPE <- NA 

QCFLAG <- 4

} else SLOPE <- WP[2]

}

# if the detected cloud-contaminated points are more than 50%, the QCFlag is 3

if(QCFLAG==1 & (length(FlagAn[FlagAn==1])/length(Flag1[Flag1==1]))<0.5) QCFLAG <- 3

# if the detected cloud-contaminated points are less than 50%, but more than 10%, the QCFlag is 2

if(QCFLAG==1 & (length(FlagAn[FlagAn==1])/length(Flag1[Flag1==1]))>=0.5 & (length(FlagAn[FlagAn==1])/length(Flag1[Flag1==1]))<0.9) QCFLAG <- 2

if(is.na(SLOPE)) WP <- NA

return(list(SLOPE, WP[1], QCFLAG, FlagAn))

}