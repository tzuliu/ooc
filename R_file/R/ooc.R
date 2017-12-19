#
# PROGRAM FORTRAN CALL FUNCTIONS
#
find_polytope <- function(nummembers,numvotes,
                         jjj,np,nrcall,ns,ndual,XMAT,LLEGERR, 
                        ZVEC,WS,MCUTS,LERROR,ltotal,mwrong, 
                        LDATA,iprint) {
   res <-.Fortran("kplegis",
     as.integer(nummembers),
     as.integer(numvotes),
     as.integer(jjj),
     as.integer(np),
     as.integer(nrcall),
     as.integer(ns),
     as.integer(ndual),
     as.double(XMAT),
     as.integer(LLEGERR),
     as.double(ZVEC),
     as.double(WS),
     as.integer(MCUTS),
     as.integer(LERROR),
     as.integer(ltotal),
     as.integer(mwrong),
     as.integer(LDATA),
     as.integer(iprint))
}
find_cutpoint <- function(nummembers,numvotes,
                        npzz,nv,np,nrcall,ns,ndual,ivot,XMAT,YSS,MVOTE,WS, 
                        LLV,LLVB,LLE,LLEB, 
                        LERROR,ZS,jch,jeh,jcl,jel,irotc,kcut,lcut, 
                        LLL,XJCH,XJEH,XJCL,XJEL){
res <-.Fortran("jan1pt",
     as.integer(nummembers),
     as.integer(numvotes),
     as.integer(npzz),
     as.integer(nv),
     as.integer(np),
     as.integer(nrcall),
     as.integer(ns),
     as.integer(ndual),
     as.integer(ivot),
     as.double(XMAT),
     as.double(YSS),
     as.integer(MVOTE),
     as.double(WS),
     as.integer(LLV),
     as.integer(LLVB),
     as.integer(LLE),
     as.integer(LLEB),
     as.integer(LERROR),
     as.double(ZS),
     as.integer(jch),
     as.integer(jeh),
     as.integer(jcl),
     as.integer(jel),
     as.integer(irotc),
     as.integer(kcut),
     as.integer(lcut),
     as.integer(LLL),
     as.double(XJCH),
     as.double(XJEH),
     as.double(XJCL),
     as.double(XJEL))
}
#
# CIRCLE SAMPLING FUNCTIONS
#
# library(hitandrun)
# hypersphere.sample(n, N)
#
#
#  PROGRAM FUNCTION "binary.balanced()" TO GENERATE REGULAR YEA/NAY BINARY MATRIX WITH CATEGORIES AS EVENLY BALANCED AS POSSIBLE
#
binary.balanced <- function(var){
#
#  NUMBER OF CATEGORIES IN THE VARIABLE (n)
#
n.categories <- length(table(factor(var, levels = min(var, na.rm=TRUE):max(var, na.rm=TRUE))))
#
#  CALCULATE HOW MANY POSSIBLE BINARY CATEGORIZATIONS (n-1)
#
n.categorizations <- n.categories - 1
#
#  CREATE A TABLE OF DIFFERENCES IN NUMBER OF RESPONDENTS IN EACH BINARY CATEGORIZATION
#  (e.g., for a 7-point scale, the difference between "1" and "2:7", "1:2" and "3:7", etc.)
#
diff.table <- rep(0,n.categorizations)
dim(diff.table) <- c(n.categorizations,1)
for (i in 1:n.categorizations){
diff.table[i,] <- abs(sum(table(factor(var, levels = min(var, na.rm=TRUE):max(var, na.rm=TRUE)))[1:i]) - sum(table(factor(var, levels = min(var, na.rm=TRUE):max(var, na.rm=TRUE)))[(i+1):n.categories]))
}
#
#  WHICH ROW HAS THE MINIMUM VALUE?
#
solution <- which.min(diff.table)
#
#  RECODE ORDINAL VARIABLE INTO MOST-BALANCED BINARY CATEGORIES
#
var.binary <- rep(9,length(var))
for (i in 1:length(var.binary)){
if (is.na(var[i])==FALSE & var[i]>=1 & var[i]<=solution) var.binary[i] <- 1
if (is.na(var[i])==FALSE & var[i]>=(solution+1) & var[i]<=(n.categories)) var.binary[i] <- 6
}
return(var.binary)
}
#
# ===================================
#
#  PROGRAM FUNCTION "binarize()" TO GENERATE ENLARGED BINARY MATRIX WITH DOMINANCE PATTERN
#
binarize <- function(var){
#
#  NUMBER OF CATEGORIES IN THE VARIABLE (n)
#
n.categories <- length(table(var))
#
#  CALCULATE NUMBER OF BINARY CATEGORIES NEEDED (n-1)
#
n.bin.categories <- length(table(var)) - 1
#
#  CALCULATE NUMBER OF RESPONDENTS
#
nresp <- length(var)
#
#  CREATE nresp*(n-1) MATRIX
#
M <- rep(NA,(nresp*n.bin.categories))
dim(M) <- c(nresp,n.bin.categories)
#
#  REPLACE "NA" VALUES WITH 999
#
var[is.na(var)] <- 999
#
#  REPLACE MATRIX
#
for (i in 1:nresp){
for (j in 1:n.bin.categories){
if (var[i] <= j & var[i]!=999) M[i,j] <- 1
if (var[i] > j & var[i]!=999) M[i,j] <- 6
}}
#
M[is.na(M)] <- 9
return(M)
}
#
# ===================================
#
#  PROGRAM FUNCTION "binarize_nondominance()" TO GENERATE ENLARGED BINARY MATRIX WITH 1 VS. 2, 2 VS. 3, ETC.
#
binarize_nondominance <- function(var){
#
#  NUMBER OF CATEGORIES IN THE VARIABLE (n)
#
n.categories <- length(table(var))
#
#  CALCULATE NUMBER OF BINARY CATEGORIES NEEDED (n-1)
#
n.bin.categories <- length(table(var)) - 1
#
#  CALCULATE NUMBER OF RESPONDENTS
#
nresp <- length(var)
#
#  CREATE nresp*(n-1) MATRIX
#
M <- rep(NA,(nresp*n.bin.categories))
dim(M) <- c(nresp,n.bin.categories)
#
#  REPLACE "NA" VALUES WITH 999
#
var[is.na(var)] <- 999
#
#  REPLACE MATRIX
#
for (i in 1:nresp){
for (j in 1:n.bin.categories){
if (var[i] == j & var[i]!=999) M[i,j] <- 1
if (var[i] == (j+1) & var[i]!=999) M[i,j] <- 6
}}
#
M[is.na(M)] <- 9
return(M)
}
#
#
# WRITE  EASY ORDERED OPTIMAL CLASSIFICATION FUNCTION (binaryeven)
#
ooc.binaryeven <- function(votemat, dims=2, minvotes=10, lop=0.001, polarity=c(1,1)){
 
	#nresp <- nrow(votemat)
	nvotes <- ncol(votemat)
	ndim <- dims

	#   !!! "votemat" is the original matrix with ordinal choices !!!
	#   !!! "votemat.even.binary" is the matrix with binary choices recoded to be as evenly-balanced as possible !!!

	votemat.even.binary <- votemat
	for (j in 1:ncol(votemat)){
	votemat.even.binary[,j] <- binary.balanced(votemat[,j])
	}
	
	ncat <- NULL
	for (j in 1:nvotes){
	ncat[j] <- length(table(factor(votemat[,j], levels = min(votemat[,j], na.rm=TRUE):max(votemat[,j], na.rm=TRUE))))
	}
	sum.ncat <- sum(ncat)

	#  ***********************************************************************************************************************
	#  NOTE THAT ncol(votemat) is the number of input columns == number of issue scales being analyzed.  Each scale can have 
	#    a different number of categories -- 3 Point Scales, 7 Point Scales, etc.  Hence, ZVEC[,s] is a set of stacked normal 
	#    vectors that are the same for each issue.  If it is a 3 Point Scale there will be two normal vectors; if it is a 7 
	#    Point scale there will be six normal vectors; etc.  
	#    Hence the length of ZVEC[,] = SUM_j=1,q{ncat_j-1}
	#
	#  **********************************************************************************************************************

	# USE OPTIMAL CLASSIFICATION ON REGULAR MATRIX TO GET Nj

	hr.even.binary <- rollcall(votemat.even.binary, yea=1, nay=6, missing=9, notInLegis=NA)

	result.even.binary <- oc(hr.even.binary, dims=dims, minvotes=minvotes, lop=lop, polarity=polarity)

	oocObject <- list(respondents = result.even.binary$legislators,
		issues = result.even.binary$rollcalls,
		fits = result.even.binary$fits,
		votemat.even.binary = votemat.even.binary)

	
oocObject


}
#
# ===================================
#
# WRITE ORDERED OPTIMAL CLASSIFICATION FUNCTION
#
ooc <- function(votemat, dims=2, minvotes=10, lop=0.001, polarity=c(1,1),
	iter=25, ols=FALSE){
 
	set.seed(1985)

	#nresp <- nrow(votemat)
	nvotes <- ncol(votemat)
	ndim <- dims

	#   !!! "votemat" is the original matrix with ordinal choices !!!
	#   !!! "votemat.even.binary" is the matrix with binary choices recoded to be as evenly-balanced as possible !!!
	#   !!! "votemat.dominance.binary" is the matrix with binary choices organized in a dominance pattern !!!
	#   !!! "votemat.stacked.binary" is the matrix with the columns of "votemat.dominance.binary" stacked for each vote !!!

	votemat.even.binary <- votemat
	for (j in 1:ncol(votemat)){
	votemat.even.binary[,j] <- binary.balanced(votemat[,j])
	}

	votemat.dominance.binary <- binarize(votemat[,1])
	for (j in 2:nvotes){
	votemat.dominance.binary <- cbind(votemat.dominance.binary, binarize(votemat[,j]))
	}

	#votemat.stacked.binary <- stack(as.data.frame(binarize(votemat[,1])))$values
	#for (j in 2:nvotes){
	#votemat.stacked.binary <- cbind(votemat.stacked.binary, stack(as.data.frame(binarize(votemat[,j])))$values)
	#}
	#colnames(votemat.stacked.binary) <- NULL
	
	ncat <- NULL
	for (j in 1:nvotes){
	ncat[j] <- length(table(factor(votemat[,j], levels = min(votemat[,j], na.rm=TRUE):max(votemat[,j], na.rm=TRUE))))
	}
	sum.ncat <- sum(ncat)

	#  ***********************************************************************************************************************
	#  NOTE THAT ncol(votemat) is the number of input columns == number of issue scales being analyzed.  Each scale can have 
	#    a different number of categories -- 3 Point Scales, 7 Point Scales, etc.  Hence, ZVEC[,s] is a set of stacked normal 
	#    vectors that are the same for each issue.  If it is a 3 Point Scale there will be two normal vectors; if it is a 7 
	#    Point scale there will be six normal vectors; etc.  
	#    Hence the length of ZVEC[,] = SUM_j=1,q{ncat_j-1}
	#
	#  **********************************************************************************************************************

	# USE OPTIMAL CLASSIFICATION ON REGULAR MATRIX TO GET Nj

	#hr.even.binary <- rollcall(votemat.even.binary, yea=1, nay=6, missing=9, notInLegis=NA)
	hr.dominance.binary <- rollcall(votemat.dominance.binary, yea=1, nay=6, missing=9, notInLegis=NA)

	#result.even.binary <- oc(hr.even.binary, dims=ndim, minvotes=minvotes, lop=lop, polarity=polarity)
	result.dominance.binary <- oc(hr.dominance.binary, dims=ndim, minvotes=minvotes, lop=lop, polarity=polarity)

	#  =======================================================================================================================  

	# THROW OUT RESPONDENTS AND ISSUES THAT DON'T MEET THRESHOLD

	use.respondents <- !is.na(result.dominance.binary$legislators[,"coord1D"])
		
	votemat.dominance.binary <- votemat.dominance.binary[use.respondents,]

	votemat <- votemat[use.respondents,]

	nresp <- nrow(votemat.dominance.binary)

	#  =======================================================================================================================  

	if (ndim==1){

		respondents <- result.dominance.binary$legislators

		issues <- result.dominance.binary$rollcalls

		fits <- result.dominance.binary$fits

		minorityvote <- apply(cbind(issues[,1] + issues[,3], issues[,2] + issues[,4]),1,min)
		classerrors <- issues[,2] + issues[,3]
		correctclass <- issues[,1] + issues[,4]

		uniqueCLASSPERC <- sum(correctclass[1:(ncat[1]-1)]) / sum(correctclass[1:(ncat[1]-1)] + classerrors[1:(ncat[1]-1)])
		for (j in 2:nvotes){
		uniqueCLASSPERC[j] <- sum(correctclass[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])]) / sum(correctclass[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])] + classerrors[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])])
		}

		uniquePRE <- sum(minorityvote[1:(ncat[1]-1)] - classerrors[1:(ncat[1]-1)]) / sum(minorityvote[1:(ncat[1]-1)])
		for (j in 2:nvotes){
		uniquePRE[j] <- sum(minorityvote[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])] - classerrors[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])]) / sum(minorityvote[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])])
		}
		
		issues.unique <- cbind(uniqueCLASSPERC,uniquePRE)
		colnames(issues.unique) <- c("correctPerc","PRE")

		nvotes.dom.binary <- nrow(issues)
		oc1 <- na.omit(respondents[,7])
		ws <- issues[,7]

		onedim.classification.errors <- matrix(NA,nrow=nresp,ncol=nvotes.dom.binary)
			for (j in 1:nvotes.dom.binary){
			ivote <- votemat.dominance.binary[,j]
			polarity <- oc1 - ws[j]
			errors1 <- ivote==1 & polarity >= 0
			errors2 <- ivote==6 & polarity <= 0
			errors3 <- ivote==1 & polarity <= 0
			errors4 <- ivote==6 & polarity >= 0
			kerrors1 <- ifelse(is.na(errors1),9,errors1)
			kerrors2 <- ifelse(is.na(errors2),9,errors2)
			kerrors3 <- ifelse(is.na(errors3),9,errors3)
			kerrors4 <- ifelse(is.na(errors4),9,errors4)
			kerrors12 <- sum(kerrors1==1)+sum(kerrors2==1)
			kerrors34 <- sum(kerrors3==1)+sum(kerrors4==1)
			if (kerrors12 >= kerrors34){
			yeaerror <- errors3
			nayerror <- errors4
			}
			if (kerrors12 < kerrors34){
			yeaerror <- errors1
			nayerror <- errors2
			}
			junk <- yeaerror
			junk[nayerror==1] <- 1
			onedim.classification.errors[,j] <- junk
			}

		correct.class.binary <- matrix(NA,nrow=nresp,ncol=nvotes.dom.binary)
	
			for (i in 1:nresp){
			for (j in 1:nvotes.dom.binary){
				if (onedim.classification.errors[i,j] == 1) correct.class.binary[i,j] <- 0
				if (onedim.classification.errors[i,j] == 0) correct.class.binary[i,j] <- 1
				if (votemat.dominance.binary[i,j] == 9) correct.class.binary[i,j] <- NA
				}}

		correct.class.scale <- matrix(NA,nrow=nresp,ncol=nvotes)		

		for (i in 1:nresp){
			if (is.na(correct.class.binary[i,1])) correct.class.scale[i,1] <- NA	
			if (!is.na(correct.class.binary[i,1]) & all(correct.class.binary[i,1:(ncat[1]-1)]==TRUE)) correct.class.scale[i,1] <- 1
			if (!is.na(correct.class.binary[i,1]) & any(correct.class.binary[i,1:(ncat[1]-1)]==FALSE)) correct.class.scale[i,1] <- 0
			}		

			for (i in 1:nresp){
			for (j in 2:nvotes){
			if (is.na(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1)])) correct.class.scale[i,j] <- NA	
			if (!is.na(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1)]) & 
				all(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])]==TRUE)) correct.class.scale[i,j] <- 1
			if (!is.na(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1)]) & 
				any(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])]==FALSE)) correct.class.scale[i,j] <- 0
			}}
	
		correct.scale.respondents <- matrix(NA,nrow=nresp,ncol=2)
	
		colnames(correct.scale.respondents) <- c("wrongScale","correctScale")

		for (i in 1:nresp){
			correct.scale.respondents[i,] <- table(correct.class.scale[i,]) 
			}

		correct.scale.issues <- matrix(NA,nrow=nvotes,ncol=6)
	
		colnames(correct.scale.issues) <- c("wrongScale","correctScale","modalChoice","nChoices","errorsNull","PREScale")

			for (j in 1:nvotes){
			correct.scale.issues[j,1:2] <- table(correct.class.scale[,j]) 
			correct.scale.issues[j,3] <- max(table(votemat[,j]))
			}
		
			correct.scale.issues[,4] <- correct.scale.issues[,1] + correct.scale.issues[,2]
			correct.scale.issues[,5] <- correct.scale.issues[,4] - correct.scale.issues[,3]
			correct.scale.issues[,6] <- (correct.scale.issues[,5] - correct.scale.issues[,1]) / correct.scale.issues[,5]

		issues.unique <- cbind(issues.unique,correct.scale.issues)
		colnames(issues.unique) <- c("correctPerc","PRE","wrongScale","correctScale","modalChoice","nChoices","errorsNull","PREScale")

		tempResp <- cbind(correct.scale.respondents,(correct.scale.respondents[,2]/(correct.scale.respondents[,1] + correct.scale.respondents[,2])))

		respondents.2 <- matrix(NA, nrow=hr.dominance.binary$n, ncol=3)
		respondents.2[use.respondents,] <- tempResp

		respondents <- cbind(respondents,respondents.2)
		colnames(respondents) <- c("rank","correctYea","wrongYea","wrongNay","correctNay","volume","coord1D","wrongScale","correctScale","percent.correctScale")

		overall.correctClassScale <- sum(issues.unique[,4]) / sum(issues.unique[,6])
		overall.APREscale <- (sum(issues.unique[,7]) - sum(issues.unique[,3])) / (sum(issues.unique[,7]))

		fits <- c(fits,overall.correctClassScale,overall.APREscale)

		oocObject <- list(respondents = respondents,
			issues = issues,
			issues.unique = issues.unique,
			fits = fits)

		}


	if (ndim>=2){

	cat("\n\tRunning Ordered Optimal Classification...\n\n")
	start <- proc.time()
	
	# XMAT = starting values for the respondent coordinates

	XMAT <- result.dominance.binary$legislators[use.respondents,grep("coord", colnames(result.dominance.binary$legislators))]
	XMAT[is.na(XMAT)] <- runif(sum(is.na(XMAT)), -1, 1)

	# ZVEC = starting values for the normal vectors

		# Random Starts for the Normal Vectors
		Nj <- hypersphere.sample(ndim, nvotes)
		for (j in 1:nvotes){
		if (Nj[j,1] < 0) Nj[j,] <- -1 * Nj[j,]
		}

		ZVEC <- matrix(rep(Nj[1,], ncat[1]-1), ncol=ndim, byrow=TRUE)
		for (j in 2:nvotes){
		ZVEC <- rbind(ZVEC, matrix(rep(Nj[j,], (ncat[j]-1)), ncol=ndim, byrow=TRUE))
		}

	# LDATA = roll call data

	LDATA <- votemat.dominance.binary
	nummembers <- nresp
	np <- nresp
	npzz <- np
	numvotes <- ncol(votemat.dominance.binary)
	nrcall <- ncol(votemat.dominance.binary)
	nv <- nrcall # nrcall = SUM_j=1,q{ncat_j-1}
	ns <- ndim
	ndual <- ndim * (nummembers + numvotes) + 111
	# ndual <- 2 * (nummembers + numvotes) + 111
	LLEGERR <- rep(0,np*ns) 
	dim(LLEGERR) <- c(np,ns)
	LERROR <- rep(0,np*nrcall)
	dim(LERROR) <- c(np,nrcall)
	MCUTS <- rep(0,nrcall*2)
	dim(MCUTS) <- c(nrcall,2)

	X3 <- rep(0,np*ns)
	dim(X3) <- c(np,ns)

	WS <- rep(0,ndual)
	dim(WS) <- c(ndual,1)
	LLV <- rep(0,ndual)
	dim(LLV) <- c(ndual,1)
	LLVB <- rep(0,ndual)
	dim(LLVB) <- c(ndual,1)
	LLE <- rep(0,ndual)
	dim(LLE) <- c(ndual,1)
	LLEB <- rep(0,ndual)
	dim(LLEB) <- c(ndual,1)
	ZS <- rep(0,ndual)
	dim(ZS) <- c(ndual,1)
	XJCH <- rep(0,25)
	dim(XJCH) <- c(25,1)
	XJEH <- rep(0,25)
	dim(XJEH) <- c(25,1)
	XJCL <- rep(0,25)
	dim(XJCL) <- c(25,1)
	XJEL <- rep(0,25)
	dim(XJEL) <- c(25,1)
	errorcount <- rep(0,iter*6)
	dim(errorcount) <- c(iter,6)
	errorcounts <- rep(0,iter*6)
	dim(errorcounts) <- c(iter,6)

	YSS <- NULL
	MM <- NULL
	XXX <- NULL
	MVOTE <- NULL
	WSSAVE <- NULL

	#  *****************************************************************
	#
	#  GLOBAL LOOP
	#
	#  *****************************************************************

	for (iiii in 1:iter){

	#  BEGIN JANICE LOOP
	kerrors <- 0
	kchoices <- 0
	kcut <- 1
	lcut <- 6
	for (jx in 1:nrcall){
	for (i in 1:np){
	sum=0.0
	for (k in 1:ns){
	sum <- sum+XMAT[i,k]*ZVEC[jx,k]
	}

	#  SAVE PROJECTION VECTORS -- RESPONDENT BY ROLL CALL MATRIX

	YSS[i]=sum
	XXX[i]=sum
	MM[i]=LDATA[i,jx]
	if(LDATA[i,jx]==0)MM[i]=9
	}

	#  END OF RESPONDENT PROJECTION LOOP

	#  SORT PROJECTION VECTOR (Y-HAT)

	LLL <- order(XXX)
	YSS <- sort(XXX)
	MVOTE <- MM[LLL]


	#  FIRST FIND MULTIPLE CUTTING POINTS CORRESPONDING TO EACH NORMAL VECTORS -- THIS STEMS FROM THE STACKING
	#  BY ISSUE EXPLAINED ABOVE
	
	ivot <- jx
	jch <- 0
	jeh <- 0
	jcl <- 0
	jel <- 0
	irotc <- 0
	rescutpoint <- find_cutpoint(nummembers,numvotes,
                        npzz,nv,np,nrcall,ns,ndual,ivot,XMAT,YSS,MVOTE,WS, 
                        LLV,LLVB,LLE,LLEB, 
                        LERROR,ZS,jch,jeh,jcl,jel,irotc,kcut,lcut, 
                        LLL,XJCH,XJEH,XJCL,XJEL)
		
	#  STORE DIRECTIONALITY OF ROLL CALL

	WSSAVE[jx] <- rescutpoint[[13]][jx]
	jch <- rescutpoint[[20]]
	jeh <- rescutpoint[[21]]
	jcl <- rescutpoint[[22]]
	jel <- rescutpoint[[23]]
	kerrors <- kerrors+jeh+jel
	kchoices <- kchoices+jch+jeh+jcl+jel
	kcut <- rescutpoint[[25]]
	lcut <- rescutpoint[[26]]
	MCUTS[jx,1] <- kcut
	MCUTS[jx,2] <- lcut
	}

	# END OF ROLL CALL LOOP FOR JANICE

	errorcount[iiii,1] <- kerrors
	WS <- WSSAVE
	jjj <- 1
	ltotal <- 0
	mwrong <- 0
	iprint <- 0
	respoly <- find_polytope(nummembers,numvotes,
			jjj,np,nrcall,ns,ndual,XMAT,LLEGERR, 
                        ZVEC,WS,MCUTS,LERROR,ltotal,mwrong, 
                        LDATA,iprint) 

	#  DATA STORED FORTRAN STYLE -- STACKED BY COLUMNS

	errorcount[iiii,2] <- respoly[[15]]
	dim(respoly[[8]]) <- c(np,ns)
	X3 <- respoly[[8]]
	r1 <- cor(X3[,1],XMAT[,1])
	r2 <- cor(X3[,2],XMAT[,2])
	errorcount[iiii,3] <- r1
	errorcount[iiii,4] <- r2
	XMAT <- X3

	cat("\t\tGetting respondent coordinates...\n")

	#  ROLL CALL LOOP (CALCULATE NORMAL VECTORS VIA OLS OR KRLS)

	XMATTEMP <- XMAT
	ZVEC_regress <- rep(0,ndim*nvotes)
	dim(ZVEC_regress) <- c(nvotes,ndim)

	for (j in 1:nvotes){

		choices <- votemat[,j]

			if (ols==FALSE){
			temp.complete <- complete.cases(XMATTEMP & choices) 
			X <- XMATTEMP[temp.complete,]
			y <- choices[temp.complete]
			select <- sample(1:length(y), 300, replace=TRUE)
			res <- krls(X=X[select,], y=y[select], print.level=0)
			denominator <- sqrt(sum(res$avgderivatives^2))
			for (q in 1:ndim){
			ZVEC_regress[j,q] <- res$avgderivatives[q] / denominator
			}}

			if (ols==TRUE){
			res <- lm(choices ~ XMATTEMP)
			denominator <- sqrt(sum(res$coef[-1]^2))
			for (q in 1:ndim){
			ZVEC_regress[j,q] <- res$coef[(q+1)] / denominator
			}}

	}

	ZVECpost <- matrix(rep(ZVEC_regress[1,], ncat[1]-1), ncol=ndim, byrow=TRUE)
	for (j in 2:nvotes){
	ZVECpost <- rbind(ZVECpost, matrix(rep(ZVEC_regress[j,], (ncat[j]-1)), ncol=ndim, byrow=TRUE))
	}
	for (j in 1:nrcall){
	if (ZVECpost[j,1] < 0) ZVECpost[j,] <- -1 * ZVECpost[j,]
	}
	r1 <- cor(ZVEC[,1],ZVECpost[,1])
	r2 <- cor(ZVEC[,2],ZVECpost[,2])
	errorcount[iiii,5] <- r1
	errorcount[iiii,6] <- r2

	ZVEC <- matrix(rep(ZVEC_regress[1,], ncat[1]-1), ncol=ndim, byrow=TRUE)
	for (j in 2:nvotes){
	ZVEC <- rbind(ZVEC, matrix(rep(ZVEC_regress[j,], (ncat[j]-1)), ncol=ndim, byrow=TRUE))
	}

	for (j in 1:nrcall){
	if (ZVEC[j,1] < 0) ZVEC[j,] <- -1 * ZVEC[j,]
	}

	cat("\t\tCalculating normal vectors...\n")

	}
	#  END OF GLOBAL LOOP

	#  RUN POLYTOPE SEARCH ONE FINAL TIME	
	respoly <- find_polytope(nummembers,numvotes,
                         jjj,np,nrcall,ns,ndual,XMAT,LLEGERR, 
                        ZVEC,WS,MCUTS,LERROR,ltotal,mwrong, 
                        LDATA,iprint) 

	dim(respoly[[8]]) <- c(np,ns)
	X3 <- respoly[[8]]
	XMAT <- X3
	totalerrors <- respoly[[15]]

	polarity <- matrix(NA,nrow=nresp,ncol=sum(ncat)-nvotes)
	correctYea <- NA
	wrongYea <- NA
	wrongNay <- NA
	correctNay <- NA
	kchoices <- NA
	kminority <- NA
	voteerrors <- NA
	PRE <- NA
	respondent.class.correctYea <- matrix(0, nrow=nresp, ncol=(sum(ncat)-nvotes))
	respondent.class.wrongYea <- matrix(0, nrow=nresp, ncol=(sum(ncat)-nvotes))
	respondent.class.wrongNay <- matrix(0, nrow=nresp, ncol=(sum(ncat)-nvotes))
	respondent.class.correctNay <- matrix(0, nrow=nresp, ncol=(sum(ncat)-nvotes))

	for (j in 1:(sum(ncat)-nvotes)){
	ivote <- votemat.dominance.binary[,j]
	#polarity <- XMAT[,1]*ZVEC[j,1] + XMAT[,2]*ZVEC[j,2] - WS[j]
	polarity <- as.vector(XMAT%*%ZVEC[j,]) - WS[j]
	errors1 <- ivote==1 & polarity >= 0
	errors2 <- ivote==6 & polarity <= 0
	errors3 <- ivote==1 & polarity <= 0
	errors4 <- ivote==6 & polarity >= 0
	kerrors1 <- ifelse(is.na(errors1),9,errors1)
	kerrors2 <- ifelse(is.na(errors2),9,errors2)
	kerrors3 <- ifelse(is.na(errors3),9,errors3)
	kerrors4 <- ifelse(is.na(errors4),9,errors4)
	kerrors12 <- kerrors34 <- 0
	kerrors12 <- sum(kerrors1==1)+sum(kerrors2==1)
	kerrors34 <- sum(kerrors3==1)+sum(kerrors4==1)
	if (kerrors12 >= kerrors34){
	yeaerror <- errors3
	nayerror <- errors4
	}
	if (kerrors12 < kerrors34){
	yeaerror <- errors1
	nayerror <- errors2
	}
	kerrorsmin <- min(kerrors12,kerrors34)

	respondent.class.correctYea[,j] <- as.numeric(ivote==1 & !yeaerror & !nayerror)
	respondent.class.wrongYea[,j] <- as.numeric(ivote==6 & nayerror)
	respondent.class.wrongNay[,j] <- as.numeric(ivote==1 & yeaerror)
	respondent.class.correctNay[,j] <- as.numeric(ivote==6 & !yeaerror & !nayerror)

	kpyea <- sum(ivote==1)
	kpnay <- sum(ivote==6)
	kpyeaerror <- sum(yeaerror)
	kpnayerror <- sum(nayerror)
	correctYea[j] <- kpyea - kpyeaerror
	wrongYea[j] <- kpnayerror
	wrongNay[j] <- kpyeaerror
	correctNay[j] <- kpnay - kpnayerror
	kchoices[j] <- kpyea+kpnay
	kminority[j] <- min(kpyea,kpnay)
	voteerrors[j] <- kerrorsmin
	PRE[j] <- (min(kpyea,kpnay) - kerrorsmin)/(min(kpyea,kpnay))
	}

	totalerrors2 <- sum(voteerrors)
	totalcorrect2 <- sum(correctYea) + sum(correctNay)
	APRE <- ((sum(kminority) - sum(voteerrors)) / sum(kminority))
	totalCorrectClass <- totalcorrect2 / (totalcorrect2 + totalerrors2)

	normVectorkD <- ZVEC
	midpoints <- WS

	issues <- cbind(correctYea,wrongYea,wrongNay,correctNay,PRE,normVectorkD,midpoints)
	colnames(issues) <- c("correctYea","wrongYea","wrongNay","correctNay","PRE",paste("normVector",1:ndim,"D",sep=""),"midpoints")

	respondent.correctYea <- rowSums(respondent.class.correctYea)
	respondent.wrongYea <- rowSums(respondent.class.wrongYea)
	respondent.wrongNay <- rowSums(respondent.class.wrongNay)
	respondent.correctNay <- rowSums(respondent.class.correctNay)

	correct.class.binary <- matrix(NA,nrow=nresp,ncol=ncol(respondent.class.correctYea))
	
		for (i in 1:nresp){
		for (j in 1:ncol(respondent.class.correctYea)){

			if (respondent.class.correctYea[i,j] == 1) correct.class.binary[i,j] <- 1
			if (respondent.class.wrongYea[i,j] == 1) correct.class.binary[i,j] <- 0
			if (respondent.class.wrongNay[i,j] == 1) correct.class.binary[i,j] <- 0
			if (respondent.class.correctNay[i,j] == 1) correct.class.binary[i,j] <- 1

			}}

	correct.class.scale <- matrix(NA,nrow=nresp,ncol=nvotes)		

		for (i in 1:nresp){
			if (is.na(correct.class.binary[i,1])) correct.class.scale[i,1] <- NA	
			if (!is.na(correct.class.binary[i,1]) & all(correct.class.binary[i,1:(ncat[1]-1)]==TRUE)) correct.class.scale[i,1] <- 1
			if (!is.na(correct.class.binary[i,1]) & any(correct.class.binary[i,1:(ncat[1]-1)]==FALSE)) correct.class.scale[i,1] <- 0
			}		

			for (i in 1:nresp){
			for (j in 2:nvotes){
			if (is.na(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1)])) correct.class.scale[i,j] <- NA	
			if (!is.na(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1)]) & 
				all(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])]==TRUE)) correct.class.scale[i,j] <- 1
			if (!is.na(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1)]) & 
				any(correct.class.binary[i,(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])]==FALSE)) correct.class.scale[i,j] <- 0
			}}
	
	correct.scale.respondents <- matrix(NA,nrow=nresp,ncol=2)
	
	colnames(correct.scale.respondents) <- c("wrongScale","correctScale")

		for (i in 1:nresp){
			correct.scale.respondents[i,] <- table(correct.class.scale[i,]) 
		}

	correct.scale.issues <- matrix(NA,nrow=nvotes,ncol=6)
	
	colnames(correct.scale.issues) <- c("wrongScale","correctScale","modalChoice","nChoices","errorsNull","PREScale")

		for (j in 1:nvotes){
			correct.scale.issues[j,1:2] <- table(correct.class.scale[,j]) 
			correct.scale.issues[j,3] <- max(table(votemat[,j]))
		}
		
			correct.scale.issues[,4] <- correct.scale.issues[,1] + correct.scale.issues[,2]
			correct.scale.issues[,5] <- correct.scale.issues[,4] - correct.scale.issues[,3]
			correct.scale.issues[,6] <- (correct.scale.issues[,5] - correct.scale.issues[,1]) / correct.scale.issues[,5]

	tempResp <- cbind(respondent.correctYea,respondent.wrongYea,respondent.wrongNay,respondent.correctNay,XMAT,
		correct.scale.respondents,(correct.scale.respondents[,2]/(correct.scale.respondents[,1] + correct.scale.respondents[,2])))

	respondents <- matrix(NA, nrow=hr.dominance.binary$n, ncol=(7+ndim))
	respondents[use.respondents,] <- tempResp
	
	colnames(respondents) <- c("correctYea","wrongYea","wrongNay","correctNay",paste("coord",1:ndim,"D",sep=""),"wrongScale","correctScale","percent.correctScale")
	
	# COMPRESS ROLL CALL MATRIX
	
	unique.normVectorkD <- unique(normVectorkD)
	unique.NV1 <- unique.normVectorkD[,1]
	unique.NV2 <- unique.normVectorkD[,2]
	minorityvote <- apply(cbind(issues[,1] + issues[,3], issues[,2] + issues[,4]),1,min)
	classerrors <- issues[,2] + issues[,3]
	correctclass <- issues[,1] + issues[,4]

	uniqueCLASSPERC <- sum(correctclass[1:(ncat[1]-1)]) / sum(correctclass[1:(ncat[1]-1)] + classerrors[1:(ncat[1]-1)])
	for (j in 2:nvotes){
	uniqueCLASSPERC[j] <- sum(correctclass[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])]) / sum(correctclass[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])] + classerrors[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])])
	}

	uniquePRE <- sum(minorityvote[1:(ncat[1]-1)] - classerrors[1:(ncat[1]-1)]) / sum(minorityvote[1:(ncat[1]-1)])
	for (j in 2:nvotes){
	uniquePRE[j] <- sum(minorityvote[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])] - classerrors[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])]) / sum(minorityvote[(cumsum(ncat-1)[j-1]+1):(cumsum(ncat-1)[j])])
	}
	
	b <- c(1.0,0)
	theta <- NULL
	for (j in 1:nvotes){
	a <- c(unique.NV1[j], unique.NV2[j])
	theta[j] <- acos_d(sum(a*b) / (sqrt(sum(a*a))*sqrt(sum(b*b))))
	if (a[2] < 0) theta[j] <- -1 * theta[j]
	}
	
	issues.unique <- cbind(uniqueCLASSPERC,uniquePRE,unique.normVectorkD,theta,correct.scale.issues)
	colnames(issues.unique) <- c("correctPerc","PRE",paste("normVector",1:ndim,"D",sep=""),"normVectorAngle2D",
		"wrongScale","correctScale","modalChoice","nChoices","errorsNull","PREScale")

	overall.correctClassScale <- sum(issues.unique[,"correctScale"]) / sum(issues.unique[,"nChoices"])
	overall.APREscale <- (sum(issues.unique[,"errorsNull"]) - sum(issues.unique[,"wrongScale"])) / (sum(issues.unique[,"errorsNull"]))

	fits <- c(totalCorrectClass,APRE,overall.correctClassScale,overall.APREscale)

	# <errorcount>
	#   First column = errors after cutpoint search
	#   Second column = errors after polytope search
	#   Third column = Correlation respondent coordinates 1st dim
	#   Fourth column = Correlation respondent coordinates 2nd dim
	#   Fifth column = Correlation normal vectors 1st dim
	#   Sixth column = Correlation normal vectors 2nd dim	

	oocObject <- list(respondents = respondents,
		issues = issues,
		issues.unique = issues.unique,
		fits = fits,
		oc.result.dominance.binary = result.dominance.binary,
		errorcount = errorcount,
		votemat.even.binary = votemat.even.binary)
		
    cat("\n\nOrdered Optimal Classification estimation completed successfully.")
    cat("\nOrdered Optimal Classification took", (proc.time() - start)[3], "seconds to execute.\n\n")

	}

	
oocObject


}