
########################
## Load  data ##
########################
seedID = 2 # seedID=1, initialize z to the truth; seedID=2, initialize one Ct per cluster
		   # 3 = initialize all Ct to one cluster
initialTruePara = FALSE
fixX_toTrue = FALSE
fixZ = FALSE
noBeta0inModel = FALSE
PrintBeta0 = FALSE
nP = 1
totalIter = 1200

set.seed(seedID)

library(mvtnorm); library(pscl); library(abind); library( plyr); library(MCMCpack)
library(doSNOW); library(foreach); library(abind)


workDir <- 'C:/Users/shirleyr/Documents/ResearchZ/'
outPath <- 'Output/Simulation/K_4/'
setwd(workDir)

###################
## load Simulated data ##
###################

load(file="Data/simuData_R10pct_K4_x0Set0_a99_intercept.RData")

covNames = c("intercept","logFinishedSquareFeet","logLotSizeSquareFeet","logBathroomCnt")
if (noBeta0inModel) {covNames = c("logFinishedSquareFeet","logLotSizeSquareFeet","logBathroomCnt")} 

ord = order(simuData$CensusTractRegionID, simuData$MonthTransaction, simuData$TransactionID)
simuData = simuData[ord,]

## CtTable ##
obs  = simuData
#Remove the Ct "306404" that has only 2 obs
#obs = obs[!obs$CensusTractRegionID == 306404,]
# 141 census tract in Seattle City
CtTable = obs[!duplicated(obs$CensusTractRegionID), 
			c("CensusTractRegionID" ,"CtRegionName")]
ord = order(CtTable$CensusTractRegionID)
CtTable = CtTable[ord,]


############################
## Gibbs Sampler Assembly ##
############################

set.seed(2014) # make sure diff chains have the same training and test split
dataAll = simuData
ord = order(dataAll$CensusTractRegionID, dataAll$MonthTransaction, dataAll$TransactionID)
dataAll = dataAll[ord,]
idxTest = NULL
nCt = nrow(CtTable)
for (i in 1:nCt){
	CtId = CtTable$CensusTractRegionID[i]
	idx.tract = which(dataAll$CensusTractRegionID==CtId)
	nObs.tract = length(idx.tract)
	nTest = round(nObs.tract*0.1)
	if (nTest>0){
		idxTest.tract = sort(sample(idx.tract, nTest, replace=FALSE))
		idxTest = c(idxTest, idxTest.tract)
	}
}	
idxTest = sort(idxTest)
dataTest  = dataAll[ idxTest,]
dataTrain = dataAll[-idxTest,]

set.seed(seedID)


##################
## Training Set ##
##################

obs = dataTrain
#obs = dataAll

scaleFactor = 200
obs$y = scaleFactor*obs$y
dataTest$y = scaleFactor*dataTest$y
obs$intercept = 1
dataTest$intercept = 1

tmLab = sort(unique(obs$MonthTransaction))
CtIds = sort(unique(obs$CensusTractRegionID))
n = length(CtIds)

if (FALSE){
	yMeanVec = NULL
	for (i in 1:n){
		idct = obs$CensusTractRegionID==CtIds[i]
		tem = mean(obs$y[idct])
		yMeanVec = c(yMeanVec, tem)
		obs$y[idct] = obs$y[idct] - mean(obs$y[idct])
	}
	names(yMeanVec) = CtIds
	yMeanVec = yMeanVec/scaleFactor
}

## Cached: count by month and tract
Lti = aggregate( y ~MonthTransaction+ CensusTractRegionID ,data=obs, FUN=length )
colnames(Lti)[3] = "count"
L = expand.grid(tmLab,CtIds)
colnames(L) = colnames(Lti)[1:2]
L = merge(L, Lti, by=colnames(L),all.x=TRUE, all.y=FALSE)
L$count[is.na(L$count)]=0
ord = order(L$CensusTractRegionID, L$MonthTransaction)
L = L[ord,]
T = length(tmLab)

## for test set
Lti = aggregate( y ~MonthTransaction+ CensusTractRegionID ,data=dataTest, FUN=length )
colnames(Lti)[3] = "count"
LTest = expand.grid(tmLab,CtIds)
colnames(LTest) = colnames(Lti)[1:2]
LTest = merge(LTest, Lti, by=colnames(LTest),all.x=TRUE, all.y=FALSE)
LTest$count[is.na(LTest$count)]=0
ord = order(LTest$CensusTractRegionID, LTest$MonthTransaction)
LTest = LTest[ord,]

## Cached: y by month
yMonth = list()
for (tt in 1:length(tmLab)){
	yMonth[[tt]] = obs[obs$MonthTransaction==tmLab[tt],
		c("CensusTractRegionID","MonthTransaction","y", covNames)]
}
names(yMonth)=tmLab

## Cached: y by tracts
yTract = list()
for (i in 1:length(CtIds)){
	yTract[[i]] = obs[obs$CensusTractRegionID==CtIds[i],
		c("CensusTractRegionID","MonthTransaction","y", covNames)]
}
names(yTract)=CtIds

## Cached [Test data]: y by month
yMonthTest = list()
for (tt in 1:length(tmLab)){
	yMonthTest[[tt]] = dataTest[dataTest$MonthTransaction==tmLab[tt],
		c("CensusTractRegionID","MonthTransaction","y", covNames)]
}
names(yMonthTest)=tmLab

## Cached [Test data]: y by tracts
yTractTest = list()
for (i in 1:length(CtIds)){
	yTractTest[[i]] = dataTest[dataTest$CensusTractRegionID==CtIds[i],
		c("CensusTractRegionID","MonthTransaction","y", covNames)]
}
names(yTractTest)=CtIds

############

CtTable = CtTable[CtTable$CensusTractRegionID %in% CtIds,]

# 0. Initialization cluster assignment

#groups = rep(1:5, each=4)
#groups = rep(1:7, each=3)[-21]
if (seedID==1) groups = z
if (seedID==2) groups = 1:length(CtIds)
if (seedID==3) groups = rep(1, length(CtIds))

# 0. Initialize parameters

if (TRUE){
	zVec = groups
	names(zVec) = CtIds
	alpha = 1
	K = length(unique(zVec))
	n = length(CtIds)
	
	if (initialTruePara){ # Initialize by the true para
		aVec = a
		names(aVec) = CtIds
		
		LambdaMat = matrix(0, n, K)
		for (k in 1:K){
			LambdaMat[,k] = lambda*scaleFactor
		}
		rownames(LambdaMat) = CtIds
		colnames(LambdaMat) = 1:K
		sigma0 = sigma0*scaleFactor
		RVec = Ri*(scaleFactor^2)
		names(RVec) = CtIds
		if (!noBeta0inModel){
			BetaMat = t(Beta)*scaleFactor
		} else {
			BetaMat = t(Beta[-1,])*scaleFactor
		}
		
		rownames(BetaMat) = CtIds
		colnames(BetaMat) = covNames
		
		muLambda = mean(lambda)*scaleFactor
		sigmaLambda = sd(lambda)*scaleFactor
		
		muA = mean(a)
		sigmaA = sd(a)
		#muBetaVec = mu_beta*scaleFactor
		muBetaVec = apply(BetaMat, MARGIN=2, FUN=mean)
		names(muBetaVec) = covNames
		#sigmaBetaVec = sigma_beta*scaleFactor
		sigmaBetaVec = rep(0.1*scaleFactor, length(covNames))
		names(sigmaBetaVec) = covNames
		if (fixX_toTrue) {
			xMat = t(x)*scaleFactor
			rownames(xMat) = CtIds
			etaMat = eta
			colnames(eta) = 1:4
			x0Vec = rep(0, length(CtIds))
			names(x0Vec) = CtIds
		}
		
	} else {
		aVec = rnorm(n, 0.95, 0.02)
		aVec[aVec>1] = 0.998
		#aVec =a
		names(aVec) = CtIds
		LambdaMat = matrix(0.015, n, K)*scaleFactor
		rownames(LambdaMat) = CtIds
		colnames(LambdaMat) = 1:K
		sigma0 = 0.005*scaleFactor
		RVec = rep(0.1,n)*(scaleFactor^2)
		names(RVec) = CtIds
		#BetaMat = matrix(0*scaleFactor, n, length(covNames) )
		# initialize beta with truth
		#BetaMat = t(Beta)*scaleFactor 
		# initialize beta with mean(y) per tract
		BetaMat = matrix(0*scaleFactor, n, length(covNames) )
		if ("intercept"%in% covNames) BetaMat[,1] = unlist(lapply(yTract, FUN=function(x){mean(x$y)}))
		rownames(BetaMat) = CtIds
		colnames(BetaMat) = covNames
		muLambda = 0.015*scaleFactor
		sigmaLambda =0.002*scaleFactor
		muA = 0.9
		#sigmaA = 0.005  
		#sigmaA = 0.05
		#sigmaA=0.2
		sigmaA = 0.5
		muBetaVec = apply(BetaMat, MARGIN=2, FUN=mean)
		names(muBetaVec) = covNames
		sigmaBetaVec = rep(0.1*scaleFactor, length(covNames))
		names(sigmaBetaVec) = covNames
	}
	# initial hyper prior parameters
		alpha.alpha0 = 2
		beta.alpha0 =4
		alpha.R0 = 3
		beta.R0 =1
		alpha.epsilon0 = 0.5
		#beta.epsilon0 = 1
		beta.epsilon0 = 0.5
		muLambda0 = 0.015*scaleFactor
		#sigmaLambda0 = 0.002*scaleFactor
		sigmaLambda0 = 0.01*scaleFactor
		#alphaLambda0 =3  # sigmaLambda^2 follows IG(alphaLambda0, betaLambda0)
		alphaLambda0 = 0.1		
		betaLambda0 = 0.1
		muA0 = 0.99       # muA follow N(muA0, sigmaA0)
		sigmaA0 = 0.1
		#alphaA0 = 3       # sigmaA^2 follow IG(alphaA0, betaA0)
		#alphaA0 = 0.5
		alphaA0 = 0.2
		betaA0 = 1/(scaleFactor^2)
		#alphaA0 = betaA0 = 0.00001
		muBetaVec0 = rep(0, length(covNames))   # muBetaVec_r follows N(muBetaVec0_r, sigmaBetaVec0_r)
		sigmaBetaVec0 = rep(1, length(covNames))*scaleFactor
		alphaBetaVec0 = rep(0.1, length(covNames))# sigmaBetaVec_r follows IG(alphaBetaVec0_r, betaBetaVec0_r)
		betaBetaVec0 = rep(0.1, length(covNames))
		mux0 = 0
		Vx0 = 0.1*(scaleFactor^2)
	
}

## plot the prior
if (FALSE){
	# e.g. for prior distribution on sigma_a, truth = 0.005
	temx = seq(0.001, 0.1, 0.001)
	#temx = seq(0.001, 1, 0.001)
	temy =  pigamma(temx^2, alphaA0,betaA0 )
	plot(temx, temy)
	head(cbind(temx, round(temy,6)))
	# e.g. for prior distribution on sigma_0, truth = 0.005
	temx = seq(0.001, 0.1, 0.001)
	temy =  pigamma((temx*scaleFactor)^2, alpha.epsilon0 , beta.epsilon0 )
	plot(temx, temy)
	head(cbind(temx, round(temy,6)))
	# e.g. for prior distribution on sigma_lambda, truth = 0.002
	temx = seq(0.001, 0.1, 0.001)
	temy =  pigamma((temx*scaleFactor)^2, alphaLambda0, betaLambda0 )
	plot(temx, temy)
	head(cbind(temx, round(temy,6)))
}


# Book keeping
zVecIter = zVec
KIter = length(unique(zVec))
alphaIter = alpha
aVecIter = aVec
RVecIter = RVec
sigma0Iter = sigma0
muLambdaIter = muLambda
sigmaLambdaIter = sigmaLambda
muAIter = muA
sigmaAIter = sigmaA
muBetaVecIter = muBetaVec
sigmaBetaVecIter = sigmaBetaVec
BetaMatIter = BetaMat
xMatIter = NULL 
x0VecIter = NULL
LambdaMatIter = LambdaMat[cbind(CtIds, as.character(zVec))] 
etaMatIter = NULL
logPosteriorIter = NULL
logPosteriorTestIter = NULL
logPosteriorIter.integrateEtaX = NULL
logPosteriorTestIter.integrateEtaX = NULL

# Sample the initial processor indicator for each cluster, piVec

kLabels = unique(zVec)
K = length(kLabels)
piVec = sample(1:nP, size=K, replace=TRUE)
names(piVec) = kLabels

# piVec maps which cluster label assigned to which processor id
# zVec maps which Ct id assigned to which cluster label

# 0. [Global] Initialize cached covariate effect
Lti = aggregate( y ~MonthTransaction+ CensusTractRegionID ,data=obs, FUN=length )
colnames(Lti)[3] = "count"
L = expand.grid(tmLab,CtIds)
colnames(L) = colnames(Lti)[1:2]
L = merge(L, Lti, by=colnames(L),all.x=TRUE, all.y=FALSE)
L$count[is.na(L$count)]=0
ord = order(L$CensusTractRegionID, L$MonthTransaction)
L = L[ord,]

source("Code/GibbsSamplerFun.R")
 
outCovEffMonth = calCovEffMonth(yMonth, BetaMat, covNames, T)
covEffMonth = outCovEffMonth$covEffMonth
yDehedonicMonthMean =  outCovEffMonth$yDehedonicMonthMean
	
outCovEffMonthTest <- calCovEffMonth(yMonthTest, BetaMat, covNames, T)
covEffMonthTest = outCovEffMonthTest$covEffMonth
yDehedonicMonthMeanTest =  outCovEffMonthTest$yDehedonicMonthMean
		
# 0. Updating para for z initial at iter = 0
iter = 0
if (TRUE){ 
for (subIter in 1:5){
	alpha = sample.alpha(alpha, n, alpha.alpha0, beta.alpha0, K)
	alphaIter = rbind(alphaIter, alpha)

	# 2. Sample the tract process x, update x. dim(xMat) = n*T
	outX = sample.x(zVec, LambdaMat, aVec, RVec, sigma0 ,L, yMonth, covEffMonth, 
					 tmLab, CtIds, covNames, mux0, Vx0)
	xMat = outX$xMat
	x0Vec = outX$x0Vec
	
	
	# 3. Sample the cluster latent factor eta_t,k. dim(etaMat) is T*KMax
	etaMat = sample.eta(zVec, LambdaMat, aVec, sigma0, xMat, tmLab, x0Vec, K)
		
	# 4. Sample the loading matrix lambda_ik. dim(LambdaMat) = n*K
	LambdaMat = sample.LambdaMat(zVec, muLambda, sigmaLambda, sigma0, etaMat, xMat, aVec, tmLab, CtIds, x0Vec, K)

	# 5. Sample AR coefficient a_i. dim(aVec) = n
	aVec <- sample.aVec(zVec, LambdaMat, xMat,sigma0, etaMat, muA, sigmaA, CtIds, tmLab, x0Vec)

	# 6. Sample the response variance R_i. length(RVec) = n
	RVec <- sample.RVec(CtIds, xMat, alpha.R0, beta.R0, yTract, tmLab, covEffMonth )

	# 7. Sample the covariate effect Beta_ir. dim(BetaMat)= n*length(covNames)
	if (length(CtIds)==1) {
		BetaMat = matrix(BetaMat,1,length(covNames))
		rownames(BetaMat) = CtIds
		colnames(BetaMat) = covNames
	}
	BetaMat <- sample.BetaMat(muBetaVec, sigmaBetaVec, RVec, xMat, BetaMat, yTract, CtIds, tmLab,covNames)

	# 9. [Global] Sample x process error, a scalar
	sigma0 <- sample.sigma0(xMat, etaMat, aVec, LambdaMat, zVec, alpha.epsilon0, beta.epsilon0, CtIds, tmLab, x0Vec,K)
	sigma0Iter= c(sigma0Iter, sigma0)

	# 10.[Global] Sample muLambda, a scalar
	muLambda <- sample.muLambda(muLambda0, sigmaLambda0, LambdaMat, zVec, sigmaLambda, K)
	muLambdaIter = c(muLambdaIter, muLambda)
	
	# 11.[Global] Sample sigmaLambda, a scalar
	sigmaLambda <- sample.sigmaLambda(zVec, LambdaMat, muLambda, alphaLambda0, betaLambda0, K)
	sigmaLambdaIter = c(sigmaLambdaIter, sigmaLambda)

	# 12.[Global] Sample muA, a scalar
	muA <- sample.muA(aVec, sigmaA, muA0, sigmaA0)
	muAIter = c(muAIter, muA)

	# 13.[Global] Sample sigmaA, a scalar
	sigmaA <- sample.sigmaA(aVec, muA, alphaA0, betaA0)
	sigmaAIter = c(sigmaAIter, sigmaA)

	# 14.[Global] Sample the mean covariate effect, muBetaVec. length(muBetaVec) = length(covNames)
	muBetaVec <- sample.muBetaVec(BetaMat, sigmaBetaVec, muBetaVec0, sigmaBetaVec0, covNames)
	muBetaVecIter = rbind(muBetaVecIter, muBetaVec)

	# 15.[Global] Sample the variability of covariate effect, sigmaBetaVec, length(sigmaBetaVec) = length(covNames)
	sigmaBetaVec <- sample.sigmaBetaVec(BetaMat, muBetaVec, alphaBetaVec0,  betaBetaVec0, covNames)
	sigmaBetaVecIter = rbind(sigmaBetaVecIter, sigmaBetaVec)
	
	# 16. [Global] Update cached covariate effect
	outCovEffMonth <- calCovEffMonth(yMonth, BetaMat, covNames, T)
	covEffMonth = outCovEffMonth$covEffMonth
	yDehedonicMonthMean =  outCovEffMonth$yDehedonicMonthMean
	
	outCovEffMonthTest <- calCovEffMonth(yMonthTest, BetaMat, covNames, T)
	covEffMonthTest = outCovEffMonthTest$covEffMonth
	yDehedonicMonthMeanTest =  outCovEffMonthTest$yDehedonicMonthMean
}
		
	# Book keeping
	xMatIter=abind(xMatIter, xMat, along=3)
	LambdaMatIter= rbind(LambdaMatIter, LambdaMat[cbind(CtIds, as.character(zVec))] )
	aVecIter = rbind(aVecIter, aVec)
	RVecIter = rbind(RVecIter, RVec)
	BetaMatIter = abind(BetaMatIter, BetaMat, along=3)
	etaCtMat = etaMat[, as.character(zVec)]
	colnames(etaCtMat) = names(zVec)
	etaMatIter = abind(etaMatIter, etaCtMat, along=3)
	x0VecIter = rbind(x0VecIter, x0Vec)
	
	#############################################
	### Calculate log posterior               ###
	#############################################
	
	logPosterior = getLogPosterior(CtIds, zVec, muA, muA0, sigmaA0, sigmaA, alphaA0,betaA0,aVec,
						muLambda,muLambda0,sigmaLambda0,sigmaLambda, alphaLambda0,betaLambda0,LambdaMat,
						etaMat, sigma0,alpha.epsilon0,beta.epsilon0,
						alpha, alpha.alpha0,beta.alpha0, RVec, alpha.R0, beta.R0,
						covNames, BetaMat, muBetaVec, sigmaBetaVec,alphaBetaVec0,betaBetaVec0,muBetaVec0,sigmaBetaVec0,
						x0Vec, xMat, mux0, Vx0,
						L, yMonth,yTract, details=FALSE)
						
	logPosteriorTest = getLogPosterior(CtIds, zVec, muA, muA0, sigmaA0, sigmaA, alphaA0,betaA0,aVec,
						muLambda,muLambda0,sigmaLambda0,sigmaLambda, alphaLambda0,betaLambda0,LambdaMat,
						etaMat, sigma0,alpha.epsilon0,beta.epsilon0,
						alpha, alpha.alpha0,beta.alpha0, RVec, alpha.R0, beta.R0,
						covNames, BetaMat, muBetaVec, sigmaBetaVec,alphaBetaVec0,betaBetaVec0,muBetaVec0,sigmaBetaVec0,
						x0Vec, xMat, mux0, Vx0,
						L, yMonthTest,yTractTest, details=FALSE)

	logPosteriorIter = c(logPosteriorIter, logPosterior)
	logPosteriorTestIter = c(logPosteriorTestIter, logPosteriorTest)
	
	logPosterior.integrateEtaX = getLogPosterior.integrateEtaX(CtIds, zVec, muA, muA0, sigmaA0, sigmaA, alphaA0,betaA0,aVec,
						muLambda,muLambda0,sigmaLambda0,sigmaLambda, alphaLambda0,betaLambda0,LambdaMat,
						sigma0,alpha.epsilon0,beta.epsilon0,
						alpha, alpha.alpha0,beta.alpha0, RVec, alpha.R0, beta.R0,
						covNames, BetaMat, muBetaVec, sigmaBetaVec,alphaBetaVec0,betaBetaVec0,muBetaVec0,sigmaBetaVec0,
						x0Vec, mux0, Vx0,
						L, yMonth, yTract, covEffMonth, details=FALSE, lambdaIsVec = FALSE,
						yDehedonicMonthMean)
	
	logPosteriorTest.integrateEtaX = getLogPosterior.integrateEtaX(CtIds, zVec, muA, muA0, sigmaA0, sigmaA, alphaA0,betaA0,aVec,
						muLambda,muLambda0,sigmaLambda0,sigmaLambda, alphaLambda0,betaLambda0,LambdaMat,
						sigma0,alpha.epsilon0,beta.epsilon0,
						alpha, alpha.alpha0,beta.alpha0, RVec, alpha.R0, beta.R0,
						covNames, BetaMat, muBetaVec, sigmaBetaVec,alphaBetaVec0,betaBetaVec0,muBetaVec0,sigmaBetaVec0,
						x0Vec, mux0, Vx0,
						LTest, yMonthTest, yTractTest, covEffMonthTest, details=FALSE, lambdaIsVec = FALSE,
						yDehedonicMonthMeanTest)

	logPosteriorIter.integrateEtaX = c(logPosteriorIter.integrateEtaX, logPosterior.integrateEtaX)
	logPosteriorTestIter.integrateEtaX = c(logPosteriorTestIter.integrateEtaX, logPosteriorTest.integrateEtaX)
			
	cat("iter", iter, ": lPost =",round(logPosterior), " lPost.Test =",round(logPosteriorTest), 
		" lPostInte =",round(logPosterior.integrateEtaX), " lPostInte.Test =",round(logPosteriorTest.integrateEtaX), 
		", K =",K, "alpha =", round(alpha,2), "\n")
	print(zVec)
	cat("\n")
	if (PrintBeta0) {
		BetaPrint = rbind(round(BetaMat[,1]/scaleFactor,2), round(Beta[1,],2) )
		rownames(BetaPrint) = c("Intercept","Inte true")
		print(BetaPrint)
	}
	cat("\n")
}

if (nP>1){
	cl<-makeCluster(nP) #specify the nP to your number of CPU cores
	registerDoSNOW(cl)
}

ptm <- proc.time()
	
for (iter in 1:totalIter){
	tm1 = proc.time()
	
	if (nP>1){
		parallelOut <- foreach(p = 1:nP) %dopar%{
			#for (p in 1:nP){
			library("mvtnorm"); library("pscl"); library("abind"); library("plyr");
			library("Matrix")

			clusterLabels.p = as.numeric(names(piVec)[piVec==p])
			zVec.p = zVec[zVec %in% clusterLabels.p ]
			CtIds.p = as.numeric(names(zVec.p))

			if (length(zVec.p)==0){
				out.p = NA
			} else {
				# 1. Sample membership zi, update zVec and K. length(zVec) = n, dim(LambdaMat)=n*K
				out = sample.z(zVec=zVec.p, alpha=alpha/nP, aVec, LambdaMat, sigma0, RVec, covEffMonth,   
								 muLambda, sigmaLambda, tmLab, CtIds=CtIds.p, L, yMonth, is.debug=FALSE, mux0, Vx0)
				#out = sample.z.suffStat(zVec=zVec.p, alpha=alpha/nP, aVec, LambdaMat, sigma0, RVec, covEffMonth,yDehedonicMonthMean,  
				#				muLambda, sigmaLambda, tmLab, CtIds=CtIds.p, L, yMonth, is.debug=FALSE, mux0, Vx0)
				zVec.p = out$zVec; K.p = out$K; LambdaMat.p = out$LambdaMat; clusterLabels.p = sort(unique(zVec.p))
				BetaMat.p = BetaMat[CtIds %in% CtIds.p,]
				if (length(CtIds.p)==1) {
					BetaMat.p = matrix(BetaMat.p,1,length(covNames))
					rownames(BetaMat.p) = CtIds.p
					colnames(BetaMat.p) = covNames
				}
				aVec.p = aVec[as.character(CtIds.p)]

				for (subIter in 1:5){
					# 2. Sample the tract process x, update x. dim(xMat) = n*T
					outX = sample.x(zVec=zVec.p, LambdaMat=LambdaMat.p, aVec, RVec, sigma0 ,L, yMonth, covEffMonth, 
								tmLab, CtIds=CtIds.p, covNames, mux0, Vx0)
					xMat.p = outX$xMat
					x0Vec.p = outX$x0Vec
					
					# 3. Sample the cluster latent factor eta_t,k. dim(etaMat) is T*K
					etaMat.p = sample.eta(zVec=zVec.p, LambdaMat=LambdaMat.p, aVec, sigma0, xMat=xMat.p, tmLab, x0Vec=x0Vec.p)

					# 4. Sample the loading matrix lambda_ik. dim(LambdaMat) = n*K
					LambdaMat.p = sample.LambdaMat(zVec=zVec.p, muLambda, sigmaLambda, sigma0, etaMat=etaMat.p, xMat=xMat.p, aVec, tmLab, CtIds=CtIds, x0Vec=x0Vec.p)

					# 5. Sample AR coefficient a_i. dim(aVec) = n
					aVec.p <- sample.aVec(zVec=zVec.p, LambdaMat=LambdaMat.p, xMat=xMat.p,sigma0, muA, sigmaA, CtIds=CtIds.p, tmLab, x0Vec=x0Vec.p)

					# 6. Sample the response variance R_i. length(RVec) = n
					RVec.p <- sample.RVec(CtIds=CtIds.p, xMat=xMat.p, alpha.R0, beta.R0, yTract, tmLab, covEffMonth )
					
					# 7. Sample the covariate effect Beta_ir. dim(BetaMat)= n*length(covNames)
					BetaMat.p = BetaMat[CtIds %in% CtIds.p,]
					if (length(CtIds.p)==1) {
						BetaMat.p = matrix(BetaMat.p,1,length(covNames))
						rownames(BetaMat.p) = CtIds.p
						colnames(BetaMat.p) = covNames
					}
					BetaMat.p <- sample.BetaMat(muBetaVec, sigmaBetaVec, RVec=RVec.p, xMat=xMat.p, BetaMat=BetaMat.p, yTract, CtIds=CtIds.p, tmLab,covNames)
				}
				out.p = list(procId=p, clusterLabels.p=clusterLabels.p, CtIds.p=CtIds.p, zVec.p = zVec.p, K.p=K.p,  
					xMat.p = xMat.p, etaMat.p = etaMat.p, LambdaMat.p = LambdaMat.p,
					aVec.p=aVec.p, RVec.p = RVec.p, BetaMat.p = BetaMat.p, x0Vec.p = x0Vec.p)
			}
			out.p
		}	
		
		###### Assign the global cluster labels to clusters generated from diff processors
		globalClusterTable = NULL
		for (p in 1:nP){
			out.p = parallelOut[[p]]
			if (!all(is.na(out.p))){
				tem = cbind(out.p$procId, out.p$clusterLabels.p)
				globalClusterTable = rbind(globalClusterTable, tem)
			}	
		}
		globalClusterTable = cbind(globalClusterTable, 1:nrow(globalClusterTable))
		colnames(globalClusterTable)=c("procId","pClusterId", "gClusterId")
		globalClusterTable=data.frame(globalClusterTable)
		
		# update piVec, using the global cluster id and updated clusters
		piVec = globalClusterTable$procId
		names(piVec) = globalClusterTable$gClusterId
		
		###### Combine the results from nP processors 
		for (p in 1:nP){
			out.p = parallelOut[[p]]
			if (!all(is.na(out.p))){
				zVec.p = out.p$zVec.p
				tab = globalClusterTable[globalClusterTable$procId==p,]
				id = match(zVec.p, tab$pClusterId)
				zVec.pg = tab$gClusterId[id]
				names(zVec.pg) = names(zVec.p)
				zVec[names(zVec.pg)] = zVec.pg
			}
		}
		K = length(unique(zVec))

		xMat = etaMat=LambdaMat=aVec=RVec=BetaMat=x0Vec=NULL
		for (p in 1:nP){
			out.p = parallelOut[[p]]
			if (!all(is.na(out.p))){
				xMat = rbind(xMat, out.p$xMat.p)
				etaMat = cbind(etaMat, out.p$etaMat.p)
				LambdaMat = cbind(LambdaMat, out.p$LambdaMat.p)
				aVec = c(aVec, out.p$aVec.p)
				RVec = c(RVec, out.p$RVec.p)
				BetaMat = rbind(BetaMat, out.p$BetaMat.p)
				x0Vec = c(x0Vec, out.p$x0Vec.p)
			}
		}
		ord = order(rownames(xMat))
		xMat = xMat[ord,]
		colnames(etaMat) = 1:K
		colnames(LambdaMat) = 1:K
		ord = order(names(aVec))
		aVec = aVec[ord]
		ord = order(names(RVec))
		RVec = RVec[ord]
		ord = order(rownames(BetaMat))
		BetaMat = BetaMat[ord,]
		x0Vec = x0Vec[ord]

	} 
	
	if (nP==1){
		if (!fixZ){
			out = sample.z.suffStat(zVec=zVec, alpha=alpha/nP, aVec, LambdaMat, sigma0, RVec, covEffMonth,yDehedonicMonthMean,  
									muLambda, sigmaLambda, tmLab, CtIds=CtIds, L, yMonth, is.debug=FALSE, mux0, Vx0)
			
			#out = sample.z(zVec, alpha, aVec, LambdaMat, sigma0, RVec, covEffMonth,   
			#		muLambda, sigmaLambda, tmLab, CtIds, L, yMonth, is.debug=TRUE, mux0, Vx0)
			
			zVec = out$zVec; K = out$K; LambdaMat = out$LambdaMat; clusterLabels = sort(unique(zVec))
		}
	}
	
	t1 = proc.time()
	
	# Book keeping
	zVecIter = rbind(zVecIter, zVec)
	KIter = c(KIter, K)
	
	for (subIter in 1:5){
		# 2. Sample the tract process x, update x. dim(xMat) = n*T
		if (T){
		outX = sample.x(zVec=zVec, LambdaMat=LambdaMat, aVec, RVec, sigma0 ,L, yMonth, covEffMonth, 
						tmLab, CtIds=CtIds, covNames, mux0, Vx0)
		xMat = outX$xMat
		x0Vec = outX$x0Vec
		}
		
		if (T){
		# 3. Sample the cluster latent factor eta_t,k. dim(etaMat) is T*K
		etaMat = sample.eta(zVec, LambdaMat, aVec, sigma0, xMat, tmLab, x0Vec)

		# 4. Sample the loading matrix lambda_ik. dim(LambdaMat) = n*K
		LambdaMat = sample.LambdaMat(zVec, muLambda, sigmaLambda, sigma0, etaMat, xMat, aVec, tmLab, CtIds, x0Vec)
		}
		# 5. Sample AR coefficient a_i. dim(aVec) = n
		if (T){
		aVec <- sample.aVec(zVec, LambdaMat, xMat,sigma0, etaMat, muA, sigmaA, CtIds, tmLab, x0Vec)
		}
		if (T){
		# 6. Sample the response variance R_i. length(RVec) = n
		RVec <- sample.RVec(CtIds, xMat, alpha.R0, beta.R0, yTract, tmLab, covEffMonth )
		
		# 7. Sample the covariate effect Beta_ir. dim(BetaMat)= n*length(covNames)
		BetaMat <- sample.BetaMat(muBetaVec, sigmaBetaVec, RVec, xMat, BetaMat, yTract, CtIds, tmLab,covNames)
		}
			
		# 8. [Global] Sample CPR alpha, update alpha. A scalar
		if (T){
		alpha = sample.alpha(alpha, n, alpha.alpha0, beta.alpha0, K)
		alphaIter = rbind(alphaIter, alpha)
		}

		xMatIter=abind(xMatIter, xMat, along=3)
		LambdaMatIter= rbind(LambdaMatIter, LambdaMat[cbind(CtIds, as.character(zVec))] )
		aVecIter = rbind(aVecIter, aVec)
		RVecIter = rbind(RVecIter, RVec)
		BetaMatIter = abind(BetaMatIter, BetaMat, along=3)
		etaCtMat = etaMat[,as.character(zVec)]
		colnames(etaCtMat) = names(zVec)
		etaMatIter = abind(etaMatIter, etaCtMat, along=3)
		x0VecIter = rbind(x0VecIter, x0Vec)

		# 9. [Global] Sample x process error, a scalar
		sigma0 <- sample.sigma0(xMat, etaMat, aVec, LambdaMat, zVec, alpha.epsilon0, beta.epsilon0, CtIds, tmLab, x0Vec)
		sigma0Iter= c(sigma0Iter, sigma0)
		
		if (T){
		# 10.[Global] Sample muLambda, a scalar
		muLambda <- sample.muLambda(muLambda0, sigmaLambda0, LambdaMat, zVec, sigmaLambda)
		muLambdaIter = c(muLambdaIter, muLambda)
		
		# 11.[Global] Sample sigmaLambda, a scalar
		sigmaLambda <- sample.sigmaLambda(zVec, LambdaMat, muLambda, alphaLambda0, betaLambda0)
		sigmaLambdaIter = c(sigmaLambdaIter, sigmaLambda)
		}
		
		if (T){
		# 12.[Global] Sample muA, a scalar
		muA <- sample.muA(aVec, sigmaA, muA0, sigmaA0)
		muAIter = c(muAIter, muA)
		
		# 13.[Global] Sample sigmaA, a scalar
		sigmaA <- sample.sigmaA(aVec, muA, alphaA0, betaA0)
		sigmaAIter = c(sigmaAIter, sigmaA)
		}
		if (T){
		# 14.[Global] Sample the mean covariate effect, muBetaVec. length(muBetaVec) = length(covNames)
		muBetaVec <- sample.muBetaVec(BetaMat, sigmaBetaVec, muBetaVec0, sigmaBetaVec0, covNames)
		muBetaVecIter = rbind(muBetaVecIter, muBetaVec)
		
		# 15.[Global] Sample the variability of covariate effect, sigmaBetaVec, length(sigmaBetaVec) = length(covNames)
		sigmaBetaVec <- sample.sigmaBetaVec(BetaMat, muBetaVec, alphaBetaVec0,  betaBetaVec0, covNames)
		sigmaBetaVecIter = rbind(sigmaBetaVecIter, sigmaBetaVec)
		
		# 16. [Global] Update cached covariate effect
		outCovEffMonth <- calCovEffMonth(yMonth, BetaMat, covNames, T)
		covEffMonth = outCovEffMonth$covEffMonth
		yDehedonicMonthMean =  outCovEffMonth$yDehedonicMonthMean
		
		outCovEffMonthTest <- calCovEffMonth(yMonthTest, BetaMat, covNames, T)
		covEffMonthTest = outCovEffMonthTest$covEffMonth
		yDehedonicMonthMeanTest =  outCovEffMonthTest$yDehedonicMonthMean
		}
	}
	
	#############################################
	### Calculate log posterior               ###
	#############################################

	logPosterior = getLogPosterior(CtIds, zVec, muA, muA0, sigmaA0, sigmaA, alphaA0,betaA0,aVec,
						muLambda,muLambda0,sigmaLambda0,sigmaLambda, alphaLambda0,betaLambda0,LambdaMat,
						etaMat, sigma0,alpha.epsilon0,beta.epsilon0,
						alpha, alpha.alpha0,beta.alpha0, RVec, alpha.R0, beta.R0,
						covNames, BetaMat, muBetaVec, sigmaBetaVec,alphaBetaVec0,betaBetaVec0,muBetaVec0,sigmaBetaVec0,
						x0Vec, xMat, mux0, Vx0,
						L, yMonth,yTract)
	logPosteriorTest = getLogPosterior(CtIds, zVec, muA, muA0, sigmaA0, sigmaA, alphaA0,betaA0,aVec,
						muLambda,muLambda0,sigmaLambda0,sigmaLambda, alphaLambda0,betaLambda0,LambdaMat,
						etaMat, sigma0,alpha.epsilon0,beta.epsilon0,
						alpha, alpha.alpha0,beta.alpha0, RVec, alpha.R0, beta.R0,
						covNames, BetaMat, muBetaVec, sigmaBetaVec,alphaBetaVec0,betaBetaVec0,muBetaVec0,sigmaBetaVec0,
						x0Vec, xMat, mux0, Vx0,
						L, yMonthTest,yTractTest)

	logPosteriorIter = c(logPosteriorIter, logPosterior)
	logPosteriorTestIter = c(logPosteriorTestIter, logPosteriorTest)

	logPosterior.integrateEtaX = getLogPosterior.integrateEtaX(CtIds, zVec, muA, muA0, sigmaA0, sigmaA, alphaA0,betaA0,aVec,
						muLambda,muLambda0,sigmaLambda0,sigmaLambda, alphaLambda0,betaLambda0,LambdaMat,
						sigma0,alpha.epsilon0,beta.epsilon0,
						alpha, alpha.alpha0,beta.alpha0, RVec, alpha.R0, beta.R0,
						covNames, BetaMat, muBetaVec, sigmaBetaVec,alphaBetaVec0,betaBetaVec0,muBetaVec0,sigmaBetaVec0,
						x0Vec, mux0, Vx0,
						L, yMonth, yTract, covEffMonth, details=FALSE, lambdaIsVec = FALSE,
						yDehedonicMonthMean)
						
	logPosteriorTest.integrateEtaX = getLogPosterior.integrateEtaX(CtIds, zVec, muA, muA0, sigmaA0, sigmaA, alphaA0,betaA0,aVec,
						muLambda,muLambda0,sigmaLambda0,sigmaLambda, alphaLambda0,betaLambda0,LambdaMat,
						sigma0,alpha.epsilon0,beta.epsilon0,
						alpha, alpha.alpha0,beta.alpha0, RVec, alpha.R0, beta.R0,
						covNames, BetaMat, muBetaVec, sigmaBetaVec,alphaBetaVec0,betaBetaVec0,muBetaVec0,sigmaBetaVec0,
						x0Vec, mux0, Vx0,
						LTest, yMonthTest, yTractTest, covEffMonthTest, details=FALSE,lambdaIsVec = FALSE,
						yDehedonicMonthMeanTest)

	logPosteriorIter.integrateEtaX = c(logPosteriorIter.integrateEtaX, logPosterior.integrateEtaX)
	logPosteriorTestIter.integrateEtaX = c(logPosteriorTestIter.integrateEtaX, logPosteriorTest.integrateEtaX)
	
	
	cat("iter", iter, ": lPost =",round(logPosterior), " lPost.Test =",round(logPosteriorTest), 
		" lPostInte =",round(logPosterior.integrateEtaX), " lPostInte.Test =",round(logPosteriorTest.integrateEtaX), 
		", K =",K, "alpha =", round(alpha,2),  "time(s) =", round((proc.time()-tm1)[3]), 
		"z sampler time =", round((t1 - tm1)[3]),"\n", "zVec =", zVec, "\n")
	if (nP>1) {print(piVec); cat("\n")}
	if (PrintBeta0){
		BetaPrint = rbind(round(BetaMat[,1]/scaleFactor,2), round(Beta[1,],2) )
		rownames(BetaPrint) = c("Intercept","Inte true")
		print(BetaPrint)
		cat("\n")	
	}
		
	# 17. [Global] Sample the processor allocation indicator piVec
	if (nP>1) piVec = sample.pi(piVec, nP, zVec, K)	
	
}

if (nP>1) stopCluster(cl)

tim = proc.time() - ptm
tim 
# 6 hr =  1200 iters

#outPath = "Output/simuParaMCMC_suffStat/x0Set0_Data_a99_intercept/alphaCondEta_Beta0InModel_recordEta/alphaCondEta_fixedOneCtPerCluster/"
#outPath = "Output/simuParaMCMC_suffStat/x0Set0_Data_a99_intercept/alphaCondEta_NoBeta0InModel_recordEta/"
#outPath = "Output/simuParaMCMC_suffStat/x0Set0_Data_a99_intercept_largeLambda/alphaCondEta_Beta0InModel_largeLambda/"
#outPath = "Output/simuParaMCMC_suffStat/x0Set0_Data_a99_intercept_largeLambda/alphaCondEta_Beta0InModel_largeLambda/alphaCondEta_fixedOneCtPerCluster/"
#outPath = "Output/simuParaMCMC_suffStat/x0Set0_Data_a99_intercept/alphaCondEta_NoBeta0InModel_recordEta/alphaCondEta_fixedOneCtPerCluster/"
#outPath = "Output/simuParaMCMC_suffStat/x0Set0_Data_a99_intercept_largeLambda_largeR/alphaCondEta_Beta0InModel/"
#outPath = "Output/simuParaMCMC_suffStat/x0Set0_Data_a99_intercept_largeLambda_largeR/alphaCondEta_Beta0InModel/alphaCondEta_fixedOneCtPerCluster/"
#outPath = "Output/simuParaMCMC_suffStat/x0Set0_Data_a60_intercept_largeLambda_largeR/alphaCondEta_Beta0InModel/"
#outPath = "Output/simuParaMCMC_suffStat/x0Set0_Data_a60_intercept_largeLambda_largeR/alphaCondEta_Beta0InModel/alphaCondEta_fixedOneCtPerCluster/"
#outPath = "Output/simuParaMCMC_suffStat/x0Set0_Data_a60_intercept_largeLambda/alphaCondEta_Beta0InModel/alphaCondEta_fixedOneCtPerCluster/"


save.image(file=paste(outPath,"unscaledResults",".RData", sep=""))



		