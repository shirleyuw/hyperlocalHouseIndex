##############################
## Cached: covariate effect ##
##############################
Ltract <- aggregate(count~CensusTractRegionID, data=L, FUN=sum)

calCovEffMonth <- function(yMonth, BetaMat, covNames, T){
	covEffMonth = list()
	yDehedonicMonthMean =  list()
	for (tt in 1:T){
		ytt = yMonth[[tt]]
		if (nrow(ytt)==0){
			covEffMonth[[tt]] = NA
			yDehedonicMonthMean[[tt]] = NA
		} else {
			covariates = ytt[,c("CensusTractRegionID" , covNames)]
			covEff = as.numeric( apply( covariates[,covNames]*BetaMat[as.character(covariates$CensusTractRegionID),covNames], 
						MARGIN=1, FUN=sum ) )
			yDehedonic = ytt$y - covEff
			covEffMonth[[tt]] = data.frame(CensusTractRegionID=covariates$CensusTractRegionID, covEff, yDehedonic)
			tem = tapply(yDehedonic, INDEX=ytt$CensusTractRegionID, FUN=mean)
			CtRegionID = sort(unique(ytt$CensusTractRegionID))
			yDehedonicMonthMean[[tt]] =data.frame(CensusTractRegionID=CtRegionID,yDehedonicMonthMean=tem )
		}
	}
	return(list(covEffMonth=covEffMonth, yDehedonicMonthMean=yDehedonicMonthMean))
}

###################
## Get the matrix (data and coefficient)ready for Kalman Filter/Smoother##
###################

get.Ct.yt <- function(L, yMonth, covEffMonth, CtIds.k, tmLab){
	CtIds.k = sort(CtIds.k)
	nk = length(CtIds.k)
	L.k = L[L$CensusTractRegionID %in% CtIds.k,]
	T = length(tmLab)
	Ct = list() # dim is varying by t, dim_t = (mt, nk), the indicator matrix from obs to tract at t
	tractMatch = list()
	yt=list()
	covEfft = list()
	
	for (tt in 1:T){
		L.k.tt = L.k[L.k$MonthTransaction==tmLab[tt],]
		mk=sum(L.k.tt$count)
		if (mk==0){
			Ct[[tt]] = NA
			tractMatch[[tt]] = NA
			yt[[tt]] = NA
			covEfft[[tt]] = NA
			
		} else {
			colMatch = rep(1:nk, L.k.tt$count)
			tem = matrix(0, mk, nk )
			for (j in 1:mk){
				tem[j,colMatch[j]] = 1
			}
			Ct[[tt]] = tem
			tractMatch[[tt]] = colMatch

			tem = yMonth[[tt]]
			idct = tem$CensusTractRegionID %in% CtIds.k		
			yt[[tt]] = tem$y[idct]

			tem = covEffMonth[[tt]]
			idct = tem$CensusTractRegionID %in% CtIds.k		
			covEfft[[tt]] = tem$covEff[idct]
		}
	}
	return(list(Ct=Ct, tractMatch=tractMatch, yt=yt, covEfft=covEfft))
}

get.Ct.yt.suffStat <- function(L, covEffMonth, yDehedonicMonthMean, CtIds.k, tmLab, tract.i){
	CtIds.k = sort(CtIds.k)
	nk = length(CtIds.k)
	L.k = L[L$CensusTractRegionID %in% CtIds.k,]
	T = length(tmLab)
	Lt = list()
	yDehedonict.i = list()
	yDehedonicMean = list()
	CtBar = list()
	
	for (tt in 1:T){
		L.k.tt = L.k[L.k$MonthTransaction==tmLab[tt],]
		Lcount = L.k.tt$count
		Lt[[tt]] = Lcount
		if (sum(Lcount)==0){
			yDehedonict.i[[tt]] = NA
			yDehedonicMean[[tt]] = NA
			CtBar[[tt]] =NA
		} else {
			tem = covEffMonth[[tt]]
			idct = tem$CensusTractRegionID == tract.i
			if (sum(idct)>0){
				yDehedonict.i[[tt]] = tem$yDehedonic[ idct]
			} else {
				yDehedonict.i[[tt]] = NA
			}
			tem = yDehedonicMonthMean[[tt]]
			yDehedonicMean[[tt]] = tem$yDehedonicMonthMean[tem$CensusTractRegionID %in% CtIds.k]
			
			if (TRUE){
			tem = matrix(0, nk, nk)
			tem[cbind(1:nk, 1:nk)] = 1
			idct =Lcount>0
			tem = tem[ idct,]
			if (sum(idct)==1) tem = matrix(tem, 1, nk)
			CtBar[[tt]] = tem 
			}
		}
	}
	return(list(Lt=Lt, yDehedonict.i=yDehedonict.i, yDehedonicMean=yDehedonicMean, CtBar=CtBar))
}

get.Ct.yt.old <- function(obs, tmLab, CtIds.k, BetaMat.k, covNames){
	CtIds.k = sort(CtIds.k)
	obs = obs[obs$CensusTractRegionID %in% CtIds.k,]
	ord = order(obs$MonthTransaction, obs$CensusTractRegionID, obs$TransactionID)
	obs = obs[ord,]
	T = length(tmLab)
	nk = length(CtIds.k)
	Ct = list() # dim is varying by t, dim_t = (mt, nk), the indicator matrix from obs to tract at t
	tractMatch = list()
	yt = list() # dim is time varying, dim_t =( mt,1)
	covEfft = list() # dim is time varying, dim_t = (mt, 1)

	for (i in 1:T){
		obs.t = obs[obs$MonthTransaction==tmLab[i],]
		
		if (nrow(obs.t)==0){
			Ct[[i]] = NA
			tractMatch[[i]] = NA
		} else {
			colMatch = match(obs.t$CensusTractRegionID, CtIds.k)
			tem = matrix(0, nrow(obs.t), nk )
			for (j in 1:nrow(obs.t)){
				tem[j,colMatch[j]] = 1
			}
			Ct[[i]] = tem
			tractMatch[[i]] = colMatch
		}
	}
	
	## Organize the response y
	for (i in 1:T){
		obs.t = obs[obs$MonthTransaction==tmLab[i],]
		if (nrow(obs.t)==0){
			yt[[i]] = NA
		} else {
			yt[[i]] = obs.t$y
		}
	}

	
	## Organize covariate effect
	for (i in 1:T){
		obs.t = obs[obs$MonthTransaction==tmLab[i],]
		if (nrow(obs.t)==0){
			covEfft[[i]] = NA
		} else {
			covEfft[[i]] = as.numeric(apply(obs.t[,covNames] * BetaMat.k[tractMatch[[i]],], MARGIN=1,FUN=sum))
		}
	}
	
	return(list(Ct=Ct, tractMatch=tractMatch, yt=yt, covEfft=covEfft))
}	

library("Matrix")
#save(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, tractMatch, Ct, yt,  covEfft, tract.i, file="inputKF.RData")

# won't use it
get.Ct.yt.SparseMat <- function(obs, tmLab, CtIds.k, BetaMat.k, covNames){
	CtIds.k = sort(CtIds.k)
	obs = obs[obs$CensusTractRegionID %in% CtIds.k,]
	ord = order(obs$MonthTransaction, obs$CensusTractRegionID, obs$TransactionID)
	obs = obs[ord,]
	T = length(tmLab)
	nk = length(CtIds.k)
	Ct = list() # dim is varying by t, dim_t = (mt, nk), the indicator matrix from obs to tract at t
	tractMatch = list()
	for (tt in 1:T){
		obs.t = obs[obs$MonthTransaction==tmLab[tt],]
		if (nrow(obs.t)==0){
			Ct[[tt]] = NA
			tractMatch[[tt]] = NA
		} else {
			colMatch = match(obs.t$CensusTractRegionID, CtIds.k)
			mt = nrow(obs.t)
			tem = sparseMatrix(i = 1:mt, j = colMatch[1:mt], x=1, dims=c(mt, nk))
			Ct[[tt]] = tem
			tractMatch[[tt]] = colMatch
		}
	}

	## Organize the response y
	yt = list() # dim is time varying, dim_t =( mt,1)
	for (tt in 1:T){
		obs.t = obs[obs$MonthTransaction==tmLab[tt],]
		if (nrow(obs.t)==0){
			yt[[tt]] = NA
		} else {
			yt[[tt]] = obs.t$y
		}
	}
	
	## Organize covariate effect
	covEfft = list() # dim is time varying, dim_t = (mt, 1)
	for (tt in 1:T){
		obs.t = obs[obs$MonthTransaction==tmLab[tt],]
		if (nrow(obs.t)==0){
			covEfft[[tt]] = NA
		} else {
			covEfft[[tt]] = as.numeric(apply(obs.t[,covNames] * BetaMat.k[tractMatch[[tt]],], MARGIN=1,FUN=sum))
		}
	}
	
	return(list(Ct=Ct, tractMatch=tractMatch, yt=yt, covEfft=covEfft))
}	

get.matrix.start <- function(obs, tmLab, CtIds,zVec, BetaMat, covNames){
	K = length(unique(zVec))
	T = length(tmLab)
	out = list()
	
	for (k in unique(zVec)){
		CtIds.k = CtIds[zVec==k]
		CtIds.k = sort(CtIds.k)
		obsk = obs[obs$CensusTractRegionID %in% CtIds.k,]
		ord = order(obsk$MonthTransaction, obsk$CensusTractRegionID, obsk$TransactionID)
		obsk = obsk[ord,]
		T = length(tmLab)
		nk = length(CtIds.k)
		mt = rep(0, T)
		names(mt) = tmLab
		mt.count=tapply(obsk$y, INDEX=obsk$MonthTransaction, FUN=length)
		mt[match(names(mt.count), as.character(tmLab))]=mt.count
		
		BetaMat.k = BetaMat[as.character(CtIds.k),]
		if (!is.matrix(BetaMat.k)) {
			BetaMat.k = matrix(BetaMat.k, nrow=1)
			rownames(BetaMat.k) = CtIds.k
		}
		
		Ct = list() # dim is varying by t, dim_t = (mt, nk), the indicator matrix from obs to tract at t
		tractMatch = list()
		yt = list() # dim is time varying, dim_t =( mt,1)
		covEfft = list() # dim is time varying, dim_t = (mt, 1)

		tm1 =proc.time()
		
		for (tt in 1:T){
			idct = obsk$MonthTransaction==tmLab[tt]
			obs.t = obsk[idct,]
			if (mt[tt]==0){
				Ct[[tt]] = NA
				tractMatch[[tt]] = NA
				yt[[tt]] = NA
				covEfft[[tt]] = NA
			} else {
				colMatch = match(obs.t$CensusTractRegionID, CtIds.k)
				tem = matrix(0, mt[tt],nk)
				tem[cbind(1:mt[tt],colMatch)] = 1
				#tem = sparseMatrix(i=1:mt[tt], j=colMatch, x=1, dims=c(mt[tt],nk))
				Ct[[tt]] = as.matrix(tem)
				tractMatch[[tt]] = colMatch
				yt[[tt]] = obs.t$y
				covEfft[[tt]] = as.numeric(apply(obs.t[,covNames] * BetaMat.k[tractMatch[[tt]],], MARGIN=1,FUN=sum))
			}
		}
		tm2 = proc.time()
		tm2-tm1
		listk = list(Ct=Ct, tractMatch=tractMatch, yt=yt, covEfft=covEfft)
		out[[as.character(k)]]=listk
	}
	return(out)
}	

get.matrix.try <- function(outMat, obs,tmLab, CtIds, i, k, zVec, BetaMat, covNames){
	kLabels = unique(zVec)
	K = length(kLabels)
	T = length(tmLab)
	
	if (zVec[i]==k){
		outMatk = outMat[[as.character(k)]]
	} else if (k=="newCluster"){
	# add tract i to cluster K+1
		obsi = obs[obs$CensusTractRegionID==CtIds[i],]
		Ct = tractMatch= yt = covEfft= list()
		for (tt in 1:T){
			idct = obsi$MonthTransaction==tmLab[tt]
			mi = sum(idct)
			if (mi==0){
				Ct[[tt]] = tractMatch[[tt]] = yt[[tt]] = covEfft[[tt]] = NA
			} else {
				Ct[[tt]] = matrix(1, mi,1)
				tractMatch[[tt]] = rep(1, mi)
				yt[[tt]] = obsi$y[idct]
				covEfft[[tt]] = as.numeric(apply(obsi[idct,covNames] * BetaMat[as.character(CtIds[i]),], MARGIN=1,FUN=sum))
			}
		}
		outMatk = list(Ct=Ct, tractMatch=tractMatch, yt=yt, covEfft=covEfft)
	} else {
	 # add tract i to cluster k
		outMatk = outMat[[as.character(k)]]
		obsi = obs[obs$CensusTractRegionID==CtIds[i],]
		Ct = outMatk$Ct
		nk = ncol(Ct[[1]])
		tractMatch = outMatk$tractMatch
		yt = outMatk$yt
		covEfft = outMatk$covEfft
		
		for (tt in 1:T){
			idct = obsi$MonthTransaction==tmLab[tt]
			mi = sum(idct)
			if (mi==0){
				Ct[[tt]] =cbind(Ct[[tt]],0)
			} else {
				add=matrix(0, mi, nk+1)
				add[,nk+1] = 1
				Ct[[tt]] =rbind(cbind(Ct[[tt]],0),add )
				tractMatch[[tt]] = c(tractMatch[[tt]], rep(nk+1, mi))
				yt[[tt]] = c(yt[[tt]], obsi$y[idct] )
				tem =as.numeric(apply(obsi[idct,covNames] * BetaMat[as.character(CtIds[i]),], MARGIN=1,FUN=sum))
				covEfft[[tt]] = c(covEfft[[tt]],tem)
			}
		}
		outMatk = list(Ct=Ct, tractMatch=tractMatch, yt=yt, covEfft=covEfft)
	}
	return(outMatk)
}	

get.matrix.update <- function(outMat, obs,tmLab, CtIds, i, zVecOld, zVecNew,newCluster=FALSE, BetaMat, covNames){
	kOld = zVecOld[i]
	kNew = zVecNew[i]
	T = length(tmLab)
	
	if (kOld!=kNew){
		# remove tract i from its old cluster 
		oldCluster = outMat[[as.character(kOld)]]
		id = match( as.character(CtIds[i] ), names(zVecOld[zVecOld==kOld]))
		Ct = tractMatch=yt=covEfft=list()
		for (tt in 1:T){		 
			tem = oldCluster$tractMatch[[tt]] 
			idct = !tem==id
			tem = tem[idct]
			tem[tem>id] = tem[tem>id]-1
			tractMatch[[tt]] = tem
			yt[[tt]] = oldCluster$yt[[tt]][idct]
			covEfft[[tt]] = oldCluster$covEfft[[tt]][idct]
			Ct[[tt]] = oldCluster$Ct[[tt]][idct,-id] 
		}
		outMat[[as.character(kOld)]] = list(Ct=Ct, tractMatch=tractMatch, yt=yt, covEfft=covEfft )
		
		if (newCluster){		
			# If tract i belongs to a brand new cluster, construc 
			obsi = obs[obs$CensusTractRegionID==CtIds[i],]
			Ct = tractMatch= yt = covEfft= list()
			for (tt in 1:T){
				idct = obsi$MonthTransaction==tmLab[tt]
				mi = sum(idct)
				if (mi==0){
					Ct[[tt]] = tractMatch[[tt]] = yt[[tt]] = covEfft[[tt]] = NA
				} else {
					Ct[[tt]] = matrix(1, mi,1)
					tractMatch[[tt]] = rep(1, mi)
					yt[[tt]] = obsi$y[idct]
					covEfft[[tt]] = as.numeric(apply(obsi[idct,covNames] * BetaMat[as.character(CtIds[i]),], MARGIN=1,FUN=sum))
				}
			}
			outMat[[as.character(kNew)]] = list(Ct=Ct, tractMatch=tractMatch, yt=yt, covEfft=covEfft)
		} else {
			# update a existing cluster k with tract i data. Sorted by the order of CtIds
			outMatk = outMat[[as.character(kNew)]]
			obsi = obs[obs$CensusTractRegionID==CtIds[i],]
			Ct = outMatk$Ct
			nk = ncol(Ct[[1]])
			tractMatch = outMatk$tractMatch
			yt = outMatk$yt
			covEfft = outMatk$covEfft
			id = match( as.character(CtIds[i] ), names(zVecNew[zVecNew==kNew]))
						
			for (tt in 1:T){
				idct = obsi$MonthTransaction==tmLab[tt]
				mi = sum(idct)
				if (mi==0){
					if (id ==1){
						Ct[[tt]] =cbind(0,Ct[[tt]])
						tractMatch[[tt]] = tractMatch[[tt]]+1
					} else if (id==sum(zVecNew==kNew)){
						Ct[[tt]] =cbind(Ct[[tt]], 0)
					} else {
						Ct[[tt]] = cbind(Ct[[tt]][,1:(id-1)],0,Ct[[tt]][,-c(1:(id-1))])
						tem = tractMatch[[tt]]
						tem[tem>=id] = tem[tem>=id]+1
						tractMatch[[tt]] = tem
					}
				} else {
					tem = tractMatch[[tt]]
					rowid = sum(tem<id)
					if (rowid==0){
						tractMatch[[tt]]=c(rep(id,mi), tem+1)
						yt[[tt]] = c(obsi$y[idct], yt[[tt]] )
						tem =as.numeric(apply(obsi[idct,covNames] * BetaMat[as.character(CtIds[i]),], MARGIN=1,FUN=sum))
						covEfft[[tt]] = c(tem, covEfft[[tt]])
						tem = Ct[[tt]]
						tem = cbind(0, tem)
						add = matrix(0, mi, ncol(tem))
						add[,1] = 1
						tem = rbind(add,tem)
						Ct[[tt]] = tem
					} else {
						tractMatch[[tt]]=c(tem[1:rowid], rep(id,mi), tem[-c(1:rowid)]+1)
						yt[[tt]] = c(yt[[tt]][1:rowid], obsi$y[idct], yt[[tt]][-c(1:rowid)])
						tem =as.numeric(apply(obsi[idct,covNames] * BetaMat[as.character(CtIds[i]),], MARGIN=1,FUN=sum))
						covEfft[[tt]] = c(covEfft[[tt]][1:rowid], tem, covEfft[[tt]][-c(1:rowid)])
						tem = Ct[[tt]]
						tem = cbind(tem[,1:(id-1)],0, tem[,-c(1:(id-1))])
						add = matrix(0, mi, ncol(tem))
						add[,id]=1
						tem = rbind(tem[1:rowid,], add ,tem[-c(1:rowid),])
						Ct[[tt]] = tem
					}
				}
			}
			outMat[[as.character(kNew)]] = list(Ct=Ct, tractMatch=tractMatch, yt=yt, covEfft=covEfft)
		}
	}
	return(outMat)
}	

# won't use it
get.matrix.updateBeta <- function(outMat, obs,tmLab, CtIds, zVec, BetaMat, covNames){
	K = length(unique(zVec))
	T = length(tmLab)
	out = outMat
	
	for (k in unique(zVec)){
		CtIds.k = CtIds[zVec==k]
		CtIds.k = sort(CtIds.k)
		obsk = obs[obs$CensusTractRegionID %in% CtIds.k,]
		ord = order(obsk$MonthTransaction, obsk$CensusTractRegionID, obsk$TransactionID)
		obsk = obsk[ord,]
		nk = length(CtIds.k)
		mt = rep(0, T)
		names(mt) = tmLab
		mt.count=tapply(obsk$y, INDEX=obsk$MonthTransaction, FUN=length)
		mt[match(names(mt.count), as.character(tmLab))]=mt.count
		
		BetaMat.k = BetaMat[as.character(CtIds.k),]
		if (!is.matrix(BetaMat.k)) {
			BetaMat.k = matrix(BetaMat.k, nrow=1)
			rownames(BetaMat.k) = CtIds.k
		}
		
		tractMatch = outMat[[as.character(k)]]$tractMatch
		covEfft = list() # dim is time varying, dim_t = (mt, 1)
	
		for (tt in 1:T){
			idct = obsk$MonthTransaction==tmLab[tt]
			obs.t = obsk[idct,]
			if (mt[tt]==0){
				covEfft[[tt]] = NA
			} else {
				covEfft[[tt]] = as.numeric(apply(obs.t[,covNames] * BetaMat.k[tractMatch[[tt]],], MARGIN=1,FUN=sum))
			}
		}

		out[[as.character(k)]]$covEfft = covEfft
	}
	return(out)
}	

###################
## Kalman Filter ##
###################

if (F){
tractMatch = out$tractMatch
Ct = out$Ct
#Lt = out$Lt
yt = out$yt
covEfft = out$covEfft
tract.i = CtIds[i] 
}

KalmanFilter <- function(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, tractMatch, Ct, yt,  covEfft, tract.i=NULL ){
	n = length(CtIds.k)
	pi0 = rep(0, n)
	V0 = 0.1*diag(rep(1,n))
	T = length(tmLab)
	A = matrix(0, n,n)
	diag(A) = aVec.k
	muPredict  =muFilter  = list()  # dim_t = (nk,1)
	VPredict = VFilter = list() # dim_t =( nk,nk)
	#Kt = list() # dim is varying by t, dim_t =(nk, mt), where mt=number of obs at time t for cluster k
	logMargLike = NULL
	logMargLike.tract.i =0
	
	for (i in 1:T){
		# assign time varying SSM parameters, Rt
		if (all(is.na(tractMatch[[i]]))){
			Rt.i = NA
		} else{
			if (length(tractMatch[[i]])==1){
				Rt.i = RVec.k[ tractMatch[[i]] ]
			} else{
				Rt.i = diag(RVec.k[ tractMatch[[i]] ])
			}
		}
		# muPredict = mu_{t|t-1}   t=1, ..., T
		if (i==1) {
			muPredict[[i]] = aVec.k*pi0
		} else {
			muPredict[[i]] = aVec.k*muFilter[[i-1]] 
		}
		
		# VPredict = V_{t|t-1}    t=1, ..., T
		if (i==1){
			VPredict[[i]] = A%*%V0 %*%A + Q.k
		} else {
			VPredict[[i]] = A%*% VFilter[[i-1]] %*% A + Q.k
		}
		
		# Kt = Kalman forward gain    t=1, ..., T
		if (all(is.na(Ct[[i]]))){
			#Kt[[i]] =  NA
			Kt =  NA
			muFilter[[i]] = muPredict[[i]]
			VFilter[[i]]  = VPredict[[i]]
		} else {
			S = Ct[[i]] %*% VPredict[[i]] %*% t(Ct[[i]]) + Rt.i
 			#Kt[[i]] = VPredict[[i]]%*%t(Ct[[i]])%*%chol2inv(chol(S))
			Kt = VPredict[[i]]%*%t(Ct[[i]])%*%chol2inv(chol(S))
			ytHat = Ct[[i]]%*%muPredict[[i]] + covEfft[[i]]
			residual = yt[[i]]-ytHat
			muFilter[[i]] = muPredict[[i]] + Kt%*% residual
			VFilter[[i]] = VPredict[[i]] - Kt%*%Ct[[i]] %*% VPredict[[i]]
		}
		
		# log marginal likelihood
		if (!all(is.na(yt[[i]]))){
			#p = length(yt[[i]])
			#log.det.S = sum(log(svd(S)$d))
			#logMargLike = logMargLike -p/2*log(2*pi) -1/2*log.det.S-1/2*t(residual) %*% chol2inv(chol(S)) %*% residual
			if (!is.null(tract.i)) {
				id = match(tract.i, CtIds.k)
				idct = tractMatch[[i]]==id
				p = sum(idct)
				if (p>0) {
					r.i = residual[idct]
					S.i = S[idct,idct]
					#log.det.S.i = ifelse(p==1, log(S.i), sum(log(svd(S.i)$d)))
					#logMargLike.tract.i = logMargLike.tract.i - p/2*log(2*pi)-1/2*log.det.S.i-1/2*t(r.i) %*% chol2inv(chol(S.i))%*%r.i
					if (p==1)
						logMargLike.tract.i = logMargLike.tract.i + dnorm(x=r.i, mean=0, sd=sqrt(S.i), log=TRUE)
					else 
						logMargLike.tract.i = logMargLike.tract.i + dmvnorm(x=r.i, mean=rep(0,length(r.i)), sigma=S.i, log=TRUE)

				}
			}
		}
	}
	return(list(logMargLike=logMargLike, logMargLike.tract.i=logMargLike.tract.i, muPredict=muPredict, 
		VPredict=VPredict, 	muFilter=muFilter,VFilter=VFilter))
}

if (FALSE){
  Lt= outSuff$Lt
  yDehedonicMean = outSuff$yDehedonicMean
  CtBar = outSuff$CtBar
}


KalmanFilter.suffStat <- function(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, Lt,  yDehedonicMean,CtBar, mux0, Vx0 ){
	n = length(CtIds.k)
	muFilter  = rep(mux0, n) # for t=0, prior
	VFilter = matrix(0, n,n)
	diag(VFilter) = Vx0
	T = length(tmLab)
	muPredictList = VPredictList = list()
	muFilterList = VFilterList = list()
	
	for (tt in 1:T){
		
		muPredict = aVec.k*muFilter
		#VPredict = A%*% VFilter %*% A + Q.k
		VPredict = (aVec.k %o% aVec.k) * VFilter + Q.k
		
		L.kt = Lt[[tt]]
		idctCtObs.t = L.kt>0
		nCtObs = sum(idctCtObs.t )
		
		if (nCtObs>0){ # cluster k has at least one obs at time t
			tem = RVec.k[idctCtObs.t]/L.kt[idctCtObs.t]
			RBar.tt = matrix(0, nCtObs, nCtObs)
			diag(RBar.tt) = tem
		
			if (nCtObs==n){
				S =  VPredict + RBar.tt
				Kt = VPredict %*%chol2inv(chol(S))

				muFilter = muPredict + Kt %*% (yDehedonicMean[[tt]]-muPredict[idctCtObs.t])
				VFilter = VPredict- Kt%*% VPredict
			} else {
				Ct = CtBar[[tt]]
				CV = Ct %*% VPredict 
				S =  CV %*% t(Ct) + RBar.tt
				Kt = t(CV) %*%chol2inv(chol(S))

				muFilter = muPredict + Kt %*% (yDehedonicMean[[tt]]-muPredict[idctCtObs.t])
				VFilter = VPredict - Kt %*% CV
			}
		} else {
			muFilter = muPredict
			VFilter  = VPredict
		}
		muPredictList[[tt]] = muPredict
		VPredictList[[tt]] = VPredict
		muFilterList[[tt]] = muFilter
		VFilterList[[tt]] = VFilter
		
	}
	return(list(muPredictList=muPredictList, VPredictList=VPredictList,
				muFilterList=muFilterList, VFilterList=VFilterList ))
}

#########
KalmanFilter.cond.eta <- function( tract.i, RVec.i, aVec.i, sigma0sq, eta.k, lambda.ik, covEffMonth, yMonth, tmLab, L, mux0, Vx0 ){
	muFilter = mux0 
	VFilter = Vx0 
	T = length(tmLab)
	A = aVec.i
	Li = L[L$CensusTractRegionID == tract.i,]
	
	muPredictList = VPredictList = rep(0, T)
	muFilterList = VFilterList =  rep(0, T)
	
	for (tt in 1:T){
		Lt = Li$count[Li$MonthTransaction==tmLab[tt]]
		
		muPredict = A*muFilter + lambda.ik*eta.k[tt]
		VPredict  = A*A * VFilter + sigma0sq
		
		if (Lt==0){
			muFilter = muPredict
			VFilter  = VPredict
		} else {
			S = VPredict + RVec.i/Lt
			Kt = VPredict/S
			
			covEff.t = covEffMonth[[tt]]
			ytHat = muPredict + covEff.t$covEff[covEff.t$CensusTractRegionID==tract.i ]
			yt = yMonth[[tt]]
			residual = yt$y[yt$CensusTractRegionID==tract.i]-ytHat
			muFilter = muPredict + Kt *mean(residual)
			VFilter = VPredict*(1 - Kt)
		}
		
		muPredictList[tt] = muPredict
		VPredictList[tt] = VPredict
		muFilterList[tt] = muFilter
		VFilterList[tt] = VFilter
	}
	return(list(muPredictList=muPredictList, VPredictList=VPredictList,
				muFilterList=muFilterList, VFilterList=VFilterList ))
}



KalmanFilterLikelihood.old <- function(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, tractMatch, Ct, yt,  covEfft, tract.i=NULL ){
	n = length(CtIds.k)
	muFilter = pi0 = rep(0, n)
	VFilter = V0 = 0.1*diag(rep(1,n))
	T = length(tmLab)
	A = matrix(0, n,n)
	diag(A) = aVec.k
	logMargLike.tract.i =0
	
	for (i in 1:T){
		# assign time varying SSM parameters, Rt
		
		if (!all(is.na(tractMatch[[i]]))){
			if (length(tractMatch[[i]])==1){
				Rt.i = RVec.k[ tractMatch[[i]] ]
			} else{
				Rt.i = diag(RVec.k[ tractMatch[[i]] ])
			}
		}
	
		muPredict = aVec.k*muFilter
		VPredict = A%*% VFilter %*% A + Q.k
		
		if (all(is.na(Ct[[i]]))){
			muFilter = muPredict
			VFilter  = VPredict
		} else {
			S = Ct[[i]] %*% VPredict %*% t(Ct[[i]]) + Rt.i
			Kt = VPredict%*%t(Ct[[i]])%*%chol2inv(chol(S))
			ytHat = Ct[[i]]%*%muPredict + covEfft[[i]]
			residual = yt[[i]]-ytHat
			muFilter = muPredict + Kt%*% residual
			VFilter = VPredict - Kt%*%Ct[[i]] %*% VPredict
		}
	
		# log marginal likelihood
		if (!all(is.na(yt[[i]]))){
			if (!is.null(tract.i)) {
				id = match(tract.i, CtIds.k)
				idct = tractMatch[[i]]==id
				p = sum(idct)
				if (p>0) {
					r.i = residual[idct]
					S.i = S[idct,idct]
					log.det.S.i = ifelse(p==1, log(S.i), sum(log(svd(S.i)$d)))
					logMargLike.tract.i = logMargLike.tract.i - p/2*log(2*pi)-1/2*log.det.S.i-1/2*t(r.i) %*% chol2inv(chol(S.i))%*%r.i
				}
			}
		}
	}
	logMargLike = NULL
	return(list(logMargLike=logMargLike, logMargLike.tract.i=logMargLike.tract.i, muPredict=muPredict, 
		VPredict=VPredict, 	muFilter=muFilter,VFilter=VFilter))
}

KalmanFilterLikelihood <- function(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, 
									tractMatch, Ct, yt, covEfft, 
									tract.i=NULL, mux0, Vx0, getYhat=FALSE){
	condLikelihood=TRUE
	n = length(CtIds.k)
	muFilter  = rep(mux0, n)
	VFilter  = matrix(0, n, n)
	diag(VFilter) = Vx0
	T = length(tmLab)
	A = matrix(0, n,n)
	diag(A) = aVec.k
	logMargLike.tract.all = 0
	logMargLike.tract.i =0
	likelihood.t = 0
	#tmSeries = rep(NA, T)
	#tmLikes = rep(NA,T)
	yHatList = list()
	
	for (tt in 1:T){
		# assign time varying SSM parameters, Rt
		if (!all(is.na(tractMatch[[tt]]))){
			if (length(tractMatch[[tt]])==1){
				Rt.tt = RVec.k[ tractMatch[[tt]] ]
				Rinv = 1/Rt.tt
			} else{
				tem = RVec.k[ tractMatch[[tt]] ]
				Rt.tt = diag(tem)
				Rinv = diag(1/tem)
			}
		}
	
		muPredict = aVec.k*muFilter
		#VPredict = A%*% VFilter %*% A + Q.k
		VPredict = (aVec.k %o% aVec.k) * VFilter + Q.k
		
		if (all(is.na(Ct[[tt]]))){
			muFilter = muPredict
			VFilter  = VPredict
		} else {
			#S = Ct[[tt]] %*% VPredict %*% t(Ct[[tt]]) + Rt.tt
			#Kt = VPredict%*%t(Ct[[tt]])%*%chol2inv(chol(S))
			# Applied woodbury formula
			Ctem = Ct[[tt]]
			VPtem = VPredict
			B = t(Ctem)%*%Rinv
			VPtem.inv = chol2inv(chol(VPtem))
			inn = chol2inv(chol(VPtem.inv + B%*%Ctem))
			Sinv = Rinv - t(B)%*%inn %*% B
			Kt = VPtem%*%t(Ctem)%*% Sinv
			S = Ctem %*% VPtem %*% t(Ctem) + Rt.tt
			
			ytHat = Ctem%*%muPredict + covEfft[[tt]]
			residual = yt[[tt]]-ytHat
			muFilter = muPredict + Kt%*% residual
			VFilter = VPredict - Kt%*%Ctem %*% VPredict
			
			if (getYhat) yHatList[[tt]] = as.numeric(ytHat)	
		}
	
		# log conditional likelihood of y_i, given y_j in the same cluster k and z_i = k
		if (!all(is.na(yt[[tt]]))){
			if (!is.null(tract.i)) { # calculate the conditional likelihood of obs for tract i cond on other tracts obs in cluster k
				id = match(tract.i, CtIds.k)
				idct = tractMatch[[tt]]==id
				p = sum(idct)
				if (p>0) {
					r.i = residual[idct]
					S.i = as.matrix(S[idct,idct])
					nObs.j = length(idct)-p
					if (condLikelihood & (nObs.j>0) ){# if other tracts have observation.
						r.j = residual[!idct]
						
						if (FALSE){
							Rj.inv = Rinv[!idct, !idct]
							Bj = B[,!idct]
							if (!is.matrix(Bj)) Bj = matrix(Bj, n, 1)
							Cj = Ctem[!idct,]
							if (!is.matrix(Cj)) Cj = matrix(Cj, 1, n)
							innj = chol2inv(chol(VPtem.inv + Bj %*% Cj))
							Sj.inv = Rj.inv - t(Bj)%*%innj %*% Bj
							Sij =S[idct, !idct]
							if (!is.matrix(Sij)) Sij = matrix(Sij, p, nObs.j)
							SijSj.inv = Sij %*% Sj.inv
							condMu = SijSj.inv %*% r.j
							condV  = S.i - SijSj.inv %*% t(Sij) 
						}
						
						## new method to get cond mu and cond V
						W = Sinv
						Wii = W[idct, idct]
						condV  = chol2inv(chol(Wii))
						Wij = W[idct, !idct]
						if ( (p==1 ) || (nObs.j==1) ) Wij = matrix(Wij, p, nObs.j)
						condMu = - condV %*% ( Wij %*% r.j)

						
						if (FALSE){
							t0 = proc.time()
							for (jj in 1:100){
								Sj = as.matrix(S[!idct, !idct])
								Sj.inv = chol2inv(chol(Sj))
							}
							(proc.time()-t0)[3]
							
							t0 = proc.time()
							for (jj in 1:100){
								Rj.inv = Rinv[!idct, !idct]
								Bj = B[,!idct]
								Cj = Ctem[!idct,]
								innj = chol2inv(chol(VPtem.inv + Bj %*% Cj))
								Sj.inv = Rj.inv - t(Bj)%*%innj %*% Bj
							}
							(proc.time()-t0)[3]
							
						}
						
						if (FALSE){
							d = 1
							Stem = Sj[1:d, 1:d]
							t0 = proc.time()
							for (jj in 1:10000){
								inv1 = solve(Stem)
							}
							(proc.time()-t0)[3]
							
							library(Matrix)
							t0 = proc.time()
							for (jj in 1:10000){
								cM = chol(Stem)
								inv2 = chol2inv(cM)
							}
							(proc.time()-t0)[3]
							
						}
						
						if (p==1){
							log.det.condV = log(condV)
						} else {
							U = chol(condV)
							ldet.U = sum(log(diag(U)))
							log.det.condV = 2*ldet.U
						}							
						rr.i = r.i - condMu
						likelihood.t = - p/2*log(2*pi) -1/2*log.det.condV -1/2*t(rr.i) %*% Wii %*%rr.i
						logMargLike.tract.i = logMargLike.tract.i  + likelihood.t
					} else { # if other tracts have no observation.
						if (p==1){
							log.det.S.i = log(S.i)
						} else {
							U = chol(S.i)
							ldet.U = sum(log(diag(U)))
							log.det.S.i = 2*ldet.U
						}	
						S.i.inv = chol2inv(chol(S.i))
						likelihood.t = - p/2*log(2*pi)-1/2*log.det.S.i-1/2*t(r.i) %*% S.i.inv %*%r.i
						logMargLike.tract.i = logMargLike.tract.i  + likelihood.t
					}
					#tmSeries[tt] = likelihood.t
					#tmLikes[tt] = logMargLike.tract.i
				}
			} else { # compute the jointly likelihood for all tracts in cluster k
				p = length(tractMatch[[tt]])
				if (p==1){
					log.det.S = log(S)
				} else {
					U = chol(S)
					ldet.U = sum(log(diag(U)))
					log.det.S = 2*ldet.U
				}
				r = residual
				likelihood.t = - p/2*log(2*pi)-1/2*log.det.S-1/2*t(r) %*% Sinv %*%r
				logMargLike.tract.all = logMargLike.tract.all  + likelihood.t
			}
		}
	}
	
	return(list(logMargLike.tract.i=as.numeric(logMargLike.tract.i),
				logMargLike.tract.all=as.numeric(logMargLike.tract.all),
				yHatList = yHatList))
}

# my solver for inversing a covariance matrix, positive, sysmetric matrix
varInverseSolve <-function(A){
	cM = chol(Stem)
	inv2 = chol2inv(cM)
} 

if (F){
	Lt = outSuff$Lt
	yDehedonict.i = outSuff$yDehedonict.i
	yDehedonicMean = outSuff$yDehedonicMean
	CtBar = outSuff$CtBar
	
	tractMatch = outCluster$tractMatch
	Ct = outCluster$Ct
	yt = outCluster$yt
	covEfft = outCluster$covEfft
	tract.i = CtIds[i] 
}

KalmanFilterLikelihood.suffStat <- function(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, 
											Lt, yDehedonict.i, yDehedonicMean, CtBar, 
											tractMatch=NULL, Ct=NULL, yt=NULL, covEfft=NULL,
											tract.i=NULL, mux0, Vx0 ){
	n = length(CtIds.k)
	muFilter  = rep(mux0, n)
	VFilter  = matrix(0, n, n)
	VFilter[cbind(1:n, 1:n)] = Vx0
	T = length(tmLab)
	logMargLike.tract.i =0
	logMargLike.tract.all =0
	#tmSeries = rep(NA, T)
	#tmLikes = rep(NA, T)
	
	for (tt in 1:T){
	
		L.kt = Lt[[tt]]
		idctCtObs.t = L.kt>0
		nCtObs = sum(idctCtObs.t )		
		
		muPredict = aVec.k*muFilter
		#VPredict = A%*% VFilter %*% A + Q.k
		VPredict = (aVec.k %o% aVec.k) * VFilter + Q.k
		
		if (nCtObs>0){ # cluster k has at least one obs at time t
			tem = RVec.k[idctCtObs.t]/L.kt[idctCtObs.t]
			RBar.tt = matrix(0, nCtObs, nCtObs)
			RBar.tt[cbind(1:nCtObs, 1:nCtObs)]= tem
		
			if (nCtObs==n){
				S.suff =  VPredict + RBar.tt
				Kt = VPredict %*% chol2inv(chol(S.suff))

				muFilter = muPredict + Kt %*% (yDehedonicMean[[tt]]-muPredict[idctCtObs.t])
				VFilter = VPredict- Kt%*% VPredict
			} else {
				Ctn = CtBar[[tt]]
				CV = Ctn %*% VPredict 
				S.suff =  CV %*% t(Ctn) + RBar.tt
				Kt = t(CV) %*% chol2inv(chol(S.suff))

				muFilter = muPredict + Kt %*% (yDehedonicMean[[tt]]-muPredict[idctCtObs.t])
				VFilter = VPredict - Kt %*% CV
			}
		} else {
			muFilter = muPredict
			VFilter  = VPredict
		}

		if (nCtObs>0){
			if (!is.null(tract.i)){
				id = match(tract.i, CtIds.k)
				p = L.kt[id]
				if (p>0) { # if tract.i has obs at time t, do likelihood calculation
					yDehedonic = yDehedonict.i[[tt]]
					r.i = yDehedonic - muPredict[id]
					v.i = VPredict[id,id]
					R.i = RVec.k[id]
					S.i = matrix(v.i, p, p )  # marginal covariance of observations in tract i
					S.i[cbind(1:p, 1:p)] = v.i + R.i
					nObs.j = sum(L.kt) - p
					if (nObs.j >0){ # if other tracts have observation.
						if (FALSE){ # conditional likelhood not using suff stats
							Ctem = Ct[[tt]]
							VPtem = VPredict
							ytHat = Ctem%*%muPredict + covEfft[[tt]]
							residual = yt[[tt]]-ytHat
							idct = tractMatch[[tt]]==id
							S = Ctem %*% VPtem %*% t(Ctem) + Rt.tt
							
							r.j = residual[!idct]
							S.j = as.matrix(S[!idct, !idct])
							S.ij =S[idct, !idct]
							if (!is.matrix(S.ij)) S.ij = matrix(S.ij, p, nObs.j)
							Sj.inv = chol2inv(chol(S.j))
							SijSj.inv = S.ij %*% Sj.inv
							condMu = SijSj.inv %*% r.j
							condV  = S.i - SijSj.inv %*% t(S.ij) 
						}
						
						### Conditional likelihood using suff stats
						id.suff = match(tract.i, CtIds.k[idctCtObs.t])
						Sjj.suff.inv = chol2inv(chol( S.suff[-id.suff, -id.suff] ))
						Sij.suff = S.suff[id.suff, -id.suff]
						condMu = Sij.suff %*% Sjj.suff.inv %*% (yDehedonicMean[[tt]]-muPredict[idctCtObs.t])[-id.suff]
						scalar = Sij.suff %*% Sjj.suff.inv %*% Sij.suff
						condV = S.i - as.numeric(scalar)
						
						rr.i = r.i - condMu
						#condV.inv = chol2inv(chol(condV))
						scalarC = (v.i-scalar)/(R.i+(v.i-scalar)*p)/R.i
						condV.inv = matrix(-as.numeric(scalarC), p, p)
						condV.inv[cbind(1:p, 1:p)] = 1/R.i - scalarC
						
						if (p==1){
							log.det.condV = log(condV)
						} else {
							U = chol(condV)
							ldet.U = sum(log(diag(U)))
							log.det.condV = 2*ldet.U
						}						

						likelihood.t = - p/2*log(2*pi) -1/2*log.det.condV -1/2*t(rr.i) %*% condV.inv %*%rr.i
						logMargLike.tract.i = logMargLike.tract.i  + likelihood.t
					} else { # # ifh other tracts have no observation, calculate log marginal likelihood for tract i
						
						#Si.inv = chol2inv(chol(S.i))
						scalarC = v.i/(R.i+v.i*p)/R.i
						Si.inv = matrix(-as.numeric(scalarC), p, p)
						Si.inv[cbind(1:p, 1:p)] = 1/R.i - scalarC

						if (p==1){
							log.det.S.i = log(S.i)
						} else {
							U = chol(S.i)
							ldet.U = sum(log(diag(U)))
							log.det.S.i = 2*ldet.U
						}	
						
						likelihood.t =  - p/2*log(2*pi)-1/2*log.det.S.i-1/2*t(r.i) %*% Si.inv %*%r.i
						logMargLike.tract.i = logMargLike.tract.i + likelihood.t
					}
					#tmSeries[tt] = likelihood.t
					#tmLikes[tt] = logMargLike.tract.i
				}
			} else { # if tract.i = NULL, calulcate the joinlty likelihood of all tracts in cluster k
				if (length(tractMatch[[tt]])==1){
					Rt.tt = RVec.k[ tractMatch[[tt]] ]
					Rinv = 1/Rt.tt
				} else{
					tem = RVec.k[ tractMatch[[tt]] ]
					Rt.tt = diag(tem)
					Rinv = diag(1/tem)
				}
				
				p = sum(L.kt)
				Ctem = Ct[[tt]]
				VPtem = VPredict
				B = t(Ctem)%*%Rinv
				VPtem.inv = chol2inv(chol(VPtem))
				inn = chol2inv(chol(VPtem.inv + B%*%Ctem))
				Sinv = Rinv - t(B)%*%inn %*% B
				S = Ctem %*% VPtem %*% t(Ctem) + Rt.tt
				ytHat = Ctem%*%muPredict + covEfft[[tt]]
				residual = yt[[tt]]-ytHat
			
				#log.det.S = ifelse(p==1, log(S), sum(log(svd(S)$d)))
				if (p==1){
					log.det.S = log(S)
				} else {
					U = chol(S)
					ldet.U = sum(log(diag(U)))
					log.det.S = 2*ldet.U
				}
				
				r = residual
				likelihood.t = - p/2*log(2*pi)-1/2*log.det.S-1/2*t(r) %*% Sinv %*%r
				
				## test suff stats (stopped here)
				rr = (yDehedonicMean[[tt]]-muPredict[idctCtObs.t])
				log.det.S.suff = log(det(S.suff))
				
				#tem = -nCtObs/2*log(2*pi) - 1/2*log.det.S.suff - 1/2*t(rr) %*% chol2inv(chol(S.suff)) %*% rr
				
				logMargLike.tract.all = logMargLike.tract.all  + likelihood.t
			}
		}
	}
	return(list(logMargLike.tract.i=as.numeric(logMargLike.tract.i) ,
				logMargLike.tract.all=as.numeric(logMargLike.tract.all) ))
}


KalmanFilterLikelihood.cond.eta <- function( tract.i , RVec.i, aVec.i, sigma0sq, eta.k, lambda.ik, covEffMonth, yMonth, tmLab, L, mux0, Vx0 ){
	muFilter = mux0 
	VFilter = Vx0 
	T = length(tmLab)
	A = aVec.i
	logMargLike.tract.i =0
	
	Li = L[L$CensusTractRegionID == tract.i,]
	
	for (tt in 1:T){
		Lt = Li$count[Li$MonthTransaction==tmLab[tt]]
		
		muPredict = A*muFilter + lambda.ik*eta.k[tt]
		VPredict  = A*A * VFilter + sigma0sq
		
		if (Lt==0){
			muFilter = muPredict
			VFilter  = VPredict
		} else {
			S = VPredict + RVec.i/Lt
			Kt = VPredict/S
			
			covEff.t = covEffMonth[[tt]]
			ytHat = muPredict + covEff.t$covEff[covEff.t$CensusTractRegionID==tract.i ]
			yt = yMonth[[tt]]
			residual = yt$y[yt$CensusTractRegionID==tract.i]-ytHat
			muFilter = muPredict + Kt *mean(residual)
			VFilter = VPredict*(1 - Kt)
			
			#log marginal likelihood
			p = Lt
			if (p>0) {
				r.i = residual
				tem = matrix(0,p,p)
				diag(tem) = RVec.i
				S.i = matrix(VPredict,p,p) + tem
				log.det.S.i = ifelse(p==1, log(S.i), sum(log(svd(S.i)$d)))
				logMargLike.tract.i = logMargLike.tract.i - p/2*log(2*pi)-1/2*log.det.S.i-1/2*t(r.i) %*% chol2inv(chol(S.i))%*%r.i
			}
		}
	}
	return(logMargLike.tract.i )
}

#####################
## Kalman smoother ##
#####################

KalmanSmoother <- function(muFilter, VFilter, muPredict, VPredict, tmLab, CtIds.k, aVec.k){ 
	T = length(tmLab)
	nk = length(CtIds.k)
	A = matrix(0, nk, nk)
	diag(A) = aVec.k
	muSmooth = list()  # dim_t = (nk,1)
	VSmooth  = list()  # dim_t =( nk,nk)
	#Jt = list()        # dim_t =(nk, nk)
	
	muSmooth[[T]] = muFilter[[T]]
	VSmooth[[T]] = VFilter[[T]]
	
	for (i in (T-1):1){
		Jt = VFilter[[i]] %*% t(A) %*% chol2inv(chol(VPredict[[i+1]]))
		muSmooth[[i]] = muFilter[[i]] + Jt%*%(muSmooth[[i+1]] - muPredict[[i+1]])
		VSmooth[[i]] = VFilter[[i]] + Jt %*%(VSmooth[[i+1]] - VPredict[[i+1]]) %*% t(Jt)
	}
	return(list(muSmooth=muSmooth, VSmooth=VSmooth))
}

########################################
## Backward Kalman Information Filter ##
########################################

BackKalmanInforFilter <- function(CtIds.k, tmLab, RVec.k, aVec.k, Q.k,  tractMatch, Ct, yt, covEfft){
	T = length(tmLab)
	precisionFilter = list()
	thetaFilter = list()
	nk = length(CtIds.k)
	A = matrix(0, nk,nk)
	diag(A) = aVec.k
	
	# Initialize filter at time T
	if (all(is.na(tractMatch[[T]]))){
		invV = matrix(0, nk, nk)
	#	diag(invV) = 1/(0.1*scaleFactor) 
		diag(invV) = 1/(0.1*scaleFactor*scaleFactor)  # should it be scaleFactor squared???
		#diag(invV) = 1/0.1
		precisionFilter[[T]] = invV
		thetaFilter[[T]] = matrix(0, nk,1)
	} else{
		if (length(tractMatch[[T]])==1){
			R = RVec.k[ tractMatch[[T]] ]
		} else{
			R = diag(RVec.k[ tractMatch[[T]] ])
		}
		C = Ct[[T]]
		precisionFilter[[T]] = t(C) %*% chol2inv(chol(R)) %*% C
		thetaFilter[[T]] = t(C) %*% chol2inv(chol(R)) %*% (yt[[T]]-covEfft[[T]])
	}
	
	# working backward
	for (tt in (T-1):1){
		# backward gain
		JtPlus1 = precisionFilter[[tt+1]] %*% chol2inv(chol(precisionFilter[[tt+1]]+chol2inv(chol(Q.k))))
		LtPlus1 = diag(nk) - JtPlus1 
		# predict backward
		precisionPredict = t(A) %*% (LtPlus1 %*% precisionFilter[[tt+1]] %*% t(LtPlus1) + 
							JtPlus1 %*% chol2inv(chol(Q.k)) %*% t(JtPlus1)) %*% A
		thetaPredict = t(A) %*% LtPlus1 %*% thetaFilter[[tt+1]]
		# update to filtered values
		if (all(is.na(tractMatch[[tt]]))){
			precisionFilter[[tt]] = precisionPredict 
			thetaFilter[[tt]] = thetaPredict
		} else{
			if (length(tractMatch[[tt]])==1){
				R = RVec.k[ tractMatch[[tt]] ]
			} else{
				R = diag(RVec.k[ tractMatch[[tt]] ])
			}
			C = Ct[[tt]]
			precisionFilter[[tt]] = precisionPredict + t(C) %*% chol2inv(chol(R)) %*% C
			thetaFilter[[tt]] = thetaPredict + t(C) %*% chol2inv(chol(R)) %*% (yt[[tt]]-covEfft[[tt]])
		}	
	}
	tt=0
	JtPlus1 = precisionFilter[[tt+1]] %*% chol2inv(chol(precisionFilter[[tt+1]]+chol2inv(chol(Q.k))))
	LtPlus1 = diag(nk) - JtPlus1 
	precisionPredict = t(A) %*% (LtPlus1 %*% precisionFilter[[tt+1]] %*% t(LtPlus1) + 
							JtPlus1 %*% chol2inv(chol(Q.k)) %*% t(JtPlus1)) %*% A
	thetaPredict = t(A) %*% LtPlus1 %*% thetaFilter[[tt+1]]
	precisionFilter0 = precisionPredict
	thetaFilter0 = thetaPredict
	
	return(list(precisionFilter0=precisionFilter0, thetaFilter0=thetaFilter0, 
			precisionFilter=precisionFilter, thetaFilter=thetaFilter))
}

################################################
## Forward sampler of x after backward filter ##
################################################

if (F){
	precisionFilter0 = outBackFilter$precisionFilter0
	thetaFilter0 = outBackFilter$thetaFilter0
	precisionFilter = outBackFilter$precisionFilter
	thetaFilter = outBackFilter$thetaFilter
}

ForwardSamplerX <- function(precisionFilter0, thetaFilter0, precisionFilter, thetaFilter, Q.k, tmLab, CtIds.k, aVec.k, x0){
	T = length(tmLab)
	nk = length(CtIds.k)
	A = matrix(0, nk,nk)
	diag(A) = aVec.k
	
	if (FALSE){
		V = chol2inv(chol(precisionFilter0))
		mu = V%*%thetaFilter0
		
		if (!isSymmetric(V)){
			# for (d in c(7:2)){
				# V = chol2inv(chol(precisionFilter0))
				# V = round(V, d)
				# if (isSymmetric(V)) break
			# }
			V = (V+t(V))/2
		}
		if (any(abs(mu)>1)) {
			idct = which(abs(mu)>1)
			mu[idct] = 0
			V[idct,idct] = 0.1
		}
		if (!isSymmetric(V)){
			x0 = matrix(0, nrow=nk, ncol=1)
		} else if (det(V)<0) {
			x0 = matrix(0, nrow=nk, ncol=1)
		} else {		
			x0 = matrix( rmvnorm(1, mu, V), ncol=1 )
		}
	}
	#x0 = matrix(0, nrow=nk, ncol=1)
	x = matrix(0, nk, T)
	rownames(x) = CtIds.k
	
	tt =1
	xtMinus1=x0
	V = chol2inv(chol( chol2inv(chol(Q.k))+ precisionFilter[[tt]] ))
	mu= V%*%(chol2inv(chol(Q.k))%*%A%*%xtMinus1 + thetaFilter[[tt]])
	x[,tt] = rmvnorm(1, mu, V)
	for (tt in 2:T){
		xtMinus1=x[,tt-1]
		V = chol2inv(chol( chol2inv(chol(Q.k))+ precisionFilter[[tt]] ))
		mu= V%*%(chol2inv(chol(Q.k))%*%A%*%xtMinus1 + thetaFilter[[tt]])
		x[,tt] = rmvnorm(1, mu, V)
	}
	return(x)
}

################################################
## Backward sampler of x after forward filter ##
################################################
if (F){
	muPredictList = outFilter$muPredictList
	VPredictList = outFilter$VPredictList
	muFilterList = outFilter$muFilterList
	VFilterList = outFilter$VFilterList
}

BackwardSamplerX <- function( muPredictList,VPredictList,muFilterList,VFilterList,  tmLab, CtIds.k, aVec.k, mux0, Vx0){
	T = length(tmLab)
	nk = length(CtIds.k)
	A = matrix(0, nk,nk)
	diag(A) = aVec.k
	x = matrix(0, nk, T)
	rownames(x) = CtIds.k
	
	V = VFilterList[[T]]
	V = (V+t(V))/2
	
	x[,T] = rmvnorm(1, muFilterList[[T]], V)
	
	for (tt in (T-1):1){
#	for (tt in (T-1):169){

		xtPlus1 = x[,tt+1]
		VFilter = VFilterList[[tt]]
		VPredict = VPredictList[[tt+1]]
		
		Jt = VFilter %*% A %*% chol2inv(chol(VPredict))
		mu = muFilterList[[tt]] + Jt %*%(xtPlus1-muPredictList[[tt+1]])
		V  = VFilter - Jt %*% VPredict %*% t(Jt)
		V = (V+t(V))/2
		tempScale = mean(V)
		if (det(V/tempScale)<0){
			Vscaled = V/tempScale + diag(rep(1, nrow(V)))
			x[,tt] = rmvnorm(1, mu, Vscaled*tempScale)
		} else {
			x[,tt] = rmvnorm(1, mu, V)
		}
	}
	
	#tt=0
	xtPlus1 = x[,1]
	VFilter = matrix(0, nk , nk)
	diag(VFilter) = Vx0
	muFilter = rep(mux0, nk)
	VPredict = VPredictList[[1]]
	
	Jt = VFilter %*% A %*% chol2inv(chol(VPredict))
	mu = muFilter + Jt %*%(xtPlus1-muPredictList[[1]])
	V  = VFilter - Jt %*% VPredict %*% t(Jt)
	V = (V+t(V))/2
	tempScale = mean(V)
	if (det(V/tempScale)<0){
		Vscaled = V/tempScale + diag(rep(1, nrow(V)))
		x0 = rmvnorm(1, mu, Vscaled*tempScale)
	} else {
		x0 = rmvnorm(1, mu, V)
	}

	x0 = as.numeric(x0)
	names(x0) = CtIds.k
	
	return(list(x=x, x0=x0))
}


BackwardSamplerX.cond.eta <- function( muPredictList,VPredictList,muFilterList,VFilterList,  tmLab, CtIds.i, aVec.i, mux0, Vx0){
	T = length(tmLab)
	A = aVec.i
	x = rep(0, T)
	
	# at tt = T
	V = VFilterList[T]
	x[T] = rnorm(1, muFilterList[T], sd=sqrt(V))
	
	# sequentially backward for tt=T-1, ..., 1
	for (tt in (T-1):1){
		xtPlus1 = x[tt+1]
		VFilter = VFilterList[tt]
		VPredict = VPredictList[tt+1]
		
		Jt = A*VFilter/VPredict
		mu = muFilterList[tt] + Jt*(xtPlus1 - muPredictList[tt+1])
		V  = VFilter - Jt * Jt * VPredict
		x[tt] = rnorm(1, mu, sd=sqrt(V))
	}
	
	# at tt=0
	xtPlus1 = x[1]
	VFilter = Vx0
	muFilter = mux0
	VPredict = VPredictList[1]
	
	Jt = A*VFilter/VPredict
	mu = muFilter + Jt*(xtPlus1-muPredictList[1])
	V  = VFilter - Jt * Jt * VPredict 
	x0 = rnorm(1, mu, sd=sqrt(V))
	names(x0) = CtIds.i
	
	return(list(x=x, x0=x0))
}


#1000000000 
if (FALSE){
	for (ii in 2:10){
		for (jj in 1:(ii-1)){
			if ( abs(V[ii,jj]-V[jj,ii])> 1* .Machine$double.eps){
				cat(ii, jj, "\n")
			}
		}
	}
}
#########################
## Gibbs Sampler Steps ##
#########################

# 1 a. Sample membership z_i, for i= 1:n
sample.z <- function(zVec, alpha, aVec, LambdaMat, sigma0, RVec, covEffMonth,   
					muLambda, sigmaLambda, tmLab, CtIds, L, yMonth, is.debug=FALSE, mux0, Vx0){
	n = length(zVec)
	
	if (n > 1){
	#set.seed(1)
	permu = sample(1:n)
	for (per in 1:n){
		i = permu[per] 
		K = length(table(zVec[-i]))
		Kwithi = length(table(zVec))
		kLabels  = as.numeric(names(table(zVec[-i])))
		part1 = c(table(zVec[-i]),alpha )/(n-1+alpha)
		part2 = rep(0,K+1)
		names(part2) = kLabels
		for (k in kLabels){
			if (k==zVec[i]) {
				CtIds.k = sort(CtIds[zVec==k])
			} else {
				CtIds.k = sort(c(CtIds[zVec==k], CtIds[i]))
			}
			sigma0diag = matrix(0, length(CtIds.k), length(CtIds.k))
			diag(sigma0diag) = sigma0^2
			Q.k = LambdaMat[as.character(CtIds.k),as.character(k)] %o% LambdaMat[as.character(CtIds.k),as.character(k)]+
					sigma0diag
			aVec.k = aVec[as.character(CtIds.k)]
			RVec.k = RVec[as.character(CtIds.k)]
			
			#startTime = proc.time()
			#for (l in 1:1000){
			out <- get.Ct.yt(L, yMonth, covEffMonth, CtIds.k, tmLab)
			# The following likelihood is evaluated in Kalman filter with Woodbury formula
			outFilter = KalmanFilterLikelihood(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, out$tractMatch,out$Ct,out$yt,out$covEfft, tract.i = CtIds[i],mux0, Vx0 )
			#}
			#totalTime = proc.time() - startTime
			
			#outSuff <- get.Ct.yt.suffStat(L, covEffMonth, yDehedonicMonthMean, CtIds.k, tmLab,  tract.i = CtIds[i])
			#outFilter = KalmanFilterLikelihood.suffStat(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, out$Lt, out$yDehedonict.i, out$yDehedonicMean,out$CtBar, tract.i = CtIds[i]  )

			
			part2[as.character(k)] = outFilter$logMargLike.tract.i
		}
		# for (K+1)-th cluster
		CtIds.k = CtIds[i]
		
		if (is.debug){
			lambda.Kplus1 = LambdaMat[as.character(CtIds[i]),1]
		} else {
			if ((Kwithi-K)==1){
				lambda.Kplus1 = LambdaMat[as.character(CtIds[i]), as.character(zVec[i])]
			} else {
				lambda.Kplus1 = rnorm(1, muLambda, sigmaLambda)
			}
		}
		Q.k = lambda.Kplus1^2+sigma0^2
		aVec.k = aVec[as.character(CtIds.k)]
		RVec.k = RVec[as.character(CtIds.k)]

		out <- get.Ct.yt(L, yMonth, covEffMonth, CtIds.k, tmLab)
		outFilter = KalmanFilterLikelihood(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, out$tractMatch,out$Ct,out$yt,out$covEfft, tract.i = CtIds[i], mux0, Vx0  )
		#out <- get.Ct.yt.suffStat(L, covEffMonth, yDehedonicMonthMean, CtIds.k, tmLab,  tract.i = CtIds[i])
		#outFilter = KalmanFilterLikelihood.suffStat(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, out$Lt, out$yDehedonict.i, out$yDehedonicMean,out$CtBar, tract.i = CtIds[i]  )

		part2[K+1] = outFilter$logMargLike.tract.i
		part2 = part2-mean(part2)
		pVec = part1*exp(part2)
		pVec = pVec/sum(pVec)
		zi = sample(1:(K+1), size=1, prob=pVec)

		if (zi==K+1){
			if (is.debug) {
				LambdaMat = cbind(LambdaMat[,as.character(kLabels)], rnorm(nrow(LambdaMat),muLambda, sigmaLambda))
			} else {
				if ((Kwithi-K)==1){
					LambdaMat = cbind(LambdaMat[,as.character(kLabels)], LambdaMat[, as.character(zVec[i])])
				} else {
					LambdaMat = cbind(LambdaMat[,as.character(kLabels)], rnorm(nrow(LambdaMat),muLambda, sigmaLambda))
				}
			}
			LambdaMat[as.character(CtIds.k),ncol(LambdaMat)] = lambda.Kplus1
			new.label = max(kLabels)+1
			zVec[i] = new.label 
			colnames(LambdaMat) = c(kLabels,new.label  )
		} else {
			zVec[i] = kLabels[zi]
			LambdaMat = as.matrix( LambdaMat[,as.character(kLabels)] )
			colnames(LambdaMat) = as.character(kLabels)
		}
		#cat("tract", per, " K =", length(unique(zVec)), "\n")
	}
	}
	
	K = length(unique(zVec))
	names(zVec) = CtIds
	
	if (n==1){
		LambdaMat = as.matrix(LambdaMat[,as.character(zVec)])
		colnames(LambdaMat) = zVec
	}
	
	return(list(zVec=zVec, K=K, LambdaMat=LambdaMat))
}


# 1 a. Sample membership z_i, for i= 1:n
sample.z.suffStat <- function(zVec, alpha, aVec, LambdaMat, sigma0, RVec, covEffMonth, yDehedonicMonthMean, 
					muLambda, sigmaLambda, tmLab, CtIds,  L, yMonth, is.debug=FALSE, mux0, Vx0){
	n = length(zVec)
	
	if (n > 1){
	#set.seed(1)
	permu = sample(1:n)
	for (per in 1:n){
		i = permu[per] # i=6
		K = length(table(zVec[-i]))
		Kwithi = length(table(zVec))
		kLabels  = as.numeric(names(table(zVec[-i])))
		part1 = c(table(zVec[-i]),alpha )/(n-1+alpha)
		part2 = rep(0,K+1)
		names(part2) = kLabels
		for (k in kLabels){
			if (k==zVec[i]) {
				CtIds.k = sort(CtIds[zVec==k])
			} else {
				CtIds.k = sort(c(CtIds[zVec==k], CtIds[i]))
			}
			sigma0diag = matrix(0, length(CtIds.k), length(CtIds.k))
			diag(sigma0diag) = sigma0*sigma0
			lambdaTem = LambdaMat[as.character(CtIds.k),as.character(k)]
			Q.k =  lambdaTem %o% lambdaTem + sigma0diag
			aVec.k = aVec[as.character(CtIds.k)]
			RVec.k = RVec[as.character(CtIds.k)]
			
			#out = get.Ct.yt.suffStat(L, covEffMonth, yDehedonicMonthMean, CtIds.k, tmLab, tract.i=CtIds[i])			
			#outFilter = KalmanFilterLikelihood.suffStat(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, out$Lt, out$yDehedonict.i, out$yDehedonicMean,out$CtBar, tract.i = CtIds[i], mux0, Vx0  )
			
			#startTime = proc.time()
			#for (l in 1:1000){
			
			outSuff = get.Ct.yt.suffStat(L, covEffMonth, yDehedonicMonthMean, CtIds.k, tmLab, tract.i=CtIds[i])
			#outCluster = get.Ct.yt(L, yMonth, covEffMonth, CtIds.k, tmLab)
			outFilter = KalmanFilterLikelihood.suffStat(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, 
								outSuff$Lt, outSuff$yDehedonict.i, outSuff$yDehedonicMean, outSuff$CtBar, 
								NULL, NULL, NULL, NULL,
								tract.i=CtIds[i], mux0, Vx0 )
			#}
			#totalTime = proc.time() - startTime
			
			part2[as.character(k)] = outFilter$logMargLike.tract.i
		}
		# for (K+1)-th cluster
		CtIds.k = CtIds[i]
		
		if (is.debug){
			lambda.Kplus1 = LambdaMat[as.character(CtIds[i]),1]
		} else {
			if ((Kwithi-K)==1){
				lambda.Kplus1 = LambdaMat[as.character(CtIds[i]), as.character(zVec[i])]
			} else {
				lambda.Kplus1 = rnorm(1, muLambda, sigmaLambda)
			}
		}
		Q.k = lambda.Kplus1*lambda.Kplus1+sigma0*sigma0
		aVec.k = aVec[as.character(CtIds.k)]
		RVec.k = RVec[as.character(CtIds.k)]

		outSuff = get.Ct.yt.suffStat(L, covEffMonth, yDehedonicMonthMean, CtIds.k, tmLab, tract.i=CtIds[i])
		outFilter = KalmanFilterLikelihood.suffStat(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, 
								outSuff$Lt, outSuff$yDehedonict.i, outSuff$yDehedonicMean, outSuff$CtBar, 
								NULL, NULL, NULL, NULL,
								tract.i=CtIds[i], mux0, Vx0 )	
		part2[K+1] = outFilter$logMargLike.tract.i
		
		part2 = part2-max(part2)
		pVec = part1*exp(part2)
		pVec = pVec/sum(pVec)
		zi = sample(1:(K+1), size=1, prob=pVec)

		if (zi==K+1){
			if (is.debug) {
				LambdaMat = cbind(LambdaMat[,as.character(kLabels)], rnorm(nrow(LambdaMat),muLambda, sigmaLambda))
			} else {
				if ((Kwithi-K)==1){
					LambdaMat = cbind(LambdaMat[,as.character(kLabels)], LambdaMat[, as.character(zVec[i])])
				} else {
					LambdaMat = cbind(LambdaMat[,as.character(kLabels)], rnorm(nrow(LambdaMat),muLambda, sigmaLambda))
				}
			}
			LambdaMat[as.character(CtIds.k),ncol(LambdaMat)] = lambda.Kplus1
			new.label = max(kLabels)+1
			zVec[i] = new.label 
			colnames(LambdaMat) = c(kLabels,new.label  )
		} else {
			zVec[i] = kLabels[zi]
			LambdaMat = as.matrix( LambdaMat[,as.character(kLabels)] )
			colnames(LambdaMat) = as.character(kLabels)
		}
		#cat("tract", per, " K =", length(unique(zVec)), "\n")
	}
	}
	
	K = length(unique(zVec))
	names(zVec) = CtIds
	
	if (n==1){
		LambdaMat = as.matrix(LambdaMat[,as.character(zVec)])
		colnames(LambdaMat) = zVec
	}
	
	return(list(zVec=zVec, K=K, LambdaMat=LambdaMat))
}

# 1 a (parallel). Sample membership z_i, for i= 1:n

sample.z.parallelK <- function(zVec, alpha, aVec, LambdaMat, sigma0, RVec, covEffMonth,   
					muLambda, sigmaLambda, tmLab, CtIds, covNames, L, yMonth, is.debug=FALSE){
	n = length(zVec)
	permu = sample(1:n)
	nCl = length(unique(zVec))
	cl<-makeCluster(nCl) #change the 2 to your number of CPU cores
	registerDoSNOW(cl)

	for (per in 1:n){
		i = permu[per]
		K = length(table(zVec[-i]))
		Kwithi = length(table(zVec))
		kLabels  = as.numeric(names(table(zVec[-i])))
		part1 = c(table(zVec[-i]),alpha )/(n-1+alpha)
		part2 = rep(0,K+1)
		names(part2) = kLabels
		
		parallelOut <- foreach(j = 1:length(kLabels)) %dopar%{
		#for (j in 1:length(kLabels)){
			k = kLabels[j]
			if (k==zVec[i]) {
				CtIds.k = sort(CtIds[zVec==k])
			} else {
				CtIds.k = sort(c(CtIds[zVec==k], CtIds[i]))
			}
			sigma0diag = matrix(0, length(CtIds.k), length(CtIds.k))
			diag(sigma0diag) = sigma0^2
			Q.k = LambdaMat[as.character(CtIds.k),as.character(k)] %o% LambdaMat[as.character(CtIds.k),as.character(k)]+
					sigma0diag
			aVec.k = aVec[as.character(CtIds.k)]
			RVec.k = RVec[as.character(CtIds.k)]
			
			#out <- get.Ct.yt(L, yMonth, covEffMonth, CtIds.k, tmLab)
			CtIds.k = sort(CtIds.k)
			nk = length(CtIds.k)
			L.k = L[L$CensusTractRegionID %in% CtIds.k,]
			T = length(tmLab)
			Ct = list() # dim is varying by t, dim_t = (mt, nk), the indicator matrix from obs to tract at t
			tractMatch = list()
			yt=list()
			covEfft = list()
			
			for (tt in 1:T){
				L.k.tt = L.k[L.k$MonthTransaction==tmLab[tt],]
				mk=sum(L.k.tt$count)
				if (mk==0){
					Ct[[tt]] = NA
					tractMatch[[tt]] = NA
					yt[[tt]] = NA
					covEfft[[tt]] = NA
				} else {
					colMatch = rep(1:nk, L.k.tt$count)
					tem = matrix(0, mk, nk )
					for (j in 1:mk){
						tem[j,colMatch[j]] = 1
					}
					Ct[[tt]] = tem
					tractMatch[[tt]] = colMatch

					tem = yMonth[[tt]]
					idct = tem$CensusTractRegionID %in% CtIds.k		
					yt[[tt]] = tem$y[idct]

					tem = covEffMonth[[tt]]
					idct = tem$CensusTractRegionID %in% CtIds.k		
					covEfft[[tt]] = tem$covEff[idct]
				}
			}
			
			##
			#outFilter = KalmanFilterLikelihood(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, out$tractMatch,out$Ct,out$yt,out$covEfft, tract.i = CtIds[i]  )
			#outFilter$logMargLike.tract.i
			tract.i = CtIds[i]
			n = length(CtIds.k)
			muFilter = pi0 = rep(0, n)
			VFilter = V0 = 0.1*diag(rep(1,n))
			T = length(tmLab)
			A = matrix(0, n,n)
			diag(A) = aVec.k
			logMargLike.tract.i =0
			
			for (tt in 1:T){
				# assign time varying SSM parameters, Rt
				if (!all(is.na(tractMatch[[tt]]))){
					if (length(tractMatch[[tt]])==1){
						Rt.i = RVec.k[ tractMatch[[tt]] ]
					} else{
						Rt.i = diag(RVec.k[ tractMatch[[tt]] ])
					}
				}
			
				muPredict = aVec.k*muFilter
				VPredict = A%*% VFilter %*% A + Q.k
				
				if (all(is.na(Ct[[tt]]))){
					muFilter = muPredict
					VFilter  = VPredict
				} else {
					S = Ct[[tt]] %*% VPredict %*% t(Ct[[tt]]) + Rt.i
					Kt = VPredict%*%t(Ct[[tt]])%*%chol2inv(chol(S))
					ytHat = Ct[[tt]]%*%muPredict + covEfft[[tt]]
					residual = yt[[tt]]-ytHat
					muFilter = muPredict + Kt%*% residual
					VFilter = VPredict - Kt%*%Ct[[tt]] %*% VPredict
				}
			
				# log marginal likelihood
				if (!all(is.na(yt[[tt]]))){
					if (!is.null(tract.i)) {
						id = match(tract.i, CtIds.k)
						idct = tractMatch[[tt]]==id
						p = sum(idct)
						if (p>0) {
							r.i = residual[idct]
							S.i = S[idct,idct]
							log.det.S.i = ifelse(p==1, log(S.i), sum(log(svd(S.i)$d)))
							logMargLike.tract.i = logMargLike.tract.i - p/2*log(2*pi)-1/2*log.det.S.i-1/2*t(r.i) %*% chol2inv(chol(S.i))%*%r.i
						}
					}
				}
			}
			logMargLike.tract.i
		}
		part2[as.character(kLabels)] = unlist(parallelOut)

		# for (K+1)-th cluster
		CtIds.k = CtIds[i]
		
		if (is.debug){
			lambda.Kplus1 = LambdaMat[i,1]
		} else {
			if ((Kwithi-K)==1){
				lambda.Kplus1 = LambdaMat[i, as.character(zVec[i])]
			} else {
				lambda.Kplus1 = rnorm(1, muLambda, sigmaLambda)
			}
		}
		Q.k = lambda.Kplus1^2+sigma0^2
		aVec.k = aVec[as.character(CtIds.k)]
		RVec.k = RVec[as.character(CtIds.k)]

		out <- get.Ct.yt(L, yMonth, covEffMonth, CtIds.k, tmLab)
		outFilter = KalmanFilterLikelihood(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, out$tractMatch,out$Ct,out$yt,out$covEfft, tract.i = CtIds[i]  )
		part2[K+1] = outFilter$logMargLike.tract.i
		part2 = part2-mean(part2)
		pVec = part1*exp(part2)
		pVec = pVec/sum(pVec)
		zi = sample(1:(K+1), size=1, prob=pVec)

		if (zi==K+1){
			if (is.debug) {
				LambdaMat = cbind(LambdaMat[,as.character(kLabels)], rnorm(n,muLambda, sigmaLambda))
			} else {
				if ((Kwithi-K)==1){
					LambdaMat = cbind(LambdaMat[,as.character(kLabels)], LambdaMat[, as.character(zVec[i])])
				} else {
					LambdaMat = cbind(LambdaMat[,as.character(kLabels)], rnorm(n,muLambda, sigmaLambda))
				}
			}
			LambdaMat[as.character(CtIds.k),ncol(LambdaMat)] = lambda.Kplus1
			new.label = max(kLabels)+1
			zVec[i] = new.label 
			colnames(LambdaMat) = c(kLabels,new.label  )
		} else {
			zVec[i] = kLabels[zi]
			LambdaMat = as.matrix( LambdaMat[,as.character(kLabels)] )
			colnames(LambdaMat) = as.character(kLabels)
		}
		#cat("tract", per, " K =", length(unique(zVec)), "\n")
	}
	K = length(unique(zVec))
	names(zVec) = CtIds
	stopCluster(cl)
	
	return(list(zVec=zVec, K=K, LambdaMat=LambdaMat))
}

# 1 b. Sample membership z_i, for i= 1:n, conditional on eta_k, integrating out x_i
# The likelihood of each tract is not correlated within a cluster, given eta_k, though the z_i are coupled via prior of CRP
sample.z.cond.eta <- function(zVec, alpha, aVec, LambdaMat, sigma0, RVec, covEffMonth,   
					muLambda, sigmaLambda, tmLab, CtIds, covNames, L, yMonth, etaMat, mux0, Vx0, is.debug=FALSE){
	n = length(zVec)
	sigma0sq = sigma0*sigma0

	if (n > 1){
	permu = sample(1:n)
	for (per in 1:n){
		i = permu[per]
		CtIds.i = CtIds[i]
		aVec.i = aVec[as.character(CtIds.i)]
		RVec.i = RVec[as.character(CtIds.i)]
		
		K = length(table(zVec[-i]))
		Kwithi = length(table(zVec))
		kLabels  = as.numeric(names(table(zVec[-i])))
		part1 = c(table(zVec[-i]),alpha )/(n-1+alpha)
		part2 = rep(0,K+1)
		names(part2) = kLabels
		
		for (k in kLabels){
			lambda.ik = LambdaMat[as.character(CtIds.i), as.character(k)]
			eta.k = etaMat[, as.character(k)]
			logMargLike.tract.i = KalmanFilterLikelihood.cond.eta(tract.i = CtIds[i], RVec.i, aVec.i, sigma0sq, 
						eta.k, lambda.ik, covEffMonth, yMonth, tmLab, L, mux0, Vx0 )
			part2[as.character(k)] = logMargLike.tract.i
		}
		
		# for (K+1)-th cluster # stop here
		
		if (is.debug){
			lambda.Kplus1 = LambdaMat[as.character(CtIds[i]),1]
		} else {
			if ((Kwithi-K)==1){
				lambda.Kplus1 = LambdaMat[as.character(CtIds[i]), as.character(zVec[i])]
			} else {
				lambda.Kplus1 = rnorm(1, muLambda, sigmaLambda)
			}
		}
		
		if (TRUE){ # the likelihood of new cluster (marginal, integrate out eta)
			Q.k = lambda.Kplus1^2+sigma0sq
			out <- get.Ct.yt(L, yMonth, covEffMonth, CtIds.i, tmLab)
			outFilter = KalmanFilterLikelihood(CtIds.i, tmLab, RVec.i, aVec.i, Q.k, 
							out$tractMatch,out$Ct,out$yt,out$covEfft, tract.i = CtIds[i], mux0, Vx0  )
			part2[K+1] = outFilter$logMargLike.tract.i
			eta.k = rnorm(T, 0, 1) # Question: propose the series for the K+1 cluster, (should we evalute likelihood of K+1 by int out eta?)
									# but it may lead to small likelihood that the new cluster is created. Check eta sampler
	
		}
		
		if (FALSE){ # the likelihood of new cluster (the time series of eta are sampled from prior, but this likelihood may be very small)
			lambda.ik = lambda.Kplus1
			eta.k = rnorm(T, 0, 1) # Question: propose the series for the K+1 cluster, (should we evalute likelihood of K+1 by int out eta?)
									# but it may lead to small likelihood that the new cluster is created. Check eta sampler
			logMargLike.tract.i = KalmanFilterLikelihood.cond.eta(tract.i = CtIds[i], RVec.i, aVec.i, sigma0sq, 
							eta.k, lambda.ik, covEffMonth, yMonth, tmLab, L )
			part2[K+1] = logMargLike.tract.i
		}
		part2 = part2-mean(part2)
		pVec = part1*exp(part2)
		pVec = pVec/sum(pVec)
		zi = sample(1:(K+1), size=1, prob=pVec)

		if (zi==K+1){
			if (is.debug) {
				LambdaMat = cbind(LambdaMat[,as.character(kLabels)], rnorm(nrow(LambdaMat),muLambda, sigmaLambda))
			} else {
				if ((Kwithi-K)==1){
					LambdaMat = cbind(LambdaMat[,as.character(kLabels)], LambdaMat[, as.character(zVec[i])])
					etaMat = cbind(etaMat[,as.character(kLabels)], etaMat[, as.character(zVec[i])])
				} else {
					LambdaMat = cbind(LambdaMat[,as.character(kLabels)], rnorm(nrow(LambdaMat),muLambda, sigmaLambda))
					etaMat = cbind(etaMat[,as.character(kLabels)], eta.k)
				}
			}
			LambdaMat[as.character(CtIds.i),ncol(LambdaMat)] = lambda.Kplus1
			new.label = max(kLabels)+1
			zVec[i] = new.label 
			colnames(LambdaMat) = c(kLabels,new.label  )
			colnames(etaMat) =  c(kLabels,new.label  )
		} else {
			zVec[i] = kLabels[zi]
			LambdaMat = as.matrix( LambdaMat[,as.character(kLabels)] )
			colnames(LambdaMat) = as.character(kLabels)
			etaMat = as.matrix( etaMat[,as.character(kLabels)] )
			colnames(etaMat) = as.character(kLabels)
		}
		# cat("tract", per, " K =", length(unique(zVec)), "\n")
	}
	}
	
	K = length(unique(zVec))
	names(zVec) = CtIds
	
	if (n==1){
		LambdaMat = as.matrix(LambdaMat[,as.character(zVec)])
		colnames(LambdaMat) = zVec
		etaMat = as.matrix(etaMat[,as.character(zVec)])
		colnames(etaMat) = zVec
	}
	return(list(zVec=zVec, K=K, LambdaMat=LambdaMat, etaMat=etaMat))
}

# 1 c. Sample membership z_i by sampling the measure (mixture weights) w, first
# then conditional on w, z_i are i.i.d multinomial(w). So z_i are sampled in parallel over tracts 
# with expensive marginal likelihood that integrate out x and eta

sample.mixtureWeights <- function(zVec, KMax, alpha){
	alphaw = data.frame(labels = 1:KMax, prior=rep(alpha/KMax, KMax))
	nks = table(zVec)
	countk = data.frame(labels = names(nks), count = as.numeric(nks))
	tab = merge(alphaw, countk, by="labels", all.x=TRUE, sort=TRUE)
	tab$count[is.na(tab$count)] =0
	w = rdirichlet(n=1, alpha=tab$prior+tab$count)
	w = as.numeric(w)
	names(w) = tab$labels
	return(w)
}

sample.zi.cond.w <- function(iVec, w, KMax,zVec, aVec, LambdaMat, sigma0, RVec, covEffMonth,   
					 tmLab, CtIds, L, yMonth, mux0, Vx0){
	kMaxLabels = as.numeric(names(w))
	
	#parallelOut <- foreach(i = 1:n) %dopar%{
	res = NULL
	for (i in iVec){
		part1 = w
		part2 = rep(0,KMax)
		names(part2) = kMaxLabels
		CtIds.i = CtIds[i]
		aVec.i = aVec[as.character(CtIds.i)]
		RVec.i = RVec[as.character(CtIds.i)]
		out.i <- get.Ct.yt(L, yMonth, covEffMonth, CtIds.i, tmLab)
		for (k in kMaxLabels ){
			if (k==zVec[i]) {
				CtIds.k = sort(CtIds[zVec==k])
			} else {
				CtIds.k = sort(c(CtIds[zVec==k], CtIds[i]))
			}
			sigma0diag = matrix(0, length(CtIds.k), length(CtIds.k))
			diag(sigma0diag) = sigma0^2
			Q.k = LambdaMat[as.character(CtIds.k),as.character(k)] %o% LambdaMat[as.character(CtIds.k),as.character(k)]+
					sigma0diag
			aVec.k = aVec[as.character(CtIds.k)]
			RVec.k = RVec[as.character(CtIds.k)]
			
			out <- get.Ct.yt(L, yMonth, covEffMonth, CtIds.k, tmLab)
			outFilter = KalmanFilterLikelihood(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, out$tractMatch,out$Ct,out$yt,out$covEfft, tract.i = CtIds[i],mux0, Vx0 )
	
			part2[as.character(k)] = outFilter$logMargLike.tract.i
		}
		part2 = part2-mean(part2)
		pVec = part1*exp(part2)
		if (any(pVec==Inf)) {
			pTem = rep(0, KMax)
			pTem[which(pVec==Inf)] = 1
			pVec = pTem
		} else {
			pVec = pVec/sum(pVec)
		}
		zi = sample(kMaxLabels, size=1, prob=pVec)
		res = c(res,zi)
	}	
	return(res)
}

sample.zi.cond.w.suffStat <- function(iVec, w, KMax,zVec, aVec, LambdaMat, sigma0, RVec, covEffMonth, yDehedonicMonthMean,  
					 tmLab, CtIds, L, mux0, Vx0){
	kMaxLabels = as.numeric(names(w))
	
	res = NULL
	for (i in iVec){
		part1 = w
		part2 = rep(0,KMax)
		names(part2) = kMaxLabels
		CtIds.i = CtIds[i]
		aVec.i = aVec[as.character(CtIds.i)]
		RVec.i = RVec[as.character(CtIds.i)]
		out.i <- get.Ct.yt(L, yMonth, covEffMonth, CtIds.i, tmLab)
		for (k in kMaxLabels ){
			if (k==zVec[i]) {
				CtIds.k = sort(CtIds[zVec==k])
			} else {
				CtIds.k = sort(c(CtIds[zVec==k], CtIds[i]))
			}
			sigma0diag = matrix(0, length(CtIds.k), length(CtIds.k))
			diag(sigma0diag) = sigma0^2
			Q.k = LambdaMat[as.character(CtIds.k),as.character(k)] %o% LambdaMat[as.character(CtIds.k),as.character(k)]+
					sigma0diag
			aVec.k = aVec[as.character(CtIds.k)]
			RVec.k = RVec[as.character(CtIds.k)]

			outSuff = get.Ct.yt.suffStat(L, covEffMonth, yDehedonicMonthMean, CtIds.k, tmLab, tract.i=CtIds[i])
			outFilter = KalmanFilterLikelihood.suffStat(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, 
								outSuff$Lt, outSuff$yDehedonict.i, outSuff$yDehedonicMean, outSuff$CtBar, 
								NULL, NULL, NULL, NULL,
								tract.i=CtIds[i], mux0, Vx0 )	
			
			part2[as.character(k)] = outFilter$logMargLike.tract.i
		}
		#part2 = part2-mean(part2)
		part2 = part2 - max(part2)
		pVec = part1*exp(part2)
		pVec = pVec/sum(pVec)
		zi = sample(kMaxLabels, size=1, prob=pVec)
		res = c(res,zi)
	}	
	return(res)
}


# 2. Sample the concentration parameter alpha, a scalar (by the method in the crime paper)
# Method of slice sampling in the nonparametric variable clustering paper doesn't work well

sample.alpha <- function(alpha, n, alpha.alpha0, beta.alpha0, K){
	kap = rbeta(1, alpha+1, n )
	above = alpha.alpha0 + K -1
	below = n*(beta.alpha0-log(kap))
	p = above/(above+below)
	u = runif(1)
	if (u<p){
		alpha = rgamma(1, shape=alpha.alpha0+K, rate=beta.alpha0-log(kap))
	} else {
		alpha = rgamma(1, shape=alpha.alpha0+K-1, rate=beta.alpha0-log(kap))
	}
	return(alpha)
}

# 3. Sample the tract process x
sample.x <- function(zVec, LambdaMat, aVec, RVec, sigma0 ,L, yMonth, covEffMonth, tmLab, CtIds, covNames, mux0, Vx0){
	K = length(table(zVec))
	n = length(CtIds)
	kLabels  = as.numeric(names(table(zVec)))
	T = length(tmLab)
	xMat = matrix(0, n, T)
	rownames(xMat) = CtIds
	x0Vec= rep(0, n)
	muTVec = rep(0,n)
	names(x0Vec) = CtIds
	names(muTVec) = CtIds

	for (k in kLabels){
		CtIds.k = sort(CtIds[zVec==k])
		nk = length(CtIds.k)
		sigma0diag = matrix(0, length(CtIds.k), length(CtIds.k))
		diag(sigma0diag) = sigma0^2
		Q.k = LambdaMat[as.character(CtIds.k),as.character(k)] %o% LambdaMat[as.character(CtIds.k),as.character(k)]+sigma0diag
		aVec.k = aVec[as.character(CtIds.k)]
		RVec.k = RVec[as.character(CtIds.k)]
		
		if (FALSE){
		out <- get.Ct.yt(L, yMonth, covEffMonth, CtIds.k, tmLab)
		outBackFilter <- BackKalmanInforFilter(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, 
												out$tractMatch,out$Ct,out$yt,out$covEfft)
		}										
		# Way 1: sample x0 directly from prior
		#x0 = matrix(rnorm(nk, mux0, sqrt(Vx0)), ncol=1)
		
		# Way 2: sample x0 from posterior = prior*filteredLikelihood(x0) , # note rmvnorm and rnorm
		if (FALSE){
			Sigma0inv = matrix(0, nk, nk)
			diag(Sigma0inv) = 1/Vx0
			post.vx0 = chol2inv(chol(Sigma0inv+outBackFilter$precisionFilter0))
			post.mux0 = post.vx0 %*% (Sigma0inv %*%rep(mux0,nk) + outBackFilter$precisionFilter0 %*%outBackFilter$thetaFilter0)
			if (nk==1) {
				x0 = matrix( rnorm(1, post.mux0, sd=sqrt(post.vx0)), ncol=1 )
			} else {
				x0 = matrix( rmvnorm(1, post.mux0, post.vx0), ncol=1 )
			}
		}
		
		# Way 3: fix x0 at 0. Most stable solution. but too much confidence put on x=0 for t=0
		if (FALSE){ 
			x0 = matrix(0, nk,1)
		}
		
		# Way 4: sample x0 from filtered results at t=0  # note rmvnorm and rnorm
		if (FALSE){
			vFilter0 = chol2inv(chol(outBackFilter$precisionFilter0))
			muFilter0 = vFilter0 %*% outBackFilter$thetaFilter0
			if (nk==1){
				x0 = matrix(rnorm(1, muFilter0, sd=sqrt(vFilter0)), ncol=1)
			} else {
				if (!isSymmetric(vFilter0)) vFilter0=(vFilter0+t(vFilter0))/2 
				x0 = matrix( rmvnorm(1, muFilter0, vFilter0), ncol=1 )
			}
			# if some Cts have extreame backward filtered x0 values (e.g. simu Id=10, too few obs), sample x0 from prior
			n.outlier = sum(abs(x0/scaleFactor)>0.1)
			x0[abs(x0/scaleFactor)>0.1] = rnorm(n.outlier, mux0, sqrt(Vx0))
		}
		
		# Way 5: The correct way to take prior of x_0 
		if (FALSE){
			mu00 = rep(mux0, nk)
			V00 = matrix(0, nk,nk)
			diag(V00) = Vx0
			V1Smooth = chol2inv(chol(outBackFilter$precisionFilter[[1]]))
			mu1Smooth = V1Smooth %*% outBackFilter$thetaFilter[[1]]
			muPredict = aVec.k*mu00
			#VPredict = A%*% VFilter %*% A + Q.k
			VPredict = (aVec.k %o% aVec.k) * V00 + Q.k
			A = matrix(0, nk, nk)
			diag(A) = aVec.k
			J0 = V00 %*% A %*% chol2inv(chol(VPredict))
			mu0SmoothPrior = mu00 +  J0 %*% (mu1Smooth - muPredict)
			V0SmoothPrior  = V00 + J0 %*% (V1Smooth - VPredict) %*% t(J0)
			if (nk==1){
				x0 = matrix( rnorm(1, mu0SmoothPrior, sd=sqrt(V0SmoothPrior)), ncol=1 )
			} else {
				x0 = matrix( rmvnorm(1, mu0SmoothPrior, V0SmoothPrior), ncol=1 )
			}
			x0[abs(x0)>1*scaleFactor] = 0
		}
		
		
		if (FALSE){
		xMat[as.character(CtIds.k),] <- ForwardSamplerX(outBackFilter$precisionFilter0, outBackFilter$thetaFilter0, 
						outBackFilter$precisionFilter, outBackFilter$thetaFilter, Q.k, tmLab, CtIds.k, aVec.k, x0)			  
		x0Vec[as.character(CtIds.k)] = x0
		}
		
		# Forward Filter and Backward Sampler
		outSuff = get.Ct.yt.suffStat(L, covEffMonth, yDehedonicMonthMean, CtIds.k, tmLab, tract.i=NULL)			
		outFilter = KalmanFilter.suffStat(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, outSuff$Lt, outSuff$yDehedonicMean, outSuff$CtBar, mux0, Vx0 )
		outSampler = BackwardSamplerX(outFilter$muPredictList, outFilter$VPredictList, outFilter$muFilterList, outFilter$VFilterList, tmLab, CtIds.k, aVec.k, mux0, Vx0)
		xMat[as.character(CtIds.k),] = outSampler$x 
		#xMat[as.character(CtIds.k),] = outSampler$x - mean(outSampler$x) # adjust for the mean, July 15 2015
		x0Vec[as.character(CtIds.k)] = outSampler$x0 
		#x0Vec[as.character(CtIds.k)] = outSampler$x0 - mean(outSampler$x) # adjust for the mean, July 15 2015
		muTVec[as.character(CtIds.k)] = as.numeric(outFilter$muFilterList[[T]])
		
		if (F){ # BFFS
			plot(x[,9], type="l", ylim=c(-0.5, 0.5))
			for (i in 1:100){
				#x0 = matrix( rmvnorm(1, mu0SmoothPrior, V0SmoothPrior), ncol=1 )
				x0 =matrix(0, nk,1)
				x0[abs(x0)>1*scaleFactor] = 0
				xMat[as.character(CtIds.k),] <- ForwardSamplerX(outBackFilter$precisionFilter0, outBackFilter$thetaFilter0, 
						outBackFilter$precisionFilter, outBackFilter$thetaFilter, Q.k, tmLab, CtIds.k, aVec.k, x0)			  
				lines(xMat[as.character(CtIds.k)[1],]/200, col="grey")
			}
			lines(x[,9])
		}
		
		if (F){ # FFBS
			plot(x[,10], type="l", ylim=c(-0.5, 1))
			out <- get.Ct.yt.suffStat(L, covEffMonth, yDehedonicMonthMean, CtIds.k, tmLab, tract.i=NULL)			
			outFilter = KalmanFilter.suffStat(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, out$Lt, out$yDehedonicMean,out$CtBar,mux0, Vx0 )
			for (i in 1:100){
				outSampler = BackwardSamplerX(outFilter$muPredictList, outFilter$VPredictList, outFilter$muFilterList, outFilter$VFilterList, tmLab, CtIds.k, aVec.k, mux0, Vx0)
				xMat[as.character(CtIds.k),] = outSampler$x
				x0Vec[as.character(CtIds.k)] = outSampler$x0
				lines(xMat[as.character(CtIds.k)[2],]/200, col="grey")
			}
			lines(x[,10])
		}
		
	}
	return(list(xMat=xMat, x0Vec=x0Vec, muTVec=muTVec))
}

# 3b. Sample the tract process x conditional on eta
sample.x.cond.eta <- function(zVec, LambdaMat, aVec, RVec, sigma0 ,L, yMonth, covEffMonth, etaMat, tmLab, CtIds, covNames, mux0, Vx0){
	K = length(table(zVec))
	n = length(CtIds)
	kLabels  = as.numeric(names(table(zVec)))
	T = length(tmLab)
	xMat = matrix(0, n, T)
	rownames(xMat) = CtIds
	x0Vec= rep(0, n)
	names(x0Vec) = CtIds
	sigma0sq = sigma0*sigma0

	for (i in 1:n){
		CtIds.i = CtIds[i]
		lambda.ik = LambdaMat[as.character(CtIds.i),as.character(zVec[i])] 
		A = aVec.i = aVec[as.character(CtIds.i)]
		RVec.i = RVec[as.character(CtIds.i)]
		Li = L[L$CensusTractRegionID == CtIds.i,]
		eta.k = etaMat[,as.character(zVec[i])]

		## FFBS
		## Forward Filter of x process
		
		outFilter = KalmanFilter.cond.eta( tract.i=CtIds[i], RVec.i, aVec.i, sigma0sq, eta.k, lambda.ik, covEffMonth, yMonth, tmLab, L, mux0, Vx0 )
		outSampler = BackwardSamplerX.cond.eta(outFilter$muPredictList, outFilter$VPredictList, outFilter$muFilterList, outFilter$VFilterList,  tmLab, CtIds.i, aVec.i, mux0, Vx0)
		xMat[i,] = outSampler$x
		x0Vec[i] = outSampler$x0
	}
	
	return(list(xMat=xMat, x0Vec=x0Vec))
}

# 3b. Sample the tract process x conditional on eta
sample.x.cond.eta.parallel <- function(zVec, LambdaMat, aVec, RVec, sigma0 ,L, yMonth, covEffMonth, etaMat, tmLab, CtIds, covNames, mux0, Vx0){
	K = length(table(zVec))
	n = length(CtIds)
	kLabels  = as.numeric(names(table(zVec)))
	T = length(tmLab)
	xMat = matrix(0, n, T)
	rownames(xMat) = CtIds
	x0Vec= rep(0, n)
	names(x0Vec) = CtIds
	sigma0sq = sigma0*sigma0

	parallelOutX <- foreach(i = 1:n) %dopar%{	
	#for (i in 1:n){
		CtIds.i = CtIds[i]
		lambda.ik = LambdaMat[as.character(CtIds.i),as.character(zVec[i])] 
		A = aVec.i = aVec[as.character(CtIds.i)]
		RVec.i = RVec[as.character(CtIds.i)]
		Li = L[L$CensusTractRegionID == CtIds.i,]
		eta.k = etaMat[,as.character(zVec[i])]

		## BFFS
		if (FALSE){
			## Backward filtering of x process
			tt = T
			precisionFilter = rep(0, T)
			thetaFilter = rep(0,T)
			Lt = Li$count[Li$MonthTransaction==tmLab[tt]]
			precisionFilter[tt] = Lt/RVec.i
			covEff.t = covEffMonth[[tt]]
			yt = yMonth[[tt]]
			residual = yt$y[yt$CensusTractRegionID==CtIds.i]-covEff.t$covEff[covEff.t$CensusTractRegionID==CtIds.i ]
			thetaFilter[tt] = sum(residual)/RVec.i
			
			for (tt in (T-1):1){
				Lt = Li$count[Li$MonthTransaction==tmLab[tt]]
				# backward gain
				JtPlus1 = precisionFilter[tt+1] /(precisionFilter[tt+1]+1/sigma0sq)
				LtPlus1 = 1 - JtPlus1 
				# predict backward
				precisionPredict = A*A * (LtPlus1*LtPlus1* precisionFilter[tt+1] + JtPlus1 *JtPlus1/ sigma0sq) 
				thetaPredict = A * LtPlus1 * (thetaFilter[tt+1] - precisionFilter[tt+1]*lambda.ik*etaMat[tt+1,as.character(zVec[i])])
				# update to filtered values
				if (Lt==0){
					precisionFilter[tt] = precisionPredict 
					thetaFilter[tt] = thetaPredict
				} else{
					precisionFilter[tt] = precisionPredict + Lt/RVec.i
					covEff.t = covEffMonth[[tt]]
					yt = yMonth[[tt]]
					residual = yt$y[yt$CensusTractRegionID==CtIds.i]-covEff.t$covEff[covEff.t$CensusTractRegionID==CtIds.i ]
					thetaFilter[tt] = thetaPredict + sum(residual)/RVec.i
				}	
			}
			
			# fix x0 at 0
			x0 = 0
			
			## Forward sample x
			xx = rep(0,  T)		
			tt =1
			xtMinus1=x0
			V = 1/( 1/sigma0sq + precisionFilter[tt] )
			mu= V*(1/sigma0sq*A*xtMinus1 + thetaFilter[tt] + lambda.ik*etaMat[tt,as.character(zVec[i])]/sigma0sq)
			xx[tt] = rnorm(1, mean=mu, sd=sqrt(V))
			for (tt in 2:T){
				xtMinus1=xx[tt-1]
				V = 1/( 1/sigma0sq + precisionFilter[tt] )
				mu= V*(1/sigma0sq*A*xtMinus1 + thetaFilter[tt] + lambda.ik*etaMat[tt,as.character(zVec[i])]/sigma0sq)
				xx[tt] = rnorm(1, mean=mu, sd=sqrt(V))
			}
			
			xx
			#xMat[as.character(CtIds.i),] = xx
			#x0Vec[as.character(CtIds.i)] = x0
		}
		
		## FFBS
		## Forward Filter of x process
		
		outFilter = KalmanFilter.cond.eta( tract.i=CtIds[i], RVec.i, aVec.i, sigma0sq, eta.k, lambda.ik, covEffMonth, yMonth, tmLab, L, mux0, Vx0 )
		outSampler = BackwardSamplerX.cond.eta(outFilter$muPredictList, outFilter$VPredictList, outFilter$muFilterList, outFilter$VFilterList,  tmLab, CtIds.i, aVec.i, mux0, Vx0)
		x = outSampler$x
		x0 = outSampler$x0
		out = list(x=x, x0=x0)
		out
	}
	
	for (i in 1:n){
		out = parallelOutX[[i]]
		xMat[i,] = out$x
		x0Vec[i] = out$x0
	}
	
	return(list(xMat=xMat, x0Vec=x0Vec))
}

if (FALSE){
plot(xMat[as.character(CtIds.k)[1],]/200, type="l", ylim=c(-1,1))
lines(xMat[as.character(CtIds.k)[2],]/200, col=2)
lines(xMat[as.character(CtIds.k)[3],]/200, col=3)
lines(xMat[as.character(CtIds.k)[4],]/200, col=4)
}

# 4. Sample the cluster latent factor eta_t,k
sample.eta <- function(zVec, LambdaMat, aVec, sigma0, xMat, tmLab, x0Vec=NULL, 	K = length(unique(zVec))){
	n = length(zVec)
	T = length(tmLab)
	CtIdsChar = names(zVec)
	kLabels  = colnames(LambdaMat)
	#LambdaMat = LambdaMat[CtIdsChar,]
	etaMat = matrix(0, T, K)
	colnames(etaMat) = kLabels

	tt = 1
	for (k in as.numeric(kLabels) ){
		idctk = zVec==k
		lambdak = LambdaMat[idctk, as.character(k)]
		Vk = 1/(1+ (lambdak %*% lambdak)/(sigma0*sigma0))
		tem = lambdak %*% (xMat[idctk,tt] - aVec[idctk]*x0Vec[idctk])
		muk = Vk/(sigma0*sigma0)*tem
		etaMat[tt,as.character(k)] = rnorm(1, muk, sd=sqrt(Vk))
	}
	
	for (tt in 2:T){
		for (k in as.numeric(kLabels) ){
			idctk = zVec==k
			lambdak = LambdaMat[idctk, as.character(k)]
			Vk = 1/(1+ (lambdak %*% lambdak)/(sigma0*sigma0))
			tem = lambdak %*% (xMat[idctk,tt] - aVec[idctk]*xMat[idctk,tt-1])
			muk = Vk/(sigma0*sigma0)*tem
			etaMat[tt,as.character(k)] = rnorm(1, muk, sd=sqrt(Vk))
		}
	}
	
	return(etaMat)
}


sample.eta.old <- function(zVec, LambdaMat, aVec, sigma0, xMat, tmLab, x0Vec=NULL, 	K = length(unique(zVec))){
	n = length(zVec)
	zMat = matrix(0, n, K)
	T = length(tmLab)
	CtIdsChar = names(zVec)
	kLabels  = colnames(LambdaMat)
	colnames(zMat) = kLabels
	rownames(zMat) = CtIdsChar
	
	LambdaMat = LambdaMat[CtIdsChar,]
	zMat[cbind(CtIdsChar, as.character(zVec))] = 1
		
	LambdaZMat = LambdaMat*zMat
	Ik = matrix(0, K,K)
	diag(Ik) = 1
	V = chol2inv(chol(1/sigma0^2*t(LambdaZMat)%*% LambdaZMat + Ik))
	A = matrix(0, n,n)
	diag(A) = aVec[names(zVec)]
	etaMat = matrix(0, T, K)
	
	tt = 1
	leftMat = 1/sigma0^2*V %*% t(LambdaZMat) 
	if (is.null(x0Vec)) {
		mu = leftMat %*% xMat[,tt]
	} else {
		mu = leftMat %*% (xMat[,tt]-A %*% matrix(x0Vec,ncol=1))
	}
	etaMat[tt, ] = rmvnorm(1, mu, V)

	for (tt in 2:T){
		mu = leftMat %*% (xMat[,tt]-A %*% xMat[,tt-1])
		etaMat[tt,] = rmvnorm(1, mu, V)
	}
	colnames(etaMat) = kLabels
	return(etaMat)
}

# 5. Sample the loading matrix lambda_ik. dim(LambdaMat) = n*K
sample.LambdaMat <- function(zVec, muLambda, sigmaLambda, sigma0, etaMat, xMat, aVec, tmLab, CtIds, x0Vec, K = length(unique(zVec))){
	n = length(zVec)
	T = length(tmLab)
	LambdaMat = matrix(0, length(CtIds), K)
	kLabels  = colnames(etaMat)
	rownames(LambdaMat) = CtIds
	colnames(LambdaMat) = kLabels
	
	for (i in 1:n){
		k = as.character(zVec[i])
		Cti = names(zVec[i])
		# for z_ik = 0
		LambdaMat[Cti,kLabels!=k] = rnorm(K-1, muLambda, sigmaLambda)
		# for z_ik = 1
		v = (1/sigmaLambda^2 + 1/sigma0^2*sum(etaMat[,k]^2))^(-1)
		epsilon = rep(0,T)
		epsilon[1] = xMat[Cti,1] - aVec[Cti]*x0Vec[Cti]
		epsilon[2:T] = xMat[Cti,2:T] - aVec[Cti]*xMat[Cti,1:(T-1)] 
		mu = v*(muLambda/sigmaLambda^2 + sum(epsilon*etaMat[,k])/sigma0^2)
		LambdaMat[names(zVec[i]), k] = rnorm(1, mu, sd=sqrt(v)) 
	}
	if (n!=length(CtIds)){
		idct = !(CtIds %in% names(zVec))
		LambdaMat[idct,] = rnorm(sum(idct)*K, muLambda, sigmaLambda)
	}
	return(LambdaMat)
}

# 6. Sample AR coefficient a_i conditional on eta. length(aVec) = n
sample.aVec <- function(zVec, LambdaMat, xMat,sigma0, etaMat, muA, sigmaA, CtIds, tmLab, x0Vec){
	n = length(zVec)
	T = length(tmLab)
	aVec = rep(0, n)
	names(aVec) = CtIds
	
	for (i in 1:n){
		x0i = x0Vec[i]
		x0_Tminus1 = c(x0i, xMat[i,-T])
		cumSumV = sum(x0_Tminus1*x0_Tminus1) / (sigma0*sigma0)
		clusteri = as.character(zVec[i])
		lambdai = LambdaMat[i, clusteri]
		cumSumMu = sum( x0_Tminus1*(xMat[i,]-lambdai*etaMat[,clusteri]) ) / (sigma0*sigma0)
		
		V = 1 / (1/(sigmaA*sigmaA) + cumSumV)
		mu = V * ( muA/(sigmaA*sigmaA) + cumSumMu )
		aVec[i] = rnorm(1, mu, sd=sqrt(V)) 
	}
	return(aVec)
}

sample.aVec.integrateEta <- function(zVec, LambdaMat, xMat,sigma0, muA, sigmaA, CtIds, tmLab, x0Vec){
	n = length(zVec)
	K = length(unique(zVec))
	T = length(tmLab)
	kLabels.cha = names(table(zVec))
	kLabels.num = as.numeric(kLabels.cha)
	aVec = rep(0, n)
	names(aVec) = CtIds
	
	for (k in kLabels.num){
		CtIds.k = CtIds[zVec==k]
		n.k = length(CtIds.k)
		sigma0diag = matrix(0, length(CtIds.k), length(CtIds.k))
		diag(sigma0diag) = sigma0^2
		Q.k = LambdaMat[as.character(CtIds.k),as.character(k)] %o% LambdaMat[as.character(CtIds.k),as.character(k)]+
					sigma0diag
		cumSum1 = cumSum2 = 0
		
		tt = 1
		if (all(x0Vec==0)){
			cumSum1 = matrix(0, n.k, n.k)
			cumSum2 = matrix(0, n.k, 1)
		} else {
			xtMinus1 = x0Vec[as.character(CtIds.k)]
			SigmaTild.t = Q.k/(xtMinus1 %o% xtMinus1)
			invSigmaTild.t=tryCatch(chol2inv(chol(SigmaTild.t)), error = function(e){matrix(0, n.k, n.k)})
			cumSum1 = cumSum1 + invSigmaTild.t
			cumSum2 = cumSum2 + invSigmaTild.t %*% (xMat[as.character(CtIds.k), tt]/xtMinus1)
		}
		
		for (tt in 2:T){
			xtMinus1 = xMat[as.character(CtIds.k),tt-1]
			SigmaTild.t = Q.k/(xtMinus1 %o% xtMinus1)
			invSigmaTild.t=tryCatch(chol2inv(chol(SigmaTild.t)), error = function(e){matrix(0, n.k, n.k)})
			cumSum1 = cumSum1 + invSigmaTild.t
			cumSum2 = cumSum2 + invSigmaTild.t %*% (xMat[as.character(CtIds.k), tt]/xtMinus1)
		}
		invSigmaAdiag = matrix(0, n.k, n.k)
		diag(invSigmaAdiag) = 1/sigmaA^2
		V =  chol2inv(chol(invSigmaAdiag + cumSum1))
		mu = V %*% (invSigmaAdiag %*% rep(muA, n.k)+ cumSum2)
		aVec[as.character(CtIds.k)] = rmvnorm(1, mu, V) 
	}
	return(aVec)
}

# 7. Sample the response variance R_i. length(RVec) = n
sample.RVec <- function(CtIds, xMat, alpha.R0, beta.R0, yTract, tmLab, covEffMonth ){
	n = length(CtIds)
	RVec=rep(0,n)
	names(RVec) = CtIds
	
	covEff <- ldply (covEffMonth, data.frame)
  
	for (i in 1:n){
		CtId.i = CtIds[i]
		obs.i = yTract[[as.character(CtId.i)]]
		covEff.i = covEff$covEff[covEff$CensusTractRegionID==CtId.i]
		mi = nrow(obs.i)
		x = xMat[as.character(CtId.i), ]
		id = match(obs.i$MonthTransaction, tmLab)
		RVec[i] = rigamma(1, alpha.R0 + mi/2, beta.R0+ sum((obs.i$y- x[id]-covEff.i)^2)/2 )
	}
	return(RVec)
}

# 7 [parallel]. Sample the response variance R_i. length(RVec) = n
sample.RVec.parallel <- function(CtIds, xMat, alpha.R0, beta.R0, yTract, tmLab, covEffMonth ){
	n = length(CtIds)
	covEff = matrix(unlist(covEffMonth),ncol=2, byrow=FALSE )
	
	covEff <- ldply (covEffMonth, data.frame)

	parallelOutR <- foreach(i = 1:n) %dopar%{		
		CtId.i = CtIds[i]
		obs.i = yTract[[as.character(CtId.i)]]
		covEff.i = covEff$covEff[covEff$CensusTractRegionID==CtId.i]
		mi = nrow(obs.i)
		x = xMat[as.character(CtId.i), ]
		id = match(obs.i$MonthTransaction, tmLab)
		Ri = rigamma(1, alpha.R0 + mi/2, beta.R0+ sum((obs.i$y- x[id]-covEff.i)^2)/2 )
		Ri
	}
	
	RVec = as.numeric(parallelOutR)
	names(RVec) = CtIds	
	return(RVec)
}

# 8. Sample the covariate effect Beta_ir. dim(BetaMat)= n*length(covNames)
sample.BetaMat <- function(muBetaVec, sigmaBetaVec, RVec, xMat, BetaMat, yTract, CtIds, tmLab,covNames){
	n = length(CtIds)
	R = length(covNames)
	
	for (i in 1:n){
		CtId.i = CtIds[i]
		obs.i = yTract[[as.character(CtId.i)]]
		xT = xMat[as.character(CtId.i), ]
		id = match(obs.i$MonthTransaction, tmLab)
		x = xT[id] 
		rSeq = sample(1:R)
		for (r in rSeq){
			U = obs.i[,covNames[r]]
			covEffect = as.matrix(obs.i[,covNames[-r]])%*% BetaMat[as.character(CtId.i),covNames[-r]] 
			residual = obs.i$y - x - as.numeric(covEffect)
			v = (1/sigmaBetaVec[r]^2 + 1/RVec[as.character(CtId.i)]*sum(U^2) )^(-1)
			mu = v*(muBetaVec[r]/sigmaBetaVec[r]^2 + 1/RVec[as.character(CtId.i)]*sum(U*residual))
			BetaMat[as.character(CtId.i),r] = rnorm(1, mu, sd=sqrt(v))
		}
	}
  
  # adjust the mean of beta0, July 16 2015
  if ("intercept" %in% covNames){
    adjMean = weighted.mean(BetaMat[,covNames[1]],Ltract$count )
    BetaMat[,covNames[1]] = BetaMat[,covNames[1]]-adjMean  
  }
	return(BetaMat)
}

# 9. Sample x process error, a scalar
sample.sigma0 <- function(xMat, etaMat, aVec, LambdaMat, zVec, alpha.epsilon0, beta.epsilon0, CtIds, tmLab,x0Vec=NULL,	K = length(unique(zVec))){
	n = length(CtIds)
	T = length(tmLab)
	kLabels  = colnames(LambdaMat)
	zMat = matrix(0, n, K)
	colnames(zMat) = kLabels
	rownames(zMat) = CtIds

	#for (i in 1:n){ zMat[i, as.character(zVec[i])] = 1 }
	zMat[cbind(CtIds, as.character(zVec))] = 1
	LambdaZMat = LambdaMat*zMat
	
	tt = 1
	if (is.null(x0Vec)) {x0Vec = matrix(0, n,1)}
	cumSum = sum( (xMat[,tt] - aVec*x0Vec-LambdaZMat %*% etaMat[tt,])^2 )
	for (tt in 2:T){
		cumSum = cumSum + sum( (xMat[,tt] - aVec*xMat[,tt-1] - LambdaZMat %*% etaMat[tt,])^2 )
	}
	sigma0square = rigamma(1, alpha.epsilon0+T*n/2, beta.epsilon0+1/2*cumSum )
	sigma0 = sqrt(sigma0square)
	return(sigma0)
}

if (F){
tem =  pigamma (seq(0.1, 5, 0.1), alpha.epsilon0+T*n/2, beta.epsilon0+1/2*cumSum )
plot(sqrt(seq(0.1, 5, 0.1))/ scaleFactor, tem)
}

# 10. Sample muLambda, a scalar
sample.muLambda <- function(muLambda0, sigmaLambda0, LambdaMat, zVec, sigmaLambda,K = length(unique(zVec))){
	n = length(zVec)
	v = (1/sigmaLambda0^2 + n/sigmaLambda^2)^(-1)
	kLabels  = colnames(LambdaMat)
	zMat = matrix(0, n, K)
	colnames(zMat) = kLabels

	for (i in 1:n){ zMat[i, as.character(zVec[i])] = 1 }
	LambdaZMat = LambdaMat*zMat
	mu = v*(muLambda0/sigmaLambda0^2 + 1/sigmaLambda^2*sum(LambdaZMat))
	muLambda = rnorm(1, mu, sqrt(v))
	return(muLambda)
}

# 11. Sample sigmaLambda, a scalar
sample.sigmaLambda <- function(zVec, LambdaMat, muLambda, alphaLambda0, betaLambda0,K = length(unique(zVec))){
	n = length(zVec)
	kLabels  = colnames(LambdaMat)
	zMat = matrix(0, n, K)
	colnames(zMat) = kLabels

	for (i in 1:n){ zMat[i, as.character(zVec[i])] = 1 }
	LambdaZMat = LambdaMat*zMat
	
	LambdaZVec = as.numeric(LambdaZMat)
	Lambdai = LambdaZVec[LambdaZVec!=0]
	
	sigmaLambdaSquare = rigamma(1, alphaLambda0+ n/2, betaLambda0+1/2*sum((Lambdai-muLambda)^2) )
	sigmaLambda = sqrt(sigmaLambdaSquare)
	return(sigmaLambda)
}

# 12. Sample muA, a scalar
sample.muA <- function(aVec, sigmaA, muA0, sigmaA0){
	n = length(aVec)
	v = (1/sigmaA0^2+n/sigmaA^2)^(-1)
	mu = v*(muA0/sigmaA0^2 + 1/sigmaA^2*sum(aVec))
	muA = rnorm(1, mu, sqrt(v))
	return(muA)
}

# 13. Sample sigmaA, a scalar
sample.sigmaA <- function(aVec, muA, alphaA0, betaA0){
	n = length(aVec)
	sigmaAsquare = rigamma(1, alphaA0+n/2, betaA0+1/2*sum((aVec-muA)^2) )
	sigmaA = sqrt(sigmaAsquare)
	return(sigmaA)
}

# 14. Sample the mean covariate effect, muBetaVec. length(muBetaVec) = length(covNames)
sample.muBetaVec <- function(BetaMat, sigmaBetaVec, muBetaVec0, sigmaBetaVec0, covNames){
	n = nrow(BetaMat)
	R = length(covNames)
	muBetaVec = rep(0, R)
	names(muBetaVec) = covNames
	rSeq = sample(1:R)
	for (r in rSeq){
		v = ( 1/sigmaBetaVec0[r]^2+n/sigmaBetaVec[r]^2 )^(-1)
		mu = v*(muBetaVec0[r]/sigmaBetaVec0[r]^2 + 1/sigmaBetaVec[r]^2*sum(BetaMat[,r]))
		muBetaVec[r] = rnorm(1, mu, sqrt(v))
	}
	return(muBetaVec)
}

# 15. Sample the variability of covariate effect, sigmaBetaVec, length(sigmaBetaVec) = length(covNames)
sample.sigmaBetaVec <- function(BetaMat, muBetaVec, alphaBetaVec0,  betaBetaVec0, covNames){
	n = nrow(BetaMat)
	R = length(covNames)
	sigmaBetaVec = rep(0, R)
	names(sigmaBetaVec) = covNames
	rSeq = sample(1:R)
	for (r in rSeq){
		 sigmaSquare = rigamma(1, alphaBetaVec0[r]+n/2, betaBetaVec0[r]+1/2*sum((BetaMat[,r]-muBetaVec[r])^2))
		 sigmaBetaVec[r] = sqrt(sigmaSquare)
	}
	return(sigmaBetaVec)
}

# 17. Sample the processor allocation indicator, pi for each cluster

sample.pi <- function(piVec, nP, zVec, K){
	piVecStar = sample(1:nP, size=K, replace=TRUE)
	names(piVecStar) = 1:K
	sizeCluster = table(zVec)
	tab = data.frame(cbind(sizeCluster, piVec, piVecStar))
	accp = 1
	for (p in 1:nP){
		size.p = tab$sizeCluster[tab$piVec==p]
		if (length(size.p)==0){
			above =1
		} else {
			a = table(size.p)
			above = prod( sapply(a, FUN=factorial) )
		}
		size.p.star = tab$sizeCluster[tab$piVecStar==p]
		if (length(size.p.star)==0){
			below=1
		} else {
			aStar = table(size.p.star)
			below = prod(sapply(aStar, FUN=factorial))
		}
		accp = accp* (above/below)
	}
	u = runif(1)
	if (u<accp) {
		piV = piVecStar
	} else {
		piV = piVec
	}
	return(piV)
}

# Sample the coefficients of NCS in the global trend
sample.g.coefNCS <- function(dataForG, sVec, RVec, wBasisVec, muWjVec, sigmaWjVec){
  RVec.rescale = RVec/(scaleFactor^2)
  # reorder RVec
  R = RVec.rescale[as.character(dataForG$CensusTractRegionID)]
  # add month 1 to sVec
  sVecExtend = c(0, sVec)
  names(sVecExtend)[1] = 1
  basisNames = names(wBasisVec)
  nBasis = length(basisNames)
  jSeq = sample(1: nBasis)
  wBasisVecCurrent = wBasisVec
  for (j in jSeq){
    Bj = dataForG[,basisNames[j]]
    precision = 1/(sigmaWjVec[basisNames[j]]^2) + sum( Bj^2/R )
    v = 1/precision
    tempV = Bj*(dataForG$residualObs - as.matrix(dataForG[,basisNames[-j]])%*%wBasisVec[basisNames[-j]] - sVecExtend[as.character(dataForG$MonthNumber)] )
    mu = sum(v*tempV/R) + v*muWjVec[basisNames[j]]/sigmaWjVec[basisNames[j]]^2
    wBasisVecCurrent[basisNames[j]] = rnorm(1, mean = mu, sd = sqrt(v))
    #cat("j =",j, "mu =", mu, "v =",v, "\n")
  }
  wBasisVec = wBasisVecCurrent
  return(wBasisVec)
}

sample.g.monthlyEffect <- function(dataForG, sVec, RVec, wBasisVec, muSjVec, sigmaSjVec, monthlyCount=NULL){
  RVec.rescale = RVec/(scaleFactor^2)
  # reorder RVec
  R = RVec.rescale[as.character(dataForG$CensusTractRegionID)]
  monthNames = names(sVec)
  nMonth = length(monthNames)
  jSeq = sample(1:nMonth)
  basisNames = names(wBasisVec)
  temp = (dataForG$residualObs -as.matrix(dataForG[,basisNames])%*%wBasisVec[basisNames] )
  for (j in jSeq){
    Ij = (dataForG$MonthNumber==as.numeric(monthNames[j]))*1
    precision = 1/(sigmaSjVec[monthNames[j]]^2) + sum(Ij/R)
    v = 1/precision
    tempV = Ij*temp
    mu = sum(v*tempV/R) + v*muSjVec[monthNames[j]]/sigmaSjVec[monthNames[j]]^2
    sVec[monthNames[j]] = rnorm(1, mean=mu, sd=sqrt(v))  
  }
  if (!is.null(monthlyCount)){
    overallMean = sVec%*%monthlyCount/sum(monthlyCount)
    sVec = sVec - overallMean
  }
  names(sVec) = names(muSjVec)
  return(sVec)
}

sample.g.monthlyEffect.m11 <- function(dataForG, sVec, RVec, wBasisVec, muSjVec, sigmaSjVec){
  RVec.rescale = RVec/(scaleFactor^2)
  # reorder RVec
  R = RVec.rescale[as.character(dataForG$CensusTractRegionID)]
  monthNames = names(sVec)
  nMonth = length(monthNames)
  jSeq = sample(1:nMonth)
  basisNames = names(wBasisVec)
  temp = (dataForG$residualObs -as.matrix(dataForG[,basisNames])%*%wBasisVec[basisNames] )
  for (j in jSeq){
    Ij = (dataForG$MonthNumber==as.numeric(monthNames[j]))*1
    precision = 1/(sigmaSjVec[monthNames[j]]^2) + sum(Ij/R)
    v = 1/precision
    tempV = Ij*temp
    mu = sum(v*tempV/R) + v*muSjVec[monthNames[j]]/sigmaSjVec[monthNames[j]]^2
    sVec[monthNames[j]] = rnorm(1, mean=mu, sd=sqrt(v))  
  }
  return(sVec)
}

#############################################
### Calculate log posterior               ###
### (one metric to monitor the sampler)   ###
#############################################

getLogPosterior <- function(CtIds, zVec, muA, muA0, sigmaA0, sigmaA, alphaA0,betaA0,aVec,
						muLambda,muLambda0,sigmaLambda0,sigmaLambda, alphaLambda0,betaLambda0,LambdaMat,
						etaMat, sigma0,alpha.epsilon0,beta.epsilon0,
						alpha, alpha.alpha0,beta.alpha0, RVec, alpha.R0, beta.R0,
						covNames, BetaMat, muBetaVec, sigmaBetaVec,alphaBetaVec0,betaBetaVec0,muBetaVec0,sigmaBetaVec0,
						x0Vec, xMat, mux0, Vx0,
						L, yMonth, yTract, details=FALSE){
	# (22)	
	n = length(CtIds)
	part1 = dnorm(muA,mean=muA0,sd=sigmaA0, log=TRUE) +
			log(densigamma(sigmaA^2, alpha=alphaA0 , beta=betaA0)) +
			sum( dnorm(aVec, mean=muA, sd = sigmaA, log=TRUE) )
    # (23)
	part2 = dnorm(muLambda,mean=muLambda0,sd=sigmaLambda0, log=TRUE)+
		log(densigamma(sigmaLambda^2, alpha=alphaLambda0 , beta=betaLambda0))+
		sum( dnorm(LambdaMat[cbind(as.character(CtIds), as.character(zVec))], mean=muLambda, sd = sigmaLambda, log=TRUE) )

	# (24)
	if (FALSE){ # marginal likelihood that integrates out eta
		part3 = sum( dnorm(as.numeric(etaMat), 0, 1, log=TRUE) ) +
			log(densigamma(sigma0^2, alpha=alpha.epsilon0, beta=beta.epsilon0))
	}
	part3 = log(densigamma(sigma0^2, alpha=alpha.epsilon0, beta=beta.epsilon0))
	
	# (25) 
	K = length(unique(zVec))
	part4 = dgamma(alpha, shape =alpha.alpha0, rate=beta.alpha0, log=TRUE) +
		lgamma(alpha) + sum(lgamma(table(zVec)))+ K*log(alpha) - lgamma(n+alpha) +
		sum( log(densigamma(RVec, alpha=alpha.R0 , beta=beta.R0)) )

	# (26)
	pBeta=0
	for (r in 1:length(covNames)){
		pBeta = pBeta + sum( dnorm(BetaMat[,r], mean=muBetaVec[r], sd=sigmaBetaVec[r], log=TRUE) )+
				log(densigamma(sigmaBetaVec[r]^2, alpha=alphaBetaVec0[r], beta=betaBetaVec0[r] )) +
				dnorm(muBetaVec[r], mean=muBetaVec0[r], sd=sigmaBetaVec0[r], log=TRUE) 
	}
	part5 = pBeta

	
	# (27)
	kLabels = as.numeric( names(table(zVec)) )
	pX = 0
	for (k in kLabels){
		CtIds.k = sort(CtIds[zVec==k])
		nk = length(CtIds.k)
		sigma0diag = matrix(0, length(CtIds.k), length(CtIds.k))
		diag(sigma0diag) = sigma0^2
		Q.k = LambdaMat[as.character(CtIds.k),as.character(k)] %o% LambdaMat[as.character(CtIds.k),as.character(k)]+sigma0diag
		aVec.k = aVec[as.character(CtIds.k)]
		x0Vec.k = x0Vec[as.character(CtIds.k)]
		xMat.k = xMat[as.character(CtIds.k),] 
		pX.k = 0
		
		if (nk==1){
			pX.k = dnorm(xMat.k[1], mean=aVec.k*x0Vec.k, sd= as.numeric(sqrt(Q.k)),log=TRUE)
			for (tt in 2:length(tmLab)){pX.k = pX.k + dnorm(xMat.k[tt], mean=aVec.k*xMat.k[tt-1], sd= as.numeric(sqrt(Q.k)), log=TRUE)}
		} else {
			pX.k = dmvnorm(xMat.k[,1], mean=aVec.k*x0Vec.k, sigma= Q.k,log=TRUE)
			for (tt in 2:length(tmLab)){pX.k = pX.k + dmvnorm(xMat.k[,tt], mean=aVec.k*xMat.k[,tt-1], sigma= Q.k, log=TRUE)}
		}
		pX = pX + pX.k
	}
	part6 = sum( dnorm(x0Vec, mean=mux0, sd=sqrt(Vx0), log=TRUE) ) + pX

	# (28)
	py = 0
	# new way by tract
	for (i in 1:n){
		CtId.i = CtIds[i]
		obs.i = yTract[[as.character(CtId.i)]]
		xT = xMat[as.character(CtId.i), ]
		id = match(obs.i$MonthTransaction, tmLab)
		xMatch = xT[id] 
		covEffect = as.matrix(obs.i[,covNames])%*% BetaMat[as.character(CtId.i),covNames] 
		residual = obs.i$y - xMatch - as.numeric(covEffect)
		Ri = RVec[as.character(CtId.i)]
		py = py + sum(dnorm(residual, mean=0, sd=sqrt(Ri), log=TRUE))
	}
	part7=py

	logPosterior = part1 + part2 + part3 + part4 + part5 + part6 + part7
	if (details){
		out = c(part1, part2, part3, part4, part5, part6, part7, logPosterior)
		return(out)
	} else {
		return(logPosterior)
	}
	
}


getLogPosterior.eta <- function(CtIds, zVec, muA, muA0, sigmaA0, sigmaA, alphaA0,betaA0,aVec,
						muLambda,muLambda0,sigmaLambda0,sigmaLambda, alphaLambda0,betaLambda0,LambdaMat,
						etaMat, sigma0,alpha.epsilon0,beta.epsilon0,
						alpha, alpha.alpha0,beta.alpha0, RVec, alpha.R0, beta.R0,
						covNames, BetaMat, muBetaVec, sigmaBetaVec,alphaBetaVec0,betaBetaVec0,muBetaVec0,sigmaBetaVec0,
						x0Vec, xMat, mux0, Vx0,
						L, yMonth, covEffMonth, yTract){
		
	n = length(CtIds)
	part1 = dnorm(muA,mean=muA0,sd=sigmaA0, log=TRUE) +
			log(densigamma(sigmaA^2, alpha=alphaA0 , beta=betaA0)) +
			sum( dnorm(aVec, mean=muA, sd = sigmaA, log=TRUE) )

	part2 = dnorm(muLambda,mean=muLambda0,sd=sigmaLambda0, log=TRUE)+
		log(densigamma(sigmaLambda^2, alpha=alphaLambda0 , beta=betaLambda0))+
		sum( dnorm(LambdaMat[cbind(as.character(CtIds), as.character(zVec))], mean=muLambda, sd = sigmaLambda, log=TRUE) )

	part3 = 0
	if (TRUE){ 
	part3 = sum( dnorm(as.numeric(etaMat), 0, 1, log=TRUE) ) +
		log(densigamma(sigma0^2, alpha=alpha.epsilon0, beta=beta.epsilon0))
	}
	
	# (25) 
	part4 = dgamma(alpha, shape =alpha.alpha0, rate=beta.alpha0, log=TRUE) +
		lgamma(alpha) + sum(lgamma(table(zVec)))+ K*log(alpha) - lgamma(n+alpha) +
		sum( log(densigamma(RVec, alpha=alpha.R0 , beta=beta.R0)) )

	# (26)
	pBeta=0
	for (r in 1:length(covNames)){
		pBeta = pBeta + sum( dnorm(BetaMat[,r], mean=muBetaVec[r], sd=sigmaBetaVec[r], log=TRUE) )+
				log(densigamma(sigmaBetaVec[r]^2, alpha=alphaBetaVec0[r], beta=betaBetaVec0[r] )) 
	}
	part5 = sum( dnorm(muBetaVec, mean=muBetaVec0, sd=sigmaBetaVec0, log=TRUE) )+pBeta

	
	# (27) conditional on eta
	pX = 0
	for (i in 1:n){
		CtIds.i = CtIds[i]
		k.i = as.character(zVec[i])
		lambda.ik = LambdaMat[as.character(CtIds.i), k.i]
		aVec.i = aVec[as.character(CtIds.i)]
		x0Vec.i = x0Vec[as.character(CtIds.i)]
		xMat.i = xMat[as.character(CtIds.i),] 
		pX.i = 0
		
		pX.i = dnorm(xMat.i[1], mean=aVec.i*x0Vec.i+lambda.ik*etaMat[1, k.i], sd= sigma0,log=TRUE)
		for (tt in 2:length(tmLab)){
			pX.i = pX.i + dnorm(xMat.i[tt], mean=aVec.i*xMat.i[tt-1]+lambda.ik*etaMat[tt, k.i], sd= sigma0, log=TRUE)
		}
		pX = pX + pX.i
	}
	part6 = sum( dnorm(x0Vec, mean=mux0, sd=sqrt(Vx0), log=TRUE) ) + pX

	# (28)
	py = 0
	# new way by tract
	for (i in 1:n){
		CtId.i = CtIds[i]
		obs.i = yTract[[as.character(CtId.i)]]
		xT = xMat[as.character(CtId.i), ]
		id = match(obs.i$MonthTransaction, tmLab)
		xMatch = xT[id] 
		covEffect = as.matrix(obs.i[,covNames])%*% BetaMat[as.character(CtId.i),covNames] 
		residual = obs.i$y - xMatch - as.numeric(covEffect)
		Ri = RVec[as.character(CtId.i)]
		py = py + sum(dnorm(residual, mean=0, sd=sqrt(Ri), log=TRUE))
	}
	part7=py

	logPosterior = part1 + part2 + part3 + part4 + part5 + part6 + part7
	
}
	
getLogPosterior.integrateEtaX <- function(CtIds, zVec, muA, muA0, sigmaA0, sigmaA, alphaA0,betaA0,aVec,
						muLambda,muLambda0,sigmaLambda0,sigmaLambda, alphaLambda0,betaLambda0,LambdaMat,
						sigma0,alpha.epsilon0,beta.epsilon0,
						alpha, alpha.alpha0,beta.alpha0, RVec, alpha.R0, beta.R0,
						covNames, BetaMat, muBetaVec, sigmaBetaVec,alphaBetaVec0,betaBetaVec0,muBetaVec0,sigmaBetaVec0,
						x0Vec, mux0, Vx0,
						L, yMonth, yTract,covEffMonth, details=FALSE, lambdaIsVec = FALSE,
						yDehedonicMonthMean = NULL, sampleMeasure=FALSE, w=NULL, KMax=NULL){
	# (22)	
	n = length(CtIds)
	part1 = dnorm(muA,mean=muA0,sd=sigmaA0, log=TRUE) +
			log(densigamma(sigmaA^2, alpha=alphaA0 , beta=betaA0)) +
			sum( dnorm(aVec, mean=muA, sd = sigmaA, log=TRUE) )
    # (23)
	if (lambdaIsVec){
		lambdaVec = LambdaMat
	} else {
		lambdaVec = LambdaMat[cbind(as.character(CtIds), as.character(zVec))]
	}
	names(lambdaVec) = CtIds

	part2 = dnorm(muLambda,mean=muLambda0,sd=sigmaLambda0, log=TRUE)+
		log(densigamma(sigmaLambda^2, alpha=alphaLambda0 , beta=betaLambda0))+
		sum( dnorm(lambdaVec, mean=muLambda, sd = sigmaLambda, log=TRUE) )
	

	part3 = log(densigamma(sigma0^2, alpha=alpha.epsilon0, beta=beta.epsilon0))
	
	# (25) 
	K = length(unique(zVec))
	if (sampleMeasure){
		alphaw = data.frame(labels = 1:KMax, prior=rep(0, KMax))
		nks = table(zVec)
		countk = data.frame(labels = names(nks), count = as.numeric(nks))
		tab = merge(alphaw, countk, by="labels", all.x=TRUE, sort=TRUE)
		tab$count[is.na(tab$count)] =0
		zCondAlpha = log(ddirichlet(w, alpha=rep(alpha/KMax, KMax))) + 
					dmultinom(x=tab$prior+tab$count, prob=w, log = TRUE)

	} else {
		zCondAlpha = lgamma(alpha) + sum(lgamma(table(zVec)))+ K*log(alpha) - lgamma(n+alpha)
	}
	part4 = dgamma(alpha, shape =alpha.alpha0, rate=beta.alpha0, log=TRUE) +
		 zCondAlpha +
		sum( log(densigamma(RVec, alpha=alpha.R0 , beta=beta.R0)) )

	# (26)
	pBeta=0
	for (r in 1:length(covNames)){
		pBeta = pBeta + sum( dnorm(BetaMat[,r], mean=muBetaVec[r], sd=sigmaBetaVec[r], log=TRUE) )+
				log(densigamma(sigmaBetaVec[r]^2, alpha=alphaBetaVec0[r], beta=betaBetaVec0[r] )) +
				dnorm(muBetaVec[r], mean=muBetaVec0[r], sd=sigmaBetaVec0[r], log=TRUE) 
	}
	part5 = pBeta

	
	# (27)
	kLabels = as.numeric( names(table(zVec)) )
	part6 = sum( dnorm(x0Vec, mean=mux0, sd=sqrt(Vx0), log=TRUE) )

	# (28)
	py = 0
	for (k in kLabels){
		
		CtIds.k = sort(CtIds[zVec==k])
		sigma0diag = matrix(0, length(CtIds.k), length(CtIds.k))
		diag(sigma0diag) = sigma0^2
		
		Q.k = lambdaVec[as.character(CtIds.k)] %o% lambdaVec[as.character(CtIds.k)] + sigma0diag
		aVec.k = aVec[as.character(CtIds.k)]
		RVec.k = RVec[as.character(CtIds.k)]
		
		out = get.Ct.yt(L, yMonth, covEffMonth, CtIds.k, tmLab)
		outFilter = KalmanFilterLikelihood(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, out$tractMatch,out$Ct,out$yt,out$covEfft, tract.i =NULL,mux0, Vx0 )
		
		
		if (FALSE){
		outSuff = get.Ct.yt.suffStat(L, covEffMonth, yDehedonicMonthMean, CtIds.k, tmLab, tract.i=CtIds.k[1])
		outCluster = get.Ct.yt(L, yMonth, covEffMonth, CtIds.k, tmLab)
		outFilter = KalmanFilterLikelihood.suffStat(CtIds.k, tmLab, RVec.k, aVec.k, Q.k, 
							outSuff$Lt, outSuff$yDehedonict.i, outSuff$yDehedonicMean, outSuff$CtBar, 
							outCluster$tractMatch, outCluster$Ct, outCluster$yt, outCluster$covEfft,
							tract.i=NULL, mux0, Vx0 )
		}
		py = py + outFilter$logMargLike.tract.all
	}
	part7=py

	logPosterior = part1 + part2 + part3 + part4 + part5 + part6 + part7
	if (details){
		out = c(part1, part2, part3, part4, part5, part6, part7, logPosterior)
		return(out)
	} else {
		return(logPosterior)
	}
	
}


################################
## Univariate Kalman Smoother ##
################################

KalmanSmoothEM <- function(CtId, toPlot=FALSE, plotDir=NULL, scaleFactor=1,
							pi0 = 0, V0 = 0.1, cleanData = obsClean, tmLab, removeOutlier=TRUE){
	#############################################################
	## Extract transactions and get sufficient statistics of y ##
	#############################################################
	
	if (is.null(CtId)){
		obs = cleanData
	} else {
		obs = cleanData[cleanData$CensusTractRegionID==CtId,]
	}
	obs$y = obs$y*scaleFactor
	obs$y = obs$y - mean(obs$y)
	if (removeOutlier){
		obs = obs[abs(obs$y)<2,]
	}
	
	# Sufficient statistics from y: count Lt, mean ybar, mean ysquare
	obs$count = 1
	Lt = aggregate(count~MonthTransaction,  FUN=sum, data=obs)
	Lt = merge(data.frame(MonthTransaction=tmLab), Lt, by="MonthTransaction", all.x=T, all.y=F, sort=T )
	Lt[is.na(Lt$count),"count"] = 0
	Lt = Lt[,2]

	ybar = aggregate(y~MonthTransaction, FUN=mean, data=obs) 
	ybar = merge(data.frame(MonthTransaction=tmLab), ybar, by="MonthTransaction", all.x=T, all.y=F, sort=T)
	ybar = ybar[,2]

	ysqBar = aggregate(y~MonthTransaction, FUN=function(x){mean(x^2)}, data=obs) 
	ysqBar = merge(data.frame(MonthTransaction=tmLab), ysqBar, by="MonthTransaction", all.x=T, all.y=F, sort=T)
	ysqBar = ysqBar[,2]
	N = sum(Lt)
	T = length(Lt)
	
	## A set of initial values for unknown parameters to start EM 
	if (T){ # 10 sets of initials
		aInitials = c(0.85, 0.95)
		QInitials = c(0.0005, 0.001, 0.005)
		RInitials = 0.2
	}
	if (F){
		aInitials = c(1)
		QInitials = c(1)
		RInitials = c(1)
	}
	logMargLike.max = -10000
	d.max = NULL
	diffInitials = expand.grid(a=aInitials, Q=QInitials, R=RInitials)
	
	for (d in 1:nrow(diffInitials)){
		a = diffInitials$a[d]
		Q = diffInitials$Q[d]
		R = diffInitials$R[d]

		stopCheck = 100
		stopCriteria = 0.00001
		OldPara = rep(0, 3)
		OldlogMargLike = 100
		iter = 1
		paraNormTrack = c()
		logMargLikeTrack = c()
		paraTrack = NULL
		para = c(a, Q, R)

		while (TRUE){
			############
			## E step ##
			############

			####################
			# E1. Kalman filtering #
			####################

			muPredict = VPredict = Kt =muFilter = VFilter = rep(NA, T)
			for (i in 1:T){
				if (i==1) {
					muPredict[i] = a*pi0
				} else {
					muPredict[i] = a*muFilter[i-1] 
				}
				
				if (i==1){
					VPredict[i] = a^2*V0 + Q
				} else {
					VPredict[i] = a^2* VFilter[i-1] + Q
				}
				
				Kt[i] = ifelse(Lt[i]==0, 0, VPredict[i]/(R/Lt[i] +VPredict[i])) 
				muFilter[i] = ifelse(Lt[i]==0, muPredict[i], muPredict[i]+Kt[i]*(ybar[i]-muPredict[i]))
				VFilter[i] = ifelse(Lt[i]==0, VPredict[i], VPredict[i]*(1-Kt[i]))
			}

			####################
			# E2. Kalman smoothing #
			####################

			Jt = muSmooth = VSmooth = rep(NA, T)
			muSmooth[T] = muFilter[T]
			VSmooth[T] = VFilter[T]
			for (i in (T-1):1){
				Jt[i] = a*VFilter[i]/VPredict[i+1]
				muSmooth[i] = muFilter[i] + Jt[i]*(muSmooth[i+1] - muPredict[i+1])
				VSmooth[i] = VFilter[i] + Jt[i]^2*(VSmooth[i+1] - VPredict[i+1])
			}

			#########################################
			# E3. backward recusion for V_{t,t-1|T} #
			#########################################

			VSmoothLag2 = rep(NA, T)
			VSmoothLag2[T] = (1-Kt[T])*a*VFilter[T-1]
			for (i in (T-1):2){
				VSmoothLag2[i] = VFilter[i]*Jt[i-1] + Jt[i]*Jt[i-1]*(VSmoothLag2[i+1]-a*VFilter[i])
			}

			################
			## Evaluation ##
			################
			mu = muSmooth              # t= 1:T
			P  = VSmooth + muSmooth^2  # t = 1:T
			Plag = rep(NA, T)          # t = 2:T
			for (i in 2:T){
				Plag[i] = VSmoothLag2[i] + muSmooth[i]*muSmooth[i-1]
			}
			
			## Oct 27: modified the likelihood, using the i-th iter para and states x.
			summation1 = sum( Lt*(ysqBar-2*mu*ybar+P) , na.rm=T )
			summation2 = sum(P[2:T])-2*a*sum(Plag[2:T])+a^2*sum(P[1:(T-1)])
			# Nov 21: Marginal likelihood of y given parameters, integrating out x
			St = VPredict+R
			#logMargLike = -sum(Lt)/2*log(2*pi) - 1/2*sum(Lt*log(St))-1/2*sum(Lt/St*(ysqBar-2*muPredict*ybar+muPredict^2))
			id.na = is.na(ybar)
			logMargLike = - 1/2*sum((Lt*log(St))[!id.na])-1/2*sum(Lt/St*(ysqBar-2*muPredict*ybar+muPredict^2), na.rm=T)

			# Oct 29. Marginal log likelihood of y
			
			paraNorm = sqrt( sum( (para-OldPara)^2 ) )
			stopCheck = abs( (logMargLike-OldlogMargLike)/OldlogMargLike )

			OldPara = para
			OldlogMargLike = logMargLike
			logMargLikeTrack[iter] = logMargLike
			paraNormTrack[iter] = paraNorm
			paraTrack = rbind(paraTrack, para)
			
			## Check the convergence
			if ( (stopCheck < stopCriteria) & (iter>10)) break
			if (iter >=200) break
			
			iter = iter +1
			
			#######################################################################
			## If not converged, propose the next optimal move of all parameters ##
			#######################################################################
			############
			## M step ##
			############
			
			# M1. update output variance R=rho^2
			R = 1/N*sum(  Lt*(ysqBar-2*mu*ybar+P) , na.rm=T  )
			# M2. update state dynamic parameter a
			a = sum(Plag[2:T]) / sum(P[1:(T-1)])
			# M3. update state noise variance Q=sigma^2
			Q = 1/(T-1)*( sum(P[2:T]) - a*sum(Plag[2:T]) )
			# M4. update initial state mean \pi_1
			#pi1 = muPredict[1] # Current code approved by Emily
			# M5. update initial state variance V_1
			#V1 = VPredict[1] # # Current code approved by Emily
			
			para = c(a, Q, R)
		}
		
		cat("Initial  :", as.numeric( round(diffInitials[d,] ,6)), "|| d =", d, "\n")
		cat("Est  para: ", round(para[1:3] ,6),  "|| d =", d, "\n")
		cat("lMargLike:", round(logMargLike,2), "|| d =", d, "\n")
		cat("Iter:", iter, "|| d =", d, "\n \n")
		
		if ( ((logMargLike> logMargLike.max)&(para[1]>0) ) ){
			d.max= d
			logMargLike.max=logMargLike
			para.max = para
			muSmooth.max = muSmooth
			VSmooth.max = VSmooth
			paraTrack.max = paraTrack 
			logMargLikeTrack.max  = logMargLikeTrack
			paraNormTrack.max = paraNormTrack
		}
	}
	if (is.null(d.max)) return(NULL)
	out = c(muSmooth.max, VSmooth.max, para.max)

	if (TRUE){
		cat("\n ")
		cat("Est  para: ", round(para.max[1:3] ,6), "\n" )
		cat("Initial  :", as.numeric( round(diffInitials[d.max,] ,6)), "|| d =", d.max, "\n")
		cat("log Marg Likelihood =", round(logMargLike.max,2), "\n"  )
		cat("----------------------------", "\n \n")
	}
	
	if (toPlot){
		para = para.max
		paraTrack = paraTrack.max
		muSmooth = muSmooth.max
		VSmooth = VSmooth.max
		paraNames = c("AR(1) coef", "State var", "Response var")
		## final parameters
		names(para) =paraNames
		round(para, 4)

		## plot of likelihood convergence
		if (is.null(plotDir)){
			pdf(paste("Output/univarKalmanSmooth/LikelihoodConverge_",CtId, ".pdf", sep=""), width=12, height=6)
		} else{
			pdf(paste(plotDir,"/LikelihoodConverge_",CtId, ".pdf", sep=""), width=12, height=6)
		}
		par(mfrow=c(1,2))
		plot(logMargLikeTrack.max, xlab="Iterations in EM", ylab="log marginal likelihood", type="l", col=2, lwd=2)
		plot(paraNormTrack.max, xlab="Iterations in EM", ylab="2-Norm of parameters", type="l", col=2, lwd=2)
		dev.off()

		## plot of para convergence
		if (is.null(plotDir)){
			pdf(paste("Output/univarKalmanSmooth/paraConverge_", CtId,".pdf",sep=""), width=10, height=8)
		} else {
			pdf(paste(plotDir,"/paraConverge_", CtId,".pdf",sep=""), width=10, height=8)
		}
		par(mfrow=c(1,1))
		colnames(paraTrack) = paraNames
		require(grDevices)
		matplot(paraTrack, type="b", pch = 15:19, bg=1:5, col=1:5, main="Parameter convergence",
			xlab="Iterations in EM", ylab="Parameter value", cex=0.7, ylim=c(0, 1.2))
		legend("topleft", pch=15:19, col=1:5, legend=paste(paraNames,"=", round(para,4) ), cex=1.2, bty="n" )
		#text(x=iter-7, y=tail(paraTrack,1)[-c(3,5)]-0.05, labels=paste(paraNames, "=", round(tail(paraTrack,1),4 ))[-c(3,5)], col=c(1,2,4))
		#text(x=iter-7, y=tail(paraTrack,1)[3]-0.05, labels=paste(paraNames,  "=", round(tail(paraTrack,1),4 ))[3], col=3)
		#text(x=iter-7, y=tail(paraTrack,1)[5]-0.1, labels=paste(paraNames, "=",  round(tail(paraTrack,1),4 ))[5], col=5)
		dev.off()

		## plot of time series
		if (is.null(plotDir)){
			pdf(paste("Output/univarKalmanSmooth/KS_",CtId,".pdf",sep=""), width=16, height=8)
		} else {
			pdf(paste(plotDir,"/KS_",CtId,".pdf",sep=""), width=16, height=8)
		}
		par(mfrow=c(1,1))
		range.y = range(c(obs$y,  muSmooth - qnorm(0.95)*sqrt(VSmooth), muSmooth + qnorm(0.95)*sqrt(VSmooth)))
		plot(obs$MonthTransaction, obs$y, xaxt="n", xlab="", xlim=range(tmLab),ylim=range.y,
			ylab=expression(paste("log(y "["t,l"], " ) - log(g "["t"]," )" )))
		axis(side=1, at= as.Date(unique(format(tmLab,"%Y-01-31" ))), labels=unique(format(tmLab,"%Y-01")))
		abline(v=as.Date(unique(format(tmLab,"%Y-01-31" ))), col="grey")
		lines(tmLab, ybar, col=3)
		lines(tmLab, muSmooth, col=2, lwd=2)
		info = CtTable[CtTable$CensusTractRegionID== CtId,]
		#title(paste("Tract ", info$CtRegionName, "in", info$CountyName, "County,", round(mean(Lt),1),"trans./month"))
		title(paste("Tract ", info$CtRegionName,":" ,round(mean(Lt),1),"trans/month", "a =", round(para[1],3), "Q =",round(para[2],4), "R =",round(para[3],4)))
		legend("topleft", col=c(3,2), lwd=1:2, legend=c("Monthly average price", "Kalman Smoothed price"), bty="n")
		lines(tmLab, muSmooth - qnorm(0.95)*sqrt(VSmooth), col="darkgrey", lwd=2, lty=1)
		lines(tmLab, muSmooth + qnorm(0.95)*sqrt(VSmooth), col="darkgrey", lwd=2, lty=1)
		
		# polygon(x=c(tmLab, rev(tmLab)), 
			# y=c(muSmooth - qnorm(0.95)*sqrt(VSmooth), rev(muSmooth + qnorm(0.95)*sqrt(VSmooth))), 
			# col="grey", border=NA)
		# lines(tmLab, ybar, col=3)
		# lines(tmLab, muSmooth, col=2, lwd=2)
		# points(obs$MonthTransaction, obs$y)
		dev.off()
	}
	
	return(out)
}


########################################################
## Univariate Kalman Smoother, with covariate effects ##
########################################################

KalmanSmoothEM.CovariateEffects <- function(CtId=NULL, toPlot=FALSE, plotDir=NULL, scaleFactor=1,
							pi0 = 0, V0 = 0.1, cleanData = obsClean, tmLab, removeOutlier=TRUE, covNames,vocal=TRUE){
	#############################################################
	## Extract transactions and get sufficient statistics of y ##
	#############################################################
	
	if (is.null(CtId)){
		obs = cleanData
	} else {
		obs = cleanData[cleanData$CensusTractRegionID==CtId,]
	}
	obs$y = obs$y*scaleFactor
	obs$y = obs$y - mean(obs$y)
	
	if (removeOutlier){
		obs = obs[abs(obs$y)<2,]
	}
	
	# Sufficient statistics from y: count Lt, mean ybar, mean ysquare
	obs$count = 1
	Lt = aggregate(count~MonthTransaction,  FUN=sum, data=obs)
	Lt = merge(data.frame(MonthTransaction=tmLab), Lt, by="MonthTransaction", all.x=T, all.y=F, sort=T )
	Lt[is.na(Lt$count),"count"] = 0
	Lt = Lt[,2]

	ybar = aggregate(y~MonthTransaction, FUN=mean, data=obs) 
	ybar = merge(data.frame(MonthTransaction=tmLab), ybar, by="MonthTransaction", all.x=T, all.y=F, sort=T)
	ybar = ybar[,2]

	ysqBar = aggregate(y~MonthTransaction, FUN=function(x){mean(x^2)}, data=obs) 
	ysqBar = merge(data.frame(MonthTransaction=tmLab), ysqBar, by="MonthTransaction", all.x=T, all.y=F, sort=T)
	ysqBar = ysqBar[,2]
	
	N = sum(Lt)
	T = length(Lt)
	
	## A set of initial values for unknown parameters to start EM 
	if (TRUE){ # 10 sets of initials
		aInitials = c(0.85, 0.95)
		QInitials = c(0.0005, 0.001, 0.005)
		RInitials = 0.2
	}
	if (FALSE){
		aInitials = c(1)
		QInitials = c(1)
		RInitials = c(1)
	}
	logMargLike.max = -10000
	d.max = NULL
	diffInitials = expand.grid(a=aInitials, Q=QInitials, R=RInitials)
	
	for (d in 1:nrow(diffInitials)){
		a = diffInitials$a[d]
		Q = diffInitials$Q[d]
		R = diffInitials$R[d]
		Betas = rep(0, length(covNames))
		names(Betas) = covNames
		
		stopCheck = 100
		stopCriteria = 0.00001
		
		OldlogMargLike = 100
		iter = 1
		paraNormTrack = c()
		logMargLikeTrack = c()
		paraTrack = NULL
		para = c(a, Q, R, Betas)
		OldPara = para

		yDehedonicsBar = ybar
		yDehedonicsSqBar = ysqBar
		
		while (TRUE){
			############
			## E step ##
			############

			####################
			# E1. Kalman filtering #
			####################

			muPredict = VPredict = Kt =muFilter = VFilter = rep(NA, T)
			for (i in 1:T){
				if (i==1) {
					muPredict[i] = a*pi0
				} else {
					muPredict[i] = a*muFilter[i-1] 
				}
				
				if (i==1){
					VPredict[i] = a^2*V0 + Q
				} else {
					VPredict[i] = a^2* VFilter[i-1] + Q
				}
				
				Kt[i] = ifelse(Lt[i]==0, 0, VPredict[i]/(R/Lt[i] +VPredict[i])) 
				muFilter[i] = ifelse(Lt[i]==0, muPredict[i], muPredict[i]+Kt[i]*(yDehedonicsBar[i]-muPredict[i]))
				VFilter[i] = ifelse(Lt[i]==0, VPredict[i], VPredict[i]*(1-Kt[i]))
			}

			####################
			# E2. Kalman smoothing #
			####################

			Jt = muSmooth = VSmooth = rep(NA, T)
			muSmooth[T] = muFilter[T]
			VSmooth[T] = VFilter[T]
			for (i in (T-1):1){
				Jt[i] = a*VFilter[i]/VPredict[i+1]
				muSmooth[i] = muFilter[i] + Jt[i]*(muSmooth[i+1] - muPredict[i+1])
				VSmooth[i] = VFilter[i] + Jt[i]^2*(VSmooth[i+1] - VPredict[i+1])
			}

			#########################################
			# E3. backward recusion for V_{t,t-1|T} #
			#########################################

			VSmoothLag2 = rep(NA, T)
			VSmoothLag2[T] = (1-Kt[T])*a*VFilter[T-1]
			for (i in (T-1):2){
				VSmoothLag2[i] = VFilter[i]*Jt[i-1] + Jt[i]*Jt[i-1]*(VSmoothLag2[i+1]-a*VFilter[i])
			}

			################
			## Evaluation ##
			################
			mu = muSmooth              # t= 1:T
			P  = VSmooth + muSmooth^2  # t = 1:T
			Plag = rep(NA, T)          # t = 2:T
			for (i in 2:T){
				Plag[i] = VSmoothLag2[i] + muSmooth[i]*muSmooth[i-1]
			}
			
			## Oct 27: modified the likelihood, using the i-th iter para and states x.
			#summation1 = sum( Lt*(ysqBar-2*mu*ybar+P) , na.rm=T )
			#summation2 = sum(P[2:T])-2*a*sum(Plag[2:T])+a^2*sum(P[1:(T-1)])
			# Nov 21: Marginal likelihood of y given parameters, integrating out x
			St = VPredict+R
			#logMargLike = -sum(Lt)/2*log(2*pi) - 1/2*sum(Lt*log(St))-1/2*sum(Lt/St*(ysqBar-2*muPredict*ybar+muPredict^2))
			id.na = is.na(ybar)
			logMargLike = - 1/2*sum((Lt*log(St))[!id.na])-1/2*sum(Lt/St*(yDehedonicsSqBar-2*muPredict*yDehedonicsBar+muPredict^2), na.rm=T)

			# Oct 29. Marginal log likelihood of y
			
			paraNorm = sqrt( sum( (para-OldPara)^2 ) )
			stopCheck = abs( (logMargLike-OldlogMargLike)/OldlogMargLike )

			OldPara = para
			OldlogMargLike = logMargLike
			logMargLikeTrack[iter] = logMargLike
			paraNormTrack[iter] = paraNorm
			paraTrack = rbind(paraTrack, para)
			
			## Check the convergence
			if ( (stopCheck < stopCriteria) & (iter>10)) break
			if (iter >=200) break
			
			iter = iter +1
			
			#######################################################################
			## If not converged, propose the next optimal move of all parameters ##
			#######################################################################
			############
			## M step ##
			############
			
			# M1. update output variance R=rho^2
			R = 1/N*sum(  Lt*(yDehedonicsSqBar - 2*mu*yDehedonicsBar + P) , na.rm=T  )
			# M2. update state dynamic parameter a
			a = sum(Plag[2:T]) / sum(P[1:(T-1)])
			# M3. update state noise variance Q=sigma^2
			Q = 1/(T-1)*( sum(P[2:T]) - a*sum(Plag[2:T]) )
			# M4. update initial state mean \pi_1
			#pi1 = muPredict[1] # Current code approved by Emily
			# M5. update initial state variance V_1
			#V1 = VPredict[1] # # Current code approved by Emily
			
			# M6. covariate effects
			id = match(obs$MonthTransaction, tmLab)
			muMatched = mu[id] 
			numCovariates = length(covNames)
			rSeq = sample(1:numCovariates)
			for (r in rSeq){
				U = obs[,covNames[r]]
				covEffect = as.matrix(obs[,covNames[-r]])%*% Betas[covNames[-r]] 
				residual = obs$y - muMatched - as.numeric(covEffect)
				answer = sum(residual*U) / sum(U*U)
				answer = ifelse(is.na(answer),0,answer)
				Betas[covNames[r]]  = answer
			}
			
			covEffect = as.matrix(obs[,covNames])%*% Betas[covNames] 
			obs$yDehedonics = obs$y - covEffect
			
			yDehedonicsBar = aggregate(yDehedonics~MonthTransaction, FUN=mean, data=obs) 
			yDehedonicsBar = merge(data.frame(MonthTransaction=tmLab), yDehedonicsBar, by="MonthTransaction", all.x=T, all.y=F, sort=T)
			yDehedonicsBar = yDehedonicsBar[,2]

			yDehedonicsSqBar = aggregate(yDehedonics~MonthTransaction, FUN=function(x){mean(x^2)}, data=obs) 
			yDehedonicsSqBar = merge(data.frame(MonthTransaction=tmLab), yDehedonicsSqBar, by="MonthTransaction", all.x=T, all.y=F, sort=T)
			yDehedonicsSqBar = yDehedonicsSqBar[,2]
			
			para = c(a, Q, R, Betas)
		}
		
		if (vocal){
			cat("Initial  :", as.numeric( round(diffInitials[d,] ,6)), "|| d =", d, "\n")
			cat("Est  para: ", round(para ,3),  "|| d =", d, "\n")
			cat("lMargLike:", round(logMargLike,2), "|| d =", d, "\n")
			cat("Iter:", iter, "|| d =", d, "\n \n")
		}
		if ( ((logMargLike> logMargLike.max)&(para[1]>0) ) ){
			d.max= d
			logMargLike.max=logMargLike
			para.max = para
			muSmooth.max = muSmooth
			VSmooth.max = VSmooth
			paraTrack.max = paraTrack 
			logMargLikeTrack.max  = logMargLikeTrack
			paraNormTrack.max = paraNormTrack
		}
	}
	
	if (is.null(d.max)) return(NULL)
	out = c(muSmooth.max, VSmooth.max, para.max)

	if (vocal){
		cat("\n ")
		cat("Est  para: ", round(para.max ,6), "\n" )
		cat("Initial  :", as.numeric( round(diffInitials[d.max,] ,6)), "|| d =", d.max, "\n")
		cat("log Marg Likelihood =", round(logMargLike.max,2), "\n"  )
		cat("----------------------------", "\n \n")
	}
	
	if (toPlot){
		para = para.max
		paraTrack = paraTrack.max
		muSmooth = muSmooth.max
		VSmooth = VSmooth.max
		paraNames = c("AR(1) coef", "State var", "Response var")
		## final parameters
		names(para) =paraNames
		round(para, 4)

		## plot of likelihood convergence
		if (is.null(plotDir)){
			pdf(paste("Output/univarKalmanSmooth/LikelihoodConverge_",CtId, ".pdf", sep=""), width=12, height=6)
		} else{
			pdf(paste(plotDir,"/LikelihoodConverge_",CtId, ".pdf", sep=""), width=12, height=6)
		}
		par(mfrow=c(1,2))
		plot(logMargLikeTrack.max, xlab="Iterations in EM", ylab="log marginal likelihood", type="l", col=2, lwd=2)
		plot(paraNormTrack.max, xlab="Iterations in EM", ylab="2-Norm of parameters", type="l", col=2, lwd=2)
		dev.off()

		## plot of para convergence
		if (is.null(plotDir)){
			pdf(paste("Output/univarKalmanSmooth/paraConverge_", CtId,".pdf",sep=""), width=10, height=8)
		} else {
			pdf(paste(plotDir,"/paraConverge_", CtId,".pdf",sep=""), width=10, height=8)
		}
		par(mfrow=c(1,1))
		colnames(paraTrack) = paraNames
		require(grDevices)
		matplot(paraTrack, type="b", pch = 15:19, bg=1:5, col=1:5, main="Parameter convergence",
			xlab="Iterations in EM", ylab="Parameter value", cex=0.7, ylim=c(0, 1.2))
		legend("topleft", pch=15:19, col=1:5, legend=paste(paraNames,"=", round(para,4) ), cex=1.2, bty="n" )
		#text(x=iter-7, y=tail(paraTrack,1)[-c(3,5)]-0.05, labels=paste(paraNames, "=", round(tail(paraTrack,1),4 ))[-c(3,5)], col=c(1,2,4))
		#text(x=iter-7, y=tail(paraTrack,1)[3]-0.05, labels=paste(paraNames,  "=", round(tail(paraTrack,1),4 ))[3], col=3)
		#text(x=iter-7, y=tail(paraTrack,1)[5]-0.1, labels=paste(paraNames, "=",  round(tail(paraTrack,1),4 ))[5], col=5)
		dev.off()

		## plot of time series
		if (is.null(plotDir)){
			pdf(paste("Output/univarKalmanSmooth/KS_",CtId,".pdf",sep=""), width=16, height=8)
		} else {
			pdf(paste(plotDir,"/KS_",CtId,".pdf",sep=""), width=16, height=8)
		}
		par(mfrow=c(1,1))
		range.y = range(c(obs$y,  muSmooth - qnorm(0.95)*sqrt(VSmooth), muSmooth + qnorm(0.95)*sqrt(VSmooth)))
		plot(obs$MonthTransaction, obs$y, xaxt="n", xlab="", xlim=range(tmLab),ylim=range.y,
			ylab=expression(paste("log(y "["t,l"], " ) - log(g "["t"]," )" )))
		axis(side=1, at= as.Date(unique(format(tmLab,"%Y-01-31" ))), labels=unique(format(tmLab,"%Y-01")))
		abline(v=as.Date(unique(format(tmLab,"%Y-01-31" ))), col="grey")
		lines(tmLab, ybar, col=3)
		lines(tmLab, muSmooth, col=2, lwd=2)
		info = CtTable[CtTable$CensusTractRegionID== CtId,]
		#title(paste("Tract ", info$CtRegionName, "in", info$CountyName, "County,", round(mean(Lt),1),"trans./month"))
		title(paste("Tract ", info$CtRegionName,":" ,round(mean(Lt),1),"trans/month", "a =", round(para[1],3), "Q =",round(para[2],4), "R =",round(para[3],4)))
		legend("topleft", col=c(3,2), lwd=1:2, legend=c("Monthly average price", "Kalman Smoothed price"), bty="n")
		lines(tmLab, muSmooth - qnorm(0.95)*sqrt(VSmooth), col="darkgrey", lwd=2, lty=1)
		lines(tmLab, muSmooth + qnorm(0.95)*sqrt(VSmooth), col="darkgrey", lwd=2, lty=1)
		
		# polygon(x=c(tmLab, rev(tmLab)), 
			# y=c(muSmooth - qnorm(0.95)*sqrt(VSmooth), rev(muSmooth + qnorm(0.95)*sqrt(VSmooth))), 
			# col="grey", border=NA)
		# lines(tmLab, ybar, col=3)
		# lines(tmLab, muSmooth, col=2, lwd=2)
		# points(obs$MonthTransaction, obs$y)
		dev.off()
	}
	
	return(out)
}



#CtIds=CtIds.k, toPlot=TRUE, n.iter=200, stopByIter=FALSE,
#plotDir=plotDir, pi1 = rep(0, length(CtIds.k)), 
#V1 = 0.1*diag(rep(1,length(CtIds.k))), cleanData = obs, scaleFactor=1/scaleFactor


MVKalmanSmoothEM <- function(CtIds, toPlot=FALSE, plotDir=NULL, n.iter=200, stopByIter = FALSE,
		pi1 = rep(0, length(CtIds)), V1 = 0.1*diag(rep(1,length(CtIds))), cleanData = obsClean, 
		simu=FALSE, scaleFactor=1){
	#############################################################
	## Extract transactions and get sufficient statistics of y ##
	#############################################################
	
	CtIds = sort(CtIds)
	obs = cleanData[cleanData$CensusTractRegionID %in% CtIds,]
	ord = order(obs$MonthTransaction, obs$CensusTractRegionID, obs$TransactionID)
	obs = obs[ord,]
	# scale the response and de-mean
	obs$y = obs$y*scaleFactor
	for (i in 1:length(CtIds)){
		idct = obs$CensusTractRegionID== CtIds[i]
		obs$y[idct] = obs$y[idct] - mean(obs$y[idct])
	}
	
	# Sufficient statistics from y: count Lt, mean ybar, mean ysquare
	obs$count = 1
	Lt = aggregate(count~CensusTractRegionID+ MonthTransaction,  FUN=sum, data=obs)
	fullGrid = expand.grid(MonthTransaction=tmLab, CensusTractRegionID=CtIds)
	Lt = merge(fullGrid, Lt, by=c("MonthTransaction", "CensusTractRegionID"), 
			all.x=TRUE, all.y=F, sort=TRUE )
	Lt[is.na(Lt$count),"count"] = 0
	LtList = list()
	nk = length(CtIds)
	for (i in 1:nk) { LtList[[i]] = Lt$count[Lt$CensusTractRegionID==CtIds[i]] }
	#Lt = Lt[,"count"]
	
	ybar = aggregate(y~MonthTransaction+CensusTractRegionID, FUN=mean, data=obs)
	ybar = merge(fullGrid, ybar, by=c("MonthTransaction", "CensusTractRegionID"),
		all.x=TRUE, all.y=F, sort=TRUE)
	ybarList = list()
	for (i in 1:nk){ ybarList[[i]] = ybar$y[ybar$CensusTractRegionID==CtIds[i]] }

	ysqBar = aggregate(y~MonthTransaction+CensusTractRegionID, FUN=function(x){mean(x^2)}, data=obs) 
	ysqBar = merge(fullGrid, ysqBar, by=c("MonthTransaction", "CensusTractRegionID"), all.x=TRUE, all.y=F, sort=TRUE)
	ysqBarList = list()
	for (i in 1:nk){ ysqBarList[[i]] = ysqBar$y[ysqBar$CensusTractRegionID==CtIds[i]] }
	
	N = aggregate(count~CensusTractRegionID,FUN=sum, data=Lt) # number of total obs in each tract
	T = length(tmLab) # number of time periods
	nk =length(CtIds) # number of tracts in cluster k

	## given unknown parameters to start EM 
	a = rep(0.99, length(CtIds))
	#a = rep(0.9, length(CtIds)) 
	A = diag(a)
	#Q = diag(rep(0.0001, length(CtIds))) # Q=  state noise variance
	Q = diag(rep(1, length(CtIds)))
	#Q = diag(rep(0.001, length(CtIds)))
	#r = rep(0.1, length(CtIds)) # R= Output variance = diag(r), for r_{i} represent tracts
	r = rep(1, length(CtIds))
	#r = rep(0.05, length(CtIds))
	
	if (F){ # start from the true para
		a=aTrue 
		A = diag(a)
		Q = QTrue
		r= rTrue
	}
	
	if (simu){ # start from near true para
		#set.seed(11)
		a = aInit = aTrue + rnorm(3,0, 0.15)
		A = diag(a)
		Q = QInit = riwish(5, QTrue)
		r = rInit = rTrue+rnorm(3,0,0.5*rTrue)
	}

	## Specify known parameter in the State space model
	Ct = list() # dim is varying by t, dim_t = (mt, nk), the indicator matrix from obs to tract at t
	tractMatch = list()
	for (i in 1:T){
		obs.t = obs[obs$MonthTransaction==tmLab[i],]
		if (nrow(obs.t)==0){
			Ct[[i]] = NULL
			tractMatch[[i]] = NULL
		} else {
			colMatch = match(obs.t$CensusTractRegionID, CtIds)
			tem = matrix(0, nrow(obs.t), nk )
			for (j in 1:nrow(obs.t)){
				tem[j,colMatch[j]] = 1
			}
			Ct[[i]] = tem
			tractMatch[[i]] = colMatch
		}
	}
	
	## Organize the response y
	yt = list() # dim is time varying, dim_t =( mt,1)
	for (i in 1:T){
		obs.t = obs[obs$MonthTransaction==tmLab[i],]
		if (nrow(obs.t)==0){
			yt[[i]] = NULL
		} else {
			yt[[i]] = obs.t$y
		}
	}
	
	stopCheck = 100
	#stopCriteria = 0.001
	#stopCriteria = 0.0001
	stopCriteria = 0.0003

	OldLike = 100
	OldPara = c(a, r, as.numeric(Q))
	iter = 1
	eLogLikeTrack = c()
	paraNormTrack = c()
	paraTrack = NULL

	while (iter<=n.iter){
		############
		## E step ##
		############

		####################
		# E1. Kalman filtering #
		####################

		muPredict  =muFilter  = list()  # dim_t = (nk,1)
		VPredict = VFilter = list() # dim_t =( nk,nk)
		Kt = list() # dim is varying by t, dim_t =(nk, mt), where mt=number of obs at time t for cluster k
		logMargLike = 0
		
		for (i in 1:T){
			# assign time varying SSM parameters
			if (is.null(tractMatch[[i]])){
				Rt.i = NULL
			} else{
				if (length(tractMatch[[i]])==1){
					Rt.i = r[ tractMatch[[i]] ]
					Rinv =  1/r[ tractMatch[[i]] ]
				} else{
					Rt.i = diag(r[ tractMatch[[i]] ])
					Rinv =  diag(1/r[ tractMatch[[i]] ])
				}
			}
			# muPredict = mu_{t|t-1}   t=1, ..., T
			if (i==1) {
				muPredict[[i]] = pi1
			} else {
				muPredict[[i]] = a*muFilter[[i-1]] 
			}
			
			# VPredict = V_{t|t-1}    t=1, ..., T
			if (i==1){
				VPredict[[i]] = V1
			} else {
				VPredict[[i]] = A%*% VFilter[[i-1]] %*% A + Q
			}
			
			# Kt = Kalman forward gain    t=1, ..., T
			#Kt[i] = ifelse(Lt[i]==0, 0, VPredict[i]/(R/Lt[i] +VPredict[i]))
			if (is.null(Ct[[i]])){
				Kt[[i]] =  NULL
				muFilter[[i]] = muPredict[[i]]
				VFilter[[i]]  = VPredict[[i]]
			} else {
				Ctem = Ct[[i]]
				VPtem = VPredict[[i]]
				B = t(Ctem)%*%Rinv
				inn = chol2inv(chol(chol2inv(chol(VPtem))+B%*%Ctem))
				Sinv = Rinv - t(B)%*%inn %*% B
				S = Ctem %*% VPtem %*% t(Ctem) + Rt.i
				Kt[[i]] = VPtem%*%t(Ctem)%*% Sinv
				# muFilter = mu_{t|t}   t=1, ..., T, note ybar[1] has not been used to get the X_1
				residual = yt[[i]]-Ct[[i]]%*%muPredict[[i]]
				muFilter[[i]] = muPredict[[i]] + Kt[[i]]%*%residual
				# VFilter = V_{t|t}   t=1, ..., T
				VFilter[[i]] = VPredict[[i]] - Kt[[i]]%*%Ct[[i]] %*% VPredict[[i]]
			}
			
			# new: log marginal likelihood of y, integrating out x
			if (!is.null(Ct[[i]])){
				resi = residual
				p = length(residual)
				log.det.S = ifelse(p==1, log(S), sum(log(svd(S)$d)))
				logMargLike = logMargLike - p/2*log(2*pi)-1/2*log.det.S-1/2*t(resi) %*% chol2inv(chol(S))%*%resi			
			}
		}

		####################
		# E2. Kalman smoothing 
		####################

		muSmooth = list()  # dim_t = (nk,1)
		VSmooth  = list()  # dim_t =( nk,nk)
		Jt = list()        # dim_t =(nk, nk)
		
		muSmooth[[T]] = muFilter[[T]]
		VSmooth[[T]] = VFilter[[T]]
		
		for (i in (T-1):1){
			# Backward gain, t = T-1, ...., 1
			Jt[[i]] = VFilter[[i]] %*% A %*% chol2inv(chol(VPredict[[i+1]]))
			# muSmooth = mu_{t|T}, t=T, ..., 1 (the next two lne are the same)
		#   muSmooth[[i]] = muFilter[[i]] + Jt[[i]]*(muSmooth[[i+1]] - muPredict[i+1])
			muSmooth[[i]] = muFilter[[i]] + Jt[[i]]%*%(muSmooth[[i+1]] - A %*% muFilter[[i]])
			# VSmooth = V_{t|T}, t=T, ..., 1
			VSmooth[[i]] = VFilter[[i]] + Jt[[i]] %*%(VSmooth[[i+1]] - VPredict[[i+1]]) %*% t(Jt[[i]])
		}

		#########################################
		# E3. backward recusion for V_{t,t-1|T} #
		#########################################

		VSmoothLag2 = list() #dim_t
		VSmoothLag2[[T]] = (diag(nk)-Kt[[T]]%*%Ct[[T]])%*%A%*%VFilter[[T-1]]
		for (i in (T-1):2){
			VSmoothLag2[[i]] = VFilter[[i]]%*%t(Jt[[i-1]]) + Jt[[i]]%*%(VSmoothLag2[[i+1]]-A%*%VFilter[[i]])%*%t(Jt[[i-1]])
		}
		
		mu = muSmooth       # t= 1:T
		P = list() 			# t = 1:T
		for (i in 1:T){ P[[i]] = VSmooth[[i]] + mu[[i]] %*% t(mu[[i]]) }
		Plag = list() 		# t = 2:T
		for (i in 2:T){ Plag[[i]] = VSmoothLag2[[i]] + mu[[i]]%*%t(mu[[i-1]]) }

		################################################################
		## Evaluate expected log complete likelihood conditional on y ## stopped here
		################################################################

		
		## old expected log likelihood
	
		sumErrSqResp = rep(NA, nk)
		L1 = rep(0, nk)
		for (i in 1:nk){
			mu.tract.i = sapply(mu, FUN=function(x){x[i,1]})
			P.tract.i  = sapply(P,  FUN=function(x){x[i,i]})
			sumErrSqResp[i] = sum( LtList[[i]]*(ysqBarList[[i]] - 2*mu.tract.i*ybarList[[i]]+P.tract.i) , na.rm=T  )
			L1[i] = sumErrSqResp[i]/(r[i])
		}
		
		sumErrSqLatent = 0
		for (i in 2:T){
			sumErrSqLatent = sumErrSqLatent + P[[i]] - A%*%t(Plag[[i]])- Plag[[i]] %*% t(A) + A%*%P[[i-1]] %*% t(A)
		}
		
		if (FALSE){
		invQ = chol2inv(chol(Q))
		L2 = rep(0, nk)
		for (i in 2:T){
			L2[i] = sum(as.numeric(P[[i]])*as.numeric(invQ))-2*sum(as.numeric(Plag[[i]])*as.numeric(invQ%*%A))+sum(as.numeric(P[[i-1]])*as.numeric(t(A)%*%invQ%*%A))
		}
		invV1 = chol2inv(chol(V1))
		L3 = sum(as.numeric(P[[1]])*as.numeric(invV1))-2*t(mu[[1]])%*%invV1%*%pi1+ pi1%*%invV1%*%pi1
		
		expectedLogLike = -sum(L1)/2- sum(N$count*log(r))/2 - sum(L2)/2 - (T-1)/2*log(det(Q)) - L3/2 - log(det(V1))/2
		expectedLogLike = as.numeric(expectedLogLike)
		}
		
		para = c(a, r, as.numeric(Q))
		paraNorm = sqrt( sum( (para-OldPara)^2 ) )
		stopCheck = abs( (logMargLike-OldLike)/OldLike )
		OldLike = logMargLike
		OldPara = para

		eLogLikeTrack[iter] = logMargLike
		paraNormTrack[iter] = paraNorm
		paraTrack = rbind(paraTrack, para)
			
		cat("iter", iter, "\n")

		## Check the convergence
		if (!stopByIter) {if ((stopCheck < stopCriteria)  ) break}
			
		iter = iter +1
		############
		## M step ##
		############

		# M1. update output variance R_i for tract i
		r = rep(NA, nk)
		for (i in 1:nk){
			r[i] = sumErrSqResp[i]/N$count[i]
		} ## updated r
		
		# M2. update state dynamic parameter a_i for tract i
		#a = sum(Plag[2:T]) / sum(P[1:(T-1)])
		invQ = chol2inv(chol(Q))
		for (j in 1:nk){
			denom = sapply(P, FUN=function(x){x[j,j]})[-T]
			denominator = sum(denom)*invQ[j,j]
			###(Stopped heare)
			numer1 = sapply(Plag, FUN=function(x){ if (is.null(x)) NA else as.numeric(t(x)[j,] %*% invQ[,j]) } )
			numerator1 = sum(numer1[2:T])
			numer2 = sapply(P, FUN=function(x){ (x[j,-j]*a[-j]) %*% invQ[-j,j] })
			numerator2 = sum(numer2[1:(T-1)])
			a[j] = 	(numerator1 - numerator2)/denominator
		}
		A = diag(a) ## updated A
		
		# M3. update state noise variance Q
		Q = sumErrSqLatent/(T-1) ## updated Q
		
		# M4. update initial state mean \pi_1
		#pi1 = muSmooth[1]
		pi1 = muPredict[[1]] # Current code approved by Emily
		#pi1 = muFilter[1]

		# M5. update initial state variance V_1
		#V1 = VSmooth[1]
		V1 = VPredict[[1]] # # Current code approved by Emily
		#V1 = VFilter[1]
	}
	
	# matrix of the smoothed mu for three tracts
	muSmoothMatMV = matrix(0, T, nk)
	colnames(muSmoothMatMV) = CtIds
	for (j in 1:nk){
		muSmoothMatMV[,j] = sapply(muSmooth, FUN=function(x){x[j,1]} )
	}
	
	# matrix of the smoothed Variance for three tracts
	VSmoothMatMV = matrix(0, T, nk)
	colnames(VSmoothMatMV) = CtIds
	for (j in 1:nk){
		VSmoothMatMV[,j] = sapply(VSmooth, FUN=function(x){x[j,j]} )
	}
	
	# matrix of the smoothed correlation among three tracts
	corrSmoothMatMV = matrix(0, T, 3)
	colnames(corrSmoothMatMV) = c("corr12","corr13","corr23")
	corrSmoothMatMV[,1] = sapply(VSmooth, FUN=function(x){x[1,2]} )/(sqrt(VSmoothMatMV[,1])*sqrt(VSmoothMatMV[,2]))
	corrSmoothMatMV[,2] = sapply(VSmooth, FUN=function(x){x[1,3]} )/(sqrt(VSmoothMatMV[,1])*sqrt(VSmoothMatMV[,3]))
	corrSmoothMatMV[,3] = sapply(VSmooth, FUN=function(x){x[2,3]} )/(sqrt(VSmoothMatMV[,2])*sqrt(VSmoothMatMV[,3]))

	out = list(muSmooth=muSmooth, VSmooth=VSmooth, muSmoothMatMV=muSmoothMatMV,
			VSmoothMatMV=VSmoothMatMV, corrSmoothMatMV=corrSmoothMatMV, para=para)
	aEst = out[['para']][1:nk]; names(aEst) = CtIds
	rEst = out[['para']][(nk+1):(nk+nk)]
	QEst = matrix( tail( out[['para']], nk^2 ), nk,nk  )

	if (toPlot){
		## plot of likelihood convergence
		if (is.null(plotDir)){
			pdf(paste("Output/multivarKalmanSmooth/LikelihoodConverge", ".pdf", sep=""), width=12, height=6)
		} else{
			pdf(paste(plotDir,"/LikelihoodConverge.pdf", sep=""), width=12, height=6)
		}
		par(mfrow=c(1,2))
		plot(eLogLikeTrack, xlab="Iterations in EM", ylab="Expected log likelihood", type="l", col=2, lwd=2)
		plot(paraNormTrack, xlab="Iterations in EM", ylab="2-Norm of parameters", type="l", col=2, lwd=2)
		dev.off()

		## plot of para convergence
		if (is.null(plotDir)){
			pdf(paste("Output/multivarKalmanSmooth/paraConverge",".pdf",sep=""), width=10, height=8)
		} else {
			pdf(paste(plotDir,"/paraConverge.pdf",sep=""), width=10, height=8)
		}
		par(mfrow=c(1,1))
		paraNames = c(paste("a",1:nk,sep=""), paste("r",1:nk,sep=""), paste("Q",c("11","12","13","22","23","33"),sep=""))
		paraTrack = paraTrack[,1:length(paraNames)] 
		colnames(paraTrack) = paraNames
		require(grDevices)
		matplot(paraTrack, type="b", pch = 15:18, col=1:5, main="Parameter convergence",
			xlab="Iterations in EM", ylab="Parameter value", cex=0.7)
		legend("topleft", pch=15:18, col=1:5, legend=paste(paraNames,"=", round(para[1:length(paraNames)],4) ), cex=1.2, bty="n" )
		dev.off()
		
		## plot of time series
		if (is.null(plotDir)){
			pdf(paste("Output/multivarKalmanSmooth/MVkalmanSmoothedSeries",".pdf",sep=""), width=6, height=8)
		} else {
			pdf(paste(plotDir,"/MVkalmanSmoothedSeries.pdf",sep=""), width=6, height=8)
		}
		r = c(-1, 1)
		par(mfrow=c(2,2), oma=c(2,2,5,1), mar=c(2,2,2,2))
		for (i in 1:length(CtIds)){
			dataOneCt = obs[obs$CensusTractRegionID==CtIds[i],]
			plot(dataOneCt$MonthTransaction, dataOneCt$y, xaxt="n", ylab="log(Price)-log(WA)", 
				xlab="", ylim=r, pch=".")
			axis(side=1, at= as.Date(unique(format(tmLab,"%Y-01-31" ))), 
				labels=unique(format(tmLab,"%Y-01")))
			abline(v=as.Date(unique(format(tmLab,"%Y-01-31" ))), col="grey")
			Ct = as.character(CtIds[i])
			arcoef = round(aEst[Ct ],2)
			nObsPerM = round( nrow(dataOneCt)/length(tmLab),1 )
			title(paste("Tract", CtIds[i], ", a =", arcoef, ", #obs/m =", nObsPerM) )
			polygon(x=c(tmLab, rev(tmLab)), 
					y=c(muSmoothMatMV[,Ct] - qnorm(0.95)*sqrt(VSmoothMatMV[,Ct]), 
						rev(muSmoothMatMV[,Ct] + qnorm(0.95)*sqrt(VSmoothMatMV[,Ct]))), 
					col="grey", border=NA)
			points(dataOneCt$MonthTransaction, dataOneCt$y, pch=".")
			lines(tmLab, muSmoothMatMV[,Ct] , col=2)
		}
		dev.off()
		
		## plot of correlation over time 
		if (is.null(plotDir)){
			pdf(paste("Output/multivarKalmanSmooth/corr.pdf",sep=""), width=6, height=8)
		} else {
			pdf(paste(plotDir,"/corr.pdf",sep=""), width=6, height=8)
		}
		rag = range(corrSmoothMatMV)
		par(mfrow=c(3,1), oma=c(2,2,5,1), mar=c(2,2,2,2))
		for (i in 1:3){
			plot(tmLab, corrSmoothMatMV[,i], xaxt="n", ylab="correlation", 
				xlab="", ylim=rag, type="l", col=2)
			axis(side=1, at= as.Date(unique(format(tmLab,"%Y-01-31" ))), 
				labels=unique(format(tmLab,"%Y-01")))
			abline(v=as.Date(unique(format(tmLab,"%Y-01-31" ))), col="grey")
		}
		dev.off()
		
		# image plot of the estimated state covariance matrix Q
		#install.packages("corrplot")
		library(corrplot)
		if (is.null(plotDir)){
			pdf(paste("Output/multivarKalmanSmooth/Q_corrplot",".pdf",sep=""), width=8, height=8)
		} else {
			pdf(paste(plotDir,"/Q_corrplot.pdf",sep=""), width=8, height=8)
		}
		par( mar = par( "mar" ) + c( 2, 4, 0, 0 ) )
		colnames(QEst) = rownames(QEst) = CtIds
		corrplot(cov2cor(QEst),method = "circle" )
		dev.off()
	}
	if (simu) {
		return(list(aEst=aEst, rEst=rEst,QEst=QEst, aInit=aInit,rInit=rInit,QInit=QInit))
	} else {return(out)}
	
}


scaleAHouse <- function(aHouse, tm1, tm2, tm3, qVar1, qVar2, qVar3){
	aHouse$logFinishedSqft = log(aHouse$FinishedSquareFeet)
	aHouse$logLotSize = log(aHouse$LotSizeSquareFeet)
	aHouse$logBathroomCnt = log(aHouse$BathroomCnt)

	aHouse$slogFinishedSqft = scale(aHouse$logFinishedSqft,attributes(tm1)$'scaled:center', attributes(tm1)$'scaled:scale' )
	aHouse$slogLotSize = scale(aHouse$logLotSize, attributes(tm2)$'scaled:center', attributes(tm2)$'scaled:scale')
	aHouse$slogBathroomCnt = scale(aHouse$logBathroomCnt, attributes(tm3)$'scaled:center', attributes(tm3)$'scaled:scale')

	#### Construct piecewise linear splines ##
	aHouse$v1p1 = pmax(0, aHouse$slogFinishedSqft - qVar1[1])
	aHouse$v1p2 = pmax(0, aHouse$slogFinishedSqft - qVar1[2])
	aHouse$v1p3 = pmax(0, aHouse$slogFinishedSqft - qVar1[3])
	aHouse$v1p4 = pmax(0, aHouse$slogFinishedSqft - qVar1[4])

	aHouse$v2p1 = pmax(0, aHouse$slogLotSize - qVar2[1])
	aHouse$v2p2 = pmax(0, aHouse$slogLotSize - qVar2[2])
	aHouse$v2p3 = pmax(0, aHouse$slogLotSize - qVar2[3])
	aHouse$v2p4 = pmax(0, aHouse$slogLotSize - qVar2[4])

	aHouse$v3p1 = pmax(0, aHouse$slogBathroomCnt - qVar3[1])
	aHouse$v3p2 = pmax(0, aHouse$slogBathroomCnt - qVar3[2])
	aHouse$v3p3 = pmax(0, aHouse$slogBathroomCnt - qVar3[3])
	aHouse$v3p4 = pmax(0, aHouse$slogBathroomCnt - qVar3[4])
	return(aHouse)
}
