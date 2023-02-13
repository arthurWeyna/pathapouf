library(PoissonMultinomial)
library(matrixStats)
library(extraDistr)

log <- c()
AddToLog <- function(string){
    cat(paste0(string, "\n"))
    log <<- c(log, string)
}

PatHapLogLikSite <- function(X, M, freq, e=0.01, probe=1){
## Computes and outputs Lp(Xi|M)
## X is a vector of allele counts per allele
## M is a vector of haplome counts per group
## freq is a GxA matrix of allele frequencies in each group
	Mtot <- sum(M)
	Nall <- length(X)
	ASLogLiks <- c()
	##probability of each allele combinations in spermathecae	
	freq2 <- t(sapply(rep(1:length(M), times=M), function(x) freq[x,]))
	ASprobs <- dpmd(freq2)
	##loop over all allele combinations and sum loglik
	for(i in 1:length(ASprobs)){
		ASprob=ASprobs[i]
		if(ASprob > 0){
			##recover what is allele combination for one prob in 
			ais <- as.vector(arrayInd(i, dim(ASprobs))-1)	
			ais[Nall] <- Mtot - sum(ais[1:(Nall-1)])
			##Compute allele frequencies 
			fis <- ais/Mtot
			##Compute lp(X|M,AS) + lp(AS|M)
			#No error
			#ASLogLiks[length(ASLogLiks)+1] <- dmultinom(x=X, prob=fis, log=T) + log(ASprob)
			#Error as coumpound multinomial
			#ASLogLiks[length(ASLogLiks)+1] <- logSumExp(c(log(1-e) + dmultinom(x=X, prob=fis, log=T), log(e) + dmultinom(x=X, prob=rep(1,Nall)/Nall, log=T))) + log(ASprob)
			#Error as Dirichlet multinomial
			fis2 <- (fis/e)+(1-fis)^e
			ASLogLiks[length(ASLogLiks)+1] <- ddirmnom(x=X, alpha=fis2, size=sum(X), log=T) + log(ASprob) + log(probe)
		}
	}
	ASLogLiks[is.na(ASLogLiks) | is.nan(ASLogLiks)] <- -Inf
	return(logSumExp(ASLogLiks))
}

PatHapLogLik <- function(Xs, Mmax, freqs, seq_e=c(0.01), gnames=NULL){
## Compute Lp(X|M) for every possible M. (1)
## Sum logliks over each group to get all Lp(X|Mi). (2)
## Computes relative likelihoods on (2). (3)
## Outputs 1 2 3
## Xs is a list of vectors of allele counts per allele
## Mmax is the maximum number of haplomess per group
## freqs is a list of GxA matrices of allele frequencies in each group
	G <- nrow(freqs[[1]]) #get number of groups from freqs, should be constant across freqs
	if(is.null(gnames)){gnames <- paste0("group", 1:G)}
	AddToLog(paste0("Computing log-likelihoods for ", ((Mmax+1)^G)*length(seq_e), " parameter combinations."))
	LogLiks <- array(0, dim=c(rep(Mmax+1, G), length(seq_e))) 
	mis<-rep(0, G)
	pb <- txtProgressBar(min = 2, max = length(LogLiks), style = 3, width = 50, char = "=")
	for(i in 1:length(LogLiks)){
	 	setTxtProgressBar(pb, i)	
		##recover what is the haplome number combination for one Loglik
		mis <- as.vector(arrayInd(i, dim(LogLiks))-1)
		if(all(mis[-length(mis)] == 0)){
			LogLiks[i] <- -Inf #All mi == 0 is impossible
		} else {
			LogLiks[i] <- sum(mapply(Xs, freqs, FUN=function(x,y) PatHapLogLikSite(X=x, M=mis[-length(mis)], freq=y, e=seq_e[mis[length(mis)]+1], probe=1/length(seq_e)), SIMPLIFY=T))
		}
	}
	cat("\n")
	#arrange results into a df
	LogLiksdf <- as.data.frame.table(LogLiks)	
	LogLiksdf[,1:(G+1)] <- arrayInd(1:length(LogLiks), dim(LogLiks))-1
	names(LogLiksdf) <- c(gnames, "e", "LogLik")
	LogLiksdf[,"e"] <- seq_e[LogLiksdf[,"e"]+1]
	RelLiksdf <- LogLiksdf
	RelLiksdf[,G+2] <- round(exp(RelLiksdf[,G+2]-max(RelLiksdf[,G+2])), 5)
	names(RelLiksdf)[G+2] <- "RelLik"
    	if(all(is.infinite(unlist(LogLiksdf["LogLik"])))){
        	AddToLog("All likelihoods are 0. Consider raising maximum haplome number or error rate, or curating the data.")
    	} else {
		bestcomb=LogLiksdf[which.max(LogLiksdf[,G+2]),1:(G+1)]
		AddToLog(paste0("Best parameter combination is ", paste0(names(bestcomb), " : ", bestcomb, collapse=", ")))
	}
	#sum logliks within groups, for each e, and compute relative liks
	LogLiksGrp <- as.data.frame(t(sapply(1:G, function(x) apply(LogLiks, x, function(y) logSumExp(y)))))
	RelLiksGrp <- as.data.frame(t(apply(LogLiksGrp, 1, function(x) round(exp(x-max(x)), 5))))
	LogLiksE <- as.data.frame(t(apply(LogLiks, G+1, function(y) logSumExp(y))))
	RelLiksE <- round(exp(LogLiksE-max(LogLiksE)), 5)
	LogLiksGrp <- cbind(gnames, LogLiksGrp)
	RelLiksGrp <- cbind(gnames, RelLiksGrp)
	names(LogLiksGrp) <- c("group", paste0("m_", 0:Mmax))
	names(RelLiksGrp) <- c("group", paste0("m_", 0:Mmax))
	names(LogLiksE) <- paste0("e_", seq_e)
	names(RelLiksE) <- paste0("e_", seq_e)
	#Print a few interesting facts about results
	bests <- c()
	for(g in 1:G){
        rl <- sort(unlist(RelLiksGrp[g,2:(Mmax+2)]))
        	if(!all(is.nan(rl))){
        		best <- names(rl)[Mmax+1]
			bests <- c(bests, as.numeric(sub("m_", "", best)))
        	    	second <- names(rl)[Mmax]
		    	AddToLog(paste0("Best number of haplomes for ", gnames[g], ": ", sub("m_", "", best), ". Second best ", sub("m_", "", second), " is ", rl[second], " times as probable."))
        	} else {
			AddToLog(paste0("No best number of haplomes was found for ", gnames[g], "..."))
        	}
	} 
	rl <- sort(unlist(RelLiksE[1,1:length(seq_e)]))
        if(!all(is.nan(rl))){
		best <- names(rl)[length(seq_e)]
		bests <- c(bests, as.numeric(sub("e_", "", best)))
		second <- names(rl)[length(seq_e)-1]
		AddToLog(paste0("Best value for e: ", sub("e_", "", best), ". Second best ", sub("e_", "", second), " is ", rl[second], " times as probable."))
	} else {
		AddToLog("No best value for e was found...")
	}
	return(list(LogLiks=LogLiksdf, RelLiks=RelLiksdf, LogLiksGrp=LogLiksGrp, RelLiksGrp=RelLiksGrp, LogLiksE=LogLiksE, RelLiksE=RelLiksE, best=bests))
}

PatHapOptimE <- function(dataset, M){
#Estimate e for each site of a dataset given one M.
        AddToLog(paste0("Running per-site e estimations with ", paste0(dataset$groups, ": ", M, collapse=", ")))
        optims <- list()
        pb <- txtProgressBar(min = 0, max = length(dataset$a), style = 3, width = 50, char = "=")
        for(i in 1:length(dataset$a)){
                setTxtProgressBar(pb, i)
                optims[[i]] <- optim(par=c(e=0.1), fn=PatHapLogLikSite, X=dataset$a[[i]], M=M, freq=dataset$f[[i]], method="Brent", lower=0, upper=1, control=list(fnscale=-1))
        }
        cat("\n")
        es <- t(sapply(optims, function(x) c(x$par, x$convergence)))
        colnames(es) <- c("e", "convergence")
        return(es=es)
}

PatHapSimul <- function(n=c(1,1,1), nallmax=4, nloc=1, covl=10){
##Simulates a dummy dataset suitable for LogLikPatHap()
## For each of nloc loci:
### samples a number of alleles from [2,nallmax]
### generate random allele frequencies within each group (independently)
### samples sum(n) alleles from allele frequencies to "fill" spermathecae, sampling ni alleles from population i with multinomial.
### Sample reads from a multinomial which frequency depend on spermathecae alleles. Number of reads is poisson distributed with lambda = covl
	lf <- list()
	las <- list()
	la <- list()
	for(i in 1:nloc){
		nall <- sample(2:nallmax, 1)
		lf[[i]] <- t(sapply(1:length(n), function(y) {f<-runif(n=nall, min=0, max=1); f/sum(f)}))
		las[[i]] <- rep(0, nall)
		for(j in 1:length(n)){
			las[[i]] <- las[[i]] + as.vector(rmultinom(size=n[j], prob=lf[[i]][j,], n=1))
		}
		la[[i]] <- as.vector(rmultinom(size=rpois(n=1, lambda=covl), prob=las[[i]]/sum(las[[i]]), n=1))
	}	
	return(list(f=lf, as=las, a=la))
}
#simu <- PatHapSimul(n=c(1,2), nallmax=4, nloc=100, cov=10)


PatHapRead <- function(file){
##Reads in PatHapOuf input 
	tab <- read.table(file, sep=" ", header=T, check.names=F)
	if(length(colnames(tab))!=2){stop("Input should have 2 whitespace-separated columns")}
	groups <- strsplit(colnames(tab)[2], ";")[[1]]
	la <- list()
	lf <- list()
	for(i in 1:nrow(tab)){
		la[[i]] <- as.numeric(strsplit(tab[i,1], ",")[[1]])
		lf[[i]] <- do.call("rbind", lapply(strsplit(tab[i,2], ";")[[1]], function(x) as.numeric(strsplit(x, ",")[[1]])))
		if(length(la[[i]]) != ncol(lf[[i]])){stop(paste0("Site ", i, " frequencies do not match number of alleles"))}
	}
	AddToLog(paste0("Data for ", i, " sites, frequencies for ",length(groups), " groups: ",paste0(groups, collapse=", "),"."))
	return(list(groups=groups, a=la, f=lf))
}

args = commandArgs(trailingOnly=TRUE)
##1: path to input
##2: max number of males per group
##3: min. error rate
##4: max. error rate
##5: error rate step
##6: output prefix

AddToLog(paste0("Parameter m_max: ", as.numeric(args[2])))
AddToLog(paste0("Parameter e: from ", as.numeric(args[3]), " to ", as.numeric(args[4]), " by steps of ", as.numeric(args[5])))

#Read input
dataset <- PatHapRead(args[1])
#Run main analysis
results <- PatHapLogLik(Xs=dataset$a, Mmax=as.numeric(args[2]), freqs=dataset$f, gnames=dataset$groups, seq_e=seq(as.numeric(args[3]), as.numeric(args[4]), as.numeric(args[5])))
#Run per-site e estimations
es <- PatHapOptimE(dataset=dataset, M=results$best[-length(results$best)])

prefix <- args[6]
write.table(results$LogLiks, file=paste0(prefix, ".LogLiks.txt"), sep=" ", quote=F, row.names=F)
write.table(results$RelLiks, file=paste0(prefix, ".RelLiks.txt"), sep=" ", quote=F, row.names=F)
write.table(results$LogLiksGrp, file=paste0(prefix, ".LogLiksGrp.txt"), sep=" ", quote=F, row.names=F)
write.table(results$RelLiksGrp, file=paste0(prefix, ".RelLiksGrp.txt"), sep=" ", quote=F, row.names=F)
write.table(results$LogLiksE, file=paste0(prefix, ".LogLiksE.txt"), sep=" ", quote=F, row.names=F)
write.table(results$RelLiksE, file=paste0(prefix, ".RelLiksE.txt"), sep=" ", quote=F, row.names=F)
write.table(es, file=paste0(prefix, ".PerSiteE.txt"), sep=" ", quote=F, row.names=F)
write.table(log, file=paste0(prefix, ".log"), sep=" ", quote=F, row.names=F, col.names=F)


