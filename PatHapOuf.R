library(PoissonMultinomial)
library(matrixStats)
library(extraDistr)

log <- c()
AddToLog <- function(string){
    cat(paste0(string, "\n"))
    log <<- c(log, string)
}

PatHapLogLikSite <- function(X, M, freq, e=0.01){
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
            if(e==0){e <- e+.Machine$double.xmin}
			fis2 <- (fis/e)+(1-fis)^e
			ASLogLiks[length(ASLogLiks)+1] <- ddirmnom(x=X, alpha=fis2, size=sum(X), log=T) + log(ASprob)
		}
	}
	return(logSumExp(ASLogLiks))
}

PatHapLogLik <- function(Xs, Mmax, freqs, e=0.01, gnames=NULL){
## Compute Lp(X|M) for every possible M. (1)
## Sum logliks over each group to get all Lp(X|Mi). (2)
## Computes relative likelihoods on (2). (3)
## Outputs 1 2 3
## Xs is a list of vectors of allele counts per allele
## Mmax is the maximum number of haplomess per group
## freqs is a list of GxA matrices of allele frequencies in each group
	G <- nrow(freqs[[1]]) #get number of groups from freqs, should be constant across freqs
	if(is.null(gnames)){gnames <- paste0("group", 1:G)}
	cat("Computing log-likelihoods.\n")
	LogLiks <- array(log(1/((Mmax+1)^G-1)), dim=rep(Mmax+1, G)) #uniform prior on number of haplomes per group between 0 and Mmax
	LogLiks[1] <- -Inf #prob 0 for no haplomes
	mis<-rep(0, G)
	pb <- txtProgressBar(min = 2, max = length(LogLiks), style = 3, width = 50, char = "=")	
	for(i in 2:length(LogLiks)){
	 	setTxtProgressBar(pb, i)	
		##recover what is the haplome number combination for one Loglik
		mis <- as.vector(arrayInd(i, dim(LogLiks))-1)
		LogLiks[i] <- LogLiks[i] + sum(mapply(Xs, freqs, FUN=function(x,y) PatHapLogLikSite(X=x, M=mis, freq=y, e=e), SIMPLIFY=T))
	}
	cat("\n")
	#arrange results into a df
	LogLiksdf <- as.data.frame.table(LogLiks)	
	LogLiksdf[,1:G] <- arrayInd(1:length(LogLiks), dim(LogLiks))-1
	names(LogLiksdf) <- c(gnames, "LogLik")
    if(all(is.infinite(unlist(LogLiksdf["LogLik"])))){
        AddToLog("The likelihoods of all male number combinations are 0. Consider raising maximum haplome number or curating the data.")
    }
	#sum logliks within groups and compute relative liks
	LogLiksGrp <- as.data.frame(t(sapply(1:G, function(x) apply(LogLiks, x, function(y) logSumExp(y)))))
	RelLiksGrp <- as.data.frame(t(apply(LogLiksGrp, 1, function(x) round(exp(x)/max(exp(x)), 5))))
	RelLiksGrp <- as.data.frame(t(apply(LogLiksGrp, 1, function(x) round(exp(x-max(x)), 5))))
	LogLiksGrp <- cbind(gnames, LogLiksGrp)
	RelLiksGrp <- cbind(gnames, RelLiksGrp)
	names(LogLiksGrp) <- c("group", paste0("n_", 0:Mmax))
	names(RelLiksGrp) <- c("group", paste0("n_", 0:Mmax))
	for(g in 1:G){
        rl <- sort(unlist(RelLiksGrp[g,2:(Mmax+2)]))
        if(!all(is.nan(rl))){
            best <- names(rl)[Mmax+1]
            second <- names(rl)[Mmax]
		    AddToLog(paste0("Best number of haplomes for ", gnames[g], " is ", sub("_", "=", best), ". Second best ", sub("_", "=", second), " is ", rl[second], " times as probable."))
        } else {
		    AddToLog(paste0("No best number of haplomes were found for ", gnames[g], "..."))
        }
	} 
	return(list(LogLiks=LogLiksdf, LogLiksGrp=LogLiksGrp, RelLiksGrp=RelLiksGrp))
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
	AddToLog(paste0("Found data for ", i, " sites, frequencies for ",length(groups), " groups."))
	return(list(groups=groups, a=la, f=lf))
}

PatHapWrite <- function(results, prefix){
	write.table(results$LogLiks, file=paste0(prefix, ".LogLiks.txt"), sep=" ", quote=F, row.names=F)
	write.table(results$LogLiksGrp, file=paste0(prefix, ".LogLiksGrp.txt"), sep=" ", quote=F, row.names=F)
	write.table(results$RelLiksGrp, file=paste0(prefix, ".RelLiksGrp.txt"), sep=" ", quote=F, row.names=F)
	write.table(log, file=paste0(prefix, ".log"), sep=" ", quote=F, row.names=F, col.names=F)
}


args = commandArgs(trailingOnly=TRUE)
##1: path to input
##2: max number of males per group
##3: error rate
##4: output prefix


dataset <- PatHapRead(args[1])
results <- PatHapLogLik(Xs=dataset$a, Mmax=as.numeric(args[2]), freqs=dataset$f, gnames=dataset$groups, e=as.numeric(args[3]))
PatHapWrite(results, args[4])

