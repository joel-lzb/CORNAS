# CORNAS: Coverage-dependent RNA-Seq analysis of gene expression data without biological replicates


######## Functions ########
# function to determine coverage and the model's slope:
calcitall <-function(sampleID,sampCount,params) {
	covID = paste(sampleID,"Coverage",sep="")

	# sum the reads and find the largest count:
	totalReads <- sum(sampCount)
	maxSample <- max(sampCount)

	# default population size prior to PCR:
	samplePopSize = 300000000

	# getting the coverage:
	if (length(params[[covID]]) > 0){ # if coverage given for Sample
		sampling_cov = as.numeric(params[[covID]][1])
	} else { # if coverage not given
		sampling_cov = totalReads/samplePopSize
	}
	
	# DEPRECATED: calculating the gradient: y=mx+c
	#m = 1.383529
	#c = -0.025983
	#mgrad = exp(m * sampling_cov + c)
	mgrad = 1/(1-sampling_cov) #hypergeometric func.
	
	# Output
	return (c(sampling_cov,mgrad,maxSample,totalReads))
}


# function that models the posterior distribution given observed count and coverage:
ModelTcD <- function(b,x) {
	# b = coverage
	# x = observed count
	
	# gamma model:
	Gm <- 1/(b -9.065e-17)
	Im <- 1/(0.001909+0.828118*b+2.209732*b^2) 
	Gs <- (1/(0.0004074+0.9553568*b+0.8773089*b^2))^2
	Is <- (1/(0.002587+0.845860*b+1.171028*b^2))^2
	
	mu <- x*Gm + Im
	sigma.sq <- x*Gs + Is

	K <- mu^2 / sigma.sq
	theta <- sigma.sq / mu
	
	return(c(K,theta))
}

# function to begin pair evaluation:
evalPair <- function(geneID,obsA,covA,htA,obsB,covB,htB,tcr,foldc) {
	# get TCR of SampleA:
	SampAout <- getOvlNow("sampleA",obsA,covA,tcr,htA)

	# get TCR of SampleB:
	SampBout <- getOvlNow("sampleB",obsB,covB,tcr,htB)

	# evaluate DEG:
	eval_ovl <- evalDEG(as.numeric(SampAout[1]),as.numeric(SampAout[2]),as.numeric(SampBout[1]),as.numeric(SampBout[2]),as.numeric(foldc))

	# Final output:
	cat(geneID,eval_ovl,obsA,obsB,SampAout,SampBout,"\n",sep="\t")
}

# main function to get TCR:
getOvlNow <- function(sampleID,obs,cov,TCR,ht) {
	# check if there is an overlap available:
	if (length(ht[[as.character(obs)]]) > 0) {
		lowerlim <- ht[[as.character(obs)]][1]
		upperlim <- ht[[as.character(obs)]][2]		
	}else {
		# PART1: Getting the probabilities
		newstuff <- ModelTcD(cov,obs)
		
		# PART2: True Count range determination
		K <- newstuff[1]
		theta <- newstuff[2]
		#find lower and upper limits such that P(l < X < u) = TCR/100
		lowcut <- (1-TCR/100)/2
		highcut <- 1-lowcut
		lowerlim <- round(qgamma(lowcut, shape=K, scale=theta))
		upperlim <- round(qgamma(highcut, shape=K, scale=theta))

		
		# PART3: make corrections:
		# lowerlimits cannot be less than the observed:
		if (lowerlim < obs) {
			lowerlim = obs
		}
		if (upperlim < lowerlim){ #precaution
			upperlim = obs
		}

		# PART4: store new in hash table:
		ht[[as.character(obs)]] <- c(lowerlim,upperlim,K,theta)		
	}
	return(c(lowerlim,upperlim))
}

# function to evaluate DEG boundery overlap and fold change:
evalDEG <- function(lowerA,upperA,lowerB,upperB,foldc) {
	if (lowerB < upperA && upperB > upperA){
		eval_ovl = c("N","-","0")
	}else if (lowerA < upperB && upperA > upperB){
		eval_ovl = c("N","-","0")
	}else if (lowerB >= lowerA && upperA >= upperB){
		eval_ovl = c("N","-","0")
	}else if (lowerA >= lowerB && upperB >= upperA){
		eval_ovl = c("N","-","0")
	}else{
		if (lowerA < upperB){
			foldcheck <- as.numeric(lowerB/upperA)
			if (foldcheck > foldc) {
				eval_ovl = c("Y","B",foldcheck)
			} else {
				eval_ovl = c("N","-",foldcheck)
			}
		}else{
			foldcheck <- as.numeric(lowerA/upperB)
			if (foldcheck > foldc) {
				eval_ovl = c("Y","A",foldcheck)
			} else {
				eval_ovl = c("N","-",foldcheck)
			}
		}
	}
	return(eval_ovl)
}


######## Main ########
cornas <- function(confFile,inFile) {
	# prepare the Tc limit hash tables:
	htA <- new.env() #for SampleA
	htB <- new.env() #for SampleB

	# read config file
	params <- new.env() # a hash for parameters
	
	confilein <- file(confFile, open= "r")

	while (length(oneLine <- readLines(confilein, n = 1, warn = FALSE)) > 0) {
		if (grepl("^#",oneLine) == 0 && grepl("^\n",oneLine) == 0) {
			myVec <- (strsplit(oneLine, ":"))
			paramType <- gsub(" ","",myVec[[1]][1])
			paramData <- gsub(",","",gsub(" ","",myVec[[1]][2]))
			params[[paramType]] <- paramData
		}
	}

	# must have these 3 columns:
	if (length(params[["SampleAcolumn"]]) == 0 || length(params[["SampleBcolumn"]]) == 0 || length(params[["GeneName"]]) == 0 ){
		stop("Please identify the column numbers for each sample and the gene name/id!\n")
	}
	geneID <- as.integer(params[["GeneName"]][1])


	# load input file as a data table
	data1 <- read.table(inFile,header=FALSE) 


	# prepare coverage and slope for Sample A
	colA <- as.integer(params[["SampleAcolumn"]][1])
	sampParamA <- calcitall("SampleA",data1[[colA]],params)
	covA <- sampParamA[[1]]
	slopeA <- sampParamA[[2]]

	cat(paste("Sample_A_Coverage:",sampParamA[[1]]),"\n")
	cat(paste("Sample_A_Slope:",sampParamA[[2]]),"\n")
	cat(paste("Sample_A_Max_Observed_Count:",sampParamA[[3]]),"\n")
	cat(paste("Sample_A_Total_Observed_Count:",sampParamA[[4]]),"\n")


	# prepare coverage and slope for Sample B
	colB <- as.integer(params[["SampleBcolumn"]][1])
	sampParamB <- calcitall("SampleB",data1[[colB]],params)
	covB <- sampParamB[[1]]
	slopeB <- sampParamB[[2]]

	cat(paste("Sample_B_Coverage:",sampParamB[[1]]),"\n")
	cat(paste("Sample_B_Slope:",sampParamB[[2]]),"\n")
	cat(paste("Sample_B_Max_Observed_Count:",sampParamB[[3]]),"\n")
	cat(paste("Sample_B_Total_Observed_Count:",sampParamB[[4]]),"\n")


	# TCR (alpha) set:
	if (length(params[["Alpha"]]) > 0){
		tcr = as.numeric(params[["Alpha"]][1])
	} else {
		tcr = as.numeric(99)
	}
	cat(paste("Alpha_pc_set_at:",tcr),"\n")


	# fold set:
	if (length(params[["Foldthreshold"]]) > 0){
		foldc = as.numeric(params[["Foldthreshold"]][1])
	} else {
		foldc = as.numeric(1.5)
	}
	cat(paste("Fold_threshold_set_at:",foldc),"\n")



	# run the comparisons for all genes:
	cat("***********************************************************\n")
	cat("Gene_Name\tDEG_call\tExpress_higher\tFold_difference\tA_O-count\tB_O-count\tA_T-lower\tA_T-upper\tB_T-lower\tB_T-upper\n")

	# eval all lines of dataset:
	apply(data1, 1, function(x) evalPair(as.character(x[[geneID]]),as.integer(x[[colA]]),as.numeric(covA),htA,as.integer(x[[colB]]),as.numeric(covB),htB,as.integer(tcr),as.numeric(foldc)))

}



######## Parameters ########
args <- commandArgs(TRUE)
arg1 <-args[1]
arg2 <-args[2]

if (is.na(arg1) || is.na(arg2)) { # if loading into R or no argument placed in Rscript
	cat("USAGE 1 (R Script): Rscript CORNAS.R <config> <datatable> \nUSAGE 2 (R Console): cornas(\"/path/to/config\" , \"/path/to/datatable\")\nNote: datatable should have no column headers.\n")
	
} else { # if using Rscript
	cornas(arg1,arg2)
}



## Created by: Joel Low Zi-Bin 20160229 ##
