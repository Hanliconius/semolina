### written by jjh, v1 completed 11Dec2018, written to automatically calculate the full width half maximum values for peaks in a heterogeneous, oscillating graph (specifically pixel grey values extracted in FIJI or a vector over an EM of the ridges or crossribs of a butterfly scale)
### updated 14Dec2018; incorporated the fourier analysis pipeline to measure avg peak periodicity, plus some tidying. 
### Updated 09Feb2019; incorpoated the 'correct' calculation of Fourier SD also rounded outputs to pixel width
### Updated 11Feb2019; FWHM changed to Full width 2/3 max Original calculation: p <- (k+n)/2 New calculation: p <- (k+n)*(2/3)
##############################
##### SET VARIABLES HERE #####
#### note to user - change only these variables, then from the command line, type "./batch.sh autopeakCaller_R.r /Users/FullPath nameofOutPut.csv" file. 
##############################
args = commandArgs(trailingOnly=TRUE)
print (args[1])
#quit()

#filename <- '001.csv' # set filename
filename <- args[1]
verbose <- TRUE # if set to TRUE, will export extra data from all analyses. 

##Fourier variables
highPassUM <- 2 # every frequency longer than this number will be ignored.
lowPassUM <- 0.7 # every frequency longer than this number will be ignored.

##FWHM variables
silverFactor <- 1 # if the scale is silver, change this to -1 to invert the whole process. 
lowPassFilter <- 0.01 # every value below this number will get filtered out at the end. Pick wisely.  
##############################

## load necessary libraries
library(quantmod)
library(ggplot2)
library(TSA)
library(data.table)

## read CSV. THIS IMPORT LINE WILL REMOVE THE HEADER
trace <- read.csv(filename, h=F, skip=1)

######## Fourier code
### plot the periodogram
p <- periodogram(trace$V2)  

## save periodgram for visual checking
pdf(file=paste(filename, '.fourier.pdf', sep=""), width = 30, height = 7)
p <- periodogram(trace$V2)  
dev.off()

### There should be one or more dominant frequency visible in that graph. Let's extract the dominant frequencies. 
# this line takes the frequency and intensity values out of the plot, and assigns to domFreq
scaleFreq <- data.frame(freq=p$freq, spec=p$spec)
order <- scaleFreq[order(-scaleFreq$spec),]
domFreq <- head(order, 25)

# get pixel width
a <- mean(diff(trace$V1)) 

# then, convert from frequency back to period, scaling by pixel width 'a'
domFreq$size <- (1/domFreq$f)*a

# filter out any with a size larger than 5um
domFreq <- domFreq[domFreq$size < highPassUM,]
domFreq2 <- domFreq[domFreq$size > lowPassUM,]

# take the remaining peak with the largest $spec value, only keeping size, to use as a filter later
q <- head(domFreq2, 1)
q <- q$size

# if verbose = TRUE, then save the table of top five dominant frequencies, post highPassUM filter
if (verbose == TRUE) { write.table(domFreq2, file=paste(filename,'.fourier.tophits.tsv', sep="")) }

######### FHMW code 
### create the modded findPeaks function, taken from user 'stas g' on stackexchange. 
### m can be set to make the procedure more stringent
find_peaks<- function (x, m = 5){
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pks <- sapply(which(shape < 0), FUN = function(i){
       z <- i - m + 1
       z <- ifelse(z > 0, z, 1)
       w <- i + m + 1
       w <- ifelse(w < length(x), w, length(x))
       if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
       pks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) < 1) + 1
    })
     pks <- unlist(pks)
     pks
}

## remove all positions where pixel grey value was 0 (prints as 11 in the fiji output). Multiply by silverFactor, set above, (default =1). The loop can't handle plateaus. 
trace <- trace[(trace$V2 != 11),]*silverFactor

## annotate the valleys in the data. Some important comments:
#### these calculations are neutral to width, we will re-convert to the correct scale at the end
#### the command find_peaks is defined above, and was modified by someone on stackexchange. the negation
#### in front of the objects means it will find the valleys, if you remove the negation it will find peaks instead. 
valleys <- find_peaks(-trace$V2)

##### FOR LOOP, CALCULATES FHMW ##### 
#set index and other dummy values
i <- 1
k <- 0
j <- 0
l <- 0
m <- 0
n <- 0
p <- 0
x1 <- 0
x2 <- 0
# make a data frame for FHMW to populate
FHMW <- data.frame(x1,x2)
for(i in 1:length(valleys)) { try({
    k <- max(trace$V2[valleys[i]:valleys[(i+1)]])
    j <- trace$V1[valleys[i]]
    m <- trace$V1[valleys[(i+1)]]
    n <- trace$V2[valleys[i]]
    l <- trace$V1[i]
    p <- (k+n)*(2/3)
	temptrace <- trace[valleys[i]:valleys[(i+1)],]
	temptrace$V2 <- temptrace$V2-p
	values <- which(sign(temptrace$V2[-1]) != sign(temptrace$V2[-length(temptrace$V2)]))
	x1 <- head(temptrace$V1[values], 1)
	x2 <- tail(values, 1) +1
	x2 <- temptrace$V1[x2]
	x1 <- silverFactor * x1
	x2 <- silverFactor * x2
    FHMW[i,] <- c(x1, x2)
    }, silent = TRUE)}
    
## output a plot to visually check that everything worked. 
pdf(file=paste(filename, '.pdf', sep=""), width = 30, height = 7)
plot(trace, type='l')
abline(v=trace$V1[valleys], col='blue')
abline(v=FHMW$x1, col='green')
abline(v=FHMW$x2, col='red')
dev.off()

## calculate the actual FHMW value, called distance
FHMW$distance <- FHMW$x2 - FHMW$x1
## Low pass filter all the values of $distance that fall under the assigned value 'lowPassFilter' (i'm using ! as a negation in the term). 
FHMW <- FHMW[!(FHMW$distance < lowPassFilter),]
## High pass filter, remove all values where the 
FHMW <- FHMW[!(FHMW$distance > q),]
FHMW <- na.omit(FHMW)
## if verbose = TRUE, output a tsv of x1, x2 and the distance between them
if (verbose == TRUE) {write.table(FHMW, file=paste(filename,'_FHMWtable.tsv', sep=""))}

## bodge calculation of deviance in inter-ridge distance (I don't know how to do stats on fourier transforms, it's very very very complex mathematics....) 
#updated Feb/9 to add rounding   
tmp <- (diff(FHMW$x1))
tmp <- tmp[tmp<2*q]
#tmp <- tmp[(0.5*q)<tmp && tmp<(1.5*q)] 


### the clean, verbose=FALSE export of data
export <- data.frame(0)
export$fourierDist <- q
export$fourierSD <- sd(tmp)
export$FWHMDist <- mean(FHMW$distance)
export$FWHMSD <- sd(FHMW$distance)
export$scaleWidth <- abs(tail(trace$V1,1))
write.table(export, file=paste(filename,'_out.tsv', sep=""))
