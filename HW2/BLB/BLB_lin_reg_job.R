#Code for problem 1

mini <- FALSE

#============================== Setup for running on Gauss... ==============================#

args <- commandArgs(TRUE)

cat("Command-line arguments:\n")
print(args)

####
# sim_start ==> Lowest possible dataset number
###

###################
sim_start <- 999
###################

if (length(args)==0){
  sim_num <- sim_start + 1
  set.seed(121231)
} else {
  # SLURM can use either 0- or 1-indexing...
  # Lets use 1-indexing here...
  sim_num <- sim_start + as.numeric(args[1])
  #sim_seed <- (762*(sim_num-1) + 121231)
}

cat(paste("\nAnalyzing dataset number ",sim_num,"...\n\n",sep=""))

# Find r and s indices:
trial = sim_num/50
s_index = floor(trial) - 19
r_index = (trial - floor(trial))*50+1
set.seed(762*s_index)  #Set seed for each unique s
#============================== Run the simulation study ==============================#

# Load packages:
library(BH)
library(bigmemory.sri)
library(bigmemory)
library(biganalytics)

# I/O specifications:
datapath <- "/home/pdbaines/data"

# mini or full?
if (mini){
  rootfilename <- "blb_lin_reg_mini"
} else {
  rootfilename <- "blb_lin_reg_data"
}


infilename <- paste0(rootfilename,".txt")
backingfilename <- paste0(rootfilename,".bin")
descriptorfilename <- paste0(rootfilename,".desc")

# Set up I/O stuff:
infile <- paste(datapath,infilename,sep="/")
backingfile <- paste(datapath,backingfilename,sep="/")
descriptorfile <- paste(datapath,descriptorfilename,sep="/")


# Attach big.matrix :
dat <- attach.big.matrix(dget(descriptorfile),backingpath=datapath)

#Number of rows and columns
n = nrow(dat); d = ncol(dat)-1;

#Specifications:
s = 5; r = 50; gamma = 0.7; b = n^gamma

#Sample b rows
samples = dat[sample(1:n,b),]

#Reset seed number for new r value
set.seed(762*s_index+r_index)

#Bootstrap n samples
bootstrap = rmultinom(1,size = n, prob = rep(1/b,b))

#Fit linear model
model = lm(samples[,d+1] ~ samples[,1:d]-1, weights = bootstrap)

#Store output
outfile = paste0("output/","coef_",sprintf("%02d",s_index),"_",sprintf("%0.0f",r_index),".txt")
write.table(model$coefficients, file = outfile, row.names = FALSE, col.names = FALSE, sep = ",")
