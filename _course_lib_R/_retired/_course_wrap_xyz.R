# source("../_course_lib_R/_course_wrap_xyz.R") ; wrap_DEoptim(seedInit=-1)

# install.packages(c("DEoptim", "hash", "tools", "utils", "bazar"))

rm(list=ls())
library(hash) 
library(tools) ; library(utils) 
source("../_course_lib_R/_course_OF_xyz.R")

# initialize globPar
globPar=hash()       ; clear(globPar)

wrap_DEoptim = function( sampleSize=1, seedInit=-1, OFname="wildN", OFpar=2, 
                        iterLmt=200, popSize=64, trace=FALSE, VTR=0, tolY=0.005) 
{
  library(DEoptim)
  
  # read the generic table of BKV (best-Known-Value) and related parameters
  fileBKV = "../_course_lib_R/_course_OF_xyz_BKV.txt"
  T = read.table(fileBKV, header=TRUE, comment.char="#", colClasses = "character") ;# print(T)
  idx      = which(T$OFname == OFname)    ;# print(idx)     
  BKV      = as.numeric(T[idx,2])         ;# print(BKV) 
  lowerBnd = rep(as.numeric(T[idx,4]))    ;# print(lowerBnd) 
  upperBnd = rep(as.numeric(T[idx,5]))    ;# print(lowerBnd)   
  
  # assign the values to global hash array for access by OF (Objective Functions)
  globPar[[paste(sep="", OFname, ".BKV")]]          = BKV
  globPar[[paste(sep="", OFname, ".tolY")]]         = tolY
  globPar[[paste(sep="", OFname, ".isValueOnly")]]  = FALSE
  print(globPar) ;# return() 
  
  if (seedInit < 0) {seedInit = round(0.5 + 9999*runif(1))} ; set.seed(seedInit) 

  nPar = length(OFpar)
  if (nPar  > 1)  {
    lowerBnd = rep(lowerBnd, nPar)
    upperBnd = rep(upperBnd, nPar)
  }
  valueTarget = 0 ;# since each OF in _course_OF_xyz.R is in the canonical form
  
  cat("\nseedInit =", seedInit,
      "\nOFname   =", OFname,
      "\nlowerBnd =", lowerBnd,
      "\nupperBnd =", upperBnd,
      "\niterLmt  =", iterLmt,
      "\npopSize  =", popSize,
      "\nBKV      =", BKV,
      "\ntolY     =", tolY,
      "\n")
  
  out = DEoptim(fn=OFname, lower=lowerBnd, upper=upperBnd,
                DEoptim.control(NP=popSize, itermax=iterLmt, trace=trace, VTR=valueTarget))
  return(out)
} # wrap_DEoptim
