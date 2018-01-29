# source("../_course_lib_R/_course_wrap_xyz_tmp.R") ; test_of_wild(c(-15.51009, -16.26599, -15.66054)) ; test_of_wild()
# source("../_course_lib_R/_course_wrap_xyz_tmp.R") ; summary(wrap.DEoptim(seedInit=-1, OFname="wild", popSize=64, iterLmt=200, OFpar=2, tolR=5e-3))

# invoke "install.packages(...) before sourcing this file!!!
# install.packages(c("DEoptim", "hash", "tools", "utils", "bazar"))

rm(list=ls())
library(hash) ; library(tools) ; library(utils) 

############## These values are "hard-wired for this file only"
OFname = "wild" ; OFpar = 2 ; BKV = 67.467735 ; tolR = 0.005
globPar=hash() ; clear(globPar)
globPar[[paste(sep="", OFname, ".BKV")]]  = BKV
globPar[[paste(sep="", OFname, ".tolR")]] = tolR

print(globPar) ; 

test_of_wild = function(x = c(-15.81559, -15.81491, -15.81529)) {
  OFname = "wild"
  cat("\n.. from test_of_wild: testing objective function with name =", OFname, "\n")
  print(globPar)
  y = wild(x) 
  print(cbind(x,y))
  return(y)
}

wrap.DEoptim = function(seedInit=-1, OFname="wild", lowerBnd=-50, upperBnd=50,  
                        popSize=64, iterLmt=100, trace=FALSE, OFpar=2, tolR=5e-3) 
{
  library(DEoptim)
  if (seedInit < 0) {seedInit = round(0.5 + 9999*runif(1))} ; set.seed(seedInit) 
  
  if (OFpar  > 1)  {
    lowerBnd = rep(lowerBnd, OFpar)
    upperBnd = rep(upperBnd, OFpar)
  }
  cat("\nseedInit =", seedInit,
      "\nlowerBnd =", lowerBnd,
      "\nupperBnd =", upperBnd,
      "\niterLmt  =", iterLmt,
      "\npopSize  =", popSize,
      "\n")
      
  valueTarget = 0 ;# see the normalization of OFname="wild" in this file or under _course_OF_xyz.R
  out = DEoptim(fn=OFname, lower=lowerBnd, upper=upperBnd, 
                DEoptim.control(NP=popSize, itermax=iterLmt, trace=trace, VTR=valueTarget)) 
  return(out) 
} ;# wrap.DEoptim

wild = function(x) { 
  # Copyright 2018: an example from the R package DEoptim-2.2-4
  # Extension by Franc Brglez: globPar=hash()
  BKV = globPar[["wild.BKV"]]  ;# Best-Known-Value
  eps = globPar[["wild.tolR"]] 
  # x = c(x1, x2, ... , xd)
  d   = length(x)
  y = sum(10 * sin(0.3 * x) * sin(1.3 * x^2) + 0.00001 * x^4 + 0.2 * x + 80 ) ;# print(y)
  if (y < BKV ) {
    cat("\n** WARNING from function wild**",
        "\n    current BKV =", BKV, "  is NOT the best-known-value",
        "\n   better value =", signif(y, 9), "  has been found for x =", signif(x, 9),
        "\n")
    y = NA
  } else if (y - d*BKV <= eps) {
    y = 0
  } else {
    y = y - d*BKV  
  }
  return(y)
  
  # VIGNETTE
  # for d=1, this function is an example from the R package DEoptim-2.2-4
  # it seems to have a unique best-known-value (BKV = 67.46773) solution at 
  # x* = -15.81515 not only within [-50, 50] but also up to [-1e5, 1e5]
  # seedInit = 5381 ; NP = 64 ; interval = -+ 50000 ; itermax = 200 
  # > wild(-15.81515537)
  # [1] 67.46773489
} # wild

