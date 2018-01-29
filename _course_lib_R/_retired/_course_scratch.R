# source("../_course_lib_R/_course_scratch.R") ; 

rm(list=ls())
source("../_course_lib_R/_course_basic.R")
source("../_course_lib_R/_course_OF_xyz.R")
library(hash) ; library(tools) ; library(utils)  



OFname = "wild" ; BKV   = 67.46773 ; OFpar = 1 ; lowerBnd=-50 ; upperBnd=50 ; coordBest = -15.81515 
signifDigits = 4  ; nBits = 8:10 ; 
tolX = (upperBnd - lowerBnd)/2^nBits ; tolY = 0.5/10^(signifDigits) 

globPar=hash() ; clear(globPar)
globPar[[paste(sep="", OFname, ".BKV")]]          = BKV
globPar[[paste(sep="", OFname, ".tolY")]]         = tolY
globPar[[paste(sep="", OFname, ".isValueOnly")]]  = FALSE
print(globPar)

x = coordBest ; y = wild(x) ; print(cbind(x,y))


return()
x = coordBest ; tolY ; y = wild(x) ; y ; is_a_hit = y <= BKV + tolY ; is_a_hit
