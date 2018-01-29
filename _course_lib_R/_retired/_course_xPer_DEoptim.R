# source("../_course_lib_R/_course_xPer_DEoptim.R") ; xPer_DEoptim()
 
rm(list=ls())
library(plotrix)
source("../_course_lib_R/_course_basic.R")
source("../_course_lib_R/_course_OF_xyz.R")

# initialize globPar
globPar=hash()       ; clear(globPar)

# exp_NMOF = function(fName="parab_17_plain", valueTarget=NA, 
#                     lowerLmt=0, upperLmt=100) {}

xPer_DEoptim = function( sampleSize=5, seedInit=-1, OFname="wildN", OFpar=1:2, 
                         iterLmt=200, popSize=64, trace=FALSE, VTR=0, tolY=0.005, asymBounds=c("") )
{
  # Copyright 2018, Franc Brglez
  # experiment with package DEoptim, asympototic variable is either 'OFpar' or 'asymBounds'
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
  #print(globPar) ;# return() 
  
  if (seedInit < 0) {seedInit = round(0.5 + 9999*runif(1))} ; set.seed(seedInit) 
  
  nPar = length(OFpar)
  if (nPar  > 1)  {
    lowerBnd = rep(lowerBnd, nPar)
    upperBnd = rep(upperBnd, nPar)
  }
  valueTarget = 0
  
  fileAsym = "fg_asym_DEoptim_wild.txt"
  header   = file_header(fileAsym)$fileHeader
  fileAsym = file_header(fileAsym)$fileStamped
  fileStamp= file_header(fileAsym)$stamp 
  header = c(
    header,
    paste("# command        = xPer_DEoptim(...)"),
    paste("# sampleSize     =", sampleSize),
    paste("# seedInit       =", seedInit),
    paste("# objectiveFunc  =", OFname),
    paste("# lowerBnd       =", lowerBnd),
    paste("# upperBnd       =", upperBnd), 
    paste("# populationSize =", popSize),
    paste("# bestKnownValue =", VTR),
    paste("# asymPar        =", OFpar),
    paste("# asymBounds     =", asymBounds), 
    "#"
  )
  fileAsymPath = paste(getwd(), "/", fileAsym, sep="")
  write(header, file=fileAsymPath, ncolumns=1, append=FALSE, sep="\t")
  # write column names for all data rows in fileAsym  
  colNames = c("sampleId", "seedInit", "asymPar", "error", "tolY", "isCensored", "iterLmt", "iterCnt", "runtime")
  write(colNames, file=fileAsymPath, ncolumns=length(colNames), append=TRUE, sep="\t") ;# return()
  
  # start the asymptotic experiments: x = specified OF_par
  sampleId = 0
  for (nDim in seq_len(length(OFpar))) {  
    # must use the same seedInit for each nDim (can't maintain reproducibility otherwise) 
    seed = seedInit 
    for (sample in seq_len(sampleSize)) {
      
      sampleId = sampleId + 1 ; set.seed(seed)
      times = system.time(
        ##!! MUST USE <- when embedding into system.time !!##
        out <- DEoptim(fn=OFname, lower=lowerBnd, upper=upperBnd,
                      DEoptim.control(NP=popSize, itermax=iterLmt, trace=trace, VTR=valueTarget))
      )      
      #print(cbind(ndim, sampleId, times))
      #print(out)
      valueBest = out$optim$bestval
      if (valueBest == 0) {
        isCensored = FALSE
      } else {
        isCensored = TRUE
      }
      error   = globPar[["wildN.error"]]
      iterCnt = out$optim$iter
      row = c(sampleId, seed, nDim, error, tolY, isCensored, iterLmt, iterCnt, format(times[["elapsed"]], digits=6)) ;# cat(row, "\n")
      #print(res) ; print(row) ; return()
      write(row, file=fileAsymPath, ncolumns=length(colNames), append=TRUE, sep="\t")
      
      # introduce a new seed for the next sample
      seed = round(1e9*runif(1, 0, 1))
    }    
  }
  cat(sep="", ".. created a table with experimental samples named as = ", fileAsymPath, "\n")
  # convert the file into easier-to-read html file:
  table2html(fileAsymPath) 
  
  # extract statistics from the data table in fileAsymPath
  file_asym_stats(fileAsymPath,  cntProbeLmt, xName="signifDigits", yNames=c("iterCnt", "runtime"))

  # VIGNETTE 
  # ############### measuring runtime
  # times = system.time(Sys.sleep(10)) 
  # times ; times[["elapsed"]] ; times[3] ; as.numeric(times[3])
  # user  system elapsed 
  # 0.012   0.003  10.014 
  # [1] 10.014
  # elapsed 
  # 10.014 
  # [1] 10.014
  # 
  
} # xPer_DEoptim

