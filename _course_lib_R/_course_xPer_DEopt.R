# source("../_course_lib_R/_course_xPer_pi.R") ; xPer_DEoptim()
 
rm(list=ls())
library(plotrix)
source("../_course_lib_R/_course_basic.R")
source("../_course_lib_R/_course_OF_xyz.R")


############## These values are "hard-wired" in this version
OFname = "wild" ; OFpar = 1 ; BKV = 67.467735 ; tolY = 0.005

globPar=hash() ; clear(globPar)
globPar[[paste(sep="", OFname, ".BKV")]]          = BKV
globPar[[paste(sep="", OFname, ".tolY")]]         = tolY
globPar[[paste(sep="", OFname, ".isValueOnly")]]  = FALSE

print(globPar) ; 

# exp_NMOF = function(fName="parab_17_plain", valueTarget=NA, 
#                     lowerLmt=0, upperLmt=100) {}

xPer_DEoptim = function( sampleSize=5, seedInit=-1, OFname="wild", OFpar=1:2, lB=-50, uB=50, 
                         iterLmt=200, NP=64, trace=FALSE, VTR=0, tolY=0.005, asymBounds=c("") )
{
  # Copyright 2018, Franc Brglez
  # experiment with package DEoptim, asympototic variable is either 'OFpar' or 'asymBounds'
  library(DEoptim)
  if (seedInit < 0) {seedInit = round(9999*runif(1, 0, 1))} ;  
  
  if (OFpar  > 1)  {
    lowerBnd = rep(lowerBnd, OFpar)
    upperBnd = rep(upperBnd, OFpar)
  }
  
  fileAsym = "fg_asym_DEoptim_wild.txt"
  header   = file_header(fileAsym)$fileHeader
  fileAsym = file_header(fileAsym)$fileStamped
  fileStamp= file_header(fileAsym)$stamp 
  header = c(
    header,
    paste("# command        = xPer_DEoptim(...)"),
    paste("# sampleSize     =", sampleSize),
    paste("# seedInit       =", seedInit),
    paste("# objectiveFunc  =", OF),
    paste("# lowerBnd       =", lB),
    paste("# upperBnd       =", uB), 
    paste("# populationSize =", NP),
    paste("# bestKnownValue =", VTR),
    paste("# asymPar        =", OFpar),
    paste("# asymBounds     =", asymBounds), 
    "#"
  )
  fileAsymPath = paste(getwd(), "/", fileAsym, sep="")
  write(header, file=fileAsymPath, ncolumns=1, append=FALSE, sep="\t")
  # write column names for all data rows in fileAsym  
  colNames = c("sampleId", "seedInit", "asymPar", "tolY", " "isCensored", "iterLmt", "cntIter", "runtime")
  write(colNames, file=fileAsymPath, ncolumns=length(colNames), append=TRUE, sep="\t")
  
  # start the asymptotic experiments: x = specified OF_par
  sampleId = 0
  for (sigD in signifDigits) {
    # must use the same seedInit for each sigD (can't maintain reproducibility otherwise) 
    seed = seedInit 
    for (sample in seq_len(sampleSize)) {
      
      sampleId = sampleId + 1 ; set.seed(seed)
      times = system.time(
        ##!! MUST USE <- when embedding into system.time !!##
        res <- pi_needles(cntProbeLmt, sigD, needleRatio) 
      )
      # pi_needles: return(list(piMC=piMC, tolR=tolR,  error=error, isCensored=isCensored, 
      #                     numerator=numerator, denominator=denominator, cntProbe=cntProbe))
      row = c(sampleId, seed, sigD,  needleRatio, signif(res$piMC, sigD),
          res$tolR, res$error, res[["isCensored"]],
          res$numerator, res$denominator, res$cntProbe, format(times[["elapsed"]], digits=6)) ;# cat(row, "\n")
      #print(res) ; print(row) ; return()
      write(row, file=fileAsymPath, ncolumns=length(colNames), append=TRUE, sep="\t")
      
      # introduce a new seed for the next sample
      seed = round(1e9*runif(1, 0, 1))
    }    
  }
  cat(sep="", ".. created a table with experimental samples named as = ", fileAsymPath, "\n")
  # convert the file into easier-to-read html file:
  table2html(fileAsymPath) ;# return()
  
  # extract statistics from the data table in fileAsymPath
  file_asym_stats(fileAsymPath,  cntProbeLmt, xName="signifDigits", yNames=c("cntProbe", "runtime"))

  
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
  
} # xPer_pi_needles

xPer_pi_darts = function(sampleSize=5, seedInit=1215, cntProbeLmt=1000, signifDigits=2:4)
{ 
  # Copyright 2018, Franc Brglez
  # experiment with function pi_darts, asympototic variable is 'signifDigits'
  
  if (seedInit < 0) {seedInit = round(9999*runif(1, 0, 1))} ; seed = seedInit 
  
  fileAsym = "fg-asym-pi_darts.txt"
  header   = file_header(fileAsym)$fileHeader
  fileAsym = file_header(fileAsym)$fileStamped
  fileStamp= file_header(fileAsym)$stamp 
  header = c(
    header,
    paste("# command   = xPer_pi_darts(...) \n#")
  )
  fileAsymPath = paste(getwd(), "/", fileAsym, sep="")
  write(header, file=fileAsymPath, ncolumns=1, append=FALSE, sep="\t")
  # write column names for all data rows in fileAsym 
  colNames = c("sampleId", "seedInit", "signifDigits", "piMC",
               "tolRadius", "error", "isCensored", "cntProbe", "runtime")
  write(colNames, file=fileAsymPath, ncolumns=length(colNames), append=TRUE, sep="\t") ;# return()
  
  # start the asymptotic experiments: x = specified significant digits
  sampleId = 0
  for (sigD in signifDigits) {
    
    for (sample in seq_len(sampleSize)) {
      
      sampleId = sampleId + 1
      times = system.time(
        ##!! MUST USE <- when embedding into system.time !!##
        res <- pi_darts(cntProbeLmt, sigD)
      )
      # pi_needles: return(list(piMC=piMC, tolR=tolR,  error=error, 
      #                         isCensored=isCensored, cntProbe=cntProbe))
      row = c(sampleId,  seed, sigD, signif(res$piMC, sigD),
              res$tolR, res$error, res[["isCensored"]],
              res$cntProbe, times[["elapsed"]]) ;# cat(row, "\n")
      write(row, file=fileAsymPath, ncolumns=9, append=TRUE, sep="\t")
      
      # introduce a new seed for the next sample
      seed = round(1e9*runif(1, 0, 1))
    }    
  }
  cat(sep="", ".. created a table with experimental samples named as = ", fileAsymPath, "\n")
  # convert the file into easier-to-read html file:
  table2html(fileAsymPath)
  
  # extract statistics from yhe data table in fileAsymPath
  file_asym_stats(fileAsymPath,  cntProbeLmt, xName="signifDigits", yNames=c("cntProbe", "runtime"))
  
} # xPer_pi_darts



