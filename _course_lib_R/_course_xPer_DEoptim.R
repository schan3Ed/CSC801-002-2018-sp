# source("../_course_lib_R/_course_xPer_DEoptim.R") ; xPer_DEoptim(sampleSize=5, seedInit=873, tolY=0.005, OFpar=1:2)
 
# install.packages(c("DEoptim", "hash", "tools", "utils", "bazar"))
rm(list=ls()) # remove all variable values!

#library(tools) ; library(utils) ; 
library(hash)  ; library(plotrix)
source("../_course_lib_R/_course_basic.R")
source("../_course_lib_R/_course_OF_xyz.R")

# initialize globPar
globPar=hash()       ; clear(globPar)

# exp_NMOF = function(fName="parab_17_plain", valueTarget=NA, 
#                     lowerLmt=0, upperLmt=100) {}

xPer_DEoptim = function( sampleSize=5, seedInit=873, OFname="wildN", OFpar=1, 
                         iterLmt=200, popSize=NA, trace=FALSE, solverTargetY=0, tolY=0.005, asymBounds=NULL)
{
  # Copyright 2018, Franc Brglez
  # experiment with package DEoptim, asympototic variable is either 'OFpar' or 'asymBounds'
  library(DEoptim) ; solver = "DEoptim" ; thisFunction = "xPer_DEoptim"
  
  # Warning in DEoptim(fn = OFname, lower = lowerBnd, upper = upperBnd, DEoptim.control(NP = popSize, ...
  # For many problems it is best to set 'NP' (in 'control') to be at least ten times the length of the parameter vector. 
  if (is.na(popSize)) {
    popSize = 20*length(OFpar)
  } else {
    if (popSize < 10*length(OFpar)) {
      cat("** WARNING from function =", thisFunction, 
         "\npopSize =", popSize, " is below minimum recommended value of 10",
         "\n... proceeding with value of popSize =", 20*length(OFpar), 
         "\n\n")
      popSize = 20*length(OFpar)
    } 
  }
  #print(popSize) ; return()
  
  # read the generic table of BKV (best-Known-Value) and related parameters
  fileBKV = "../_course_lib_R/_course_OF_xyz_BKV.txt"
  T = read.table(fileBKV, header=TRUE, comment.char="#", colClasses = "character") ;# print(T)
  idx      = which(T$OFname == OFname)    ;# print(idx)     
  BKV      = as.numeric(T[idx,2])         ;# print(BKV) 
  lowerBnd = as.numeric(T[idx,4])    ;# print(lowerBnd) 
  upperBnd = as.numeric(T[idx,5])    ;# print(lowerBnd)  
  BKCoord  = as.numeric(T[idx,6])    ;# print(BKCoord) ; print(rep(BKCoord, 2)) ; return()
  rangeBnd = upperBnd - lowerBnd
  
  # assign the values to global hash array for access by OF (Objective Functions)
  globPar[[paste(sep="", OFname, ".BKV")]]          = BKV
  globPar[[paste(sep="", OFname, ".tolY")]]         = tolY
  globPar[[paste(sep="", OFname, ".isValueOnly")]]  = FALSE
  print(globPar) ;# return() 
  
  # initialize RNG with seedInit
  if (seedInit < 0) {seedInit = round(0.5 + 9999*runif(1))} ; set.seed(seedInit) 
 
  if (length(asymBounds) == 0) {
    # run the asymptotic experiment with OFpar  and write results to fileAsym
    # NOTE: valueTarget = 0 since each OF in _course_OF_xyz.R is in the canonical form!!!
    
    if (seedInit < 0) {seedInit = round(0.5 + 9999*runif(1))} ; set.seed(seedInit) 
    
    # create header  
    fileAsym = paste("fg_asym_", solver, "_", OFname, "_nPar.txt", sep="")
    header   = file_header(fileAsym)$fileHeader
    fileAsym = file_header(fileAsym)$fileStamped
    fileStamp= file_header(fileAsym)$stamp 
    
    sizeOFpar   = length(OFpar)
    BKCoord_all = rep(BKCoord, sizeOFpar)
    header = c(
      header,
      paste("# command        = xPer_DEoptim(...)"),
      paste("# sampleSize     =", sampleSize),
      paste("# seedInit       =", seedInit),
      paste("# objectiveFunc  =", OFname),
      paste("# OFpar          =", paste(OFpar, collapse=",")),
      paste("# lowerBnd       =", paste(rep(lowerBnd, sizeOFpar), collapse=",")),
      paste("# upperBnd       =", paste(rep(upperBnd, sizeOFpar), collapse=",")),
      paste("# populationSize =", popSize),
      paste("# iterLmt        =", iterLmt),
      paste("# tolY           =", formatC(tolY, format="e", digits=3)),
      paste("# BKV            =", BKV),
      paste("# BKCoord        =", BKCoord),
      paste("# BKV_all        =", BKV*sizeOFpar),
      paste("# BKCoord_all    =", paste(BKCoord_all, collapse=",")),
      "#")
    cat(paste(header, collapse="\n")) ;# return()

    fileAsymPath = paste(getwd(), "/", fileAsym, sep="")
    write(header, file=fileAsymPath, ncolumns=1, append=FALSE, sep="\t")
    # write column names for all data rows in fileAsym  
    colNames = c("sampleId", "seedInit", "xAymp", "errorY", "tolY", "isCensored", "popSize", 
                 "iterLmt", "iterCnt", "runtime", "errorX", "coordBest")
    write(colNames, file=fileAsymPath, ncolumns=length(colNames), append=TRUE, sep="\t") ;# return()
    
    # start the asymptotic experiment with OFpar and write it to fileAsym
    sampleId = 0 ; xAsymp = OFpar ; cat("\nxAsymp =", xAsymp, "\n\n")
    
    for (nPar in seq_len(sizeOFpar)) {  
      # must use the same seedInit for each nPar (can't maintain reproducibility otherwise) 
      seed = seedInit 
      if (nPar > 1) {
        lowerBnd = c(lowerBnd, lowerBnd)
        upperBnd = c(upperBnd, upperBnd)
        BKCoord = c(BKCoord, BKCoord) 
        popSize = c(popSize, popSize) 
      }
    
      for (sample in seq_len(sampleSize)) {
        
        sampleId = sampleId + 1 ; set.seed(seed) ; cat("\nnPar,sampleId = ", nPar, sampleId, "\n")
        #cat("\nupperBnd,popSize = ", upperBnd, popSize, "\n")
        times = system.time(
          ##!! MUST USE <- when embedding into system.time !!##
          out <- DEoptim(fn=OFname, lower=lowerBnd, upper=upperBnd,
                         DEoptim.control(NP=popSize, itermax=iterLmt, trace=trace, VTR=solverTargetY))
        )      
        # !!VERIFY the solution!! 
        globPar[["wildN.isValueOnly"]] = TRUE
        valueTarget = BKV*nPar
        coordBest   = as.numeric(out$optim$bestmem)
        valueBest    = wildN(coordBest)
        # we are computing relative error since valueTarget >= 1
        errorY      = formatC( abs(1 - valueBest/valueTarget), format="e", digits=3)
        errorsX     = formatC( BKCoord - coordBest, format="e", digits=5)
        errorX      = formatC( distance_xyz(BKCoord, coordBest), format="e", digits=5)
        #...
        xBits = ceiling(log2(rangeBnd/as.numeric(errorX))) ; # do NOT use xBits = log2(rangeBnd/errorX)
        
        
        solverBest = out$optim$bestval
        if (solverBest == solverTargetY ) {
          isCensored = FALSE
        } else {
          isCensored = TRUE
        }
        iterCnt = out$optim$iter
        cat("valueTarget =", valueTarget, "\nvalueBest =", valueBest,  "\noordBest =", coordBest,
            "\nerrorsX =", errorsX, "\nerrorX =", errorX, "\nerrorY =", errorY, "\nxBits =",  xBits, "\nisCensored =", isCensored, "\n")
        row = c(sampleId, seed, nPar, errorY, tolY, isCensored, iterLmt, iterCnt, 
                format(times[["elapsed"]], digits=6),  errorX, paste(signif(coordBest, 7), collapse=",")) ;  cat("row=", row, "\n")
        write(row, file=fileAsymPath, ncolumns=length(colNames), append=TRUE, sep="\t")
        
        # create a new seed for the next sample and unset globPar[["wildN.isValueOnly"]]
        seed = round(1e9*runif(1, 0, 1))
        globPar[["wildN.isValueOnly"]] = FALSE
      }    
    }    
  } 
  cat(sep="", ".. created a table with experimental samples named as = ", fileAsymPath, "\n")
  # convert the file into easier-to-read html file:
  table2html(fileAsymPath) 
  # 
  # # extract statistics from the data table in fileAsymPath
  # file_asym_stats(fileAsymPath,  cntProbeLmt, xName="signifDigits", yNames=c("iterCnt", "runtime"))

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

