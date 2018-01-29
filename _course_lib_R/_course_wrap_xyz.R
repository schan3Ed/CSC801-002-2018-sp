# source("../_course_lib_R/_course_wrap_xyz.R") ; out = wrap_DEoptim(seedInit=-1) ; summary(out) ;# print(out)

# install.packages(c("DEoptim", "hash", "tools", "utils", "bazar"))

rm(list=ls())
library(hash)  ; library(tools) ; library(utils) ; library(plotrix)
source("../_course_lib_R/_course_basic.R")
source("../_course_lib_R/_course_OF_xyz.R")


# initialize globPar
globPar=hash()       ; clear(globPar)

wrap_DEoptim_plain = function( sampleSize=1, seedInit=-1, OFname="wild", nPar=1, lowerBnd=-50, upperBnd=50,
                         iterLmt=2000, popSize=128, trace=FALSE, reltol=1e-8, tolY=1e-3, asymBounds=c("") ) 
{
  # package DEoptim
  library(DEoptim) ; solver = "DEoptim"
  
  # initialize RNG with seedInit
  if (seedInit < 0) {seedInit = round(0.5 + 9999*runif(1))} ; set.seed(seedInit)  
  
  nDim = nPar
  if (OFname == "trefethen") {lowerBnd = c(-10, -10) ; upperBnd = c(10, 10) }
  if (OFname == "wild") {
    lowerBnd = -50 ; upperBnd = 50
    if (nDim  > 1)  {
      lowerBnd = rep(lowerBnd, nDim)
      upperBnd = rep(upperBnd, nDim)
    }    
  }
  times = system.time(
    ##!! MUST USE <- when embedding into system.time !!##
    out <- DEoptim(OFname, lower=lowerBnd, upper=upperBnd,
                   DEoptim.control(NP=popSize, itermax=iterLmt, trace=trace, VTR=0))
  )
  runtime   = times[["elapsed"]]
  coordBest = as.numeric(out$optim$bestmem)
  valueBest = out$optim$bestval
  iterCnt   = out$optim$iter
  
  errorY = valueBest
  if (abs(errorY) < tolY) {
    isCensored = FALSE
  } else {
    isCensored = TRUE
  }
  
  cat("\nseedInit   =", seedInit,
      "\nsolver     =", solver,
      "\nOFname     =", OFname,
      "\nnPar       =", nPar,
      "\nnDim       =", nDim,
      "\nlowerBnd   =", lowerBnd,
      "\nupperBnd   =", upperBnd,
      "\npopSize    =", popSize,
      "\niterLmt    =", iterLmt,
      "\nruntime    =", runtime,
      "\niterCnt    =", iterCnt,
      "\ncntProbe   =", NA,
      "\nisCensored =", isCensored,
      "\ntolY       =", tolY,
      "\nerrorY     =", errorY,
      "\nreltol     =", reltol,
      "\nvalueBest  =", valueBest,
      "\ncoordBest  =", signif(coordBest, 7),
      "\n")
  
  return(out)
} # wrap_DEoptim_plain

 
wrap_DEoptim = function( sampleSize=1, seedInit=-1, OFname="wildN", nPar=1, 
                        iterLmt=200, popSize=64, trace=FALSE, valueTarget=0, tolY=0.005, asymBounds=c("") ) 
{
  # package DEoptim
  library(DEoptim) ; solver = "DEoptim"
  
  # read the generic table of BKV (best-Known-Value) and related parameters
  fileBKV = "../_course_lib_R/_course_OF_xyz_BKV.txt"
  T = read.table(fileBKV, header=TRUE, comment.char="#", colClasses = "character") ;# print(T)
  idx      = which(T$OFname == OFname)    ;# print(idx)     
  BKV      = as.numeric(T[idx,2])         ;# print(BKV) 
  lowerBnd = as.numeric(T[idx,4])    ;# print(lowerBnd) 
  upperBnd = as.numeric(T[idx,5])    ;# print(upperBnd)  
  BKcoord  = as.numeric(T[idx,6])    ;# print(BKcoord) ; print(rep(BKcoord, 2)) ; return()
  
  # assign the values to global hash array for access by OF (Objective Functions)
  globPar[[paste(sep="", OFname, ".BKV")]]          = BKV
  globPar[[paste(sep="", OFname, ".tolY")]]         = tolY
  globPar[[paste(sep="", OFname, ".isValueOnly")]]  = FALSE
  print(globPar) ;# return() 
  
  # initialize RNG with seedInit
  if (seedInit < 0) {seedInit = round(0.5 + 9999*runif(1))} ; set.seed(seedInit) 

  nDim = nPar
  if (nDim  > 1)  {
    lowerBnd = rep(lowerBnd, nDim)
    upperBnd = rep(upperBnd, nDim)
  }
  BKcoord     = rep(BKcoord, nDim)
  valueTarget = 0 ;# since each OF in _course_OF_xyz.R is in the canonical form
  
  times = system.time(
    ##!! MUST USE <- when embedding into system.time !!##
    out <- DEoptim(fn=OFname, lower=lowerBnd, upper=upperBnd,
                   DEoptim.control(NP=popSize, itermax=iterLmt, trace=trace, VTR=valueTarget))
  )
  # !!VERIFY the solution!! 
  globPar[["wildN.isValueOnly"]] = TRUE
  coordBest = as.numeric(out$optim$bestmem)
  errorY    = formatC(abs(1 - wildN(coordBest)/(nDim*BKV)), format="e", digits=3)
  valueBest = out$optim$bestval
  if (valueBest == valueTarget) {
    isCensored = FALSE
  } else {
    isCensored = TRUE
  }
  iterCnt = out$optim$iter
  
  cat("\nseedInit   =", seedInit,
      "\nsolver     =", solver,
      "\nOFname     =", OFname,
      "\nnPar       =", nPar,
      "\nnDim       =", nDim,
      "\nlowerBnd   =", lowerBnd,
      "\nupperBnd   =", upperBnd,
      "\npopSize    =", popSize,
      "\niterLmt    =", iterLmt,
      "\niterCnt    =", iterCnt,
      "\ntolY       =", formatC(tolY, format="e", digits=3),
      "\nerrorY     =", errorY,
      "\nisCensored =", isCensored,
      "\nBKV        =", BKV,
      "\nBKVcomb    =", BKV*nDim,
      "\nBKcoord    =", BKcoord,
      "\nvalueBest  =", valueBest,
      "\ncoordBest  =", signif(coordBest, 7),
      "\n")
  
  return(out)
} # wrap_DEoptim

wrap_DEopt_plain = function( sampleSize=1, seedInit=-1, OFname="wild", nPar=2, lowerBnd=-50, upperBnd=50,
                               iterLmt=1000, popSize=64, trace=FALSE, reltol=1e-6, tolY=0.005, asymBounds=c("") ) 
{
  # within the package NMOF
  library(NMOF) ; solver = "DEopt"
  
  # initialize RNG with seedInit
  if (seedInit < 0) {seedInit = round(0.5 + 9999*runif(1))} ; set.seed(seedInit)  

  nDim = nPar
  if (nDim  > 1)  {
    lowerBnd = rep(lowerBnd, nDim)
    upperBnd = rep(upperBnd, nDim)
  }
  algoPars = list(min=lowerBnd, max=upperBnd)
  
  times = system.time(
    ##!! MUST USE <- when embedding into system.time !!##
    out <- DEopt(match.fun(OFname), algo=algoPars) #!!! match.fun(OFname) WORKS!!!
    ##out <- DEopt(OFname, algo=algoPars)          #!!! plain OFname DOES NOT WORK!!!
  )
  runtime   = times[["elapsed"]]
  coordBest = as.numeric(out$xbest)
  valueBest = as.numeric(out$OFvalue)
  iterCnt   = out$optim$iter
  
  if (OFname == "trefethen")     {BKV = -3.306868}
  if (OFname == "wild")          {BKV = nDim*67.46773}
  valueTarget = BKV
  errorY = 1 - abs(valueBest/valueTarget) 
  if (abs(errorY) < tolY) {
    isCensored = FALSE
  } else {
    isCensored = TRUE
  }

  cat("\nseedInit   =", seedInit,
      "\nsolver     =", solver,
      "\nOFname     =", OFname,
      "\nnPar       =", nPar,
      "\nnDim       =", nDim,
      "\nlowerBnd   =", lowerBnd,
      "\nupperBnd   =", upperBnd,
      "\npopSize    =", NA,
      "\niterLmt    =", NA,
      "\nruntime    =", runtime,
      "\niterCnt    =", NA,
      "\ncntProbe   =", NA,
      "\nisCensored =", isCensored,
      "\ntolY       =", tolY,
      "\nerrorY     =", errorY,
      "\nreltol     =", reltol,
      "\nvalueBest  =", valueBest,
      "\ncoordBest  =", signif(coordBest, 7),
      "\n")
  
  return(out)
} # wrap_DEopt_plain


wrap_DEopt = function( sampleSize=1, seedInit=-1, OFname="trefethenN", nPar=2, 
                         iterLmt=200, popSize=64, trace=FALSE, tolX=0.0001, tolY=0.005, asymBounds=c("") ) 
{ 
  # within the package NMOF
  library(NMOF) ; solver = "DEopt"
  
  # read the generic table of BKV (best-Known-Value) and related parameters
  fileBKV = "../_course_lib_R/_course_OF_xyz_BKV.txt"
  T = read.table(fileBKV, header=TRUE, comment.char="#", colClasses = "character") ;# print(T)
  idx      = which(T$OFname == OFname)    ;# print(idx)     
  BKV      = as.numeric(T[idx,2])         ;# print(BKV) 
  lowerBnd = as.numeric(T[idx,4])    ;# print(lowerBnd) 
  upperBnd = as.numeric(T[idx,5])    ;# print(upperBnd)  
  BKcoord  = as.numeric(T[idx,6])    ;# print(BKcoord) ; print(rep(BKcoord, 2)) ; return()
  
  # assign the values to global hash array for access by OF (Objective Functions)
  globPar[[paste(sep="", OFname, ".BKV")]]          = BKV
  globPar[[paste(sep="", OFname, ".tolY")]]         = tolY
  globPar[[paste(sep="", OFname, ".isValueOnly")]]  = FALSE
  print(globPar) ;# return() 
  
  # initialize RNG with seedInit
  if (seedInit < 0) {seedInit = round(0.5 + 9999*runif(1))} ; set.seed(seedInit) 
  
  nDim = nPar
  if (nDim  > 1)  {
    lowerBnd = rep(lowerBnd, nDim)
    upperBnd = rep(upperBnd, nDim)
  }
  # BKcoord     = rep(BKcoord, nDim)
  # valueTarget = 0 ;# since each OF in _course_OF_xyz.R is in the canonical form
  algoPar = list(min=lowerBnd, max=upperBnd)
  times = system.time(
    ##!! MUST USE <- when embedding into system.time !!##
    out <- DEopt(OF=trefethenN, algo=algoPar)
  )
  print(times)
  print(out)
  # # under construction 
  # > x = c(-0.0244, 0.2106)
  # > trefethenN(x)
  # [1] -3.306868
  # > globPar[["trefethenN.isValueOnly"]] = FALSE
  # > globPar
  # <hash> containing 3 key-value pair(s).
  # trefethenN.BKV : -3.306868
  # trefethenN.isValueOnly : FALSE
  # trefethenN.tolY : 0.005
  # > x
  # [1] -0.0244  0.2106
  # > trefethenN(x)
  # y valueTarget       errorY  tolY
  # [1,] -3.306868   -3.306868 4.277945e-08 0.005
  # [1] 0
  # > x = c(-0.1244, 0.2106)
  # > trefethenN(x)
  # y valueTarget    errorY  tolY
  # [1,] -2.113282   -3.306868 0.3609416 0.005
  # [1] 0.3609416
  # > lowerLmt = c(-3, -3)
  # > upperLmt = c( 3,  3) 
  # > algoPar = list(min=lowerLmt, max=upperLmt)
  # > out = DEopt(OF=trefethenN, algo=algoPar)
  # |==============================================                                                                                           |  34%            y valueTarget   errorY  tolY                                                                                                          
  # [1,] 3.239778   -3.306868 1.979712 0.005
  # y valueTarget   errorY  tolY
  
} # wrap_DEopt


OF_target_verify = function(OFname="functionNames") 
{
  # Copyright 2018, Franc Brglez
  # OF_target_verify: verify Objective Function paramerers in table _course_OF_xyz_BKV.txt
  thisFunction = "OF_target_verify"
  
  functionNames = c(
    "aParab13_2", "aParab13_4", "aParab13_8",
    "rastriginN", "schwefelN", "wildN"
  )
  
  # read the generic table of BKV (best-Known-Value) and related parameters
  fileBKV = "../_course_lib_R/_course_OF_xyz_BKV.txt"
  T = read.table(fileBKV, header=TRUE, comment.char="#", colClasses = "character") ;# print(T)
  
  if (OFname == "functionNames") {
    
    # proceed with verification of each function
    for (i in seq_len(length(functionNames))) { 

      OFname  = functionNames[i]
      idx     = which(T$OFname == OFname)   ;#  print(idx)
      # > length(which(1:3 == 2))
      # [1] 1
      # > length(which(1:3 == 4))
      # [1] 0
      # https://stackoverflow.com/questions/8343509/better-error-message-for-stopifnot
      if (length(idx) == 0) {
        stop(paste("!! the function name =", OFname, "has not been tabulated !! \n"))
      } else {
        BKV     = as.numeric(T[idx,2])       ;#  print(BKV) 
        BKCoord = as.numeric(T[idx,6])           
      }
      
      # assign the values to global hash array for access by OF (Objective Functions)
      globPar[[paste(sep="", OFname, ".BKV")]]          = BKV
      globPar[[paste(sep="", OFname, ".tolY")]]         = 0
      globPar[[paste(sep="", OFname, ".isValueOnly")]]  = TRUE ;# exit the OF-code BEFORE normalization!!
      #print(globPar) ;# return() 
      
      #### verify that BKV tabulated for this function is correct
      x = BKCoord 
      y = match.fun(OFname)(x)
      valueTarget = BKV ;# print(cbind(x, y, valueTarget))
      if (abs(valueTarget - y) > 2e-5) {
        cat(sep="", "\n** ERROR about function " , OFname, " **",
            "\n    tabulated BKV = ", BKV, "  is NOT the best-known-value",
            "\n       better BKV = ", y,   "  has been found for BKcoord = ", x,
            "\n .. please enter the new BKV  and BKCoord into the table _course_OF_xyz_BKV.txt",
            "\n")
      } else {
        cat(sep="", "\n** Successful verification about function " , OFname, " **",
            "\n    the best-known-value,      BKV     = ", BKV, 
            "\n    the best-known-coordinate, BKCoord = ", BKCoord, 
            "\n")
      }
    }
  } else {
    # argument OFname represents a name defined in the table _course_OF_xyz_BKV.txt"
    idx     = which(T$OFname == OFname)   ;#  print(idx)
    # > length(which(1:3 == 2))
    # [1] 1
    # > length(which(1:3 == 4))
    # [1] 0
    # https://stackoverflow.com/questions/8343509/better-error-message-for-stopifnot
    if (length(idx) == 0) {
      stop(paste("!! the function name =", OFname, "has not been tabulated !! \n"))
    } else {
      BKV     = as.numeric(T[idx,2])       ;#  print(BKV) 
      BKCoord = as.numeric(T[idx,6])           
    }
    
    # assign the values to global hash array for access by OF (Objective Functions)
    globPar[[paste(sep="", OFname, ".BKV")]]          = BKV
    globPar[[paste(sep="", OFname, ".tolY")]]         = 0
    globPar[[paste(sep="", OFname, ".isValueOnly")]]  = TRUE ;# exit the OF-code BEFORE normalization!!
    #print(globPar) ;# return() 
    
    #### verify that BKV tabulated for this function is correct
    x = BKCoord 
    y = match.fun(OFname)(x)
    valueTarget = BKV ;# print(cbind(x, y, valueTarget))
    if (abs(valueTarget - y) > 2e-5) {
      cat(sep="", "\n** ERROR about function " , OFname, " **",
          "\n    tabulated BKV = ", BKV, "  is NOT the best-known-value",
          "\n       better BKV = ", y,   "  has been found for BKcoord = ", x,
          "\n .. please enter the new BKV  and BKCoord into the table _course_OF_xyz_BKV.txt",
          "\n")
    } else {
      cat(sep="", "\n** Successful verification about function " , OFname, " **",
          "\n    the best-known-value,      BKV     = ", BKV, 
          "\n    the best-known-coordinate, BKCoord = ", BKCoord, 
          "\n") 
    }     
  }
  cat("\n.. values of global parameters created:\n") ; print(globPar)
      
  # VIGNETTE
  # > source("../_course_lib_R/_course_wrap_xyz.R") ; OF_target_verify(OFname="wildN")
  # <hash> containing 3 key-value pair(s).
  # wildN.BKV : 67.46773
  # wildN.isValueOnly : TRUE
  # wildN.tolY : 0
  # x        y valueTarget
  # [1,] -15.81515 67.46773    67.46773
  # 
  # ** ERROR about function wildN **
  #   tabulated BKV = 67.46773  is NOT the best-known-value
  # better BKV = 67.46773  has been found for BKcoord = -15.81515
  # .. please enter the new BKV  and BKCoord into the table _course_OF_xyz_BKV.txt
  # > wildN(-15.81515)
  # [1] 67.46773
  # > wildN(-15.81515) -  67.46773
  # [1] 4.752264e-06
  # > signif(wildN(-15.81515), 7) -  67.46773
  # [1] 0
  # > 
} # OF_target_verify 



OF_quant_exh = function(OFname="wildN", xLO=-18, xHI=18,  xBits=11)
{
  # Copyright 2018, Franc Brglez
  # OF_quant_exh: exhaustively evaluate Objective Function on a quantized grid (single dimension only)
  thisFunction = "OF_quant_exh"
  # this code extends the original scratch file fg_motiv-wild-11-bits.R
  
  # read the generic table of BKV (best-Known-Value) and related parameters
  fileBKV = "../_course_lib_R/_course_OF_xyz_BKV.txt"
  T = read.table(fileBKV, header=TRUE, comment.char="#", colClasses = "character") ;# print(T)
  idx     = which(T$OFname == OFname)    ;# print(idx)     
  BKV     = as.numeric(T[idx,2])         ;# print(BKV) 
  BKCoord = as.numeric(T[idx,6])  
  
  # assign the values to global hash array for access by OF (Objective Functions)
  globPar[[paste(sep="", OFname, ".BKV")]]          = BKV
  globPar[[paste(sep="", OFname, ".tolY")]]         = 0
  globPar[[paste(sep="", OFname, ".isValueOnly")]]  = TRUE
  globPar[[paste(sep="", OFname, ".BKCoord")]]      = BKCoord
  print(globPar) ;# return() 
  
  xRange = xHI - xLO  
  xSize = 2^xBits ; qx = 1/xSize ; dx = xRange/xSize
  xVec = NULL ; yVec = NULL ; cntMin = 0
  xMin = NULL ; idxMin = NULL
  idxVec = NULL ; coordBVec = NULL
  for (i in seq_len(xSize)) {
    coordIdx = i - 1 ; idxVec = c(idxVec, coordIdx)
    coordB = int2bin(coordIdx, xBits)
    coordIdx = bin2dec(coordB)
    coordD   = xLO + coordIdx*dx
    x      = coordD  ; xVec = c(xVec, x)
    y      = match.fun(OFname)(x) - BKV ; yVec = c(yVec, y) 
    coordBVec = c(coordBVec, coordB)
    if (y == 0) {
      cntMin = cntMin + 1
      xMin   = c(xMin, x)
      idxMin = c(idxMin, coordIdx)
    }
    #cat(coordIdx, max(xVec), length(xVec), length(yVec),  "\n")
    #cat(coordIdx, coordB, coordD, y, "\n")
  }
  #cat(coordIdx, max(xVec), length(xVec), length(yVec),  "\n")
  if (cntMin > 0) {
    cat(sep="", 
        "\n.. under xLO = ", xLO, ", xHI = ", xHI, ", and xBits = ", xBits, " (stepSize = ", dx, ")",
        "\n   least one normalized value of y = 0 has been reached\n")
    print(cbind(idxMin, xMin))
  } else {
    cat(sep="", 
        "\n.. under xLO = ", xLO, ", xHI = ", xHI, ", and xBits = ", xBits, " (stepSize = ", dx, ")",
        "\n   NO normalized value of y = 0 has been reached\n")
  }
  if (length(xVec) < 100) {print(cbind(idxVec, coordBVec, xVec, yVec))}
  
  xlimits = c(xLO-0.1, xHI+0.1)  ; ylimits = c(-0.05* max(yVec), max(yVec))  
  
  
  fileOF_quant_exh = paste(getwd(), "/fg_OF_", OFname, "_quant_exh.pdf", sep="")
  pdf(file = fileOF_quant_exh, height=7.0,  width=14.0, pagecentre=TRUE)
  
  plot(xlimits, ylimits, type="n", xaxt='n',  cex.lab=1.2,
       xlab = paste(sep="", "x:  domain of function ", OFname, ", quantized with ", xBits, " bits  on the range [",
                    xLO, ",", xHI, "] ... (BKCoord =", signif(BKCoord, 7), ")"),
       ylab = paste(sep="", OFname,"(x): normalized  with-respect-to  BKV =", BKV))
  
  # https://stackoverflow.com/questions/11775692/how-to-specify-the-actual-x-axis-values-to-plot-as-x-axis-ticks-in-r
  axis(1, at = seq(xLO, xHI, xRange/10), las=1)
  
  # add title and subtitle
  #title("Function 'wild'")  #  shifted to range [-4, 4] under 5-bit grid")
  # add lines 
  #lines(xVec, yVec, type="b", lwd=1.5, lty=3, col="blue", pch=18)
  points(xVec, yVec, col="blue", pch=18, cex=1.)
  abline(h=0.0, lty="dashed", col="red")
  abline(v=xLO, lty="dotted", lwd=1, col="gray64")
  abline(v=xHI, lty="dotted", lwd=1, col="gray64")
  
  points(BKCoord, 0,  pch=8, cex=1.4, col="brown")
  legend("topleft",  inset=0.01, ## keep a SINGLE line for each paste!!  
         paste("xLO =", xLO, " ; xHI =", xHI, " ; xRange =", xRange, " ; xBits =", xBits, 
               "\n1/2^xBits =", qx, " ; stepSize =", dx),
         cex = 1.2, 
         text.col = "black", 
         horiz = F , 
         bg = "white") 
  
  cat(sep="", "\n.. created the R-plot with file name = ", basename(fileOF_quant_exh), "\n")
  dev.off()
  
} # OF_quant_exh

OF_quant_boxes = function(OFname="aParab13_2", tolXmin=0.001, tolXmax=0.05)
{
  # Copyright 2018, Franc Brglez
  # OF_quant_boxes: evaluate Objective Function on a dynamically quantized grid, centered
  # at BVCoord. This function plots tolerance boxes of permissible errors.
  # Errors reported within these boxes are not reported as "censored".
  thisFunction = "OF_quant_boxes"
  
  # read the generic table of BKV (best-Known-Value) and related parameters
  fileBKV = "../_course_lib_R/_course_OF_xyz_BKV.txt"
  T = read.table(fileBKV, header=TRUE, comment.char="#", colClasses = "character") ;# print(T)
  
  idx      = which(T$OFname == OFname)   ;# print(idx) 
  # https://stackoverflow.com/questions/8343509/better-error-message-for-stopifnot
  if (length(idx) == 0) {
    stop(paste("!! the function name =", OFname, "has not been tabulated !! \n"))
  } else {
    lowerBnd = as.numeric(T[idx,4])        #; print(lowerBnd) 
    upperBnd = as.numeric(T[idx,5])        ;# print(upperBnd) 
    BKV      = as.numeric(T[idx,2])        ;# print(BKV) 
    BKCoord  = as.numeric(T[idx,6])
  }
  # assign the values to global hash array for access by OF (Objective Functions)
  globPar[[paste(sep="", OFname, ".BKV")]]          = BKV
  globPar[[paste(sep="", OFname, ".tolY")]]         = 0
  globPar[[paste(sep="", OFname, ".isValueOnly")]]  = FALSE
  globPar[[paste(sep="", OFname, ".BKCoord")]]      = BKCoord
  #print(globPar) ;# return() 
  
  xRange = upperBnd - lowerBnd             ;# print(xRange)
  xBitsMin = ceiling(log2(xRange/tolXmax)) ;# print(xBitsMin) 
  xBitsMax = ceiling(log2(xRange/tolXmin)) ;# print(xBitsMax) ; return()
  cat(sep=" ", "\nxBitsMin:xBitsMax = ", xBitsMin:xBitsMax, 
              " ... length = ", length(xBitsMin:xBitsMax), "\n")
  if (length(xBitsMin:xBitsMax) != 7) {
    stop(paste("!! length(xBitsMin:xBitsMax) != 7) !! 
               \n   to remove this error, adjust values of tolXmin, tolXmax ...\n"))
  } 
  
  xGrid = NULL ; yGrid = NULL ; xBitSeq = NULL
  # to keep xVec as an increasing sequence, we use one loop for xL,  one loop for xR 
  for (xBits in xBitsMin:xBitsMax) {
    tolX = xRange/2^xBits ; xBitSeq = c(xBitSeq, xBits)
    xL = BKCoord -tolX 
    yL =  match.fun(OFname)(xL) ; xGrid = c(xGrid, xL) ; yGrid = c(yGrid, yL)
  }
  for (xBits in xBitsMax:xBitsMin) {
    tolX = xRange/2^xBits; xBitSeq = c(xBitSeq, xBits)
    xR = BKCoord +tolX ; yR =  match.fun(OFname)(xR) ; xGrid = c(xGrid, xR) ; yGrid = c(yGrid, yR)  
  }
  print(cbind(xGrid, yGrid, xBitSeq)) ;# return()

  fileOF_quant_box = paste(getwd(), "/fg_OF_", OFname, "_quant_box.pdf", sep="")
  pdf(file = fileOF_quant_box, height=7.0,  width=7.0, pagecentre=TRUE)
  
  xlimits = c(min(xGrid), max(xGrid))  ; ylimits = c(min(yGrid), max(yGrid))  
  plot(xlimits, ylimits, log="y", type="n", axes=F,
       xlab = paste(sep="", "BKCoord +- xRange/2^xBits (BKCoord =", signif(BKCoord, 7), ")",
                    "\ndomain of function ", OFname, "(x): xRange = [", lowerBnd, ", ", upperBnd, "]"),
       ylab = paste(sep="", "normalized errorY, with-respect-to  BKV = ", BKV, ", from  " , OFname, "(x)"))
  
 axis(side = 1, at=formatC(xGrid, 5))
 axis(side = 2, at=formatC(yGrid , format="e", digits=1))   
  
  points(xGrid, yGrid)
  abline(h=yGrid[1:7], lty="dotted", col="blue")
  
  # draw horizontal lines that define the target box
  lines(c(xGrid[1], xGrid[14]) , c( yGrid[1], yGrid[1]), col="red", lty="solid", lwd=2 )
  lines(c(xGrid[2], xGrid[13]) , c( yGrid[2], yGrid[2]), col="red", lty="solid", lwd=2 )
  lines(c(xGrid[3], xGrid[12]) , c( yGrid[3], yGrid[3]), col="red", lty="solid", lwd=2 )
  lines(c(xGrid[4], xGrid[11]) , c( yGrid[4], yGrid[4]), col="red", lty="solid", lwd=2 )
  lines(c(xGrid[5], xGrid[10]) , c( yGrid[5], yGrid[5]), col="red", lty="solid", lwd=2 )
  lines(c(xGrid[6], xGrid[9])  , c( yGrid[6], yGrid[6]), col="red", lty="solid", lwd=2 )
  lines(c(xGrid[7], xGrid[8])  , c( yGrid[7], yGrid[7]), col="red", lty="solid", lwd=2 )
  # draw vertical lines that define the target box
  lines(c(xGrid[1],  xGrid[1]),  c(yGrid[1],  0.9*min(yGrid)), col="red", lty="solid", lwd=2)
  lines(c(xGrid[14], xGrid[14]), c(yGrid[14], 0.9*min(yGrid)), col="red", lty="solid", lwd=2)
  lines(c(xGrid[2],  xGrid[2]),  c(yGrid[2],  0.9*min(yGrid)), col="red", lty="solid", lwd=2)
  lines(c(xGrid[13], xGrid[13]), c(yGrid[13], 0.9*min(yGrid)), col="red", lty="solid", lwd=2)
  lines(c(xGrid[3],  xGrid[3]),   c(yGrid[3], 0.9*min(yGrid)), col="red", lty="solid", lwd=2)
  lines(c(xGrid[12], xGrid[12]), c(yGrid[12], 0.9*min(yGrid)), col="red", lty="solid", lwd=2)
  lines(c(xGrid[4],  xGrid[4]),  c(yGrid[4],  0.9*min(yGrid)), col="red", lty="solid", lwd=2)
  lines(c(xGrid[11], xGrid[11]), c(yGrid[11], 0.9*min(yGrid)), col="red", lty="solid", lwd=2)
  lines(c(xGrid[5],  xGrid[5]),  c(yGrid[5],  0.9*min(yGrid)), col="red", lty="solid", lwd=2)
  lines(c(xGrid[10], xGrid[10]), c(yGrid[10], 0.9*min(yGrid)), col="red", lty="solid", lwd=2)
  lines(c(xGrid[6],  xGrid[6]),  c(yGrid[6],  0.9*min(yGrid)), col="red", lty="solid", lwd=2)
  lines(c(xGrid[9],  xGrid[9]),  c(yGrid[9],  0.9*min(yGrid)), col="red", lty="solid", lwd=2)
  lines(c(xGrid[7],  xGrid[7]),  c(yGrid[7],  0.9*min(yGrid)), col="red", lty="solid", lwd=2)
  lines(c(xGrid[8],  xGrid[8 ]), c(yGrid[8],  0.9*min(yGrid)), col="red", lty="solid", lwd=2)
  
  # template for this table: fg_asym_pi_digits.pdf
  # addtable2plot(7.9, 0.12, myTable, xjust=0,yjust=0, bty="n", hlines=FALSE, vlines=FALSE, 
  #               title=NULL, cex=0.7 ,xpad=0.15, ypad=0.7)
  # yLine=0.050 ; lines(c(7.8, 13.8), c(yLine,yLine), col="black", lwd=0.5)
  errorY = formatC(yGrid[1:7] , format="e", digits=2)
  xBits = xBitSeq[1:7]
  myTable = cbind(errorY, xBits) ; print(myTable)

  addtable2plot(xGrid[1]*(1+0.000), yGrid[1]*(1-0.000), myTable, xjust=0, yjust=0, bty="n", hlines=FALSE, vlines=FALSE, 
                title=NULL, cex=0.9 , bg="white")
  #yLine=0.050 ; lines(c(xGrid[1], yGrid[1]), c(yLine,yLine), col="black", lwd=0.5)
  
  
  cat(sep="", "\n.. created the R-plot with file name = ", basename(fileOF_quant_box), "\n")
  dev.off()
  
} # OF_quant_boxes

OF_target_asym_notYet = function(OFname="wildN", lowerBnd=-18, upperBnd=-14, xBits=11, tolX=5e-2, tolY=0.05, OFpar=1:1) 
{
  # Copyright 2018, Franc Brglez
  # OF_target_verify: verify the Objective Function represented in a normalized form
  thisFunction = "OF_target_asym"
  
  # read the generic table of BKV (best-Known-Value) and related parameters
  fileBKV = "../_course_lib_R/_course_OF_xyz_BKV.txt"
  T = read.table(fileBKV, header=TRUE, comment.char="#", colClasses = "character") ;# print(T)
  idx      = which(T$OFname == OFname)    ;# print(idx)     
  BKV0     = as.numeric(T[idx,2])         ;# print(BKV0) 
  if (length(OFpar) == 0) {
    OFpar = as.numeric(T[idx,3])
  }
  if (length(lowerBnd) == 0) {lowerBnd = as.numeric(T[idx,4]) }
  if (length(upperBnd) == 0) {upperBnd = as.numeric(T[idx,5]) }
  BKCoord0 = as.numeric(T[idx,6])    
  xRange = upperBnd - lowerBnd
  print(OFpar)
  if (OFname == "clusterLJ") {
    # agrep for approximate matching ...
    nPar   = length(OFpar)
    nDim   = 3 
    nAtoms = as.integer(nPar/nDim)
    cat("OFname =", OFname, ", nPar =", nPar, ", nDim =", nDim, ", nAtoms =", nAtoms, 
        "\n.. coding for this OFname is in progress ...\n")
    return()
  } else {
    # proceed with remainder of OFname in normal form ....
    print(OFpar)
    if (length(OFpar) == 0) {
      cat("ERROR .. an integer vector expected for OFname =", OFname, "\n") ; return()
    }
    nParMax = max(OFpar)
    if (nParMax == 1) {
      valueTarget = BKV0
      BKCoord     = rep(BKCoord0, 1)
      #print(cbind(valueTarget, BKCoord)) 
    }
    if (nParMax == 2) {
      valueTarget = 2*BKV0
      BKCoord     = rep(BKCoord0, 2)
      #print(cbind(valueTarget, BKCoord))    
    }
    if (nParMax == 3) {
      valueTarget = 3*BKV0
      BKCoord    = rep(BKCoord0, 3)
      #print(cbind(valueTarget, BKCoord)) 
    }
  }
  print(cbind(valueTarget, BKCoord)) ;# return()
  
  # assign the values to global hash array for access by OF (Objective Functions)
  globPar[[paste(sep="", OFname, ".BKV")]]          = BKV0
  globPar[[paste(sep="", OFname, ".tolY")]]         = tolY
  globPar[[paste(sep="", OFname, ".isValueOnly")]]  = TRUE ;# exit the OF-code BEFORE normalization!!
  print(globPar) ;# return() 
  
  if (nParMax == 1) { 
    valueTarget = 1*BKV0
    BKCoord     = rep(BKCoord0, 1)
    notCensoredCnt = 0
    # given the value of xBit, create a bounded uniform grid for x-parameters x
    xSize = 2^xBits ; qx = 1/xSize ;  dx = xRange/xSize
    xGrid = seq(lowerBnd, upperBnd, dx) ;# print(xGrid) ;# return()
    for (x in xGrid) {
      y      = match.fun(OFname)(x) 
      # since BKV >= 1, we evaluate relative error 
      errorY = abs( 1  - y/valueTarget )  
      errorX = errorX = distance_xyz(x ,BKCoord)
      if (errorX < tolX  && errorY < tolY) { 
        isCensored = FALSE ; notCensoredCnt = notCensoredCnt + 1
        cat("\nx =", x, "y =", y, 
            "\nerrorX =", formatC(errorX, format="e", digits=3), 
            "\nerrorY =", formatC(errorY, format="e", digits=3), "isCensored =", FALSE, "\n")
      }
    }
    cat("notCensoredCnt =", notCensoredCnt, "\n")
    return()    
  }
  
  if (nParMax == 2) { 
    valueTarget = 2*BKV0
    BKCoord     = rep(BKCoord0, 2)
    notCensoredCnt = 0
    xSize = 2^xBits ; qx = 1/xSize ;  dx = xRange/xSize
    xGrid = seq(lowerBnd, upperBnd, dx)  
    for (x1 in xGrid) {
      for (x2 in xGrid) {
        x = c(x1, x2)
        y = match.fun(OFname)(x) 
        # since BKV >= 1, we evaluate relative error 
        errorY = abs( 1  - y/valueTarget )  
        errorX = distance_xyz(x ,BKCoord)
        if (errorX < tolX  && errorY < tolY) { 
          isCensored = FALSE ; notCensoredCnt = notCensoredCnt + 1
          cat("\nx =", x, "y =", y, 
              "\nerrorX =", formatC(errorX, format="e", digits=3), 
              "\nerrorY =", formatC(errorY, format="e", digits=3), "isCensored =", FALSE, "\n")
        }
      }
    }
    cat("notCensoredCnt =", notCensoredCnt, "\n")
    return() 
  }
  
  if (nParMax == 3) {
    valueTarget = 3*BKV0
    BKCoord     = rep(BKCoord0, 3)
    notCensoredCnt = 0
    xSize = 2^xBits ; qx = 1/xSize ;  dx = xRange/xSize
    xGrid = seq(lowerBnd, upperBnd, dx)  
    for (x1 in xGrid) {
      for (x2 in xGrid) {
        for (x3 in xGrid) {
          x = c(x1, x2, x3)
          y = match.fun(OFname)(x) 
          # since BKV >= 1, we evaluate relative error 
          errorY = abs( 1  - y/valueTarget )  
          errorX = distance_xyz(x ,BKCoord)
          if (errorX < tolX  && errorY < tolY) { 
            isCensored = FALSE ; notCensoredCnt = notCensoredCnt + 1
            cat("\nx =", x, "y =", y, 
                "\nerrorX =", formatC(errorX, format="e", digits=3), 
                "\nerrorY =", formatC(errorY, format="e", digits=3), "isCensored =", FALSE, "\n")
          }           
        }
      }
    }
    cat("notCensoredCnt =", notCensoredCnt, "\n")
    return()   
  }
  return()
  
} # OF_target_asym_notYet

OF_target_box_retired = function(OFname="wildN", lowerBnd=-18, upperBnd=-14, xBits = 9:12) {
  # Copyright 2018, Franc Brglez
  # OF_target_box_retired: plot rectangular bounds for tolX*tolY region
  thisFunction = " OF_target_box_retired"
  
  # https://www.statmethods.net/advgraphs/axes.html
  # Math Annotations
  # You can add mathematically formulas to a graph using TEX-like rules. 
  # Text Annotations
  # Text can be added to graphs using the text( ) and mtext( ) functions. 
  # text( ) places text within the graph while 
  # mtext( ) places text in one of the four margins.
  # 
  # text(location, "text to place", pos, ...)
  # mtext("text to place", side, line=n, ...)
  # 
  # option	   description
  # location	location can be an x,y coordinate. 
  #           Alternatively, the text can be placed interactively 
  #           via mouse by specifying location as locator(1).
  # pos	      position relative to location. 1=below, 2=left, 3=above, 4=right. 
  #           If you specify pos, you can specify offset= in percent of character width.
  # side	    which margin to place text. 1=bottom, 2=left, 3=top, 4=right. 
  #           you can specify line= to indicate the line in the margin starting 
  #           with 0 and moving out. you can also specify adj=0 for left/bottom 
  #           alignment or adj=1 for top/right alignment.
  #           
  #           Other common options are cex, col, and font 
  #           (for size, color, and font style respectively).
  # 
  # See help(plotmath) for details and examples.
  # https://stackoverflow.com/questions/1154242/remove-plot-axis-values
  
  
  # read the generic table of BKV (best-Known-Value) and related parameters
  fileBKV = "../_course_lib_R/_course_OF_xyz_BKV.txt"
  T = read.table(fileBKV, header=TRUE, comment.char="#", colClasses = "character") ;# print(T)
  idx      = which(T$OFname == OFname)    ;# print(idx)     
  BKV      = as.numeric(T[idx,2])         ;# print(BKV) 
  BKCoord  = as.numeric(T[idx,6])    
  xRange = upperBnd - lowerBnd 
  cat("lowerBnd =", lowerBnd , "  upperBnd =", upperBnd, "  xRange =", xRange, "\n")
  BKC   = BKCoord
  
  # assign the values to global hash array for access by OF (Objective Functions)
  globPar[[paste(sep="", OFname, ".BKV")]]          = BKV
  globPar[[paste(sep="", OFname, ".tolY")]]         = 0
  globPar[[paste(sep="", OFname, ".isValueOnly")]]  = TRUE
  globPar[[paste(sep="", OFname, ".BKCoord")]]      = BKCoord
  print(globPar) ;# return() 
  
  dx    = xRange/2^xBits ; cat("  dx =", dx , "\n2*dx =", 2*dx, "\n")
  xGrid = sort(c(BKC - dx, BKC + dx)) ; xGrid 
  yGrid = NULL
  for (x in xGrid) { yGrid = c(yGrid, match.fun(OFname)(x)) }
  errorY  = formatC( abs(BKV - yGrid) ,  format="e", digits=3)
  errorYR = formatC( abs(1 - yGrid/BKV), format="e", digits=3)
  
  fileOF_box = paste(getwd(), "/fg_OF_box.pdf", sep="")
  pdf(file = fileOF_box, width=7.0, height=7.0, pagecentre=TRUE)
  
  print(cbind(xGrid, yGrid, errorY, errorYR)) ;# return()
  
  plot(xGrid, errorYR, log="y", axes=F, ylim=c(1e-4, 1e-2),
       xlab = "x = BKCoord +- rangeX/2^xBits",
       ylab = "error = |1 - wildN(x)/BKV|" 
  )
  axis(side = 1, formatC(xGrid, 6))
  axis(side = 2, at=c(1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2))
  abline(h=errorYR, col="blue", lty="dotted", lwd=1)
  lines(c(xGrid[1], xGrid[8]) , c( errorYR[1], errorYR[1]), col="red", lty="solid", lwd=2 )
  lines(c(xGrid[2], xGrid[7]) , c( errorYR[2], errorYR[2]), col="red", lty="solid", lwd=2  )
  lines(c(xGrid[3], xGrid[6]) , c( errorYR[3], errorYR[3]), col="red", lty="solid", lwd=2  )
  lines(c(xGrid[4], xGrid[5]) , c( errorYR[4], errorYR[4]), col="red", lty="solid", lwd=2  )
  lines(c(xGrid[1], xGrid[1]), c(errorYR[1], 9e-5), col="red", lty="solid", lwd=2)
  lines(c(xGrid[8], xGrid[8]), c(errorYR[1], 9e-5), col="red", lty="solid", lwd=2)
  lines(c(xGrid[2], xGrid[2]), c(errorYR[2], 9e-5), col="red", lty="solid", lwd=2)
  lines(c(xGrid[7], xGrid[7]), c(errorYR[2], 9e-5), col="red", lty="solid", lwd=2)
  lines(c(xGrid[3], xGrid[3]), c(errorYR[3], 9e-5), col="red", lty="solid", lwd=2)
  lines(c(xGrid[6], xGrid[6]), c(errorYR[3], 9e-5), col="red", lty="solid", lwd=2)
  lines(c(xGrid[4], xGrid[4]), c(errorYR[4], 9e-5), col="red", lty="solid", lwd=2)
  lines(c(xGrid[5], xGrid[5]), c(errorYR[4], 9e-5), col="red", lty="solid", lwd=2)
  
  text(-15.815, 4e-3, "2*tolX = rangeX/2^xBits\nrangeX = -18 - (-16) ; xBits=c(9:12)")
  text(-15.815, 8.8e-3, "2*tolX=0.015625,  xBits=9")
  text(-15.815, 2.2e-3, "2*tolX=0.0078125,  xBits=10")
  text(-15.815, 5.6e-4, "2*tolX=0.00390625,  xBits=11")
  text(-15.815, 1.75e-4, "2*tolX=\n0.001953125,\n xBits=12")
  
  cat(sep="", "\n.. created the R-plot with file name = ", basename(fileOF_box), "\n")
  dev.off()
  
} # OF_target_box_retired