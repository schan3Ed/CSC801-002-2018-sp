# source("_course_xPer_pi.R") ; xPer_pi_BPP()
# source("_course_xPer_pi.R") ; xPer_pi_darts()
# require(stats)
rm(list=ls())
library(plotrix)
source("../_course_lib_R/_course_basic.R")
source("../_course_lib_R/_course_OF_xyz.R")

xPer_pi_BPP = function(signifDigits=10)
{
  # Copyright 2018, Franc Brglez
  # experiment with function pi_BPP, asympototic variable is 'signifDigits'
  
  fileAsym = "fg_asym_pi_BPP.txt"
  header   = file_header(fileAsym)$fileHeader
  #cat(header,"\n") ;# return()
  header = c(
    header,
    paste("# command   = xPer_pi_BPP(...) \n#")
  )
  cat(header,"\n")
  fileAsymPath = paste(getwd(), "/", fileAsym, sep="")
  write(header, file=fileAsym, ncolumns=1, append=FALSE, sep="\t")
  # write column names for all data rows in fileAsym 
  colNames = c("signifDigits", "numTerms", "piHat", "tolRad",
               "error", "isCensored")
  write(colNames, file=fileAsymPath, ncolumns=6, append=TRUE, sep="\t") ;# return()
  
  for (sigD in signifDigits) {
    res = pi_BPP(sigD)              
    # pi_BPP : return(list(numTerms=numTerms, piHat=piHat, tolR=tolR,  
    #                 error=error, isCensored=isCensored))
    row =c(sigD, res$numTerms, res$piHat, res$tolR, signif(res$error,4), res$isCensored)
    cat(row, "\n")
    write(row, file=fileAsymPath, ncolumns=length(colNames), append=TRUE, sep="\t")
  }
  # convert the file into easier-to-read html file:
  table2html(fileAsymPath)
  
  # read data from file and use it for plotting
  df = read.table(fileAsymPath, header=TRUE, colClasses="character")
  #print(df) ; print(df$piHat) ; return()
  
  # plot signifDigits vs numTerms
  fileAsymPdf_numTerms = paste(getwd(), "/fg_asym_pi_BPP_numTerms.pdf", sep="")
  pdf(file = fileAsymPdf_numTerms, width=7.0, height=7.0, pagecentre=TRUE)
  xAsym = signifDigits ; yAsym = as.numeric(df$numTerms)
  print(cbind(xAsym, yAsym))
  fit    = lm(yAsym ~ xAsym) ;# print(fit)
  co     = coef(fit)         ; print(co)
  a      = signif(co[1], 4)  ; b = signif(co[2], 4)
  yModel = a + b*xAsym  
  dif = yAsym - yModel
  print(cbind(xAsym, yAsym, yModel, dif))
  
  plot(3, 2, type="n", xlim=c(3,15), ylim=c(2,12), xlab="significant digits", 
       ylab="number of required terms in pi_BPP formula" )
  points(xAsym, yAsym)
  abline(a=a, b=b, col="red")
  text(8.5,11.50, "pi_BPP formula, invented in 1997, converges rapidly to value of pi:")
  text(8.5,10.75, "First 11 terms return pi = 3.14159265358979 (15 significant digits).")
  text(7.5,10.00, "Asymptotic model to predict the number of terms:.")
  text(7.5, 9.25, "    -0.01648 + 0.7198*significantDigits")
  
  cat(sep="", ".. created the R-plot with file name = ", fileAsymPdf_numTerms, "\n")
  dev.off()   ;# return()
  
  # plot signifDigits vs tolR (toleranceRadius)
  fileAsymPdf_tolR = paste(getwd(), "/fg_asym_pi_BPP_tolR.pdf", sep="")
  pdf(file = fileAsymPdf_tolR, width=7.0, height=7.0, pagecentre=TRUE) 
  xAsym = signifDigits ; yAsym = as.numeric(df$tolR)
  print(cbind(xAsym, yAsym))
  fit    = lm(log10(yAsym) ~ xAsym) ;# print(fit)
  co     = coef(fit)         ; print(co)
  a      = 10^(co[1])    ; b = co[2]/log10(10) ;# c=log10(2)/co[2]
  a      = signif(a, 4)  ; b = signif(b, 4)
  base   = 10^b          ; print(cbind(a,b,base))
  yModel = a*base^xAsym 
  dif = yAsym - yModel
  print(cbind(xAsym, yAsym, yModel, dif)) 
  #plot(xAsym, yAsym, log="y")
  plot(3, 1e-15, type="n",xlim=c(3,15), ylim=c( 1e-15, 1e-3), log="y",
       xlab="significant digits", ylab="maximum error associated with pi_BPP formula" )
  points(xAsym, yAsym)
  abline(fit,col="red")
  text(7.5,1e-13, "Asymptotic model to predict the maximum error:")
  text(7.5,1e-14, "    (0.5)*(0.1)^(significantDigits -1)")
  text(8.5,1e-15, "pi_BPP formula, invented in 1997, converges rapidly to value of pi")
  
  cat(sep="", ".. created the R-plot with file name = ", fileAsymPdf_tolR, "\n")
  dev.off()  
} # xPer_pi_BPP

xPer_pi_needles = function(sampleSize=10, seedInit=1215, 
                          cntProbeLmt=1000, needleRatio=1.0, signifDigits=2:4)
{
  # Copyright 2018, Franc Brglez
  # experiment with function pi_needles, asympototic variable is 'signifDigits'
  
  if (seedInit < 0) {seedInit = round(9999*runif(1, 0, 1))} ; seed = seedInit 
  
  fileAsym = "fg_asym_pi_needles.txt"
  header   = file_header(fileAsym)$fileHeader
  fileAsym = file_header(fileAsym)$fileStamped
  fileStamp= file_header(fileAsym)$stamp 
  header = c(
    header,
    paste("# command   = xPer_pi_needles(...) \n#")
  )
  fileAsymPath = paste(getwd(), "/", fileAsym, sep="")
  write(header, file=fileAsymPath, ncolumns=1, append=FALSE, sep="\t")
  # write column names for all data rows in fileAsym  
  colNames = c("sampleId", "seedInit", "signifDigits", "needleRatio", "piMC",
               "tolRadius", "error", "isCensored", "numerator", "denominator", "cntProbe", "runtime")
  write(colNames, file=fileAsymPath, ncolumns=length(colNames), append=TRUE, sep="\t")
  
  # start the asymptotic experiments: x = specified significant digits
  sampleId = 0
  for (sigD in signifDigits) {
    
    for (sample in seq_len(sampleSize)) {
      
      sampleId = sampleId + 1
      times = system.time(
        ##!! MUST USE <- when embedding into system.time !!##
        res <- pi_needles(cntProbeLmt, sigD, needleRatio) 
      )
      # pi_needles: return(list(piMC=piMC, tolR=tolR,  error=error, isCensored=isCensored, 
      #                     numerator=numerator, denominator=denominator, cntProbe=cntProbe))
      row = c(sampleId, seed, sigD,  needleRatio, signif(res$piMC, sigD),
          res$tolR, res$error, res[["isCensored"]],
          res$numerator, res$denominator, res$cntProbe, times[["elapsed"]]) ;# cat(row, "\n")
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
  # convert the file into easier-to-read html file:
  table2html(fileAsymPath)
  
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



