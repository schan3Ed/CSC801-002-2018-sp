# source("_course_xPer_pi.r") ; xPer_pi_BPP()
# source("_course_xPer_pi.r") ; xPer_pi_darts()
# require(stats)
rm(list=ls())
library(plotrix)
source("../_course_lib_r/_course_basic.r")
source("../_course_lib_r/_course_OF_xyz.r")
# (https://www.csc.ncsu.edu/research/colloquia/seminar-post.php?id=779)

xPer_pi_BPP = function(signifDigits=10)
{
  # Copyright 2018, Franc Brglez
  # experiment with function pi_BPP, asympototic variable is 'signifDigits'
  
  fileAsym = "fg-asym-pi_BPP.txt"
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
    res = pi_BPP(sigD) # saving results in 'res'
    # return(list(numTerms=numTerms, piHat=piHat, tolR=tolR,  
    #             error=error, isCensored=isCensored))
    row =c(sigD, res$numTerms, res$piHat, res$tolR, signif(res$error,4), res$isCensored)
    cat(row, "\n")
    write(row, file=fileAsymPath, ncolumns=6, append=TRUE, sep="\t")
  }
  
  # read data from file and use it for plotting
  df = read.table(fileAsymPath, header=TRUE, colClasses="character")
  #print(df) ; print(df$piHat) ; return()
  
  # plot signifDigits vs numTerms
  fileAsymPdf_numTerms = paste(getwd(), "/fg-asym-pi_BPP_numTerms.pdf", sep="")
  pdf(file = fileAsymPdf_numTerms, width=7.0, height=7.0, pagecentre=TRUE)
  x = signifDigits ; y = df$numTerms
  print(cbind(x, y))
  plot(x, y)
  
  cat(sep="", ".. created the R-plot with file name = ", fileAsymPdf_numTerms, "\n")
  dev.off()  
  
  # plot signifDigits vs tolR (toleranceRadius)
  fileAsymPdf_tolR = paste(getwd(), "/fg-asym-pi_BPP_tolR.pdf", sep="")
  pdf(file = fileAsymPdf_tolR, width=7.0, height=7.0, pagecentre=TRUE) 
  x = signifDigits ; y = df$tolR
  print(cbind(x, y))
  plot(x, y, log="y")
  
  cat(sep="", ".. created the R-plot with file name = ", fileAsymPdf_tolR, "\n")
  dev.off()  
} # xPer_pi_BPP

xPer_pi_needles = function(sampleSize=5, seedInit=1215, 
                          cntProbeLmt=500, needleRatio=1.0, signifDigits=2:4)
{
  # Copyright 2018, Franc Brglez
  # experiment with function pi_needles, asympototic variable is 'signifDigits'
  
  if (seedInit < 0) {seedInit = round(9999*runif(1, 0, 1))} ; seed = seedInit 
  
  fileAsym = "fg-asym-pi_needles.txt"
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
               "tolRadius", "error", "isCensored", "cntProbe", "runtime")
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
      # return(list(piMC=piMC, tolR=tolR,  error=error, 
      #             isCensored=isCensored, cntProbe=cntProbe))
      row = c(sampleId, seed, sigD, needleRatio, res$piMC,
          res$tolR, res$error, res[["isCensored"]],
          res$cntProbe, times[["elapsed"]]) ;# cat(row, "\n")
      write(row, file=fileAsymPath, ncolumns=length(colNames), append=TRUE, sep="\t")
      
      # introduce a new seed for the next sample
      seed = round(1e9*runif(1, 0, 1))
    }    
  }
  cat(sep="", ".. created a table with experimental samples named as = ", fileAsymPath, "\n")
  
  # extract statistics from yhe data table in fileAsymPath
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

xPer_pi_darts = function(sampleSize=5, seedInit=1215, cntProbeLmt=1000, signifDigits=3:5)
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
  write(colNames, file=fileAsymPath, ncolumns=9, append=TRUE, sep="\t") ;# return()
  
  # start the asymptotic experiments: x = specified significant digits
  sampleId = 0
  for (sigD in signifDigits) {
    
    for (sample in seq_len(sampleSize)) {
      
      sampleId = sampleId + 1
      times = system.time(
        ##!! MUST USE <- when embedding into system.time !!##
        res <- pi_darts(cntProbeLmt, sigD)
      )
      # return(list( piMC=piMC, tolR=tolR,  error=error, 
      #             isCensored=isCensored, cntProbe=cntProbe))
      row = c(sampleId,  seed, sigD, signif(res$piMC, sigD),
              res$tolR, res$error, res[["isCensored"]],
              res$cntProbe, times[["elapsed"]]) ;# cat(row, "\n")
      write(row, file=fileAsymPath, ncolumns=9, append=TRUE, sep="\t")
      
      # introduce a new seed for the next sample
      seed = round(1e9*runif(1, 0, 1))
    }    
  }
  cat(sep="", ".. created a table with experimental samples named as = ", fileAsymPath, "\n")
  
  # extract statistics from yhe data table in fileAsymPath
  file_asym_stats(fileAsymPath,  cntProbeLmt, xName="signifDigits", yNames=c("cntProbe", "runtime"))
  
} # xPer_pi_darts



