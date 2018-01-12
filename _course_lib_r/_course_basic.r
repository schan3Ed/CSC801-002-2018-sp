# source("./_course_basic.r")

signifDigits = function(digit) 
{
  # Copyright 2018, Franc Brglez
  # AND contributors to
  # https://stackoverflow.com/questions/27767841/how-do-i-determine-the-number-of-significant-figures-in-data-in-r
  #
  # > xChar <- c("4", "4.0", "4.00", "28.382", "120",
  #              +         "82.3", "100", "100.0", "30.0003")       ; xChar
  # [1] "4"       "4.0"     "4.00"    "28.382"  "120"     "82.3"    "100"     "100.0"   "30.0003"
  # > xNum <- xChar
  # > decimal <- grep(".", xChar, fixed=TRUE)          ; decimal
  # [1] 2 3 4 6 8 9
  # > xNum[decimal] <- gsub(".$", "1", xChar[decimal]) ; xNum
  # [1] "4"       "4.1"     "4.01"    "28.381"  "120"     "82.1"    "100"     "100.1"   "30.0001"
  # > xNum <- as.numeric(xNum)                         ; xNum
  # [1]   4.00000000000000   4.10000000000000   4.01000000000000  28.38100000000000 120.00000000000000
  # [6]  82.09999999999999 100.00000000000000 100.09999999999999  30.00010000000000
  # > n =  30  
  # > # true counts
  #   > xSigDig <- vector(length=length(xNum))
  #   > for (i in n:1) xSigDig[xNum == signif(xNum, i)] <- i
  #   > xSigDig
  #   [1] 1 2 3 5 2 3 1 4 6
  #   >   
  
  xChar         = as.character(digit)
  xNum          = xChar
  decimal       = grep(".", xChar, fixed=TRUE) 
  xNum[decimal] = gsub(".$", "1", xChar[decimal])
  xNum          = as.numeric(xNum)
  n       = 30  # maxiumum number of significant digits
  xSigDig = vector(length=length(xNum))
  for (i in n:1) {xSigDig[xNum == signif(xNum, i)] = i}
  
  return(xSigDig)
  # VIGNETTES
  # > signifDigits(4)
  # [1] 1
  # > signifDigits(4.1)
  # [1] 2
  # > signifDigits(4.01)
  # [1] 3
  # > signifDigits(28.382)
  # [1] 5
  # > signifDigits(120)
  # [1] 2
  # > signifDigits(82.3)
  # [1] 3
  # > signifDigits(100)
  # [1] 1
  # > signifDigits(100.1)
  # [1] 4
  # > signifDigits(30.0001)
  # [1] 6
  
} # signifDigits

file_ext = function(x="fileTest.txt") 
{
  # Copyright 2018   
  # Contributors to
  # https://stackoverflow.com/questions/29627833
  # /extract-extension-of-file-from-working-directory-and-check-condition
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
  # for explanation about '1L' and '1' see 
  # https://stackoverflow.com/questions/7014387/whats-the-difference-between-1l-and-1
} # file_ext

file_root = function(path) 
{
  # Copyright 2018   
  # A code snippet from Brian Ripley
  # http://r.789695.n4.nabble.com
  # /Extracting-File-Basename-without-Extension-td878817.html  
  # NOTE: the function file_path_sans_ext(path), from  fileutils, has issues with DEoptim 
  fileBase = basename(path)
  fileRoot = sub("\\.[^.]*$", "", fileBase)
  return(fileRoot)
} # file_root

file_header = function(fileName="fileTest.txt") 
{
  # Copyright 2018, Franc Brglez
  userId      = Sys.info()[["user"]] 
  nodeName    = Sys.info()[["nodename"]] 
  sysName     = Sys.info()[["sysname"]] 
  date        = format(Sys.time(), "%Y%m%d")
  timeStamp   = format(Sys.time(), "%d%H%M%S")
  fileExtension = file_ext(fileName)
  fileRoot      = file_root(fileName)
  fileStamped   = paste(sep="", fileRoot, "_", timeStamp, ".", fileExtension)
  
  fileHeader  = paste( "# fileName  =", fileName,
                     "\n# userId    =", userId,
                     "\n# nodeName  =", nodeName,
                     "\n# sysName   =", sysName,
                     "\n# date      =", date,
                     "\n# timeStamp =", timeStamp,
                     "\n#")
  
  return(list(fileHeader=fileHeader, fileStamped=fileStamped, timeStamp=timeStamp))
} # file_header


file_asym_stats = function(fileName="fg-asym-pi_darts_test.txt", cntProbeLmt,
                           xName="signifDigits", yNames=c("cntProbe", "runtime")) 
{
  # Copyright 2018, Franc Brglez
  # read data from file and compute statistics
  df = read.table(fileName, header=TRUE, colClasses="character")
  # x=df$signifDigits ; print(x) ; x=df[[xName]] ; print(x) ;# is OK but x=df$xName is NOT OK
  # print(df) ; x=df[[xName]] ; y1=df[[yNames[1]]] ; y2=df[[yNames[2]]] ; print(cbind(x,y1,y2)) ;  return() 
  
  # get name of fileStats and prepare the header
  fileExtension = file_ext(fileName)
  fileRoot      = file_root(fileName)
  fileStats     = paste(fileRoot, "_stats.", fileExtension, sep="")
  header = c(
    file_header(fileStats)$fileHeader,
    paste("# command   = file_asym_stats(...) \n#")
  )
  cat(header,"\n")
  fileStatsPath = paste(getwd(), "/", fileStats, sep="")
  write(header, file=fileStatsPath, ncolumns=1, append=FALSE, sep="\t")
  # write column names for all data rows in fileStatsPath 
  colNames = c(xName, "sampleSize", "cntHitsPerc", "cntProbeLmt")  
  for (yIdx in seq_len(length(yNames))) {
    colNames = c(colNames, 
                 paste(yNames[yIdx],"_med", sep=""),
                 paste(yNames[yIdx],"_mean", sep=""),
                 paste(yNames[yIdx],"_sd", sep=""),
                 paste(yNames[yIdx],"_min", sep=""),
                 paste(yNames[yIdx],"_max", sep="")
    )
  }
  write(colNames, file=fileStatsPath, ncolumns=length(colNames), append=TRUE, sep="\t") ;#  return()
  
  # subset data from table by 'xName' and 'isCensored'  
  xVals = sort(unique(df[[xName]]))
  cat(sep="", "values associated with the asymptotic variable(", xName, ") = ",xVals, "\n")
  for (xIdx in seq_len(length(xVals))) {
    row = xVals[xIdx]
    dfs = df[df[[xName]] == xVals[xIdx] ,]     ;# print(dfs)
    sampleSizeFull = length(dfs[[xName]])      ;# print(sampleSizeFull)
    for (yIdx in seq_len(length(yNames))) {
      cntCensored = sum(as.integer(dfs$isCensored))
      cntHitsPerc = format(100*(1 - cntCensored/sampleSizeFull), digits=4, nsmall=2)
      # print(yIdx) ; print(yNames[yIdx])
      if (cntCensored == sampleSizeFull) {
        if (yIdx == 1) {
          row = c(row, sampleSizeFull, cntHitsPerc, cntProbeLmt, rep(NA, 5))
        } else {
          row = c(row, sampleSizeFull, rep(NA, 5))
        }
        
      } else {
        tmp = dfs[dfs$isCensored == 0 ,][[yNames[yIdx]]] ;# print(tmp)
        if (yIdx == 1) {
          row = c(row,
                  sampleSizeFull,
                  cntHitsPerc, cntProbeLmt,
                  median(as.numeric(tmp)),
                  mean(as.numeric(tmp)),
                  sd(as.numeric(tmp)),
                  min(as.numeric(tmp)),
                  max(as.numeric(tmp))
          )          
        } else {
          row = c(row,
                  median(as.numeric(tmp)),
                  mean(as.numeric(tmp)),
                  sd(as.numeric(tmp)),
                  min(as.numeric(tmp)),
                  max(as.numeric(tmp))
          )
        }
      }
    }
    cat("row =", row,"\n")
    write(row, file=fileStatsPath, ncolumns=length(colNames), append=TRUE, sep="\t")
  }
  cat(sep="", ".. created an asymptotic table of statistics named as = ", fileStatsPath, "\n")
  
  # VIGNETTES
  # > df[df[xName] == xVals[2]  & df$isCensored == 0,]$cntProbe
  # [1] "452" "699" "438" "233"
  # > df[df[xName] == xVals[2]  & df$isCensored == 0,][[yNames[1]]]
  # [1] "452" "699" "438" "233"
  # 
  # > foo = NULL
  # > foo = c(a=1)
  # > foo = c(foo, b=2)
  # > foo = c(foo, c=3)
  # > foo
  # a b c 
  # 1 2 3 
  # > foo$a
  # Error in foo$a : $ operator is invalid for atomic vectors
  # > foo[1]
  # a 
  # 1 
  # > as.numeric(foo[1])
  # [1] 1
  # > foo[["a"]]
  # [1] 1
  # > foo[["b"]]
  # [1] 2
  # > foo[["c"]]
  # [1] 3
  # > 
} ;# file_asym_stats


distance.xyz = function(vec1=c(0,0,0), vec2=c(1,1,1)) 
{
  # Copyright 2018, Franc Brglez
  # euclidian distance of two vectors
  d = length(vec1) ; dist = 0.
  for (i in seq_len(d)) {
    dif = vec1[i] - vec2[i]
    dist = dist + dif*dif
  }
  return(sqrt(dist))
  # VIGNETTE
  # > distance.xyz()
  # [1] 1.732051 # = sqrt(1+1+1)
} # distance.xyz

