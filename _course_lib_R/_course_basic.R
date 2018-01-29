# source("../_course_lib_R/_course_basic.R")

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

printVecAsis <- function(x) {
  # Copyright 2018
  # https://stackoverflow.com/questions/10860020/output-a-vector-in-r-in-the-same-format-used-for-inputting-it-into-r
  
  ifelse(length(x) == 1, x, 
    ifelse(is.character(x), 
    paste0("c(", paste(sapply(x, function(a) paste0("\'",a,"\'")), collapse=", "), ")"),
    paste0("c(", paste(x, collapse=", "), ")"))
  )
  # VIGNETTE
  # > x=1:5 ; x
  # [1] 1 2 3 4 5
  # > printVecAsis(x)
  # [1] "c(1, 2, 3, 4, 5)"
} ;# printVecAsis


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
  #print(df) ; x=df[[xName]] ; y1=df[[yNames[1]]] ; y2=df[[yNames[2]]] ; print(cbind(x,y1,y2)) ;#  return() 
  
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
  colNames = c(xName, "sampleSize", "cntHits", "cntHitsPerc", "cntProbeLmt")  
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
  cat("values associated with the asymptotic variable(", xName, ") = ",xVals, "\n")
  for (xIdx in seq_len(length(xVals))) {
    row = xVals[xIdx]
    dfs = df[df[[xName]] == xVals[xIdx] ,]     ;# print(dfs)
    sampleSizeFull = length(dfs[[xName]])      ;# print(sampleSizeFull)
    cat("\nxIdx =", xIdx, "xName =", xName, "sampleSizeFull =", sampleSizeFull, "\n"); print(row) ; print(dfs) 
    for (yIdx in seq_len(length(yNames))) {
      cntCensored = sum(as.logical(dfs$isCensored))
      cntHits     = sampleSizeFull - cntCensored
      cntHitsPerc = format(100*(cntHits/sampleSizeFull), digits=4, nsmall=2)
      # print(yIdx) ; print(yNames[yIdx])
      if (cntCensored == sampleSizeFull) {
        if (yIdx == 1) {
          row = c(row, sampleSizeFull, cntHits, cntHitsPerc, cntProbeLmt, rep(NA, 5))
        } else {
          row = c(row, sampleSizeFull, rep(NA, 5))
        }
        
      } else {
        tmp = dfs[dfs$isCensored == FALSE ,][[yNames[yIdx]]] ;# print(tmp)
        if (yIdx == 1) {
          row = c(row,
                  sampleSizeFull,
                  cntHits, cntHitsPerc, cntProbeLmt,
                  format(median(as.numeric(tmp)), digits=2, nsmall=2), # since tmp=cntProbe
                  format(mean(as.numeric(tmp)), digits=2, nsmall=2),
                  format(sd(as.numeric(tmp)),  digits=2, nsmall=2),
                  format(min(as.numeric(tmp)), digits=2, nsmall=2),
                  format(max(as.numeric(tmp)), digits=2, nsmall=2)
          )          
        } else {
          row = c(row,
                  format(median(as.numeric(tmp)), digits=4, format="e"),  # since tmp=runtime_seconds
                  format(mean(as.numeric(tmp)), digits=4, format="e"),
                  format(sd(as.numeric(tmp)), digits=4, format="e"),
                  format(min(as.numeric(tmp)), digits=4, format="e"),
                  format(max(as.numeric(tmp)), digits=4, format="e")
          )
        }
      }
    }
    cat("row =", row,"\n")
    write(row, file=fileStatsPath, ncolumns=length(colNames), append=TRUE, sep="\t")
  }
  cat(sep="", ".. created an asymptotic table of statistics named as = ", basename(fileStatsPath), "\n")
  # convert the file into easier-to-read html file:
  table2html(fileStatsPath)
  
  # VIGNETTE (FB)
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


table2html = function(filepath) 
{
  # Copyright 2018, George Mathew, george2@ncsu.edu
  
  make_table_cell = function(x) {paste("<td>",x,"</td>", sep="")};
  make_table_header = function(x) {paste("<th>",x,"</th>", sep="")};
  style = "<style> td,th {border: 1px solid black;} </style>";
  file_root      = strsplit(filepath, split = ".", fixed = TRUE)[[1]][1];
  file_extension = strsplit(filepath, split = ".", fixed = TRUE)[[1]][2];
  file_html      = paste(file_root,".", "html", sep="")
  body = "";
  rows = "";
  first_row_seen = FALSE;
  con = file(filepath, "r");
  while ( TRUE ) {
    line = readLines(con, n = 1);
    if ( length(line) == 0 ) {
      html = paste0("<html>", style,"<div>", body, "</div><br/>", "<table>", rows, "</table>" ,"</html>", sep="\n");
      fileConn<-file(paste(file_root, ".html", sep=""));
      writeLines(html, fileConn);
      close(fileConn);
      break;
    }
    first_char = substr(line, 1, 1);
    if (first_char == "#") {
      body = paste(body, line, sep = "<br/>");
    } else {
      if (first_row_seen) {
        row = lapply(strsplit(line, "\t"), make_table_cell);
      } else {
        row = lapply(strsplit(line, "\t"), make_table_header);
        first_row_seen = TRUE;
      }
      row = paste(row[[1]], collapse = " ");
      row = paste("<tr>", row, "</tr>", sep="");
      rows = paste(rows, row, sep="");
    }
  }
  close(con);
  cat(sep="", ".. table2html: created a file named as = ", basename(file_html), "\n");
} # table2html
# table2html('fg_asym_pi_darts_unfinished.txt')

distance_xyz = function(vec1=c(0,0,0), vec2=c(1,1,1)) 
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
} # distance_xyz

bin2dec = function(x="10111")
{
  dec = strtoi(x, base=2)
  return(dec)
  # VIGNETTE (FB)
  # > bin2dec("11111")
  # [1] 31
  # > bin2dec("10111")
  # [1] 23
  # > bin2dec("10001")
  # [1] 17
} # bin2dec

int2bin<- function(x, ndigits)
{
  # Copyright, 2011, https://stat.ethz.ch/pipermail/r-help/2011-December/297262.html
  # FB to check for speed with int2bin from rnn-package
  base <- array(NA, dim=c(length(x), ndigits))
  for (q in 1:ndigits) {
    base[, ndigits-q+1] <- (x %% 2)
    x <- (x %/% 2)
  }
  bin <- apply(base,1,paste,collapse="")
  return(bin)
  # VIGNETTE (FB)
  # > int2bin(7, 4)
  # [1] "0111"
  # > int2bin(15, 4)
  # [1] "1111"
  # > int2bin(3, 5)
  # [1] "00011"
  # > int2bin(4, 5)
  # [1] "00100"
  # > 
} # int2bin

# bin2dec = function(coordB="10111", xLB=-4, dx=0.25) 
# {
#   coordIdx = strtoi(coordB, base=2)
#   coordD   = xLB + coordIdx*dx
#   #print(cbind(coordB, coordD, coordIdx, dx))
#   return(coordD)
# } # bin2dec




