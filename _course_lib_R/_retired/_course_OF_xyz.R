# source("../_course_lib_R/_course_OF_xyz.R") ; pi_digits()
# source("_course_OF_xyz.R")  ; aParab_17(17)
# source("_course_OF_xyz.Rr") ; aParab_17s(17)
# source("_course_OF_xyz.R)   ; _course_OF_xyz.R_table()

# OF_xyz_table = function(file="_course_OF_xyz.R_table.txt") {
#   table = read.table(file, header=TRUE, comment.char="#", colClasses = "character")
#  
#   print(table) ; funcName = table[,1] ; print(funcName)
#   return(table)
# }

require(plotrix)  ## https://cran.r-project.org/web/packages/plotrix/plotrix.pdf

pi_digits =function()
{
  # Copyright 2018, Franc Brglez
  # Large number of pi digits
  # .........................
  # from R help about print.default {base}
  #   Note that for large values of digits, currently for digits >= 16, the calculation of 
  #   the number of significant digits will depend on the platform's internal (C library) 
  #   implementation of sprintf() functionality.
  
  options(digits=15) ;# default is options(digits=7)
  signifDigits = 1:16 ; piHat = signif(pi, signifDigits) ; 
  error = as.numeric(format(piHat -pi, digits=4, nsmall=4))
  cat("\npi =", pi ,"  (15 significant digits)\n\n")
  myTable = cbind(signifDigits, piHat, error) ; print(myTable)

  # reduce signifDigits since error = 0 for signifDigits[15]
  xAsym = as.numeric(signifDigits[1:14]) ; yAsym = abs(as.numeric(error[1:14]))  
  # cat(length(xAsym), xAsym, "\n") ; cat(length(yAsym), yAsym, "\n") 

  # plot signifDigits vs abs(error)
  # with a minor hint from http://www.theanalysisfactor.com/r-tutorial-5/
  fileAsymPdf_digits = paste(getwd(), "/fg_asym_pi_digits.pdf", sep="")
  pdf(file = fileAsymPdf_digits, width=7.0, height=7.0, pagecentre=TRUE)
  #print(fileAsymPdf_digits) 
  
  # reduce signifDigits since error reduces to 0 (exactly!!) for signifDigits = 15
  xAsym = signifDigits[1:14] ; yAsym = abs(as.numeric(error[1:14]))
  #cat(length(xAsym), xAsym, "\n") ; cat(length(yAsym), yAsym, "\n")  

  plot(xAsym, yAsym, log="y", 
       xlab="significant digits", ylab="error controlled by significant digits" ) ;# dev.off() ; return()
  
  # get the (exponential) asymptotic model
  modelFit = lm(log10(abs(yAsym)) ~ xAsym) ; print(modelFit)
  co     = coef(modelFit) ; summary(modelFit)
  a      = 10^(co[1])      ; b = co[2]/log10(10)  
  a      = signif(a, 4)   ; b = signif(b, 4)
  base   = 10^b   
  cat("modelAsym =  a*base^xAsym, where: \n")
  print(cbind(a, b, base)) ; cat("\n")
  modelAsym = a*base^xAsym 
  dif = yAsym - modelAsym
  print(cbind(xAsym, yAsym, modelAsym, dif)) ;# return()
  # get a fine grid for ploting modelAsym using "lines"
  xGrid = seq(1,14, 0.1)
  yGrid = a*base^xGrid ;#print(cbind(xGrid, yGrid))
  lines(xGrid, yGrid, lwd=2, col="red") 
  text(6.5,1e-12, "Asymptotic model to predict the maximum error on the")
  text(6.5,1e-13, "estimate of pi, given the number of significant digits")
  text(7.5,1e-14, "   (3.288)*(0.09268298)^significantDigits")

  # template for this table: fg_labs_almostSurely.R
  addtable2plot(7.9, 0.12, myTable, xjust=0,yjust=0, bty="n", hlines=FALSE, vlines=FALSE, 
                title=NULL, cex=0.7 ,xpad=0.15, ypad=0.7)
  yLine=0.050 ; lines(c(7.8, 13.8), c(yLine,yLine), col="black", lwd=0.5)
 
  cat(sep="", "\n.. created the R-plot with file name = ", basename(fileAsymPdf_digits), "\n")
  dev.off()

} ;# pi_digits


pi_newton = function(signifDigits=3) 
{
  # Copyright 2018, Franc Brglez
  # http://mathworld.wolfram.com/PiFormulas.html, formula (20)
  t    = rep(NA, 6)
  t[1] =  (3/4)*sqrt(3) 
  t[2] =  24*(1/12)
  t[3] = -24/(5*2^5) - 1/(28*2^7) - 1/(72*2^9)
  t[4] = -24/(28*2^7) 
  t[5] = -24/(72*2^9) 
  t[6] = -(24*5)/(704*2^11)
  
  piHat = 0 ; numTerms = 0 ; isCensored = TRUE
  tolR = 5/10^(signifDigits+1)      # toleranceRadius
  LB = pi - tolR 
  UB = pi + tolR 
  for (numTerms in seq_len(length(t))) {
    piHat = piHat + t[numTerms]
    if (piHat > LB && piHat < UB) {
      isCensored = FALSE ; break
    }
  }
  error = piHat - pi
  options(digits=signifDigits)
  piHat =signif(piHat, signifDigits)
  return(list(numTerms=numTerms, piHat=piHat, tolR=tolR,  error=error, isCensored=isCensored))
} # pi_newton

pi_BPP = function(signifDigits=3) 
{
  # Copyright 2018, Franc Brglez
  # http://mathworld.wolfram.com/PiFormulas.html, formula (30)
  t = NULL ; tMax = 15
  for (i in seq_len(tMax)) {
    n = i - 1
    t[i] = (4/(8*n + 1) - 2/(8*n + 4) - 1/(8*n + 5) - 1/(8*n + 6))*(1/16)^n
  }
  piHat = 0 ; numTerms = 0 ; isCensored = TRUE
  tolR = 5/10^(signifDigits+1)      # toleranceRadius
  LB = pi - tolR 
  UB = pi + tolR 
  for (numTerms in seq_len(length(t))) {
    piHat = piHat + t[numTerms]
    if (piHat > LB && piHat < UB) {
      isCensored = FALSE ; break
    }
  }
  error = piHat - pi
  options(digits=signifDigits)
  piHat =signif(piHat, signifDigits)
  return(list(numTerms=numTerms, piHat=piHat, tolR=tolR,  error=error, isCensored=isCensored))
  
  # VIGNETTE
  # source("../_course_lib_r/_course_OF_xyz.r") ; pi_BPP(15)
  # $numTerms
  # [1] 11
  # $piHat
  # [1] 3.14159265358979
  # $tolR
  # [1] 5e-16
  # $error
  # [1] 0
  # $isCensored
  # [1] FALSE
  
} # pi_BPP

pi_needles_ClayFord = function(n=100) {
  
  # FROM: http://www.clayford.net/statistics/buffons-needle-problem-or-how-to-use-probability-to-estimate-pi/
  a <- 1 # length of needle
  L <- 1 # distance between lines
  #n <- 100000 # number of dropped needles
  hit <- 0
  for(i in 1:n) {
    x <- runif(1,0,1)
    y <- runif(1,0,1)
    while(x^2 + y^2 > 1) { # no points outside of unit circle (SEE slides to appreciate this advantage)
      x <- runif(1,0,1)
      y <- runif(1,0,1)
    }	
    theta <- atan(y/x) # the random angle
    d <- runif(1,0,(L/2)) # distance of needle midpoint to nearest line
    if(d <= (a/2)*sin(theta)) {
      hit <- hit + 1
    } 
  }
  pi.est <- (n*2)/(hit)
  return(list(needlesDropped=n, piHat=pi.est))
  # VIGNETTES
  # >  source("../_course_lib_r/_course_OF_xyz.r") ; pi_needles_ClayFord(1000)
  # $needlesDropped
  # [1] 1000
  # $piHat
  # [1] 3.062787
  # 
  # > source("../_course_lib_r/_course_OF_xyz.r") ; pi_needles_ClayFord(10000)
  # $needlesDropped
  # [1] 10000
  # $piHat
  # [1] 3.139717
  # 
  # > source("../_course_lib_r/_course_OF_xyz.r") ; pi_needles_ClayFord(100000)
  # $needlesDropped
  # [1] 1e+05
  # $piHat
  # [1] 3.131655
  
} # pi_needles_ClayFord 

pi_needles = function(cntProbeLmt=1000, signifDigits=3, needleRatio=1.0) 
{
  # Copyright 2018, Franc Brglez
  # Laplace's MC experiments with Buffon's needles to approximate pi: horizontal lines only 
  tolR = 5/10^(signifDigits+1)  ;# tolerance radius
  piLB = pi - tolR  
  piUB = pi + tolR  
  cntProbe  = 0 ; isCensored  = TRUE ; error = NA 
  
  while (cntProbe < cntProbeLmt) {
    cntProbe = cntProbe + 1
    hits = 0 ; numerator = NA ; denominator = NA
    for(i in 1:cntProbe) {
      x = runif(1,0,1)
      y = runif(1,0,1)
      while(x^2 + y^2 > 1) { # no points outside of unit circle
        x = runif(1,0,1)
        y = runif(1,0,1)
      }	
      theta <- atan(y/x) # the random angle
      d = runif(1,0,0.5) # distance of needle midpoint to nearest line
      if(d <= (0.5*needleRatio)*sin(theta)) {
        hits = hits + 1
      }
    }
    piMC = (2*cntProbe)/hits
    not_a_hit = piMC < piLB || piMC > piUB
    if ( !not_a_hit ) {
      isCensored  = FALSE
      piMC  = signif(piMC, signifDigits)
      error = formatC(piMC - pi, format="e", digits=3)
      numerator = 2*cntProbe ; denominator = hits
      break
    } else {
      piMC  = signif(piMC, signifDigits+2)
      error = formatC(piMC - pi, format="e", digits=3)
    }
    if (cntProbe >= cntProbeLmt) {isCensored  = TRUE}
  }
  return(list(piMC=piMC, tolR=tolR,  error=error, 
              isCensored=isCensored, numerator=numerator, denominator=denominator, cntProbe=cntProbe))
} # pi_needles

pi_darts = function(cntProbeLmt=10, signifDigits=3) 
{
  # Copyright 2018, Franc Brglez
  # MC experiments with darts to approximate pi: a unit diameter circle is embedded into a unit square
  tolR = 5/10^(signifDigits+1)
  piLB = pi - tolR  
  piUB = pi + tolR  
  cntProbe = 0   
  isCensored  = TRUE
  
  while (cntProbe < cntProbeLmt) {
    cntProbe = cntProbe + 1
    hits = 0
    for (i in seq_len(cntProbe)) {
      x = runif(1,0,1)
      y = runif(1,0,1)
      if (x^2 + y^2 <= 1) {
        hits = hits + 1  
      }	      
    }
    piMC = (4*hits)/cntProbe
    not_a_hit = piMC < piLB || piMC > piUB
    if ( !not_a_hit ) {
      isCensored  = FALSE
      piMC  = signif(piMC, signifDigits)
      error = formatC(piMC - pi, format="e", digits=3)
      numerator = 4*hits ; denominator = cntProbe
      break
    }
    if (cntProbe >= cntProbeLmt) {isCensored  = TRUE}
  }
  return(list(piMC=piMC, tolR=tolR,  error=error, 
              isCensored=isCensored, numerator=numerator, denominator=denominator, cntProbe=cntProbe))
} # pi_darts

aNeedle = function(x) {
  delta   = globPar[["delta"]]               ;# delta
  xTarget = globPar[["aNeedle.upperLmt"]] /2 ;# upperLmt
  if (abs(x - xTarget) < delta) {
    return(0)
  } else {
    return(1)
  }
}

aParab_17_plain = function(x) {
  y   = 1000/17 + (x - 17)^2
  return(y)
} 

aParab_17 = function(x) {

  BKV = globPar[["aParab_17.BKV"]]  ;# Best-Known-Value
  tolY = globPar[["aParab_17.tolY"]]  ;# NOTE: globPar=hash()
  #BKV = BKV + tolY  ;# use for testing only, see comments below return(y)!!
  
  y   = 1000/17 + (x - 17)^2      ;# print(cbind(BKV,y))
  if (y < BKV ) {
    cat("\n** WARNING from function aParab_17**",
    "\n    current BKV =", BKV, "  is NOT the best-known-value",
    "\n   better value =", signif(y, 13), "  has been found for x =", signif(x, 10),
    "\n")
    y = NA  
  } else if (y - BKV < tolY) {
    y = 0
  } else {
    y = y - BKV  
  }
  return(y)
  
  # > signif(1.*2^32, digits=10) - 2^32
  # [1] 0
  # For signif the recognized values of digits are 1...22, and non-missing values are rounded 
  # to the nearest integer in that range. Each element of the vector is rounded individually, 
  # unlike printing.
  
  # > options(digits=10)
  # > source("lib.OF_xyz.r") ; aParab_17(17)
  # ** WARNING from function aParab_17** 
  #   current BKV = 58.82852941   is NOT the best-known-value 
  # better value = 58.82352941   has been found for x = 17 
  # [1] NA
  # > 
  
  # > source("lib.OF_xyz.r") ; aParab_17(17+0.07)
  # [1] 0
  # > source("lib.OF_xyz.r") ; aParab_17(17+0.08)
  # [1] 0.006400000005
  # > 
} # aParab_17

aParab_17s = function(x) {
  
  # NOTE: globPar=hash()
  BKV = globPar[["parab17s.BKV"]]  ;# Best-Known-Value
  tolY = globPar[["aParab_17.tolY"]]       
  #BKV = BKV + tolY  ;# uncomment for testing only, see comments below return(y)!!
  
  y   = 1 + (x - 17)^2 
  if (y < BKV ) {
    cat("\n** WARNING from function aParab_17s**",
        "\n    current BKV =", BKV, "  is NOT the best-known-value",
        "\n   better value =", signif(y, 13), "  has been found for x =", signif(x, 10),
        "\n")
    y = NA
  } else if (y - BKV < tolY) {
    y = 0
  } else {
    y = y - BKV  
  }
  return(y)
  
  # > aParab_17s(17)
  # ** WARNING from function aParab_17s** 
  #   current BKV = 1.005   is NOT the best-known-value 
  # better value = 1   has been found for x = 17 
  # [1] NA
  # > 
  
  # > aParab_17s(17+0.07)
  # [1] 0
  # > aParab_17s(17+0.08)
  # [1] 0.0064
  # > 
} # aParab_17s

wildN = function(x) { 
  # NOTE: globPar=hash()
  BKV         = globPar[["wildN.BKV"]]  ;# Best-Known-Value
  tolY        = globPar[["wildN.tolY"]]
  isValueOnly = globPar[["wildN.isValueOnly"]]
  
  # x = c(x1, x2, ... , xd)
  d   = length(x)
  y = sum(10 * sin(0.3 * x) * sin(1.3 * x^2) + 0.00001 * x^4 + 0.2 * x + 80 ) 
  y = signif(y, 7)
  if (isValueOnly) {
    print(cbind(x, y, BKV, tolY))
    return(y)
  }
  yBest = d*BKV
  if (y < d*BKV ) {
    cat("\n** WARNING from function wild**",
        "\n    current BKV =", BKV, "  is NOT the best-known-value",
        "\n   better value =", signif(y, 7), "  has been found for x =", signif(x, 7),
        "\n")
    y = NA
  } 
  # since BKV != 0, we use relative error 
  error = abs((y - yBest)/yBest)  ; globPar[["wildN.error"]] = error
  if (error < tolY) {
    y = 0
  } else {
    y = error  
  }
  return(y)
  
  # for d=1, this function is an example from the R package DEoptim-2.2-4
  # it seems to have a unique best-known-value (BKV = 67.46773) solution at 
  # x* = -15.81515 not only within [-50, 50] but also up to [-1e5, 1e5]
  # seedInit = 5381 ; NP = 64 ; interval = -+ 50000 ; itermax = 200 
  # > wild(-15.81515537)
  # [1] 67.46773489
} # wildN

rastriginN = function(x) 
{
  # NOTE: globPar=hash()
  BKV  = globPar[["rastrigin.BKV"]]  ;# Best-Known-Value
  tolY = globPar[["rastrigin.tolY"]]
  isValueOnly = globPar[["rastrigin.isValueOnly"]]
  # x = c(x1, x2, ... , xd)
  d   = length(x)
  y    = 10*d + sum(x^2 - 10*cos(2*pi*x))
  y = signif(y, 7)
  if(isValueOnly) {
    print(cbind(x, y, BKV, tolY))
    return(y)
  }
  # since BKV == 0, we use absolute error 
  error = abs(y - d*BKV)
  if (error < tolY) {
    y = 0
  } else {
    y = error  
  }
  return(y)
} # rastriginN

ModRosenbrockN = function(x) 
{ 
  y = 100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2 
  return(y)
} # ModRosenbrock

trefethenN = function(x) 
{
  # copied from testFunctions in {NMOF}
  # Trefethen’s function takes a two-dimensional x (here written as x and y); it is defined as
  # exp(sin(50x))+sin(60exp(y)+sin(70sin(x))+sin(sin(80y))−sin(10(x+y))+ (x2 +y2)/4
  # The minimum function value is -3.3069; reached at c(-0.0244, 0.2106).
  # > x[1] = -0.0244 ; x[2] = 0.2106
  # > exp(sin(50 * x[1])) + sin(60 * exp(x[2])) + sin(70 * sin(x[1])) + sin(sin(80 * x[2])) - sin(10 * (x[1] + x[2])) + (x[1]^2 + x[2]^2)/4
  # [1] -3.306868
  # > options(digits=11)
  # > exp(sin(50 * x[1])) + sin(60 * exp(x[2])) + sin(70 * sin(x[1])) + sin(sin(80 * x[2])) - sin(10 * (x[1] + x[2])) + (x[1]^2 + x[2]^2)/4
  # [1] -3.3068678585
  # http://infinity77.net/globPar_optimization/test_functions_nd_T.html#go_benchmark.Trefethen
  y = exp(sin(50 * x[1])) + sin(60 * exp(x[2])) + sin(70 * sin(x[1])) + sin(sin(80 * x[2])) - sin(10 * (x[1] + x[2])) + (x[1]^2 + x[2]^2)/4
  return(y)
} # trefethenN

clusterLJ = function(x) 
{
  # computing the normalized value of Lennard-Jones potential
  # https://en.wikipedia.org/wiki/Lennard-Jones_potential
  nPar   = length(x)
  nDim   = 3 
  nAtoms = as.integer(nPar/nDim) 
  coords = matrix(x, ncol=nDim, byrow=T) ;# print(coords) 
  
  pot = 0 # initializing the value of LJ potential
  for (i in seq_len(nAtoms - 1)) {
    for (j in (i+1):nAtoms) {
      v_i = coords[i,] ; v_j = coords[j,]   
      r   = distance.xyz(v_i, v_j)
      pot = pot + 4*(r^(-12) - r^(-6))
      # cat(sep=" ", "\nv_", i, " = ", v_i, "\nv_", j, " = ", v_j, 
      #     "\nr =", r, " ; pot =", pot, "\n" )
    }
  }
  # create a canonical form by normalization wrt BKV and the aParab_17.tolY threshold
  y = pot - clusterLJ.bkv(nAtoms)
  if (y < aParab_17.tolY()) {y = 0}
  return(y)
  # tolY = 0.01
  # diff = pot - clusterLJ.bkv[nAtoms]
  # # if (diff < 0) {cat("*********** diff =", diff, " ***********\n")}
  # # stopif(diff < 0) # found a solution better than BKV
  # if (diff <= tolY) {
  #   pot =  BKV[nAtoms]
  # }
  # return(pot)
} # clusterLJ


schwefel26N = function(x) 
{ 
  # x = c(x1, x2, ..., xd)
  d = length(x) 
  BKV = 418.9829 
  sum  = sum(x*sin(sqrt(abs(x)))) 
  y    = sum  + BKV*d  
  if (y < aParab_17.tolY()) {y = 0}
  return(y)
  # y* = 0 at x* = c(420.9687, ..., 420.9687) for xi = [-500, 500] for all i =1, 2, ..., d
} # schwefel26N
