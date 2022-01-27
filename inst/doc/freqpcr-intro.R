## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE------------------------------------------------------------
#  library(remotes)
#  install_github("sudoms/freqpcr")
#  
#  library(freqpcr)
#  packageVersion("freqpcr")

## ---- eval = FALSE------------------------------------------------------------
#  ** byte-compile and prepare package for lazy loading
#  Error: (converted from warning) package 'cubature' was built under R version 3.6.3

## ---- eval = FALSE------------------------------------------------------------
#  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
#  install_github("sudoms/freqpcr")

## -----------------------------------------------------------------------------
library(freqpcr)

## ---- results = "hide"--------------------------------------------------------
# Example: R and S are mixed at the exact ratios 1:9, 1:3, 1:1, and 1:0.
# Four dummy bulk samples are generated for each combination. 
# Template DNA amounts follows the gamma distribution.
# K:2, scaleDNA:1e-11, targetScale:1.5, baseChange:0.3, zeroAmount:1e-3,
# sdMeasure:0.3, and EPCR:0.95. 

A <- rep(1, 16)
trueY <- c(rep(0.1, 4), rep(0.25, 4), rep(0.5, 4), rep(1, 4))
housek0 <- c( 19.39, 19.78, 19.28, 19.58,  18.95, 19.91, 19.66, 19.96,
              20.05, 19.86, 19.55, 19.61,  19.86, 19.27, 19.59, 20.21 )
target0 <- c( 19.16, 19.08, 19.28, 19.03,  19.17, 19.67, 18.68, 19.52,
              18.92, 18.79, 18.8, 19.28,   19.57, 19.21, 19.05, 19.15 )
housek1 <- c( 21.61, 21.78, 21.25, 21.07,  22.04, 21.45, 20.72, 21.6,
              21.51, 21.27, 21.08, 21.7,   21.44, 21.46, 21.5, 21.8 )
target1 <- c( 24.3, 24.22, 24.13, 24.13,   22.74, 23.14, 23.02, 23.14,
              21.65, 22.62, 22.28, 21.65,  20.83, 20.82, 20.76, 21.3 )
d.cmp <- data.frame(A, trueY, housek0, target0, housek1, target1)
print(d.cmp)

## -----------------------------------------------------------------------------
# housek0: control sample (intact), housekeeping gene
# target0: control sample (intact), target gene
# housek1: test sample (restriction enzyme), housekeeping gene
# target1: test sample (restriction enzyme), target gene
# trueY  : exact frequency of R allele
# A      : relative concentration of template DNA (not necessary)

result <- knownqpcr(housek0=d.cmp$housek0, target0=d.cmp$target0,
                    housek1=d.cmp$housek1, target1=d.cmp$target1,
                    trueY=d.cmp$trueY, A=d.cmp$A, verbose=FALSE)
print(result)

## -----------------------------------------------------------------------------
# housek1: sample with known ratio, housekeeping gene
# target1: sample with known ratio, target gene, allele specific primer set

result <- knownqpcr(housek1=d.cmp$housek1, target1=d.cmp$target1,
                    trueY=d.cmp$trueY, A=d.cmp$A, verbose=FALSE)
print(result)

## ---- results = "hide"--------------------------------------------------------
A <- rep(1, 16)
trueY <- c(rep(0.1, 4), rep(0.25, 4), rep(0.5, 4), rep(1, 4))
housek0 <- c( 19.39, 19.78, 19.28, 19.58,  18.95, 19.91, 19.66, 19.96,
              20.05, 19.86, 19.55, 19.61,  19.86, 19.27, 19.59, 20.21 )
target0 <- c( 19.16, 19.08, 19.28, 19.03,  19.17, 19.67, 18.68, 19.52,
              18.92, 18.79, 18.8, 19.28,   19.57, 19.21, 19.05, 19.15 )
housek1 <- c( 21.61, 21.78, 21.25, 21.07,  22.04, 21.45, 20.72, 21.6,
              21.51, 21.27, 21.08, 21.7,   21.44, 21.46, 21.5, 21.8 )
target1 <- c( 24.3, 24.22, 24.13, 24.13,   22.74, 23.14, 23.02, 23.14,
              21.65, 22.62, 22.28, 21.65,  20.83, 20.82, 20.76, 21.3 )

# When the combination of (with or without restriction enzyme) is incomplete,
# the R data frame can be prepared in "long" format.
# First, make a complete data frame and then extract a subset.
d.long.all <- data.frame(
    trueY=rep(trueY, 4), 
    Digest=c(rep(0, 16 + 16), rep(1, 16 + 16)),
    Gene=c(rep(0, 16), rep(1, 16), rep(0, 16), rep(1, 16)),
    A=rep(1, 16*4), 
    Cq=c(housek0, target0, housek1, target1) )

# For example, samples without restriction enzyme (Digest == 0) are only available if trueY == 1.
d.long <- d.long.all[d.long.all$Digest == 1 | d.long.all$trueY == 1, ]
print(d.long)

result <- knownqpcr_unpaired(   Digest=d.long$Digest, Gene=d.long$Gene,
                                trueY=d.long$trueY, Cq=d.long$Cq, A=d.long$A   )

## ---- eval = FALSE------------------------------------------------------------
#  > print(d.long)
#     trueY Digest Gene A    Cq
#  13  1.00      0    0 1 19.86
#  14  1.00      0    0 1 19.27
#  15  1.00      0    0 1 19.59
#  16  1.00      0    0 1 20.21
#  29  1.00      0    1 1 19.57
#  30  1.00      0    1 1 19.21
#  31  1.00      0    1 1 19.05
#  32  1.00      0    1 1 19.15
#  33  0.10      1    0 1 21.61
#  34  0.10      1    0 1 21.78
#  35  0.10      1    0 1 21.25
#  36  0.10      1    0 1 21.07
#  37  0.25      1    0 1 22.04
#  38  0.25      1    0 1 21.45
#  39  0.25      1    0 1 20.72
#  40  0.25      1    0 1 21.60
#  41  0.50      1    0 1 21.51
#  ...
#  all following rows are Digest = 1, i.e., qPCR with restricted enzyme.

## -----------------------------------------------------------------------------
targetScale <- 1.2
sdMeasure <- 0.2
scaleDNA <- 1e-06
baseChange <- 0.2

# For the following parameters, we input arbitrary value
P <- 0.15
K <- 4
EPCR <- 0.97
zeroAmount <- 0.0016

## -----------------------------------------------------------------------------
dmy_cq <- make_dummy(   rand.seed=71, P=P, K=K, ntrap=4, npertrap=8,
                        scaleDNA=scaleDNA, 
                        targetScale=targetScale, 
                        baseChange=baseChange,
                        EPCR=EPCR, 
                        zeroAmount=zeroAmount,
                        sdMeasure=sdMeasure, 
                        diploid=FALSE   )
print(dmy_cq)

## ---- results = "hide"--------------------------------------------------------
# Access the keys of an S4 object with @ instead of $. They are extracted as vector.
N <- dmy_cq@N
housek0 <- dmy_cq@housek0
target0 <- dmy_cq@target0
housek1 <- dmy_cq@housek1
target1 <- dmy_cq@target1

# freqpcr() with beta assumption, assuming haploidy, print every optimization step
result <- freqpcr(  N=N, housek0=housek0, target0=target0,
                    housek1=housek1, target1=target1,
                    EPCR=EPCR, zeroAmount=zeroAmount, 
                    beta=TRUE, diploid=FALSE, print.level=2  )
print(result)

## ---- eval = FALSE------------------------------------------------------------
#  > print(result)
#  An object of class "CqFreq"
#  Slot "report":
#                                      Estimate Fixed  (scaled) (scaled.SE)       2.5%       97.5%
#  P (R-allele frequency)            0.09773189     0 -2.222684  0.60914923 0.03178086   0.2633220
#  K (gamma shape parameter)        20.92728471     0  3.041054  1.83522515 0.57354355 763.5884712
#  targetScale (rel. DNA quantity)   1.11922896     0  0.112640  0.08911953 0.93985370   1.3328388
#  sdMeasure (Cq measurement error)  0.20973065     0 -1.561931  0.32845070 0.11017528   0.3992451
#  EPCR (amplification per cycle)    0.97000000     1        NA          NA         NA          NA
#  
#  Slot "obj":
#  $minimum
#  [1] 6.094915
#  
#  $estimate
#  [1] -2.222684  3.041054  0.112640 -1.561931
#  
#  $gradient
#  [1] -3.400973e-05 -8.275889e-05 -5.170087e-05  8.878422e-05
#  
#  $hessian
#              [,1]        [,2]       [,3]       [,4]
#  [1,]  2.71023719  0.05094362   1.168535 -0.1766756
#  [2,]  0.05094362  0.37630056   2.045198 -0.6539723
#  [3,]  1.16853469  2.04519835 147.389558  6.4578190
#  [4,] -0.17667556 -0.65397228   6.457819 11.1504636
#  
#  $code
#  [1] 1
#  
#  $iterations
#  [1] 12
#  
#  
#  Slot "cal.time":
#     user  system elapsed
#     0.72    0.15    0.88

## ---- eval = FALSE------------------------------------------------------------
#  result <- freqpcr(  N=N, housek0=housek0, target0=target0,
#                      housek1=housek1, target1=target1,
#                      K = 1,
#                      EPCR=EPCR, zeroAmount=zeroAmount, beta=TRUE, print.level=1  )

## ---- eval = FALSE------------------------------------------------------------
#  result <- freqpcr(  N=N, housek1=housek1, target1=target1,
#                      targetScale=1.2, sdMeasure=0.2,
#                      EPCR=EPCR, zeroAmount=zeroAmount, beta=TRUE, print.level=1  )

## -----------------------------------------------------------------------------
# Insert the expected size for target1 
x <- 16.2
N <- c(20, 20, 20)
housek0 <- c(20.5, 25.0, 25.6)
target0 <- c(20.5, 25.5, 26.5)
housek1 <- c(25.0, 26.0, 25.5)
target1 <- c(35.3, 26.0 + x, 34.9)
result <- freqpcr(  N=N, housek0=housek0, target0=target0,
                    housek1=housek1, target1=target1,
                    EPCR=0.9, zeroAmount=0.00005, targetScale=0.6,
                    sdMeasure=0.4, beta=TRUE, diploid=TRUE, print.level=0  )

# Try the expected size + 10
x <- 16.2 + 10
target1 <- c(35.3, 26.0 + x, 34.9)
result <- freqpcr(  N=N, housek0=housek0, target0=target0,
                    housek1=housek1, target1=target1,
                    EPCR=0.9, zeroAmount=0.00005, targetScale=0.6, 
                    sdMeasure=0.4, beta=TRUE, diploid=TRUE, print.level=0  )

## -----------------------------------------------------------------------------
citation("freqpcr")

