################################################################################
# Moore's test for a common distribution for paired samples
# (Code taken from Andrews et al.)
################################################################################

MooreRStats <- function(ldat1, ldat2) { 
  x <- cos(ldat1)-cos(ldat2); y <- sin(ldat1)-sin(ldat2)
  r <- sqrt((x*x)+(y*y)) ; Ranks <- rank(r)
  cosphi <- x/r; sinphi <- y/r 
  return(list(cosphi, sinphi, Ranks))
}

MooreRTestStat <- function(cosphi, sinphi, Ranks){ 
  n <- length(cosphi)
  RbarC <- (1/n)*sum(Ranks*cosphi); RbarS <- (1/n)*sum(Ranks*sinphi)
  Rval <- sqrt(((RbarC*RbarC)+(RbarS*RbarS))/n) ; return(Rval)
}

MooreRTestRand <- function(cosphi, sinphi, Ranks, NR) {
  RObs <- MooreRTestStat(cosphi, sinphi, Ranks) ; nxtrm <- 1
  n <- length(cosphi)
  for (r in 1:NR) {
    cosphirand <- 0 ; sinphirand <- 0
    for (j in 1:n) {
      if (runif(1) < 0.5) {
        cosphirand[j] <- cosphi[j] ; sinphirand[j] <- sinphi[j] }
      else {
        cosphirand[j] <- -cosphi[j] ; sinphirand[j] <- -sinphi[j] } }
    RRand <- MooreRTestStat(cosphirand, sinphirand, Ranks)
    if (RRand >= RObs) { nxtrm <- nxtrm+1 }
  }
  pval <- nxtrm/(NR+1) ; return(c(RObs, pval))
}


################################################################################
# Wallraff's test for a common concentration
################################################################################

WalraffTest <- function(cdat, ndat, g) {
  N <- length(cdat) ; ndatcsum <- cumsum(ndat) ; tbar <- circular(0) ; distdat <- 0
  for (k in 1:g) {
    dist <- 0 ; sample <- circular(0)  
    if (k==1) {low <- 0} else
      if (k > 1) {low <- ndatcsum[k-1]}
    for (j in 1:ndat[k]) { sample[j] <- cdat[j+low] }
    tm1 <- trigonometric.moment(sample, p=1) ; tbar[k] <- tm1$mu
    for (j in 1:ndat[k]) { dist[j] <- pi-abs(pi-abs(sample[j]-tbar[k])) }
    distdat <- c(distdat, dist)
  }
  distdat <- distdat[-1]
  gID <- c(rep(1,39), rep(2,39))
  TestRes <- kruskal.test(distdat, g=gID)
  return(TestRes)
} 

