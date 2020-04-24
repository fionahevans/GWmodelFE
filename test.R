

library(microbenchmark)

###
### Bugs in GWmodel function gwr.multiscale when only one predictor ----
###

data(LondonHP)
EUDM <- gw.dist(coordinates(londonhp))

res2 <- GWmodel::gwr.multiscale(PURCHASE ~ FLOORSZ + PROF, 
                                data = londonhp, 
                                criterion="dCVR",
                                kernel = "gaussian", 
                                adaptive = T, 
                                bws0 = c(100, 100, 100),
                                bw.seled = rep(T, 3), 
                                dMats = list(EUDM,EUDM,EUDM))

res3 <- GWmodelFE::gwr.multiscale(PURCHASE ~ FLOORSZ + PROF, 
                                  data = londonhp, 
                                  criterion="dCVR",
                                  kernel = "gaussian", 
                                  adaptive = T, 
                                  bws0 = c(100, 100, 100),
                                  bw.seled = rep(T, 3), 
                                  dMats = list(EUDM,EUDM,EUDM))

# These are the same - my version has fixed the bug!!!

###
### GWmodel function gwr.mixed is excuciatingly slow ----
###
### Goal is to re-code into C++ to improve speed
###
### New functions:
###  gwr.q.fast
###  gwr.mixed.trace.fast
###  gwr.mixed.2.fast
###  gwr.mixed.dMat
###
### Note: To ease coding overheads, only works with dMat supplied  and
###       no regression points.
###

###
### Load data for testing----
###

load("testData.RData")

dat <- data.linear[1:200,]
DMdat <- gw.dist(dp.locat=coordinates(dat)) 


###
### Code adaptive weighting functions into C++ and test ----
###

# Gaussian
x0 <- GWmodel::gw.weight.gau.ad(DM[,1], 10)
x1 <- GWmodelFE::gau.wt.vec.ad(DM[,1], 10) 
plot(x0, x1)
abline(0, 1, col="red") 
sum(x0 == x1) / length(x0)*100 # 100% correct

microbenchmark(
  GWmodel::gw.weight.gau.ad(DM[,1], 10),
  GWmodelFE::gau.wt.vec.ad(DM[,1], 10) 
)  # same speed

# Bisquare
x0 <- GWmodel::gw.weight.bis.ad(DM[,1], 10) 
x1 <- GWmodelFE::bis.wt.vec.ad(DM[,1], 10) 
plot(x0, x1)
abline(0, 1, col="red") 
sum(x0 == x1) / length(x0)*100 # 100% correct

# Tricube
x0 <- GWmodel::gw.weight.tri.ad(DM[,1], 10) 
x1 <- GWmodelFE::tri.wt.vec.ad(DM[,1], 10) 
plot(x0, x1)
abline(0, 1, col="red") 
sum(x0 == x1) / length(x0)*100 # 100% correct

# Exponential
x0 <- GWmodelFE::gw.weight.exp.ad(DM[,1], 10) # wasn't exported from GWmodel 
x1 <- GWmodelFE::exp.wt.vec.ad(DM[,1], 10) 
plot(x0, x1)
abline(0, 1, col="red") 
sum(x0 == x1) / length(x0)*100 # 100% correct

# Boxcar
x0 <- GWmodel::gw.weight.box.ad(DM[,1], 10) 
x1 <- GWmodelFE::box.wt.vec.ad(DM[,1], 10) 
plot(x0, x1)
abline(0, 1, col="red") 
sum(x0 == x1) / length(x0)*100 # 100% correct

microbenchmark(
  GWmodel::gw.weight.box.ad(DM[,1], 10),
  GWmodelFE::box.wt.vec.ad(DM[,1], 10) 
)  # not much faster

###
### Code gwr.q into C++ and test ----
###

# Test gwr.q
test <- function(kernel, bw, adaptive) {
  t0 <- system.time(z0 <- GWmodel::gwr.q(matrix(data.linear$rate, ncol=1), 
               matrix(data.linear$yield, ncol=1), 
               loc=coordinates(data.linear), dMat=DM, bw=bw, kernel=kernel, adaptive=adaptive))
  t1 <- system.time(z1 <- GWmodelFE::gwr.q.fast(matrix(data.linear$rate, ncol=1), 
               matrix(data.linear$yield, ncol=1), 
               loc=coordinates(data.linear), dMat=DM, bw=bw, kernel=kernel, adaptive=adaptive))
  plot(z0, z1)
  abline(0, 1, col="red")
  list(times = c(t0[3], t1[3]), percentCorrect = sum(z0 == z1) / length(z0)*100,
       cor=cor(z0, z1))
}

test("gaussian", 4.5, F)      # 100% correct
test("exponential", 4.5,  F)  # 100% correct
test("bisquare", 4.5,  F)     # 100% correct
test("tricube", 4.5,  F)      # 100% correct
test("boxcar", 4.5,  F)       # 100% correct

test("gaussian", 20, T)       # 100% correct
test("exponential", 20,  T)   # 100% correct
test("bisquare", 20,  T)      # 100% correct
test("tricube", 20,  T)       # 100% correct
test("boxcar", 20, T)         # 100% correct

dat <- data.linear[1:200,]
DMdat <- gw.dist(dp.locat=coordinates(dat)) 

microbenchmark(
  GWmodel::gwr.q(matrix(dat$rate, ncol=1), 
                 matrix(dat$yield, ncol=1), 
                 loc=coordinates(dat), dMat=DMdat, bw=4.5, kernel="gaussian", adaptive=F),
  GWmodelFE::gwr.q.fast(matrix(dat$rate, ncol=1), 
                 matrix(dat$yield, ncol=1), 
                 dMat=DMdat, bw=4.5, kernel="gaussian", adaptive=F)
) 
# min       lq     mean   median       uq       max neval cld
# 4.682761 4.939308 6.125318 5.323308 6.148745 26.944927   100   b
# 2.458624 2.586441 2.978486 2.691467 3.027695  7.787583   100  a 
# Around twice as fast

microbenchmark(
  GWmodel::gwr.q(matrix(dat$rate, ncol=1), 
                 matrix(dat$yield, ncol=1), 
                 loc=coordinates(dat), dMat=DMdat, bw=4.5, kernel="boxcar", adaptive=T),
  GWmodelFE::gwr.q.fast(matrix(dat$rate, ncol=1), 
                   matrix(dat$yield, ncol=1), 
                   dMat=DMdat, bw=4.5, kernel="boxcar", adaptive=T)
) 
# min        lq      mean    median        uq       max neval cld
# 17.624991 18.431464 21.840058 19.068000 23.569153 45.046544   100   b
# 3.044652  3.215501  3.428138  3.282053  3.444515  6.192141   100  a 
# Around 5 times as fast




###
### Code gwr.mixed.trace into C++ and test ----
###

# Check results are the same
GWmodel::gwr.mixed.trace(x1 = matrix(rep(1, nrow(dat))), 
                         x2 = matrix(dat$rate, ncol=1),
                         y = matrix(dat$yield, ncol=1), 
                         loc = coordinates(dat), dMat = DMdat, bw = 4.5, 
                         kernel = "gaussian", adaptive = F)


GWmodelFE::gwr.mixed.trace.fast(matrix(rep(1, nrow(dat))), 
                           x2 = matrix(dat$rate, ncol=1),
                           y = matrix(dat$yield, ncol=1), 
                           dMat = DMdat, bw = 4.5, 
                           kernel = "gaussian", adaptive = F)


# Benchmark
microbenchmark(
  GWmodel::gwr.mixed.trace(x1 = matrix(rep(1, nrow(dat))), 
                           x2 = matrix(dat$rate, ncol=1),
                           y = matrix(dat$yield, ncol=1), 
                           loc = coordinates(dat), dMat = DMdat, bw = 4.5, 
                           kernel = "gaussian", adaptive = F),
  GWmodelFE::gwr.mixed.trace.fast(matrix(rep(1, nrow(dat))), 
                             x2 = matrix(dat$rate, ncol=1),
                             y = matrix(dat$yield, ncol=1), 
                             dMat = DMdat, bw = 4.5, 
                             kernel = "gaussian", adaptive = F),
  times = 50
) 
# min       lq     mean   median       uq      max neval cld
# 5.124067 5.506247 5.849987 5.809478 6.049470 7.577778    50   b
# 1.147785 1.190027 1.247828 1.215158 1.274302 1.686457    50  a 
# Around 5 times faster


###
### Code gwr.mixed.2 into C++ and test ----
###

# Check results are the same
GWmodel::gwr.mixed.2(x1 = matrix(rep(1, nrow(dat))), 
                     x2 = matrix(dat$rate, ncol=1),
                     y = matrix(dat$yield, ncol=1), 
                     loc = coordinates(dat), dMat = DMdat, bw = 4.5, 
                     kernel = "gaussian", adaptive = F)


GWmodelFE::gwr.mixed.2.fast(x1 = matrix(rep(1, nrow(dat))), 
                       x2 = matrix(dat$rate, ncol=1),
                       y = matrix(dat$yield, ncol=1), 
                       dMat = DMdat, bw = 4.5, 
                       kernel = "gaussian", adaptive = F)

# Benchmark
microbenchmark(
  GWmodel::gwr.mixed.2(x1 = matrix(rep(1, nrow(dat))), 
                           x2 = matrix(dat$rate, ncol=1),
                           y = matrix(dat$yield, ncol=1), 
                           loc = coordinates(dat), dMat = DMdat, bw = 4.5, 
                           kernel = "gaussian", adaptive = F),
  GWmodelFE::gwr.mixed.2.fast(matrix(rep(1, nrow(dat))), 
                                  x2 = matrix(dat$rate, ncol=1),
                                  y = matrix(dat$yield, ncol=1), 
                                  dMat = DMdat, bw = 4.5, 
                                  kernel = "gaussian", adaptive = F),
  times = 50
) 

# min       lq     mean   median       uq       max neval cld
# 50.44114 53.18531 63.76629 59.91224 67.58841 106.28309    50   b
# 10.42125 11.09516 12.00038 11.26182 12.64612  25.91655    50  a
# Around 5 times faster

###
### Time improvement aftering re-coding gwr.mixed.2 ----
###

# Check results are the same
res0 <- GWmodel::gwr.mixed(yield ~ rate, 
                           data = dat, 
                           bw = 4.5,
                           fixed.vars = "rate",
                           kernel = "gaussian",
                           adaptive = FALSE,
                           diagnostic = T,
                           dMat = DMdat)

res1 <- GWmodelFE::gwr.mixed.dMat(yield ~ rate, 
                             data = dat, 
                             bw = 4.5,
                             fixed.vars = "rate",
                             kernel = "gaussian",
                             adaptive = FALSE,
                             dMat = DMdat)

# Check results are the same (with regression.points)
DMdatR <- GWmodel::gw.dist(dp.locat=coordinates(dat), rp.locat=coordinates(dat)[1:100,])

res0 <- GWmodel::gwr.mixed(yield ~ rate, 
                           data = dat, 
                           regression.points = coordinates(dat)[1:100,],
                           bw = 4.5,
                           fixed.vars = "rate",
                           kernel = "gaussian",
                           adaptive = FALSE,
                           diagnostic = T,
                           dMat = DMdatR)

res1 <- GWmodelFE::gwr.mixed.dMat(yield ~ rate, 
                                  data = dat, 
                                  regression.points = coordinates(dat)[1:100,],
                                  bw = 4.5,
                                  fixed.vars = "rate",
                                  kernel = "gaussian",
                                  adaptive = FALSE,
                                  dMat = DMdat,
                                  dMat.rp = DMdatR)

GWmodelFE::gwr.mixed.AICc(yield ~ rate, 
               data = dat,
               bw = 4.5,
               fixed.vars = "rate",
               kernel = "gaussian",
               adaptive = FALSE,
               dMat = DMdat)

# Benchmark
microbenchmark(
  GWmodel::gwr.mixed(yield ~ rate, data = dat,bw = 4.5, fixed.vars = "rate",
                     kernel = "gaussian", adaptive = FALSE, diagnostic = T,
                     dMat = DMdat),
  GWmodelFE::gwr.mixed.dMat(yield ~ rate, data = dat,bw = 4.5, fixed.vars = "rate",
                       kernel = "gaussian", adaptive = FALSE, diagnostic = T,
                       dMat = DMdat),
  times = 50
) 
# min       lq     mean   median       uq      max neval cld
# 5.327581 5.629210 6.013973 5.968463 6.193683 7.955761    50   b
# 1.168743 1.214306 1.285892 1.250883 1.323195 1.596110    50  a 

# Around 5 times faster!
