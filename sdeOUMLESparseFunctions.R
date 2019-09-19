# Aim: calulate likelihood and estimate MLE using optim
# of OU process
# dX = - gamma X dt + sigma dB


# - 1D time series
# - breaks
# - spline (a, b, c) values for each break
# - gamma0
# - sigma0
# Output: 
# - gamma and sigma estimates

logLikelihoodOU <- function(X, t, gamma, sigma) {
        n <- length(X)
        deltaT <- t[2] - t[1] # assume equal spacing
        termOne <- (n / 2) * (log(gamma) - 
                              log(pi * sigma^2 * (1 - exp(-2 * gamma * deltaT))))
        termTwo <- gamma / (sigma^2 * (1 - exp(-2 * gamma * deltaT))) *
                    sum((X[2:n] - X[1:(n-1)] * exp(-gamma * deltaT))^2)
        logL <- termOne - termTwo
        logL
}

.dldgamma <- function(X, t, gamma, sigma) {
        n <- length(X)
        deltaT <- t[2] - t[1] # assume equal spacing
        termOne <- (n / 2) * (1 - exp(-2 * gamma * deltaT) - 
                              2 * gamma * deltaT * exp(-2 * gamma * deltaT)) / 
                                (gamma * (1 - exp(-2 * gamma * deltaT)))
        termTwo1 <- 2 * gamma * deltaT * exp(-gamma * deltaT) * 
                        sum(X[1:(n - 1)] * (X[2:n] - X[1:(n - 1)] * exp(-gamma * deltaT)))
        termTwo2 <- (1 - exp(-2 * gamma * deltaT) * (1 + 2 * sigma^2 * deltaT)) / 
                        (1 - exp(-2 * gamma * deltaT))
        termTwo <- (1 / (sigma^2 * (1 - exp(-2 * gamma * deltaT)))) * 
                        (termTwo1 - termTwo2)
        res <- termOne - termTwo
        res
}

.dldsigma <- function(X, t, gamma, sigma) {
        n <- length(X)
        deltaT <- t[2] - t[1] # assume equal spacing
        termOne <- -n / 2 * pi
        termTwo1 <- 2 * gamma * 
                        sum((X[2:n] - X[1:(n - 1)] * exp(-gamma * deltaT))^2)
        termTwo2 <- sigma^3 * (1 - exp(-2 * gamma * deltaT))
        termTwo <- termTwo1 / termTwo2
        res <- termOne + termTwo
        res
}

mleOptim <- function(X, t, gamma0, sigma0) {
        n <- length(X)
        deltaT <- t[2] - t[1] # assume equal spacing
        if(length(gamma0) == 1 & length(sigma0) == 1) {
                opt <- optim(par = c(gamma0, sigma0), function(gs) {
                                     gamma <- gs[1]
                                     sigma <- gs[2]
                                     resTmp <- logLikelihoodOU(X = X, t = t, 
                                                               gamma = gamma, sigma = sigma)
                                     -resTmp
                                      }, method = "BFGS", control = list(maxit = 500))
        } else {
                lenG <- length(gamma0)
                lenS <- length(sigma0)
                gamma0 <- rep(gamma0, each = lenS)
                sigma0 <- rep(sigma0, length = lenG * lenS)
                gs0 <- cbind(gamma0, sigma0)
                optAll <- apply(gs0, 1, function(par) {
                        opt <- optim(par = par, function(gs) {
                                     gamma <- gs[1]
                                     sigma <- gs[2]
                                     resTmp <- logLikelihoodOU(X = X, t = t, 
                                                               gamma = gamma, sigma = sigma)
                                     -resTmp
                                      }, method = "BFGS", control = list(maxit = 500))
                      })
                optValue <- sapply(optAll, function(o) o$value)
                whichMinValue <- which.min(optValue)
                opt <- optAll[[whichMinValue]]
        }
        names(opt$par) <- c("gammaHat", "sigmaHat")
        opt
}

.muT <- function(t, spline, breaks) {
        N <- length(t) - 1 
        numBreaks <- nrow(spline)
        maxBreak <- breaks[numBreaks + 1]
        deltaT <- t[2:(N+1)] - t[1:N]
        maxT <- max(t)
        # split t up into relevant breaks
        breakMultiplier <- max(t) / maxBreak
        fullBreaks1 <- rep(breaks, length = (length(breaks) + 1) * breakMultiplier)
        fullBreaks2 <- rep(0:ceiling(breakMultiplier), each = length(breaks), 
                          length = (length(breaks) + 1) * breakMultiplier) * 
                                (maxBreak + 1)
        fullBreaks <- fullBreaks1 + fullBreaks2 
        fullBreaks <- fullBreaks[1:min(which(fullBreaks > maxT))]
        whichTBreaks <- lapply(seq_len(length(fullBreaks) - 1), function(i) {
                                  bF <- fullBreaks[i+1]
                                  bB <- fullBreaks[i]
                                  whichT <- which(t < bF & t >= bB)
                                  whichT
                        })
        degree <- c(2, 1, 0)
        muTBreaks <- lapply(seq_len(length(fullBreaks) - 1), function(i) {
                                  w <- whichTBreaks[[i]]
                                  tmpSeq <- t[w]
                                  tmpSeqMod <- tmpSeq %% maxBreak
                                  xQuadratic <- outer(tmpSeqMod, degree, "^")
                                  iMod <- i %% maxBreak
                                  if (iMod == 0) iMod <- maxBreak
                                  splineTmp <- spline[iMod,]
                                  muTmp <- as.vector(tcrossprod(splineTmp, xQuadratic))
                                  list(t = tmpSeq, mu = muTmp)
                        })
        muTBreaks
}

mleSplineOptim <- function(Y, t, spline, breaks, gamma0, sigma0) {
        n <- length(Y)
        muT <- .muT(t = t, spline = spline, breaks = breaks)
        muT <- unlist(lapply(seq_len(length(muT)), function(i) {
                                     muT[[i]]$mu
                        }))
        X <- Y - muT
        opt <- mleOptim(X = X, t = t, gamma0 = gamma0, sigma0 = sigma0) 
        opt
}
