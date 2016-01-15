# MULTIPLICATIVE UPDATES FOR NMF ---------------------------
# nmf.R
#
# Tue Jan 12 2016
# From Lee and Seung, "Algorithms for NMF", NIPS 2001
#

### Set up parameters
pp <- 5
nn <- 10
rr <- 2
W0 <- matrix(rexp(pp*rr), nrow = pp, ncol = rr)
H0 <- matrix(rexp(nn*rr), nrow = rr, ncol = nn)
V <- W0%*%H0

### Initial guesses
W <- matrix(rexp(pp*rr), nrow = pp, ncol = rr)
H <- matrix(rexp(nn*rr), nrow = rr, ncol = nn)

### Check stationary case
# W <- W0
# H <- H0

### Set up error tracing
res <- matrix(NA, ncol = 3, nrow = iters)
res <- data.frame(res)
names(res) <- c("L2", "L1", "max")

### Run the main algorithm
iters <- 1000
for(i in 1:iters){
    res[i, ] <- c(sqrt(sum((V-W%*%H)^2)), sum(abs(V-W%*%H)), max(V-W%*%H))
    H <- H * (t(W)%*%V) / (t(W)%*%W%*%H)
    W <- W * (V%*%t(H)) / (W%*%H%*%t(H))
}

matplot(res, log="y", type="l")
