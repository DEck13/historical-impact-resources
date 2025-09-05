

########################################
## Full House Model for Pantheon data
########################################

## load in data set
persons2020 <- read_csv("person_2020_update.csv")
head(persons2020)
colnames(persons2020)

## get a subset variables
persons2020_small <- persons2020 %>% dplyr::select(name, gender, age, birthyear, 
                                            hpi, l_, non_en_page_views, coefficient_of_variation) %>% 
  mutate(A = log(age, base = 4) - (70-age)/7 *ifelse(age < 70,1,0))

## investigate birthyear points 
quantile(persons2020_small %>% 
           ## birthyear 1973 chosen so that a person has to be 50+ to be eligible
           filter(birthyear <= 1973) %>% 
           dplyr::select(A), 
         probs = c(0.01,0.10,0.25,0.50,0.75,0.90,0.99),
         na.rm = TRUE)



library(parallel)
nCores <- detectCores() - 2

## calculate number of bios in a window around a target year
n_calc <- function(y, x, gen){
  z <- x[abs(x - y) <= gen]
  length(z[!is.na(z)])
}
gen <- 50
range <- seq(from = -3450, to = 1973, by = 99)
nSystem <- cbind(range, unlist(mclapply(range, mc.cores = nCores, 
  function(y) n_calc(y = y, x = persons2020_small$birthyear, gen = gen))))
colnames(nSystem) <- c("birthyear", "nSystem")
nSystem <- as.data.frame(nSystem)
tail(nSystem)
persons2020_small <- persons2020_small %>% left_join(nSystem)

## calculate eligible population size around a target year
n_pop <- function(y, x, gen){
  diff(range(x$cumulative_population[abs(x$year - y) <= gen]))
}
dat_HPI_eligible <- read_csv("dat_HPI_eligible.csv")
nEligible <- cbind(range, unlist(mclapply(range, mc.cores = nCores, 
  function(y) n_pop(y = y, x = dat_HPI_eligible, gen = gen))))
colnames(nEligible) <- c("birthyear", "nEligible")
nEligible <- as.data.frame(nEligible)
tail(nEligible)
persons2020_small <- persons2020_small %>% left_join(nEligible)
persons2020_small %>% dplyr::select(name, gender, hpi, birthyear, nEligible, nSystem) %>% 
  as.data.frame()

## extract relevant variables
hpi_extract <- function(y, x, gen){
  z <- x[abs(x$birthyear - y) <= gen, ] %>% dplyr::select(name, gender, birthyear, hpi)
  z[!is.na(z$hpi), ]
}


## consider generations with increasing year ranges 
# gen - (2*gen+1) * 0:ceiling((1973+3500)/(2*gen))
yr <- 25
gen <- yr + 9 * 0:20
generations <- 1973 - c(yr, 2*cumsum(gen)+1+yr)
generations
hpi_generations <- mclapply(seq_along(generations), mc.cores = nCores, 
                            function(j) hpi_extract(y = generations[j], 
                                                    x = persons2020_small, 
                                                    gen = c(yr,gen)[j]))
Epop_generations <- mclapply(seq_along(generations), mc.cores = nCores, 
                             function(j) n_pop(y = generations[j], 
                                               x = dat_HPI_eligible, 
                                               gen = c(yr,gen)[j]))
Spop_generations <- mclapply(seq_along(generations), mc.cores = nCores, 
                             function(j) n_calc(y = generations[j], 
                                                x = persons2020_small$birthyear, 
                                                gen = c(yr,gen)[j]))


## load in software
library(tidyverse)
library(orderstats)
library(Pareto)
library(parallel)
library(doParallel)
library(xtable)

###############################
## parametric specification
###############################

# we first create some functions and demonstrate 
# them with a simple example

## create era adjustment functions: pareto talent/normal outcome
set.seed(13)
n <- 100
Z <- sort(rnorm(n))
Z[n] <- 1.5*Z[n]

## order_pnorm function
# converts normal order stats to their percentiles
order_pnorm <- function(q = 0, mean = 0, sd = 1, k = 1, n = 1e4){
  p <- pnorm(q, mean = mean, sd = sd)
  pbinom(k - 1, prob = p, size = n, lower.tail = FALSE)
}

## order_pnorm_vec function
# this function converts a vector of normal order stats 
# to their percentiles. This vector should be the entire 
# sample sorted in increasing order
order_pnorm_vec <- function(q, mean = 0, sd = 1){
  q <- sort(q) # just in case
  n <- length(q)
  unlist(lapply(1:n, function(j){
    order_pnorm(q[j], k = j, n = n, mean = mean, sd = sd)
  }))
}

# try it out
u <- order_pnorm_vec(Z)


## order_Pareto_vec function 
# this function transform percentiles from order stats (in increasing order)
# to Pareto values corresponding to the general population 
# of a greater than or equal to size
# default alpha is that of the Pareto principle 80-20
order_Pareto_vec <- function(u, t = 1, alpha = 1.16, npop = 1e4){
  n <- length(u)
  if(length(npop) == 1) npop <- rep(npop, n)
  unlist(lapply(1:n, function(j){
    qPareto(qbeta(u[j], j + npop[j]-n, n + 1 - j), t = t, alpha = alpha)
  }))
}

# try it out on the percentiles obtained from order_norm_vec
X <- order_Pareto_vec(u, npop = 1e6)
#X <- order_Pareto_vec(b, npop = 11907616)


## map_Pareto_vals_vec function 
# this function transforms ordered Pareto values corresponding to 
# the general population to percentiles from order stats 
# (in increasing order)
map_Pareto_vals_vec <- function(x, t = 1, alpha = 1.16, npop = 1e4){
  n <- length(x)
  if(length(npop) == 1) npop <- rep(npop, n)  
  unlist(lapply(1:n, function(j){
    pbeta(pPareto(x[j], t = t, alpha = alpha), j + npop[j]-n, n + 1 - j)
  }))
}

# try it out on the percentiles obtained from order_Pareto_vec
u2 <- map_Pareto_vals_vec(X, npop = 1e6)
sqrt(crossprod(u2 - u)) # check


## oder_qnorm function
# take the ordered percentiles and convert them to order statistics 
# from a normal distribution
order_qnorm <- function(u, mean = 0, sd = 1){
  n <- length(u)
  qnorm(qbeta(u, shape1 = 1:n, shape2 = n:1), mean = mean, sd = sd)
}


# try it out on the percentiles obtained from map_Pareto_vals_vec
Z2 <- order_qnorm(u2)
sqrt(crossprod(Z - Z2)) # check
#cbind(Z,Z2)[which.max(abs(Z - Z2)) + c(-7:7), ]
#tail(cbind(Z,Z2), 10)


## runs the method with normal distribution for hpi and 
## Pareto for talent distribution
test <- do.call(rbind, lapply(1:20, function(j){
  foo <- hpi_generations[[j]]
  bar <- order_pnorm_vec(foo$hpi, 
                         mean = mean(foo$hpi), 
                         sd = sd(foo$hpi))
  
  baz <- order_Pareto_vec(bar, npop = Epop_generations[[j]])
  qux <- cbind(which(baz != Inf), baz[baz != Inf])
  baz[baz == Inf] <- approx(x = qux[, 1], 
                            y = qux[, 2], 
                            xout = which(baz == Inf), 
                            yleft = min(qux[, 2]))$y
  foo$hpi_talent <- sort(baz, decreasing = TRUE)
  foo  
})) %>% arrange(desc(hpi_talent))


## check results
test %>% as.data.frame() %>% head(50)
test %>% as.data.frame() %>% head(50) %>% arrange(desc(birthyear))
test %>% as.data.frame() %>% arrange(desc(hpi)) %>% head(50) %>% arrange(desc(birthyear))


## play around
foo <- hpi_generations[[2]]
print(nrow(foo))
bar <- order_pnorm_vec(foo$hpi, 
                       mean = mean(foo$hpi), 
                       sd = sd(foo$hpi))

baz <- order_Pareto_vec(bar, npop = Epop_generations[[2]])
qux <- cbind(which(baz != Inf), baz[baz != Inf])
print(nrow(qux))
baz[baz == Inf] <- approx(x = qux[, 1], 
                          y = qux[, 2], 
                          xout = which(baz == Inf), 
                          yleft = min(qux[, 2]) )$y

#print(length(baz))
#print(which(is.na(baz)))
#print( length(sort(baz, decreasing = TRUE)) )

foo$hpi_talent <- sort(baz, decreasing = TRUE)
foo  





####################################
## nonparametric specification
####################################

# These functions follow from the era-adjustment baseball paper

Ftilde <- function(y, t, ystar){
  y <- sort(y)
  n <- length(y)
  ytilde <- rep(0, n + 1)
  
  ytilde[n+1] <- y[n] + ystar
  ytilde[2:n] <- unlist(lapply(2:n, function(j){
    (y[j]+y[j-1])/2 
  }))
  
  if (t >= ytilde[n+1]) {
    1 - 0.1^7
  } else if (t <= ytilde[1]) {
    0
  } else {
    j <- length(which(ytilde < t))
    (j - 1) / n + (t - ytilde[j]) / (n*(ytilde[j+1] - ytilde[j]))
  }
  
}


## This function computes the y** value in the era-adjusment baseball paper.
## This is an upper bound so that F(Y_(n)) <= 1 where the distribution function 
##   estimate is similar to the standard empirical cdf estimator
## A large y** value means that F(Y_(n)) is relatively far away from 1. This 
##   corresponds to a scenario in which Y_(n) does not stand out from the sample.
## A small y** value means that F(Y_(n)) is relatively close to 1. This 
##   corresponds to a scenario in which Y_(n) stands out from the sample.

#component = hpi_generations[[15]] %>% pull(hpi)
#stab = 0.0001

bat_thres <- function(component, stab = 0.0001){
  ## stab means we add a small number to the largest value to avoid numerical stability problems. 
  ## cutoff means we select a certain percentage of systems that include maximal possible components in the tail
  ## rather than include k possible components in the tail based on the adjusted R adjusted square. 
  ## If the distance between the largest value and second largest value in the right tail is relatively large,
  ## adjusted R square may not be a good quantity to represent how well the fit is. 
  #component = hpi_generations[[17]] %>% pull(hpi)

  # obtain initial quantities for linear approximation
  Y <- sort(as.matrix(component))
  n <- length(Y)
  Y[n] <- Y[n] + stab # for stability
  pi <- 1 - (n:1 - 1/3)/(n + 1/3)
  W <- log(pi/(1-pi))
  K1 = max(5, floor(1.3*sqrt(n))); K2 = 2*floor(log10(n)*sqrt(n))
  k <- 5
  
  # use arguments from Scholz section 3 for estimating k
  #
  # this argument is based on model fit and not longest stretch of 
  # contiguous I0
  ind <- NULL
  try({
    k_selector <- do.call(rbind, lapply(5:min(K2,500), function(k){
      
      Ytil <- Y - median(Y)
      Ztil <- tail(Ytil, k)
      M1k <- 1/(k-1) * sum( log(Ztil[2:k]/Ztil[1]) )
      M2k <- 1/(k-1) * sum( log(Ztil[2:k]/Ztil[1])^2 )
      ck <- M1k + 1 - 0.5*(1 - M1k^2/M2k)^{-1}
      fck <- ((-n*log(pi))^{-ck} - 1)/ck
      
      Sigma <- matrix(0, k, k)
      for(i in 1:k){
        for(j in 1:i){
          Sigma[i,j] <- i^{-ck-1} * j^{-ck}
        } 
      }
      for(j in 1:k){
        for(i in 1:(j-1)){
          Sigma[i,j] <- j^{-ck-1} * i^{-ck}
        } 
      }
      
      rotate <- function(x) t(apply(x, 2, rev))
      Sigma <- rotate(rotate(Sigma))
      #eigen(Sigma)
      #Sigma.inv <-  solve(Sigma)
      eig <- eigen(Sigma)
      C <- eig$vec %*% diag(1/sqrt(eig$val)) %*% t(eig$vec)
      Zk <- C %*% tail(Y, k)
      Xk <- cbind(1, tail(fck, k))
      Wk <-  C %*% Xk
      # try linear and quadratic model
      m1 <- lm(tail(Y, k) ~ tail(fck, k))
      m2 <- lm(tail(Y, k) ~ tail(fck, k) + I(tail(fck, k)^2))
      m3 <- lm(Zk ~ -1 + Wk)
      delta.sq <- summary(m3)$sigma^2
      Tk <- coef(m3)[2] / summary(m3)$sigma
      
      kappa.sq <- solve(crossprod(Wk))[2,2]
      kappa <- sqrt(kappa.sq)
      I0 <- c(kappa * qt(0.25, df = k - 2, ncp = 1/kappa),
              kappa * qt(0.75, df = k - 2, ncp = 1/kappa))
      I1 <- c(kappa * qt(0.05, df = k - 2, ncp = 1/kappa), 
              kappa * qt(0.95, df = k - 2, ncp = 1/kappa))
      I0int <- ifelse(I0[1] <= Tk && Tk <= I0[2], 1, 0)
      I1int <- ifelse(I1[1] <= Tk && Tk <= I1[2], 1, 0)
      c(k, Tk, I0int, I1int, summary(m1)$adj.r.squared, 
        summary(m2)$adj.r.squared)
      
    }))
    
    #k <- k_selector[max(which(k_selector[, 3] == 1)), 1]
    #k <- k_selector[which.max(k_selector[, 5]), 1]
    k_selector <- as.data.frame(k_selector)
    colnames(k_selector) <- c("k", "Tk", "I0", "I1", "R.sq", "Rquad.sq")
    k_selector_I0 <- k_selector %>% filter(I0 == 1)
    a <- which.max(k_selector_I0$R.sq)
    b <- which.max(k_selector_I0$Rquad.sq)
    ind <- which.max(c(k_selector_I0[a, ]$R.sq, 
                       k_selector_I0[b, ]$Rquad.sq))
    k <- k_selector_I0[c(a,b)[ind] , 1]
    #if(diff(Y)[n-1] > cutoff){ 
    #  k <- max(k_selector_I0$k)
    #  if(k < 0) k <- K2
    #}
    
  }, silent = TRUE)
  
  if(length(k) == 0) k <- round(mean(K1,K2))
  if(is.na(k)) k <- round(mean(K1,K2))
  if(k == 0) k <- round(mean(K1,K2))
  if(k >= n) k <- K2
  
  
  # find probability value using linear tail behavior
  Z <- tail(Y, k)
  m1 <- lm(tail(Y, k) ~ tail(pi, k))
  beta <- m1$coefficients
  ystar <- ub <- 0
  f <- function(x) beta[1] + beta[2] * x - max(Y)
  try({
    foo <- uniroot(f, c(0.0001, 5), tol = 1e-10)
    ub <- foo$root        
  })
  #plot(tail(pi, k), tail(Y, k))
  #lines(tail(pi, k), predict(m1), col = "red")
  if(ub >= 1){
    m1 <- lm(tail(Y, k) ~ tail(pi, k) + tail(pi^2, k))
    beta <- m1$coefficients
    ystar <- ub <- 0
    f <- function(x){
      beta[1] + beta[2] * x + beta[3] * x - max(Y)
    } 
    try({
      foo <- uniroot(f, c(0.0001, 5), tol = 1e-10)
      ub <- foo$root        
    })
  }
  
  if(ub >= 1){
    m1 <- loess(tail(Y, k) ~ tail(pi, k), degree = 2)
    f <- function(w){
      max(Y) - predict(m1, newdata =  w)
    }
    try({
      foo <- uniroot(f, c(tail(pi, k)[1], max(pi)), tol = 1e-10)
      ub <- foo$root        
    }, silent = TRUE)
    #plot(tail(pi, k), tail(Y, k))
    #lines(tail(pi, k), predict(m1))    
  }

  # find probability value using logistic tail behavior
  if(ub >= 1){
    m1 <- lm(tail(Y,k) ~ tail(W, k))
    beta <- m1$coefficients
    f <- function(x) beta[1] + beta[2] * log(x/(1-x)) - max(Y)
    try({
      foo <- uniroot(f, c(0.000001, 0.999999), tol = 1e-10)
      ub <- foo$root        
    })  
  }
  
  # if possible, find ystar 
  if(ub >= Ftilde(y = Y, t = max(Y), ystar = 10)){
    try({
      g <- function(ystar) ub - Ftilde(y = Y, t = max(Y), ystar = ystar)
      bar <- uniroot(g, c(0, 10), tol = 1e-10)
      ystar <- bar$root
    })
  }
  
  # try higher order polynomial fits on the logistic scale
  flag1 <- flag2 <- flag3 <- NULL
  if(ub < Ftilde(y = Y, t = max(Y), ystar = 10)){
    m1 <- lm(tail(Y,k) ~ tail(W, k) + I(tail(W, k)^2))
    beta <- m1$coefficients
    f <- function(x){
      beta[1] + beta[2] * log(x/(1-x)) + 
        beta[3] * log(x/(1-x))^2 - max(Y)
    }
    x <- seq(1e-6, 1-1e-6, length = 3e4)
    y <- sapply(x, f)
    flag1 <- try({
      x_upper <- x[max(which(y > 0))]
      y_lower_vec <- y[1:max(which(y > 0))]
      x_lower <- x[max(which(y_lower_vec < 0))]
    }, silent = TRUE)
    flag2 <- try({
      foo <- uniroot(f, c(x_lower, x_upper), tol = 1e-10)
      ub <- foo$root        
    }, silent = TRUE)  
    flag3 <- try({
      g <- function(ystar) ub - Ftilde(y = Y, t = max(Y), ystar = ystar)
      bar <- uniroot(g, c(0, 10), tol = 1e-10)
      ystar <- bar$root
    }, silent = TRUE)
  }
  
  if(class(flag1) == "try-error" | class(flag2) == "try-error" | class(flag3) == "try-error" ){
    if(ub < Ftilde(y = Y, t = max(Y), ystar = 10)){
      m1 <- lm(tail(Y,k) ~ tail(W, k) + I(tail(W, k)^2) + I(tail(W, k)^3))
      beta <- m1$coefficients
      f <- function(x){
        beta[1] + beta[2] * log(x/(1-x)) + 
          beta[3] * log(x/(1-x))^2 + 
          beta[4] * log(x/(1-x))^3 - max(Y)
      }
      x <- seq(1e-6, 1-1e-6, length = 3e4)
      y <- sapply(x, f)
      flag1 <- try({
        x_upper <- x[max(which(y > 0))]
        y_lower_vec <- y[1:max(which(y > 0))]
        x_lower <- x[max(which(y_lower_vec < 0))]
      }, silent = TRUE)  
      flag2 <- try({
        foo <- uniroot(f, c(x_lower, x_upper), tol = 1e-10)
        ub <- foo$root        
      }, silent = TRUE)  
      flag3 <- try({
        g <- function(ystar) ub - Ftilde(y = Y, t = max(Y), ystar = ystar)
        bar <- uniroot(g, c(0, 10), tol = 1e-10)
        ystar <- bar$root
      }, silent = TRUE)
    }
  }
  
  ## try loess fit on logistic scale
  if(class(flag1) == "try-error" | class(flag2) == "try-error" | class(flag3) == "try-error" ){
    m2 <- loess(tail(Y, k) ~ tail(W, k), degree = 2)
    f <- function(w){
      max(Y) - predict(m2, newdata =  w)
    }
    flag1 <- try({
      foo <- uniroot(f, c(tail(W, k)[1], max(W)), tol = 1e-10)
      ub_w <- foo$root        
      ub <- 1/(1+exp(-ub_w))
    }, silent = TRUE)  
    flag2 <- try({
      g <- function(ystar) ub - Ftilde(y = Y, t = max(Y), ystar = ystar)
      bar <- uniroot(g, c(0, 10), tol = 1e-10)
      ystar <- bar$root
    }, silent = TRUE)
  }
  #plot(tail(W, k), tail(Y, k))
  #lines(tail(W, k), predict(m2), col = "red")

  
  ## try the loess method on fewer extreme points and iterate the number 
  ## of extreme points until a ystar value is obtained
  ystar2 <- ystar3 <- 0 #ystar
  #if(ystar == 0){
    
    ## K1 down 
    k1 <- K1    
    while(ystar2 == 0 & k1 >= 5){
      m1 <- loess(tail(Y, k1) ~ tail(W, k1), degree = 2)
      f <- function(w){
        max(Y) - predict(m1, newdata =  w)
      }
      try({
        foo <- uniroot(f, c(tail(W, k1)[1], max(W)), tol = 1e-10)
        ub_w <- foo$root        
        ub <- 1/(1+exp(-ub_w))
      }, silent = TRUE)
      #plot(tail(W, k), tail(Y, k))
      #lines(tail(W, k), predict(m1))
      try({
        g <- function(ystar) ub - Ftilde(y = Y, t = max(Y), ystar = ystar)
        bar <- uniroot(g, c(0, 10), tol = 1e-10)
        ystar2 <- bar$root
      }, silent = TRUE)
      
      if(ystar2 == 0){
        m1 <- lm(tail(Y,k1) ~ tail(W, k1) + I(tail(W, k1)^2))# + I(tail(W, k)^3))
        beta <- m1$coefficients
        #plot(tail(W, k1+1), tail(Y, k1+1))
        #lines(tail(W, k1+1), predict(m1))
        f <- function(x){
          beta[1] + beta[2] * log(x/(1-x)) + 
            beta[3] * log(x/(1-x))^2 -
            #beta[4] * log(x/(1-x))^3 
            max(Y)
        }
        x <- seq(1e-6, 1-1e-6, length = 3e4)
        y <- sapply(x, f)
        try({
          x_upper <- x[max(which(y > 0))]
          y_lower_vec <- y[1:max(which(y > 0))]
          x_lower <- x[max(which(y_lower_vec < 0))]
        }, silent = TRUE)  
        try({
          foo <- uniroot(f, c(x_lower, x_upper), tol = 1e-10)
          ub <- foo$root        
        }, silent = TRUE)  
        try({
          g <- function(ystar) ub - Ftilde(y = Y, t = max(Y), ystar = ystar)
          bar <- uniroot(g, c(0, 10), tol = 1e-10)
          ystar2 <- bar$root
        }, silent = TRUE)
      }
      k1 <- k1 - 1
    }    
    
    ## try cubic method
    #k1 <- k
    k1 <- K1
    while(ystar2 == 0 & k1 >= 5){
      
      m1 <- lm(tail(Y,k1) ~ tail(W, k1) + I(tail(W, k1)^2) + I(tail(W, k1)^3))
      beta <- m1$coefficients
      #plot(tail(W, k1+1), tail(Y, k1+1))
      #lines(tail(W, k1+1), predict(m1))
      f <- function(x){
        beta[1] + beta[2] * log(x/(1-x)) + 
          beta[3] * log(x/(1-x))^2 +
          beta[4] * log(x/(1-x))^3 - 
          max(Y)
      }
      x <- seq(1e-6, 1-1e-6, length = 3e4)
      y <- sapply(x, f)
      try({
        x_upper <- x[max(which(y > 0))]
        y_lower_vec <- y[1:max(which(y > 0))]
        x_lower <- x[max(which(y_lower_vec < 0))]
      }, silent = TRUE)  
      try({
        foo <- uniroot(f, c(x_lower, x_upper), tol = 1e-10)
        ub <- foo$root        
      }, silent = TRUE)  
      try({
        g <- function(ystar) ub - Ftilde(y = Y, t = max(Y), ystar = ystar)
        bar <- uniroot(g, c(0, 10), tol = 1e-10)
        ystar2 <- bar$root
      }, silent = TRUE)
      k1 <- k1 - 1
    }
    
    
    ## k 5 up
    k1 <- 5    
    while(ystar3 == 0 & k1 <= K2){
      m1 <- loess(tail(Y, k1) ~ tail(W, k1), degree = 2)
      f <- function(w){
        max(Y) - predict(m1, newdata =  w)
      }
      try({
        foo <- uniroot(f, c(tail(W, k1)[1], max(W)), tol = 1e-10)
        ub_w <- foo$root        
        ub <- 1/(1+exp(-ub_w))
      }, silent = TRUE)
      #plot(tail(W, k), tail(Y, k))
      #lines(tail(W, k), predict(m1))
      try({
        g <- function(ystar) ub - Ftilde(y = Y, t = max(Y), ystar = ystar)
        bar <- uniroot(g, c(0, 10), tol = 1e-10)
        ystar3 <- bar$root
      }, silent = TRUE)
      
      if(ystar3 == 0){
        m1 <- lm(tail(Y,k1) ~ tail(W, k1) + I(tail(W, k1)^2))# + I(tail(W, k)^3))
        beta <- m1$coefficients
        #plot(tail(W, k), tail(Y, k))
        #lines(tail(W, k), predict(m1))
        f <- function(x){
          beta[1] + beta[2] * log(x/(1-x)) + 
            beta[3] * log(x/(1-x))^2 -
            #beta[4] * log(x/(1-x))^3 
            max(Y)
        }
        x <- seq(1e-6, 1-1e-6, length = 3e4)
        y <- sapply(x, f)
        try({
          x_upper <- x[max(which(y > 0))]
          y_lower_vec <- y[1:max(which(y > 0))]
          x_lower <- x[max(which(y_lower_vec < 0))]
        }, silent = TRUE)  
        try({
          foo <- uniroot(f, c(x_lower, x_upper), tol = 1e-10)
          ub <- foo$root        
        }, silent = TRUE)  
        try({
          g <- function(ystar) ub - Ftilde(y = Y, t = max(Y), ystar = ystar)
          bar <- uniroot(g, c(0, 10), tol = 1e-10)
          ystar3 <- bar$root
        }, silent = TRUE)
      }
      k1 <- k1 + 1
    }    
    
    
    ## try cubic method
    #k1 <- k
    k1 <- 5
    while(ystar3 == 0 & k1 <= K2){
      
      m1 <- lm(tail(Y,k1) ~ tail(W, k1) + I(tail(W, k1)^2) + I(tail(W, k1)^3))
      beta <- m1$coefficients
      #plot(tail(W, k), tail(Y, k))
      #lines(tail(W, k), predict(m1))
      f <- function(x){
        beta[1] + beta[2] * log(x/(1-x)) + 
          beta[3] * log(x/(1-x))^2 +
          beta[4] * log(x/(1-x))^3 - 
          max(Y)
      }
      x <- seq(1e-6, 1-1e-6, length = 3e4)
      y <- sapply(x, f)
      try({
        x_upper <- x[max(which(y > 0))]
        y_lower_vec <- y[1:max(which(y > 0))]
        x_lower <- x[max(which(y_lower_vec < 0))]
      }, silent = TRUE)  
      try({
        foo <- uniroot(f, c(x_lower, x_upper), tol = 1e-10)
        ub <- foo$root        
      }, silent = TRUE)  
      try({
        g <- function(ystar) ub - Ftilde(y = Y, t = max(Y), ystar = ystar)
        bar <- uniroot(g, c(0, 10), tol = 1e-10)
        ystar3 <- bar$root
      }, silent = TRUE)
      k1 <- k1 + 1
    }
    
    
  #}
  #plot(tail(W, k+1), tail(Y, k+1))
  #lines(tail(W, k+1), predict(m2))
  
  c(ystar, ystar2, ystar3)
  
}
 






## compute y**
ystar_mat <- do.call(rbind, mclapply(1:length(hpi_generations), function(j){
  bat_thres(component = hpi_generations[[j]] %>% pull(hpi))
}, mc.cores = nCores))
ystar_mat
ystar_mat2 <- apply(ystar_mat, 1, function(x) mean(x[x > 0]))
#ystar_mat2[17] <- 5

## This function estimates underlying talent values
talent_computing_nonpara <- function(ystar, data, npop, alpha = 1.16){
  y <- data$hpi
  
  Ftilde <- function(y, t, ystar){
    y <- sort(y)
    n <- length(y)
    ytilde <- rep(0, n + 1)
    
    ytilde[n+1] <- y[n] + ystar
    ytilde[2:n] <- unlist(lapply(2:n, function(j){
      (y[j]+y[j-1])/2 
    }))
    
    if (t >= ytilde[n+1]) {
      1 - 0.1^7
    } else if (t <= ytilde[1]) {
      0
    } else {
      j <- length(which(ytilde < t))
      (j - 1) / n + (t - ytilde[j]) / (n*(ytilde[j+1] - ytilde[j]))
    }
    
  }
  
  Aptitude_nonpara <- function(p, alpha = 1.16, npop){
    
    # converts order stats to their percentiles
    order_pbino <- function(p = 0, k = 1, n = 1e4){
      pbinom(k - 1, prob = p, size = n, lower.tail = FALSE)
    }
    
    # converts a vector of order stats 
    # to their percentiles. This vector should be the entire 
    # sample sorted in increasing order
    p <- sort(p) # just in case
    n <- length(p)
    u = unlist(lapply(1:n, function(j){
      order_pbino(p[j], k = j, n = n)
    }))
    
    # transforms percentiles from order stats (in increasing order)
    # to Pareto values corresponding to the general population 
    # of a greater than or equal to size
    # default alpha is that of the Pareto principle 80-20
    n <- length(u)
    #if(length(npop) == 1) npop <- rep(npop, n)
    unlist(lapply(1:n, function(j){
      qPareto(qbeta(u[j], j + npop -n , n + 1 - j), t = 1, alpha = alpha)
    }))
  }
  
  ## hpi talent
  hpi_talent <- Aptitude_nonpara(p = unlist(lapply(y, function(xx) 
    Ftilde(y = y, t = xx, ystar = ystar))), npop = npop)
  
  #max_hpi_talent <- max(hpi_talent) - 1
  
  ## using the distribution from full time players
  #bar <- rbind(bar, do.call(rbind, lapply(range, function(j){
  #  rbind(bar %>% dplyr::select(-WAR_talent), foo[j, ]) %>% arrange(comp) %>%
  #    mutate(WAR_talent = Aptitude_nonpara(p = unlist(lapply(comp, function(xx) 
  #      Ftilde(y = full_comp, t = xx, ystar = ystar, component_name = component_name))), npop = pops)) %>%
  #    filter(PA < thres) %>% 
  #    mutate(WAR_talent = ifelse(WAR_talent > max_WAR_talent+1, max_WAR_talent, WAR_talent))
  #})))
  data %>% mutate(hpi_talent = sort(hpi_talent, decreasing = TRUE))
}


## now compute the talent values for all people in the Pantheon list
# ystar_mat <- apply(ystar_mat, 2, as.numeric)
# ystar_mat[is.na(ystar_mat), ] <- 6
# era <- 2
# talent_computing_nonpara(ystar = as.numeric(ystar_mat[[era]]), 
#                          data = hpi_generations[[era]], 
#                          npop = Epop_generations[[era]])
### length(hpi_generations)
hpi_adjust3 <- do.call(rbind, mclapply(1:length(hpi_generations), 
                                      mc.cores = nCores, function(j){
                                        talent_computing_nonpara(ystar = as.numeric(ystar_mat2[j]), 
                                                                 data = hpi_generations[[j]], 
                                                                 npop = Epop_generations[[j]])
                                      })) %>% arrange(desc(hpi_talent))

hpi_adjust3 %>% as.data.frame() %>% head(50)
## write_csv(hpi_adjust, "~/research/pantheon/hpi_adjust_22eras.csv")
## write_csv(hpi_adjust, "~/research/pantheon/hpi_adjust_22eras_kdown.csv")
## write_csv(hpi_adjust2, "~/research/pantheon/hpi_adjust_22eras_k6up.csv")
## write_csv(hpi_adjust3, "~/research/pantheon/hpi_adjust_22eras_mean.csv")
## write_csv(ystar_mat, "~/research/pantheon/ystar_mat_22eras.csv")
hpi_adjust %>% pull(gender) %>% table()
View(hpi_adjust3)

dat <- read_csv("~/research/pantheon/hpi_adjust_22eras.csv")
View(dat)


## Construct quantities based on eras defined in the link below: 
## https://phi.history.ucla.edu/nchs/preface/developing-standards/

# Era 1: The Beginnings of Human Society
# Era 2: Early Civilizations and the Emergence of Pastoral Peoples, 4000-1000 BCE
# Era 3: Classical Traditions, Major Religions, and Giant Empires, 1000 BCE-300 CE
# Era 4: Expanding Zones of Exchange and Encounter, 300-1000 CE
# Era 5: Intensified Hemispheric Interactions, 1000-1500 CE
# Era 6: Emergence of the First Global Age, 1450-1770
# Era 7: An Age of Revolutions, 1750-1914
# Era 8: A Half-Century of Crisis and Achievement, 1900-1945
# Era 9: The 20th Century Since 1945: Promises and Paradoxes

era_assigner <- function(yr, lag = 30){
  yr <- yr + lag
  era <- 2 + ifelse(yr > -1000,1,0) + 
    ifelse(yr > 300,1,0) + 
    ifelse(yr > 1000,1,0) + 
    ifelse(yr > 1500,1,0) + 
    ifelse(yr > 1770,1,0) + 
    ifelse(yr > 1900,1,0) + 
    ifelse(yr > 1945,1,0)
  era
}

persons2020_small_era <- persons2020_small %>% 
  mutate(era = era_assigner(birthyear))
persons2020_small_era %>% group_by(era) %>% summarise(N = n())

#persons2020_small
hpi_generations <- lapply(split(persons2020_small_era, f = as.factor(persons2020_small_era$era)), 
                          function(xx) xx %>% dplyr::select(name, gender, birthyear, hpi))


Epop_generations <- dat_HPI_eligible %>% 
  filter(year %in% c(-3500, -1000, 300, 1000, 1500, 1770, 1900, 1945, 1973)) %>% 
  pull(cumulative_population) %>% diff()


len <- length(hpi_generations)
foo <- do.call(rbind, lapply(1:len, function(group){
  cbind(hpi_generations[[group]], Epop_generations[[group]], group)
}))
head(foo)
colnames(foo)[5] <- "npop"


ystar_mat <- do.call(rbind, mclapply(1:length(hpi_generations), function(j){
  bat_thres(component = hpi_generations[[j]] %>% pull(hpi))
}, mc.cores = nCores))
ystar_mat
ystar_mat2 <- apply(ystar_mat, 1, function(x) mean(x[x > 0]))


# talent_computing_nonpara(ystar = as.numeric(ystar_mat[[era]]), 
#                          data = hpi_generations[[era]], 
#                          npop = Epop_generations[[era]])
hpi_adjust_hist <- do.call(rbind, mclapply(1:length(hpi_generations), 
                                      mc.cores = nCores, function(j){
                                        talent_computing_nonpara(ystar = as.numeric(ystar_mat[[j]]), 
                                                                 data = hpi_generations[[j]], 
                                                                 npop = Epop_generations[[j]])
                                      })) %>% arrange(desc(hpi_talent))

hpi_adjust_hist %>% as.data.frame() %>% head(50)
## write_csv(hpi_adjust_hist, "~/research/pantheon/hpi_adjust_hist_eras.csv")
## write_csv(hpi_adjust_hist, "~/research/pantheon/hpi_adjust_hist_eras_mean.csv")
## write_csv(ystar_mat, "~/research/pantheon/ystar_mat_hist_eras.csv")


## investigate agreement
length(which((hpi_adjust_hist %>% as.data.frame() %>% head(50) %>% pull(name)) %in% (persons2020_small %>% head(50) %>% pull(name))))
length(which((hpi_adjust3 %>% as.data.frame() %>% head(50) %>% pull(name)) %in% (persons2020_small %>% head(50) %>% pull(name))))

length(which((hpi_adjust_hist %>% as.data.frame() %>% head(25) %>% pull(name)) %in% (persons2020_small %>% head(25) %>% pull(name))))
length(which((hpi_adjust3 %>% as.data.frame() %>% head(25) %>% pull(name)) %in% (persons2020_small %>% head(25) %>% pull(name))))

length(which((hpi_adjust_hist %>% as.data.frame() %>% head(10) %>% pull(name)) %in% (persons2020_small %>% head(10) %>% pull(name))))
length(which((hpi_adjust3 %>% as.data.frame() %>% head(10) %>% pull(name)) %in% (persons2020_small %>% head(10) %>% pull(name))))

