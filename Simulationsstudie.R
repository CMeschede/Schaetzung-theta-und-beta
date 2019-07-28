## Simulation Study ####

# packages

#install.packages("extRemes")       for Frechet distribution
#install.packages("MittagLeffleR")
#install.packages("CTRE")
#install.packages("ReIns")          for Pareto distribution

# functions:

# (1) function to get the k highest events and the associated interarrival times
fct_thin <- function(data , k){
  JJ <- data[,1]
  WW <- data[,2]
  n <- length(JJ)
  if (k > n) 
    stop("Can't threshold to ", k, " observations if I only have ", 
         n)
  idxJ <- sort(order(JJ, decreasing = TRUE)[1:k])
  b <- rep(0 , times = n+1)
  b[idxJ+1] <- 1
  a <- (1+cumsum(b==1))[-length(b)]
  newJJ <- JJ[idxJ]
  newWW <- aggregate(WW , list(a) , sum)$x[1:k]
  out <- cbind(newJJ,newWW)
  return(out)
}

# (2) Data generating:
# MAR-Process (events)
fct_MARdata <- function(n ,EI){
  X <- rep(0, n)
  eps <- extRemes::revd(n , scale = 1 , shape = 1 , loc = 1)
  X[1] <- eps[1]/EI
  for (i in 2:n){
    X[i] <- max((1 - EI) * X[i - 1], eps[i])
  }
  return(X)
}
# Pareto distributed waiting times
fct_Paretodata <- function( n , tail){
    sigma <- gamma(1-tail)^{-1/tail}
    X <- ReIns::rpareto(n , shape = tail , scale = sigma)
    return(X)
}

# (3) functions for estimating the parameters
fct_empfracmom <- function( q , data ){       # emp. frac. moment
  r <- 1/(length(data))*sum(data^(2*q))
  s <- (1/(length(data))*sum(data^q))^2
  return(s/r)       }
fct_theta <- function( beta , q , data ){     # estimator for theta depending on beta
  r <- fct_empfracmom(q,data)*2*beta/(q*pi)*
    (gamma(1-q)^2)/gamma(1-2*q)*(sin(q*pi/beta)^2)/sin(2*q*pi/beta)
  return(r)}
fct_beta <- function(q1 , q2 , data){         # estimator for beta
  fct_root <- function( beta , q1 , q2 , data){
        r <- fct_theta( beta , q1 , data) - fct_theta( beta , q2 , data)
        return(r)                  }
  r <- uniroot(fct_root , c(2*max(q1,q2)+0.001,1) , q1 = q1 , 
               q2 = q2 , data = data)
  return(r$root)      }

# alternativen:
fct_empfracmom2 <- function( q , data ){       # emp. frac. moment
  r <- 1/(length(data))*sum(data^q*(pmax(data-1,0)^q))
  s <- (1/(length(data))*sum(data^q))^2
  return(s/r)       }
fct_theta2 <- function( beta , q , data ){     # estimator for theta depending on beta
  r <- fct_empfracmom2(q,data)*2*beta/(q*pi)*
    (gamma(1-q)^2)/gamma(1-2*q)*(sin(q*pi/beta)^2)/sin(2*q*pi/beta)
  return(r)}
fct_beta2 <- function(q1 , q2 , data){         # estimator for beta
  fct_root <- function( beta , q1 , q2 , data){
    r <- fct_theta2( beta , q1 , data) - fct_theta2( beta , q2 , data)
    return(r)                  }
  r <- uniroot(fct_root , c(2*max(q1,q2)+0.001,1) , q1 = q1 , 
               q2 = q2 , data = data)
  return(r$root)      }
   



# parameter choice
n <- 10000  # sample size
k <- 500      # number of exceedances
tail <- 0.8   # tail parameter beta
EI <- 0.6     # extremal index theta
q1 <- 0.01
q2 <- 0.02

daten <- cbind(fct_MARdata(n = n , EI = EI) , fct_Paretodata(n = n , tail = tail))

daten_k <- fct_thin(data = daten , k = 500)
daten_k[,2] <- floor(daten_k[,2])
head(table(daten_k[,2]))

beta_dach <- fct_beta(q1 = q1 , q2 = q2 , data = daten_k[,2]); beta_dach

fct_theta(beta_dach , q1 , data = daten_k[,2])
fct_theta(beta_dach , q2 , data = daten_k[,2])
fct_theta(tail, 0.01, data = daten_k[,2])
# works so far

beta_dach2 <- fct_beta2(q1 = q1 , q2 = q2 , data = daten_k[,2]); beta_dach2

fct_theta2(beta_dach2 , q1 , data = daten_k[,2])
fct_theta2(beta_dach2 , q2 , data = daten_k[,2])
fct_theta2(tail, 0.01, data = daten_k[,2])
# not as good as first try




# test
u <- sort(daten[, 1], decreasing = TRUE)[k]
b <- (1-extRemes::pevd(EI * u , scale = 1 , shape = 1 , loc = 1))^(-1/tail)
rmisch <- function(n, theta, beta, b) {
  p <- runif(n)
  ml <- MittagLeffleR::rml(n, beta, b * theta^(-1/beta))
  ifelse(p < 1-theta, 0, ml)
}
plot(ecdf(daten_k[,2]))
lines(ecdf(rmisch(100000, .6, .8, b)), col = "red")

fct_beta(q1, q2, rmisch(500, .6, .8, b))
fct_theta(.8 , q1 , data = rmisch(500, .6, .8, b))

summary(daten_k)


