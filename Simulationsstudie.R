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
  new_magnitudes <- JJ[idxJ[1:k]]
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

# (3) functions for estiamting the parameters
fct_empfracmom <- function( q , data ){       # emp. frac. moment
  r <- 1/(length(data))*sum(data^q)
  return(r)       }
fct_theta <- function( beta , q , data ){     # estimator for theta depending on beta
  r <- 2*(fct_empfracmom(q,data)^2)/fct_empfracmom(2*q , data)*beta/(q*pi)*
    (gamma(1-q)^2)/gamma(1-2*q)*(sin(q*pi/beta)^2)/sin(2*q*pi/beta)
  return(r)               }

fct_beta <- function(q1 , q2 , data){         # estimator for beta
  fct_root <- function( beta , q1 , q2 , data){
        r <- fct_theta( beta , q1 , data) - fct_theta( beta , q2 , data)
        return(r)                  }
  r <- uniroot(fct_root , c(2*max(q1,q2)+0.001,1) , q1 = q1 , 
               q2 = q2 , data = data)
  return(r$root)      }   



# parameter choice
n <- 1000000    # sample size
k <- 500      # number of exceedances
tail <- 0.8   # tail parameter beta
EI <- 0.6     # extremal index theta
q1 <- 0.2
q2 <- 0.3

daten <- cbind(fct_MARdata(n = n , EI = EI) , fct_Paretodata(n = n , tail = tail))

daten_k <- fct_thin(data = daten , k = 500)

beta_dach <- fct_beta(q1 = q1 , q2 = q2 , data = daten_k)


