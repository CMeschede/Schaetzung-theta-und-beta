# R-Code zu Mischverteilung.pdf
## functions: ####
fct_empfracmom <- function( q , data ){       # emp. frac. moment
  r <- 1/(length(data))*sum(data^q)
  return(r)       }
fct_theta <- function( beta , q , data ){     # estimator for theta
  r <- ( ( beta*gamma(1-q)*sin(pi*q/beta) )/( q*pi )*
           fct_empfracmom(q , data) )^(beta/(beta-q))
  return(r)               }
fct_root <- function( beta , q1 , q2 , data){
  r <- fct_theta( beta , q1 , data) - fct_theta( beta , q2 , data)
  return(r)                  }
fct_beta <- function(q1 , q2 , data){         # estimator for beta
  r <- uniroot(fct_root , c(max(q1,q2)+0.001,1) , q1 = q1 , 
               q2 = q2 , data = data)
  return(r$root)      }   

## plot: ####
set.seed(1234)
n <- 10000
beta <- 0.5
theta <- 0.5
# simulated sample:
prob <- sample( c(0,1) , size = n , replace = T , prob = c(1-theta , theta))
# install.packages("MittagLeffleR") 
# (we need the R-package 'MittagLeffleR' to sample ML-distributed values)
daten <- sapply(prob ,                       # sample
                function(x){ ifelse( x == 0 , return(0) , 
                                     return( MittagLeffleR::rml(1 , tail = beta , 
                                                                scale = theta^(-1/beta)) 
                                     ) ) } )


b <- seq(0.001,1,0.001)
q <- c(0.01,0.1,0.4,0.5)
par(mar = c(4,4,0.5,0.5))
plot(b[b>q[1]] , fct_theta(b[b>q[1]], q=q[1] , data=daten) , xlim = c(0,1) , ylim = c(0,1) , type = "l" , 
     xlab = "beta" , ylab = "estimated theta")
lines(b[b>q[2]] , fct_theta(b[b>q[2]], q=q[2] , data=daten) , col="blue" )
lines(b[b>q[3]] , fct_theta(b[b>q[3]], q=q[3] , data=daten) , col="darkgreen" )
lines(b[b>q[4]] , fct_theta(b[b>q[4]], q=q[4] , data=daten) , col="violet" )
points(beta , theta , col = "red" , pch = 19)
legend("bottomright" , c("q=0.01" , "q=0.1" , "q=0.4" , "q=0.5" , "true parameters") , 
       col = c("black" , "blue" , "darkgreen" , "violet" , "red") , pch = c(16,16,16,16,19))

## parameter estimation ####
set.seed(2345)
## chosen parameters:
# 1) true parameters
beta <- 0.5
theta <- 0.5
# 2) random sample
n <- 10000                                   # sample size
# 3) fractions
q1 <- 0.01
q2 <- 0.05
## sample:
prob <- sample( c(0,1) , size = n , replace = T , prob = c(1-theta , theta))
# we need the R-package 'MittagLeffleR' to sample ML-distributed values:
# install.packages("MittagLeffleR") 
daten <- sapply(prob ,                       # sample
                function(x){ ifelse( x == 0 , return(0) , 
                                     return( MittagLeffleR::rml(1 , tail = beta , 
                                                                scale = theta^(-1/beta)) 
                                     ) ) } )
## functions:
fct_empfracmom <- function( q , data ){       # emp. frac. moment
  r <- 1/(length(data))*sum(data^q)
  return(r)       }
fct_theta <- function( beta , q , data ){     # estimator for theta
  r <- ( ( beta*gamma(1-q)*sin(pi*q/beta) )/( q*pi )*
           fct_empfracmom(q , data) )^(beta/(beta-q))
  return(r)               }
fct_root <- function( beta , q1 , q2 , data){
  r <- fct_theta( beta , q1 , data) - fct_theta( beta , q2 , data)
  return(r)                  }
fct_beta <- function(q1 , q2 , data){         # estimator for beta
  r <- uniroot(fct_root , c(max(q1,q2)+0.001,1) , q1 = q1 , 
               q2 = q2 , data = data)
  return(r$root)      }
## estimating beta and theta
# estimating beta by calculating the root:
beta_hat1 <- fct_beta(q1=q1 , q2=q2 , data = daten); beta_hat1
# estimating theta by using beta_hat1
theta_hat11 <- fct_theta(beta = beta_hat1 , q=q1 , data = daten); theta_hat11
theta_hat12 <- fct_theta(beta = beta_hat1 , q=q2 , data = daten); theta_hat12
# error:
abs(beta-beta_hat1); abs(theta-theta_hat11); abs(theta-theta_hat12)
## -> it works

# estimating theta by calculating the q-th emp. frac. moment with q very small: 
theta_hat2 <- fct_empfracmom(10^{-6} , data=daten); theta_hat2
beta_hat2 <- uniroot(function(x){fct_theta(x , q=q1 , data=daten)-theta_hat2} , 
                     interval = c(q1,1))$root; beta_hat2
abs(beta-beta_hat2); abs(theta-theta_hat2)
## -> works as well; 
## -> unsolved problem: which fraction q is small enough 
##    to estimate theta reliably?