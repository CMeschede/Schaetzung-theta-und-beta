# Example:  Solar Flare Data
# aus dem CTRE-Paket

## Daten:
# install.packages("CTRE")
daten1 <- CTRE::flares
daten2 <- CTRE::ctre(daten1)

magnitudes <- CTRE::magnitudes(daten2)
interarrival <- CTRE::interarrival(daten2)

## functions: ####
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
  r <- uniroot(fct_root , c(max(q1,q2)+0.001,1) , q1 = q1 , 
               q2 = q2 , data = data)
  return(r$root)      }   

k <- length(CTRE::interarrival(daten2)):10
q1 <- 0.01
q2 <- 0.05
est <- sapply( k , function(x){
  daten <- CTRE::interarrival(CTRE::thin(daten2, k = x))
  beta_hat1 <- fct_beta(q1 , q2 , data = daten)
  theta_hat11 <- fct_theta(beta = beta_hat1 , q = q1 , data = daten)
  theta_hat12 <- fct_theta(beta = beta_hat1 , q = q2 , data = daten)
  return(c(beta_hat1,theta_hat11,theta_hat12))
} , simplify = "array")
