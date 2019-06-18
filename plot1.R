## first plot

# true parameters 
beta <- c(0.3,0.6,0.9)
theta <- c(0.3,0.6,0.9)
gitter <- expand.grid(beta,theta)
b <- seq(0.001,1,0.001)
q <- c(0.01 , 0.1 , 0.4 , 0.8)

# functions
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

# generating data
n <- 10000

daten <- array(0 , dim = c(3,3,n))
for(i in 1:3){
  for(j in 1:3){
    prob <- sample( c(0,1) , size = n , replace = T , prob = c(1-theta[j] , theta[j]))
    daten[i,j,] <- sapply(prob , 
                          function(x){ ifelse( x == 0 , return(0) , 
                                               return( MittagLeffleR::rml(1 , tail = beta[i] , 
                                                                    scale = theta[j]^(-1/beta[i])) 
                                         ) ) } )
  }
}

theta_hat <- array(0 , dim = c(length(beta) , length(theta) , length(q) , length(b)))

for(i in 1:length(beta)){
  for(j in 1:length(theta)){
    theta_hat[i,j,,] <-sapply(b , function(x){
      sapply(q , function(y){
        fct_theta(beta = x , q = y , data = daten[i,j,])
      } , simplify = "array")
    } , simplify = "array")
  }
}

par( mfrow = c(3,3) , mar = c(4.1 , 4.1 , 1.5 , 1) ,  las = 1 , oma = c(2, 1, 2, 1) )
for(j in 1:length(theta)){
  for(i in 1:length(beta)){
    plot(b[b>q[1]] , theta_hat[i,j,1,][b>q[1]] , xlim = c(0,1) , ylim = c(0,1) , 
         type = "l" , col = "blue" , xlab = "beta" , ylab = "theta_hat(beta)")
    lines(b[b>q[2]] , theta_hat[i,j,2,][b>q[2]] , col ="violet")
    lines(b[b>q[3]] , theta_hat[i,j,3,][b>q[3]] , col = "darkorange")
    points(b[b>q[4]] , theta_hat[i,j,4,][b>q[4]] , col = "darkgreen")
    points(beta[i] , theta[j] , pch = 16 , col = "red")
  }
}

par( mfrow = c(1,1))
