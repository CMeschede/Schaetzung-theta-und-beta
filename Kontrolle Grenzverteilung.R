# Untersuchung, ob die Grenzverteilung von Hees & Fried 2019 passt:

# implemeniterung der mixed distribution:
rmisch <- function(n, theta, beta, b) {
  p <- runif(n)
  ml <- MittagLeffleR::rml(n, beta, b * theta^(-1/beta))
  ifelse(p < 1-theta, 0, ml)
}
pmisch <- function(q , theta , beta , b){
  x <- (1-theta)+theta*MittagLeffleR::pml(q , beta , b*theta^(-1/beta))
  x[q<0] <- 0
  return(x)
}
qmisch <- function(p ,theta , beta , b){
  ifelse(p <= (1-theta) , 0 , MittagLeffleR::qml((p-1+theta)/theta , beta , b*theta^(-1/beta)))
}

# function to get the k highest events and the associated interarrival times
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
# Data generating:
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

# parameter choice
n <- 100000000    # sample size
k <- 500      # number of exceedances
tail <- 0.9   # tail parameter beta
EI <- 0.7   

# data:
daten <- cbind(fct_MARdata(n = n , EI = EI) , fct_Paretodata(n = n , tail = tail))

daten_k <- fct_thin(data = daten , k = k) # thinned data/Exceedances


u <- sort(daten[,1], decreasing = TRUE)[k]; u # threshold
b <- (1-extRemes::pevd(EI * u , scale = 1 , shape = 1 , loc = 1))^(-1/tail); b # 1-F(u)


# plotted empirical cdf and theoretical cdf
plot(ecdf(daten_k[,2]))
lines(seq(0,60000000,5000),pmisch(seq(0,60000000,5000) , EI , tail , b) , col = "blue")
## sieht ganz gut aus

# qq-plot
z1 <- sapply((1:k)/(k+1) , qmisch , theta = EI , beta = tail , b = b)
plot(sort(z1) , sort(daten_k[,2]) , xlab = "Quantile der Mischverteilung" , 
     ylab = "Exzess waiting times u0" #, ylim = c(0,50000) , xlim = c(0,50000)
)
points(c(0,500000000) , c(0,500000000) , type = "l" , col = "red")
## der rechte Rand passt nicht so gut


# Simulation:
## hier ziehe ich N mal k Zufallszahlen der Mischverteilung und schaue mir so die Verteilung
## der empirischen Quantile an. Wir sehen, dass gerade bei den hohen Quantilen eine hohe Variabilitaet
## herrscht und dies den Rand des obrigen QQ-Plot erklärt

N <- 10000 # number of iterations

sim <- replicate(N , rmisch(k , EI , tail , b) , simplify = "array")
sim_sort <- apply(sim , 2 , sort , decreasing  = T)
sim_sort_sort <- t(apply(sim_sort , 1 , sort , decreasing = T))

plot((1:k)/(k+1) , rev(sim_sort_sort[,1]) , type = "l" , col = "blue" , log = "y") # max
lines((1:k)/(k+1) , rev(sim_sort_sort[,10000]) , col = "blue") # min
lines((1:k)/(k+1) , rev(sim_sort_sort[,250]) , col = "red") # 0.975
lines((1:k)/(k+1) , rev(sim_sort_sort[,9750]) , col = "red") # 0.025
lines((1:k)/(k+1) , z1 , col = "black") # theoretical
lines((1:k)/(k+1) , sort(daten_k[,2]) , col = "orange") # sample

plot((1:k)/(k+1) , sort(daten_k[,2]) , col = "orange" , log = "y") # sample

## Es gibt keinen Anhaltspunkt, dass die Grenzverteilung nicht stimmt. 