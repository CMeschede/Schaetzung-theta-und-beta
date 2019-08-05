# Example:  Solar Flare Data
# aus dem CTRE-Paket

## Daten:
# install.packages("MittagLeffleR)
# install.packages("CTRE")
# install.packages("extRemes")
# install.packages("fExtremes")
daten <- CTRE::flares
daten <- CTRE::ctre(daten)

magnitudes <- CTRE::magnitudes(daten)
interarrival <- CTRE::interarrival(daten)
data <- cbind(magnitudes,interarrival)
summary(magnitudes)
summary(interarrival)
plot(cumsum(sort(interarrival,decreasing = F)))
plot(sort(interarrival,decreasing = F) , log = "y")
plot(sort(interarrival,decreasing = F))
## functions ####
fct_thin <- function(data , k=NULL , u=NULL){
  JJ <- data[,1]
  WW <- data[,2]
  n <- length(JJ)
  if(is.null(k) & is.null(u)){
    stop("Es fehlt die Eingabe des Argumentes k oder u. Eins von beiden muss gegeben werden")
  }
  if(!is.null(k) & !is.null(u)){
    k_1 <- sum(JJ>u)
    if(k_1 != k){
      stop("Die Anzahl der Ueberschreitungen ", k, " und der Threshold ", u, "passen nicht zusammen. 
         Geben Sie nur eins von beiden an")
    }
  }
  if(is.null(k)){
    k <- sum(JJ>u)
  }
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

## Verteilung der Wartezeit implementieren 
## T_EI 

## Mischverteilung:    T_EI ~ (1-EI)*e_0 + EI*Exp(EI)
# EI \in (0,1) , 
# de_o  = Dirac-Verteilung am Punkt 0
# Exp(EI) = Exponentialverteilung mit Rate EI
## Verteilungsfunktion
# F(x) = 1-EI*exp(-EI*x) , x>=0
## Quantilsfunktion
# Q(y) = - (log(1-y)-log(EI))/(EI*(1-F(u)))

## Realisierung ziehen:

rT_EI <- function(n , EI){
  s <- sample(c(0,1) , size = n , replace = T , prob = c(1-EI,EI))
  rT_theta <- sapply(s , function(x){
    if(x == 1){
      rexp(1 , rate = EI) 
    }
    else {0}
  } )
  return(rT_theta)
}

pT_EI <- function(q , EI){
  x <- 1-EI*exp(-EI*q)
  x[q<0] <- 0
  return(x)
}

qT_EI <- function(p , EI){
  x <- -(log(1-p)-log(EI))/EI
  x[p<(1-EI)] <- 0
  return(x)
}

## T_u0 ~ T_EI/(1-F(u0))

## F(u0) wird durch die GP-Verteilung geschaetzt. Es wird daher fuer
## einen gewaehlten Threshold 'u < u0' der scale- und shape-Parameter geschaetzt
## F(u) selbst wird durch die empirische Vtlgsfkt geschaetzt.


rT_u0 <- function(n , EI , u0 , u , Fu , scale , shape){
  Tail <- (1-((1-Fu)*extRemes::pevd(u0-u, scale = scale , 
                          shape = shape , type = "GP")+Fu))
  rT_EI(n , EI)/Tail
}

pT_u0 <- function(q , EI , u0 , u , Fu , scale , shape){
  Tail <- (1-((1-Fu)*extRemes::pevd(u0-u, scale = scale , 
                          shape = shape , type = "GP")+Fu))
  pT_EI(q*Tail , EI)
}

qT_u0 <- function(p , EI , u0 , u , Fu , scale , shape){
  Tail <- (1-((1-Fu)*extRemes::pevd(u0-u, scale = scale , 
                          shape = shape , type = "GP")+Fu))
  qT_EI(p , EI)/Tail
}

fct_EI_Interval <- function(wait){
  if(max(wait)>2){
    min(2*(sum(wait[-1]-1))^2/( length(wait[-1])*sum((wait[-1]-1)*(wait[-1]-2)) ),1)
  }
  else{
    min(2*sum(wait[-1])^2/(length(wait[-1])*sum(wait[-1]^2)),1)
  }
}

fct_myextremalindex <- function(data , u0 = NULL , q0 = seq(0.8,0.999,0.001) ){
  magnitudes <- data[,1]
  if( any(u0 > max(magnitudes)) ){
    stop( "die thresholds duerfen nicht groesser als das maximum der magnitudes sein ")
  }

  if(is.null(u0)){
    u0 <- quantile(magnitudes , q0 , names = F)
  }
  data_u0 <- sapply(u0 , function(x){fct_thin(data , u=x)})
  x <- sapply(data_u0 , function(y){
    fct_EI_Interval(y[,2])
  })
  return(EI = x)
}

rmse <- function(error){
  sqrt(mean(error^2))
}

## ####

# Model 1: ####
# X1,...,Xn uiv magnitudes with deterministic, equidistant waiting times
# -> we have to estimate the tailfunction with the extreme value theory

q <- round(seq(0.80, 0.999 , 0.001),3) 
u <- quantile(data[,1], q , names = F); u # vector
k <- sapply(u , function(x){sum(data[,1]>x)}); k # vector
CTREmodel <- sapply(k , function(x){fct_thin(data , k = x)}) # list

fit_GP <- sapply( u , function(y){
  extRemes::fevd( data[,1] , threshold = y , type = "GP")$results$par
  #z <- extRemes::fevd( x[,1] , threshold = y , type = "GP")
  #extRemes::ci.fevd(z , type = "parameter" , method = "boot")
})

# stability plots
plot(q*100 , fit_GP[2,] , ylab = "shape" , xlab = "Threshold Quantile" , type = "o")
plot(q*100 , fit_GP[1,]-u , ylab = "shape" , xlab = "Threshold Quantile" , 
     type = "o")
# keine Stabilisierung zu sehen


shape1 <- extRemes::fevd( data[,1] , threshold = quantile(data[,1], 0.95 , names = F) , 
                         type = "GP")$results$par[2]; shape
scale1 <- extRemes::fevd( data[,1] , threshold = quantile(data[,1], 0.95 , names = F) , 
                         type = "GP")$results$par[1]; scale

u0 <- u[q==0.95]

# qqplot der der angepassten GP-Verteilung an die Daten
data_u0 <- fct_thin(data , u = u0)
excess_u0 <- data_u0[,1]-u0

x0 <- sapply(1:length(excess_u0) , function(x){extRemes::qevd(x/(length(excess_u0)+1) , scale = scale1 , 
                                             shape = shape1 , type = "GP")})
plot(sort(x0) , sort(excess_u0) , xlab = "Quantile der GP Verteilung" , 
     ylab = "Exzesse zum Threshold u0" , xlim = c(0,300000) , ylim = c(0,300000)
     )
qqline(excess_u0 , distribution = function(p) extRemes::qevd(p , scale = scale1 , shape = shape1 , 
                                            type = "GP") ,
       col = "red")
plot(sort(x0) , sort(excess_u0) , xlab = "Quantile der GP Verteilung" , 
     ylab = "Exzesse zum Threshold u0" , xlim = c(0,5000) , ylim = c(0,5000)
)
qqline(excess_u0 , distribution = function(p) extRemes::qevd(p , scale = scale1 , shape = shape1 , 
                                                             type = "GP") ,
       col = "red")
## sieht okay aus

# qqplot zur Wartezeitverteilung:
Fu0 <- 1-sum(data[,1]>u0)/sum(data[,2])
z0 <- sapply(1:length(data_u0[,2]) , function(x){qexp(x/(length(data_u0[,2])+1))/(1-Fu0)
  })
plot(sort(z0) , sort(data_u0[,2]) , xlab = "Quantile der Exponential Verteilung" , 
     ylab = "Exzess waiting times u0" #, ylim = c(0,20000) , xlim = c(0,20000)
     )
points(c(0,5000000) , c(0,5000000) , type = "l" , col = "red")

rmse(abs(sort(z0)-sort(data_u0[,2])))
# 378474.9

u1 <- u[q==0.99]

data_u1 <- fct_thin(data , u = u1)

Fu1 <- (1-Fu0)*extRemes::pevd(u1-u0 , scale = scale1 , shape = shape1 , type = "GP")+Fu0
z1 <- sapply(1:length(data_u1[,2]) , function(x){qexp(x/(length(data_u1[,2])+1))/(1-Fu1)})
plot(sort(z1) , sort(data_u1[,2]) , xlab = "Quantile der Exponential Verteilung" , 
     ylab = "Exzess waiting times u0" #, ylim = c(0,10000)
)
points(c(0,5000000) , c(0,5000000) , type = "l" , col = "red")

rmse(abs(sort(z1)-sort(data_u1[,2])))
# 1818417
## long waiting times would be underestimated -->  baaad

# Model 2: ####
# X1,...,Xn strict stationary sequence of magnitudes with deterministic, equidistant waiting times
# -> we have to estimate extremal index and the tailfunction with the extreme value theory

q <- round(seq(0.80, 0.999 , 0.001),3) 
u <- quantile(data[,1], q , names = F); u # vector
k <- sapply(u , function(x){sum(data[,1]>x)}); k # vector
CTREmodel <- sapply(k , function(x){fct_thin(data , k = x)}) # list

ei <- fct_myextremalindex(data , u0 = u)
summary(ei)
plot(q, ei , ylim = c(0,1) , ylab = "Interval-Schaetzer")
abline(h = 0.33 , col = "red")

EI <- 0.33 # relative große Cluster
Nc <- ceiling(k*EI)
q <- q[Nc>4]
u <- u[Nc>4]
k <- k[Nc>4]
Nc <- Nc[Nc>4]

cluster_stichprobe <- sapply(u , function(x){
  k <- sum(data[,1]>x)
  nc <- ceiling(k*EI)
  data_u <- fct_thin(data , u = x)
  mag <- data_u[,1]
  wait <- data_u[-1,2]
  owait <- sort(wait , decreasing = T)
  if(nc == 1){
    w <- 1
    r <-max(wait)+1
    } else if(nc == (length(wait)+1)){
    w <- nc
    r <- 0
    } else if(nc > (length(wait)+1)){
    w <- length(wait)+1
    r <- 0
    } else{
  r <- owait[nc-1]
  w <- nc # anzahl cluster
    if(owait[nc-1] == owait[nc]){
      if( any(owait > owait[nc]) ){
        r <- min(owait[which(owait>owait[nc])])
        w <- max(which(wait == r))+1
      } else{ 
        w <- 1
        r <- max(owait)+1
        }
    } else{
      w <- w # Anzahl Cluster
      r <- r # Länge der kuerzesten Interclusterzeit
      }
  }
  c <- rep(1,length(wait)+1)
  c[2:length(c)] <- 1 + cumsum(wait >= r)
  max_mag <- aggregate(mag , list(c) , max)$x
  max_wait <- aggregate(data_u[,2] , list(c) , sum)$x
  return(cbind(max_mag,max_wait))
})

fit_GP_max <- mapply( function(x,y){
  extRemes::fevd( x[,1] , threshold = y , type = "GP")$results$par
} , cluster_stichprobe , u)

# stability plots
plot(q*100 , fit_GP_max[2,] , ylab = "shape" , xlab = "Threshold Quantile" , type = "o")
points(q*100 , fit_GP[2,1:length(q)] , col = "red" , type = "o")
plot(q*100 , fit_GP_max[1,]-u , ylab = "shape" , xlab = "Threshold Quantile" , 
     type = "o" , log = "y")
points(q*100 , fit_GP[1,1:length(q)] , col = "red" , type = "o")

u2 <- u[q==0.90]
shape2 <- fit_GP_max[2,q==0.90]; shape2
scale2 <- fit_GP_max[1,q==0.90]; scale2

# qqplot der der angepassten GP-Verteilung an die Daten

excess_u2 <- cluster_stichprobe[[which(q==0.9)]][,1]-u2

x2 <- sapply(1:length(excess_u2) , function(x){extRemes::qevd(x/(length(excess_u2)+1) , scale = scale2 , 
                                                              shape = shape2 , type = "GP")})
plot(sort(x2) , sort(excess_u2) , xlab = "Quantile der GP Verteilung" , 
     ylab = "Exzesse zum Threshold u0" #, xlim = c(0,300000) , ylim = c(0,300000)
)
qqline(excess_u2 , distribution = function(p) extRemes::qevd(p , scale = scale2 , shape = shape2 , 
                                                             type = "GP") ,
       col = "red")
plot(sort(x2) , sort(excess_u2) , xlab = "Quantile der GP Verteilung" , 
     ylab = "Exzesse zum Threshold u0" , xlim = c(0,5000) , ylim = c(0,5000)
)
qqline(excess_u2 , distribution = function(p) extRemes::qevd(p , scale = scale2 , shape = shape2 , 
                                                             type = "GP") ,
       col = "red")
## naja

# qqplot zur Wartezeitverteilung:
Fu2 <- 1-sum(data[,1]>u2)/sum(data[,2])
data_u2 <- fct_thin(data , u = u2)
z2 <- sapply(1:length(data_u2[,2]) , function(x){qT_EI(x/(length(data_u2[,2])+1) , EI = EI )/(1-Fu2)})
plot(sort(z2) , sort(data_u2[,2]) , xlab = "Quantile der Exponential Verteilung" , 
     ylab = "Exzess waiting times u0" #, ylim = c(0,10000)
)
points(c(0,5000000) , c(0,5000000) , type = "l" , col = "red")

rmse(abs(sort(z2)-sort(data_u2[,2])))
# 77388.06

Fu01 <- (1-Fu2)*extRemes::pevd(u0-u2 , scale = scale2 , shape = shape2 , type = "GP")+Fu2
z01 <- sapply(1:length(data_u0[,2]) , function(x){qT_EI(x/(length(data_u0[,2])+1) , EI = EI )/(1-Fu01)})
plot(sort(z01) , sort(data_u0[,2]) , xlab = "Quantile der Exponential Verteilung" , 
     ylab = "Exzess waiting times u0" #, ylim = c(0,10000)
)
points(c(0,5000000) , c(0,5000000) , type = "l" , col = "red")

rmse(abs(sort(z01)-sort(data_u0[,2])))  
# 257202.1  

u3 <- u[q==0.99]
data_u3 <- fct_thin(data , u = u3)

Fu3 <- (1-Fu2)*extRemes::pevd(u3-u2 , scale = scale2 , shape = shape2 , type = "GP")+Fu2
z3 <- sapply(1:length(data_u3[,2]) , function(x){qT_EI(x/(length(data_u3[,2])+1) , EI = EI )/(1-Fu3)})
plot(sort(z3) , sort(data_u3[,2]) , xlab = "Quantile der Exponential Verteilung" , 
     ylab = "Exzess waiting times u0" #, ylim = c(0,10000)
)
points(c(0,5000000) , c(0,5000000) , type = "l" , col = "red")

rmse(abs(sort(z3)-sort(data_u3[,2]))) 
# 1724185

## das scheint es ganz sicher nicht zu sein!!

# Model 3: ####
# X1,...,Xn uiv magnitudes with heavy-tailed distributed waiting times
# -> Mittag-Leffler

q <- round(seq(0.80, 0.999 , 0.001),3) 
u <- quantile(data[,1], q , names = F); u # vector
k <- sapply(u , function(x){sum(data[,1]>x)}); k # vector
CTREmodel <- sapply(k , function(x){fct_thin(data , k = x)}) # list

fit_ML <- sapply(CTREmodel , function(x){
  MittagLeffleR::logMomentEstimator(x[,2])
})

fit_ML_MLE <- sapply(CTREmodel , function(x){
  MittagLeffleR::mlmle(x[,2])$par
})

plot(q*100,fit_ML[1,] , type = "l" , ylim = c(0,1.5))
abline(h = 0.9 , col = "red")

plot(q*100,fit_ML_MLE[1,] , type = "l" , ylim = c(0,1.5))
abline(h = 0.83 , col = "red")

shape_ML <- 0.83

plot(q*100 , fit_ML_MLE[2,]*(k/sum(data[,2]))^{1/shape_ML} , type = "l")
abline(h = 0.05 , col = "red")
scale0_ML <- 0.3

## qqplot fuer die angepasste ML-Verteilung

u0
k0 <- sum(data[,1]>u0)
data_u0 <- fct_thin(data , k = k0)

y1 <- sapply(1:length(data_u0[,2]) , function(p){
  MittagLeffleR::qml(p/(length(data[,2])+1) , tail = shape_ML , 
                     scale = scale0_ML*(k0/sum(data[,2]))^{-1/shape_ML})
})

plot(sort(y1) , sort(data_u0[,2]) , ylab = "empirische Quantile" , 
     xlab = "theoretische Quantile" #, ylim = c(0,10000) 
     , cex = 0.5
)
points(c(0,2e+07),c(0,2e+07), type ="l" , col = "red")

rmse(abs(sort(y1)-sort(data_u0[,2])))
# 614963.7

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

k <- length(CTRE::interarrival(daten)):10
q1 <- 0.01
q2 <- 0.05
est <- sapply( k , function(x){
  daten <- CTRE::interarrival(CTRE::thin(daten2, k = x))
  daten <- floor(daten)
  beta_hat1 <- fct_beta(q1 , q2 , data = daten)
  theta_hat11 <- fct_theta(beta = beta_hat1 , q = q1 , data = daten)
  theta_hat12 <- fct_theta(beta = beta_hat1 , q = q2 , data = daten)
  return(c(beta_hat1,theta_hat11,theta_hat12))
} , simplify = "array")
plot(k,est[1,])
plot(k,est[2,])
plot(k,est[3,])
