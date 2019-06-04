## Schaetzung der Parameter beta in (0,1] und theta in (0,1] der Mischverteilung 
## des Dirac-Mass am Pkt 0 und der ML(beta,theta^(-1/beta))
## (1-theta) * Dirac-Mass + theta * ML-Mass
## mithilfe von fraktionierten Momenten
## test test



# theoretisches fraktioniertes Moment (quelle: Cahoy, 2013. Estimation of Mittag-Leffler parameters.):
# Voraussetzung an das frak. Moment: q < beta (sonst existiert es nicht)
frac_mom <- function( beta , theta , q ){
  r <- ( q*pi*theta^(-1/beta*q+1) )/( beta*gamma(1-q)*sin((pi*q)/beta) )  
  return(r)
}

# empirisches fraktioniertes Moment:
frac_mom_emp <- function( q , data ){
  r <- 1/(length(data))*sum(data^q)
  return(r)
}

## nun versuchen wir mithilfe der Subtraktion der theoretischen und empirischen Momente
## die beiden Parameter zu bestimmen, indem wir die Nullstellen einer von frac_mom und frac_mom_emp abhaengigen Fkt. suchen

# die nach (beta,theta) zu optimierende Funktionen (2 Varianten):
fct <- function( beta , theta , q , data ){                     
  r <- frac_mom( beta , theta , q ) - frac_mom_emp( q , data )
  return(r)
}

fct2 <- function( beta , theta , q1 , q2 , data){
  r <- frac_mom( beta , theta , q1 )/frac_mom( beta , theta , q2) - 
    frac_mom_emp( q1 , data )/frac_mom_emp( q2 , data )
  return(r)
}

## zunächst plotten wir die Funktionen auf einem Gitter (0,1]x(0,1]
## um ein Gefuehl dafuer zu bekommen:

## Simulation von passenden Daten

theta <- 0.7
beta <- 0.9
n <- 10000
prob <- sample( c(0,1) , size = n , replace = T , prob = c(1-theta , theta))
daten <- sapply(prob , function(x){ ifelse( x == 0 , return(0) , 
                                           return(MittagLeffleR::rml( 1 , tail = beta , scale = theta^(-1/beta) )) )})

## Plots erstellen:

g1 <- seq(0.01 , 1 , 0.01)
g2 <- seq(0.01 , 1 , 0.01)
gitter <- expand.grid(g1,g2)

q <- 0.6
q1 <- 0.6
q2 <- 0.3
# Werte von fct
out <- mapply(function(x,y){ fct(x , y , q = q , data = daten) } , gitter$Var1 , gitter$Var2)

# Werte von fct2
out2 <- mapply(function(x,y){ fct2(x , y , q1 = q1 , q2 = q2 , data = daten) } , gitter$Var1 , gitter$Var2)


# ggplot
data <- data.frame(cbind(gitter$Var1,gitter$Var2,out))
# Nur Daten die nahe Null sind
data2 <- subset( abs(data) , out < 0.1 )

ggplot(data2 , aes(x = V1 , y = V2 , color = out)) +geom_point() 


data <- data.frame(cbind(gitter$Var1,gitter$Var2,out2))
data2 <- subset( abs(data) , out2 < 0.1 )

ggplot(data2 , aes(x = V1 , y = V2 , color = out2)) + geom_point() 



# 3D-Plots (nicht sehr aussagekraeftig)
par(mfrow=c(1,2))
persp(g1, g2, matrix(out , col = 100), col="lightgreen" , main="fct" , zlim = c(-25,50) , 
      ticktype = "detailed",nticks=2 , theta = 30 , phi = 30)

persp(g1, g2, matrix(out2 , col  = 100), col="lightgreen" , main="fct2" , zlim = c(-25,50) , 
      ticktype = "detailed",nticks=2 , theta = 30 , phi = 30)
