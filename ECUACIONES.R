pacman::p_load(ggplot2,reshape2,ggpubr,dplyr)

k = 1.3807E-16 #cm2 g s-2 K-1

#k = 1.3807E-23 # #m2 kg s-2 K-1
# # Temp = 298
# # D = 2.93 
# n = 4 # CONSTANTE SUTHERLAND
# visc = 2.90
# M=4
# Mm =4
# ###1.16


func_1_16 <- function(k,temp,Dens,n,visc,M,Mw){
  A1 = (k*temp*(Dens^(1/3)))
  B1a = ((3/4)^(1/3))*(n*(pi^(2/3))*visc*(M^(1/3)))
  B2a = ((3/4)^(1/3))*(n*(pi^(2/3))*visc*(M^(4/3)))
  D_C= ((4/3)*(A1/B1a))-((1/3*(A1/B2a))*Mw)

  return(D_C)
}


###1.17
func_1_17 <- function(k,temp,D,n,visc,M){

A1 = (k*Temp*(D^(1/3)))
B1a = ((3/4)^(1/3))*(n*(pi^(2/3))*visc*(M^(1/3)))
B2a = ((3/4)^(1/3))*(n*(pi^(2/3))*visc*(M^(4/3)))
Mm0 = ((4/3)*(A1/B1a))/((1/3*A1/B2a))

return(Mm0)

}  


###1.18
func_1_18 <- function(k,temp,D,n,visc,M){
  
  A1 = (k*Temp*(D^(1/3)))
  B1a = ((3/4)^(1/3))*(n*(pi^(2/3))*visc*(M^(1/3)))
  a = ((4/3)*(A1/B1a))
  
  return(a)
  
}    


###1.19
func_1_19 <- function(k,temp,D,n,visc,M){
  
  A1 = (k*Temp*(D^(1/3)))
  B1a = ((3/4)^(1/3))*(n*(pi^(2/3))*visc*(M^(1/3)))
  B2a = ((3/4)^(1/3))*(n*(pi^(2/3))*visc*(M^(4/3)))
  m = ((4/3)*(A1/B1a))/((1/3*A1/B2a))
  return(m)
  
}    

#a=(4 /3)((K T ??^0.333333))/((0.9085603) n ??^0.666667 ??(M^0.333333)))
n_a <- function(k,temp,D,visc,M,a){
  
a1 <- 4*k*temp*(D^(1/3))
a2 <- ((3/4)^(1/3))*(pi^(2/3)) *3*a* visc*(M^(1/3))
n_a <- a1/a2
return(n_a)
}    
n_a_alt <- function(k,temp,D,visc,M,a){
  a1 <- 4*(2^(2/3))*(D^(1/3))*k*temp
  a2 <- 3*(3^(1/3))*(pi^(2/3))*a*visc*(M^(1/3))
  # a1 <- 0.68415*(D^(1/3))*k*temp #4*k*temp*(D^(1/3))
  # a2 <- a*visc*(M^(1/3))#((3/4)^(1/3))*(pi^(2/3)) *3*a* visc*(M^(1/3))
   n_a <- a1/a2
  return(n_a)
} 



n_m <- function(k,temp,D,visc,M,m){
  
  a1 <- k*temp*(D^(1/3))
  a2 <- 3*((3/4)^(1/3))*m*(pi^(2/3))* visc*(M^(4/3))
  n_m <- a1/a2
  return(n_m)
}


n_m_alt <- function(k,temp,D,visc,M,m){
  
  # a1 <- 0.171038*(D^(1/3))*k*temp#k*temp*(D^(1/3))
  # a2 <- m*visc*(M^(4/3))#3*((3/4)^(1/3))*m*(pi^(2/3))* visc*(M^(4/3))
  a1 <- (2^(2/3))*(D^(1/3))*k*temp
  a2 <- 3*(3^(1/3))*(pi^(2/3))*(M^(4/3))*m*visc
  n_m <- a1/a2
  return(n_m)
}

