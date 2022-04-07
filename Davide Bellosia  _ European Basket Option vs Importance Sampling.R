############### European Basket Call Option ###############
####### Monte Carlo Method vs Importance Sampling #########

# Packages:
library(matrixStats)

# Input:
T     <- 1                  # Maturity date
steps <- 1                  # Time steps
m     <- 10000              # Number of simulations

I.v <- c(1, 1)              # Initial Value of the Assets
K   <- 2.1                  # Strike Price
r   <- 0.1                  # Interest rate
mu_vect  <- c(0.1, 0.1)     # Interest rate (No Arbitrage) = mu vector

n_assets <- 2               # Number of assets 
qt_a     <- c(1, 1)         # Quantity of each asset

sigma.1 <- 0.2              # SD Asset 1
sigma.2 <- 0.15             # SD Asset 2
sigma   <- c(0.2, 0.15)     # SD Vector
rho     <- -0.5             # Correlation Coefficent 

# Correlation between asset 1 and asset 2
corr_mat = matrix(c(1, -0.5, -0.5, 1), nrow = 2)


## Functions needed ##

# Arithmetic Brownian Motion N(O,I)
ABM <- function(T, m, n, corr_mat, mu, sigma, I.v){
  M <- matrix(rnorm(m*n), ncol = n)
  A <- t(chol(corr_mat))
  Z <- M%*%A
 dt <- T
dX1 <- mu*dt + sigma*sqrt(dt)*Z[,1]
dX2 <- mu*dt + sigma*sqrt(dt)*Z[,2]
 X1 <- rowCumsums(cbind(rep(I.v [1], m), dX1))
 X2 <- rowCumsums(cbind(rep(I.v [2], m), dX2))
 output <- list (X1,X2)
 return (output)}

# Geometric Brownian Motion N(O,I)
 GBM <- function(T, m, n, corr_mat, mu, sigma, I.v){
   X <- ABM (T, m, n, corr_mat, mu = mu_vect-(sigma^2)/2, sigma= sigma, I.v=log(I.v))
GBM1 <- exp(X[[1]])
GBM2 <- exp(X[[2]])
output <- list (GBM1, GBM2)
return (output)}

# Arithmetic Brownian Motion N(B,I)
ABM_B <- function(Beta, T, m, n, corr_mat, mu, sigma, I.v){
    M <- matrix(rnorm(m*n), ncol = n)
    N <- matrix((rep(t(Beta), m)), ncol = 2, byrow = T)
    G <- M + N
    A <- t(chol(corr_mat))
    Z <- G %*% A
   dt <- T
  dX1 <- mu*dt + sigma*sqrt(dt)*Z[,1]
  dX2 <- mu*dt + sigma*sqrt(dt)*Z[,2]
   X1 <- rowCumsums(cbind(rep(I.v [1], m), dX1))
   X2 <- rowCumsums(cbind(rep(I.v [2], m), dX2))
  output <- list (X1, X2, Z)
  return (output)}

# Geometric Brownian Motion N(B,I)
GBM_B <- function (Beta, T, m, n, corr_mat, mu, sigma, I.v){
    X <- ABM_B (Beta, T, m, n, corr_mat, mu= mu_vect-(sigma^2)/2, sigma=sigma, I.v= log(I.v))
 GBM1 <- exp(X[[1]])
 GBM2 <- exp(X[[2]])
    Z <- X[[3]]
  output <- list (GBM1, GBM2, Z)
  return (output)}


## Further Importance Sampling Input ##

# wi
weigth <- function(initial_value,r,sigma,qt_a){
     w <- rep(0,2)
   den <- 0
  for(i in 1:2){den  <- den+(qt_a[i]*initial_value[i]*exp(r-(sigma[i]^2/2))*(T))}
  for(i in 1:2){w[i] <- qt_a[i]*initial_value[i]*exp(r-(sigma[i]^2/2))*(T)/den}
  return (w)}
w <- weigth (I.v, r, sigma, qt_a)

# x0
x0 <- 0 # x0 = 0 (Only to initialize the variable)
for(i in 1:2){ x0 = x0 + (qt_a[i]*I.v[i]*exp(r-(sigma[i]^2/2))*(T))}

# ci
c <- rep(0,n_assets)
for (i in 1:n_assets) {c[i] = w[i]*sigma[i]}

## Lambda  
#to obtain Lambda we need to solve a non-linear equation, 
#to do this we use the build-in R-function uniroot
function_lambda= function(lambda){
  sqrt(T) * exp(lambda * sqrt(T)* t(sigma) %*% corr_mat %*% sigma)/
    (exp(lambda * sqrt(T)* t(sigma) %*% corr_mat %*% sigma)-K/x0) - lambda}

 lower <- (log(K / x0) )/(sqrt(T)* t(sigma) %*% corr_mat %*% sigma)
 upper <- 6 
lambda <- uniroot(function_lambda, c(lower,upper))$root

## Beta*
# to compute beta we use the formula: B* = lambda*b'*c
B = lambda*t(chol(corr_mat))%*%c


      ### Option pricing plain Monte Carlo ###
GBM_mc <- GBM(T, m, n_assets, corr_mat, r, sigma, I.v)
    V0 <- exp(-r*T)*mean(pmax((GBM_mc[[1]][,2] + GBM_mc[[2]][,2]) - K,0)) # Option Pricing with Monte Carlo Method
var_mc <- var(pmax((GBM_mc[[1]][,2] + GBM_mc[[2]][,2]) - K,0))            # Var Monte Carlo Method

      ### Option pricing Importance Sampling ###
GBM_beta <- GBM_B (B, T, m, n_assets, corr_mat, mu, sigma, I.v)
Z <- GBM_beta[[3]]
second_term <- rep(0,m)
for (i in 1:m) {second_term[i] <- exp(0.5*(t(B)%*%B)-t(B)%*%Z[i,])}

V0_beta <- exp(-r*T)*mean((pmax((GBM_beta[[1]][,2]+
               GBM_beta[[2]][,2]) - K,0)*second_term))             # Option Pricing with Importance Sampling

var_IS  <- var((pmax((GBM_beta[[1]][,2]  + GBM_beta[[2]][,2]) 
                 - K,0)*second_term))                              # Var Importance Sampling

       ########## RESULTS #############
RESULTS = matrix(c(V0, var_mc, V0_beta, var_IS), ncol = 2, byrow = T)
colnames(RESULTS) <- c("Option Price","Variance")
rownames(RESULTS) <- c("Monte Carlo Method", "Importance Sampling")
RESULTS

