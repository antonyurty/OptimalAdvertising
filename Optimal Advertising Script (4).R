# require(MittagLeffleR)
require(latex2exp)

K = function(t){
  t^2
}


########## Function that calculates kappas in Bernstein polynomials
calculate_kappa = function(K, n, Th = 1){
  kappa = numeric(n+1)
  
  for (k in 0:n){
    aux.vector = numeric(k+1)
    for (i in 0:k){
      aux.vector[i+1] = (-1)^(k - i) * K((Th*i)/n)*choose(n, i)*choose(n-i, k-i)
    }
    kappa[k+1] = 1/(Th^k)*sum(aux.vector)
  }
  kappa
}

########## Function that calculates Bernstein polynomial of order n that corresponds to K
#### Do not put n > 35, because calculations start being unstable
K_Bernstein = function(t, K, n, Th){
  kappa = calculate_kappa(K, n, Th)
  
  x = numeric(length(t))
  
  for(k in 0:n){
    x = x+ kappa[k+1]*t^k
  }
  
  x
}

################# TEST #################

plot((0:1000)/1000, K_Bernstein((0:1000)/1000, K, 10, 1), type = "l")
lines((0:1000)/1000, K((0:1000)/1000), col = "blue")

########################################


########## Function that calculates matrix of gammas; n is the degree of Bernstein polynomial
calculate_gamma = function(K, M, n, Th, beta, kappa = calculate_kappa(K, n, Th)){
  
  coeff = numeric(n+1)
  for(j in 1:(n+1)){
    coeff[j] = -beta*factorial(j)*kappa[j]
  }
  
  
  gamma = matrix(0, M+1, M+1)
  gamma[1,1] = 1
  
  for(k in 2:(M+1)){
    gamma[2:(M+1),k] = gamma[1:M, k-1]
    gamma[1,k] = sum(coeff[1:min(M+1, n+1)]*gamma[1:min(M+1, n+1),k-1])
  }
  
  gamma
}

################# TEST #################

calculate_gamma(K, 10, 5, 2, 1)

########################################


########## Function that calculates approximated optimal control

u_approx = function(t, K, M, n, Th=1, beta=1, alpha=1, a1=1, a2=1){
  
  kappa = calculate_kappa(K, n, Th)
  gamma = calculate_gamma(K, M, n, Th, beta, kappa)
  
  coeff = numeric(M+1)
  kappa_factorial = factorial(0:n)*kappa
  
  for(k in 1:(M+1)){
    coeff[k] = sum( gamma[1:min(k,n+1),k]*kappa_factorial[1:min(k,n+1)] )/factorial(k-1)
  }
  
  u = numeric(length(t))
  for(k in 1:(M+1)){
    u = u + coeff[k]*(Th - t)^(k-1)
  }
  
  u*(0.5*a2*alpha/a1)
}

################# TEST #################

ML = function(t, a, b, N){
  v = factorial(a*(0:N) + b - 1)
  
  res = numeric(length(t))
  for(i in 0:N){
    res = res + t^i/v[i+1]
  }
  res
}

plot((0:1000)/1000, ML(-2*(1-(0:1000)/1000)^3, 3, 3, 50)*(1-(0:1000)/1000)^2, type = "l", main = "Optimal control of Volterra Ornstein-Uhlenbeck process with monomial kernel, n=2", xlab = "t", ylab = "u(t)")
lines((0:1000)/1000, u_approx((0:1000)/1000, K, 100, 30), col = "red")
legend("topright", inset = 0.01, legend = c(TeX(r'(Theoretical $\widehat{u}(t)$)'), TeX(r'(Approximated $\widehat{u}_{n,M}(t)$, $n=30$, $M=100$)')), lty = c(1,1), col = c("black", "red"), box.lty=0)


#################################################################################
############################### Fractional kernel ###############################
#################################################################################

K1 = function(t){
  t^0.3
}

############################### Fractional kernel, different alphas
plot((0:1000)/1000, u_approx((0:1000)/1000, K1, 100, 30, alpha = 2), type = "l", col = "black", main = TeX(r"(Optimal control of fractional Volterra Ornstein-Uhlenbeck process, $K(t) = t^{0.3}$)"), xlab = "t", ylab = "u(t)")
lines((0:1000)/1000, u_approx((0:1000)/1000, K1, 100, 30, alpha = 1.5), col = "blue")
lines((0:1000)/1000, u_approx((0:1000)/1000, K1, 100, 30, alpha = 1), col = "red")
legend("bottomleft", inset = 0.01, legend = c( TeX(r'($\alpha = 2$)'), TeX(r'($\alpha = 1.5$)'), TeX(r'($\alpha = 1$)') ), lty = c(1,1, 1), col = c("black", "blue", "red"), box.lty=0)

############################### Fractional kernel, different betas
plot((0:1000)/1000, u_approx((0:1000)/1000, K1, 100, 30, beta = 1), type = "l", col = "black", main = TeX(r"(Optimal control of fractional Volterra Ornstein-Uhlenbeck process, $K(t) = t^{0.3}$)"), xlab = "t", ylab = "u(t)")
lines((0:1000)/1000, u_approx((0:1000)/1000, K1, 100, 30, beta = 1.5), col = "blue")
lines((0:1000)/1000, u_approx((0:1000)/1000, K1, 100, 30, beta = 2), col = "red")
legend("bottomleft", inset = 0.01, legend = c( TeX(r'($\beta = 1$)'), TeX(r'($\beta = 1.5$)'), TeX(r'($\beta = 2$)') ), lty = c(1,1, 1), col = c("black", "blue", "red"), box.lty=0)

##################################################################################
################################## Gamma kernel ##################################
##################################################################################

K2 = function(t){
  exp(-t)*(t)^0.3
}
plot((0:1000)/1000, K_Bernstein((0:1000)/1000, K2, 30, 1), type = "l")
lines((0:1000)/1000, K2((0:1000)/1000), col = "blue")

############################### Gamma kernel, different alphas
plot((0:1000)/1000, u_approx((0:1000)/1000, K2, 100, 30, alpha = 2), type = "l", col = "black", main = TeX(r"(Optimal control of Gamma Volterra Ornstein-Uhlenbeck process, $K(t) = e^t t^{0.3}$)"), xlab = "t", ylab = "u(t)")
lines((0:1000)/1000, u_approx((0:1000)/1000, K2, 100, 30, alpha = 1.5), col = "blue")
lines((0:1000)/1000, u_approx((0:1000)/1000, K2, 100, 30, alpha = 1), col = "red")
legend("bottomleft", inset = 0.01, legend = c( TeX(r'($\alpha = 2$)'), TeX(r'($\alpha = 1.5$)'), TeX(r'($\alpha = 1$)') ), lty = c(1,1, 1), col = c("black", "blue", "red"), box.lty=0)

############################### Gamma kernel, different betas
plot((0:1000)/1000, u_approx((0:1000)/1000, K2, 100, 30, beta = 1), type = "l", col = "black", main = TeX(r"(Optimal control of Gamma Volterra Ornstein-Uhlenbeck process, $K(t) = e^t t^{0.3}$)"), xlab = "t", ylab = "u(t)")
lines((0:1000)/1000, u_approx((0:1000)/1000, K2, 100, 30, beta = 1.5), col = "blue")
lines((0:1000)/1000, u_approx((0:1000)/1000, K2, 100, 30, beta = 2), col = "red")
legend("bottomleft", inset = 0.01, legend = c( TeX(r'($\beta = 1$)'), TeX(r'($\beta = 1.5$)'), TeX(r'($\beta = 2$)') ), lty = c(1,1, 1), col = c("black", "blue", "red"), box.lty=0)

##################################################################################
################################## Modified Gamma kernel ##################################
##################################################################################

K3 = function(t){
  exp(-t)*(t+1)^0.3
}
plot((0:1000)/1000, K_Bernstein((0:1000)/1000, K3, 30, 1), type = "l")
lines((0:1000)/1000, K3((0:1000)/1000), col = "blue")

############################### Modified Gamma kernel, different alphas
plot((0:1000)/1000, u_approx((0:1000)/1000, K3, 100, 30, alpha = 2), type = "l", col = "black", main = TeX(r"(Optimal control of modified Gamma Volterra Ornstein-Uhlenbeck process, $K(t) = e^t (t+1)^{0.3}$)"), xlab = "t", ylab = "u(t)", ylim = c(0,1))
lines((0:1000)/1000, u_approx((0:1000)/1000, K3, 100, 30, alpha = 1.5), col = "blue")
lines((0:1000)/1000, u_approx((0:1000)/1000, K3, 100, 30, alpha = 1), col = "red")
legend("topleft", inset = 0.01, legend = c( TeX(r'($\alpha = 2$)'), TeX(r'($\alpha = 1.5$)'), TeX(r'($\alpha = 1$)') ), lty = c(1,1, 1), col = c("black", "blue", "red"), box.lty=0)

############################### Modified Gamma kernel, different betas
plot((0:1000)/1000, u_approx((0:1000)/1000, K3, 100, 30, beta = 1), type = "l", col = "black", main = TeX(r"(Optimal control of modified Gamma Volterra Ornstein-Uhlenbeck process, $K(t) = e^t (t+1)^{0.3}$)"), xlab = "t", ylab = "u(t)", ylim = c(0,0.5))
lines((0:1000)/1000, u_approx((0:1000)/1000, K3, 100, 30, beta = 1.5), col = "blue")
lines((0:1000)/1000, u_approx((0:1000)/1000, K3, 100, 30, beta = 2), col = "red")
legend("topleft", inset = 0.01, legend = c( TeX(r'($\beta = 1$)'), TeX(r'($\beta = 1.5$)'), TeX(r'($\beta = 2$)') ), lty = c(1,1, 1), col = c("black", "blue", "red"), box.lty=0)

