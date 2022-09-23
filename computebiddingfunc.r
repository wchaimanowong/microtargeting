
# Distribution (CDF) of contextual value.
G <- function(c) {
  return(c)
}

# Derivative of G.
dG <- function(c) {
  return(1)
}

Gtilde <- function(v) {
  return(p*G(max((v-kappa)/(1-kappa),0))+(1-p)*G(min(v/(1-kappa),1)))
}

dGtilde <- function(v) {
  if (v < 1-kappa) {
    return((1-p)*dG(v/(1-kappa))/(1-kappa))
  }
  if ((v >= 1-kappa) && (v <= kappa)) {
    return(0)
  }
  if (v > kappa) {
    return(p*dG((v-kappa)/(1-kappa))/(1-kappa))
  }
}

interpolate_vec <- function(c,vec) {
  if ((c < 0) || (c > 1)) {
    return(0)
  }
  x <- (length(vec)-1)*c+1
  xi <- floor(x)
  xip <- xi + 1
  if (xip > length(vec)) {
    xip <- length(vec)
  }
  xr <- x %% 1
  interpolate_val <- (1-xr)*vec[xi] + xr*vec[xip]
  interpolate_val <- max(min(interpolate_val, 1),0)
  return(interpolate_val)
}

p <- 0.7
kappa <- 0.6
n1 <- 5
n2 <- 5

N <- 20 # Number steps in [0,1] interval
alpha <- 3.5 # Convergence speed.

dc <- 1/N

# The analytical equilibrium if beta < 1-kappa.
beta_low <- (1-kappa)*(0:N)/N

# Turning vector index into c location in [0,1] 
cval <- function(k) {
  return((k-1)*dc)
}

# Inverse a monotonically increasing function (return max beta^{-1}[0,f)).
Inverse_Fn <- function(f, vec) {
  k <- 1
  while ((vec[k] <= f) && (k <= length(vec))) {
    k <- k+1
  }
  k <- k-1
  c <- cval(k)
  return(c)
}

u <- function(beta, beta_, k) {
  
  Int1 <- 0
  j <- 1
  while (cval(j) < max((beta[k]-kappa)/(1-kappa),0)) {
    Int1 <- Int1 + dc*p*(1-kappa)*(cval(k)-cval(j))*(G(Inverse_Fn(kappa + (1-kappa)*cval(j), beta_))^(n2-1))*(n1*(G(cval(j))^(n1-1))*dG(cval(j)))
    j <- j+1
  }
  if (is.na(Int1)){
    Int1 <- 0
  }
  
  Int2 <- 0
  j <- 1
  while (cval(j) < Inverse_Fn(beta[k], beta_)) {
    Int2 <- Int2 + dc*p*(kappa + (1-kappa)*cval(k) - beta_[j])*(G(max((beta_[j]-kappa)/(1-kappa),0))^n1)*((n2-1)*(G(cval(j))^(n2-2))*dG(cval(j)))
    j <- j+1
  }
  if (is.na(Int2)){
    Int2 <- 0
  }
  
  Int3 <- 0
  j <- 1
  while (cval(j) < min(beta[k]/(1-kappa),1)) {
    Int3 <- Int3 + dc*(1-p)*(1-kappa)*(cval(k)-cval(j))*(G(Inverse_Fn((1-kappa)*cval(j), beta_))^(n2-1))*(n1*(G(cval(j))^(n1-1))*dG(cval(j)))
    j <- j+1
  }
  if (is.na(Int3)){
    Int3 <- 0
  }
  
  Int4 <- 0
  j <- 1
  while (cval(j) < Inverse_Fn(beta[k], beta_)) {
    Int4 <- Int4 + dc*(1-p)*((1-kappa)*cval(k) - beta_[j])*(G(min(beta_[j]/(1-kappa),1))^n1)*((n2-1)*(G(cval(j))^(n2-2))*dG(cval(j)))
    j <- j+1
  }
  if (is.na(Int4)){
    Int4 <- 0
  }
  
  return(Int1 + Int2 + Int3 + Int4)
}

beta <- kappa + (1-kappa)*(0:N)/N # Initial guess for beta > kappa.
max_iter <- 1000

ii <- 0
while (ii < max_iter) {
  beta_ <- beta
  for (i in 1:length(beta)) {
    dbeta <- rep(0,length(beta))
    dbeta[i] <- alpha*dc
    u_ <- u(beta, beta, i)
    up <- u(beta + dbeta, beta, i)
    um <- u(beta - dbeta, beta, i)
    
    if ((u_ < up) || (u_ < um)) {
      if (up >= um) {
        beta[i] <- min(beta[i] + (up-u_),1)
      }
      if (um > up) {
        beta[i] <- max(beta[i] - (um-u_), kappa)
      }
    }
  }
  
  # The final bidding function (also will be denoted by beta) is beta_low[i] at cval(i) if u(beta_low, beta_low, i) > u(beta, beta, i) and beta[i] otherwise.
  if (ii > max_iter/2) {
    for (i in 1:length(beta)) {
      u1 <- u(beta_low, beta_low, i)
      u2 <- u(beta, beta, i)
      if (u1 > u2) {
        beta[i] <- beta_low[i]
      }
    }
  }

  Dbeta <- sum((beta-beta_)^2)
  ii <- ii + 1
  print(paste("iter = ", ii, "Dbeta = ", Dbeta))
}

# for (i in 1:length(beta)) {
#   u1 <- u(beta_low, beta_low, i)
#   u2 <- u(beta, beta, i)
#   if (u1 > u2) {
#     beta[i] <- beta_low[i]
#   }
# }

plot((0:N)/N, beta, type='o', xlab="c", ylab="beta", main="Bidding function")

deviation_test <- function(beta) {
  pass <- TRUE
  dup <- rep(0, length(beta))
  dum <- rep(0, length(beta))
  for (i in 1:length(beta)) {
    dbeta <- rep(0, length(beta))
    dbeta[i] <- alpha*dc
    u_ <- u(beta, beta, i)
    dup[i] <- u(beta + dbeta, beta, i) - u_
    dum[i] <- u(beta - dbeta, beta, i) - u_
    
    pass <- pass && ((dup[i] <= 0) && (dum[i] <= 0))
  }
  return(list(pass, dup, dum))
}

global_deviation_test <- function(beta) {
  pass <- TRUE
  for (i in 1:length(beta)) {
    dbeta <- rep(0, length(beta))
    dbeta[i] <- 1
    u_ <- u(beta, beta, i)
    for (j in 1:length(beta)) {
      pass <- pass && (u_ >= u((j-1)*dc*dbeta, beta, i))
    }
  }
  return(pass)
}

deviation_test(beta)

################################################################################
# Compute common--value publisher revenues

w_pin = 0
for (i in 1:N) {
  w_pin <- w_pin + dc*(n1+n2)*(n1+n2-1)*(kappa*p+(1-kappa)*cval(i))*(1-G(cval(i)))*(G(cval(i))^(n1+n2-2))*dG(cval(i))
}

w_mtia = 0
for (i in 1:N) {
  w_mtia <- w_mtia + dc*p*n1*n2*(kappa + (1-kappa)*cval(i))*(1 - G(Inverse_Fn(kappa + (1-kappa)*cval(i), beta)))*(G(cval(i))^(n1-1))*(G(Inverse_Fn(kappa + (1-kappa)*cval(i), beta))^(n2-1))*dG(cval(i))
  w_mtia <- w_mtia + dc*p*n1*n2*beta[i]*(1-G(max((beta[i]-kappa)/(1-kappa),0)))*(G(max((beta[i]-kappa)/(1-kappa),0))^(n1-1))*(G(cval(i))^(n2-1))*dG(cval(i))
  if (n1 > 1) {
    w_mtia <- w_mtia + dc*p*n1*(n1-1)*(kappa + (1-kappa)*cval(i))*(1 - G(cval(i)))*(G(cval(i))^(n1-2))*(G(Inverse_Fn(kappa + (1-kappa)*cval(i), beta))^n2)*dG(cval(i))
  }
  if (n2 > 1) {
    w_mtia <- w_mtia + dc*p*n2*(n2-1)*beta[i]*(1-G(cval(i)))*(G(max((beta[i]-kappa)/(1-kappa),0))^n1)*(G(cval(i))^(n2-2))*dG(cval(i))
  }
  
  w_mtia <- w_mtia + dc*(1-p)*n1*n2*(1-kappa)*cval(i)*(1 - G(Inverse_Fn((1-kappa)*cval(i), beta)))*(G(cval(i))^(n1-1))*(G(Inverse_Fn((1-kappa)*cval(i), beta))^(n2-1))*dG(cval(i))
  w_mtia <- w_mtia + dc*(1-p)*n1*n2*beta[i]*(1-G(min(beta[i]/(1-kappa),1)))*(G(min((beta[i]-kappa)/(1-kappa),1))^(n1-1))*(G(cval(i))^(n2-1))*dG(cval(i))
  if (n1 > 1) {
    w_mtia <- w_mtia + dc*(1-p)*n1*(n1-1)*(1-kappa)*cval(i)*(1 - G(cval(i)))*(G(cval(i))^(n1-2))*(G(Inverse_Fn((1-kappa)*cval(i), beta))^n2)*dG(cval(i))
  }
  if (n2 > 1) {
    w_mtia <- w_mtia + dc*(1-p)*n2*(n2-1)*beta[i]*(1-G(cval(i)))*(G(min(beta[i]/(1-kappa),1))^n1)*(G(cval(i))^(n2-2))*dG(cval(i))
  }
}

w_ct = w_pin

print("Compute common--value publisher revenues:")
print(paste("w_pin =", w_pin, "w_mtia =", w_mtia, "w_ct =", w_ct))

################################################################################
# Compute common--value conversion rates

v_pin = 0
for (i in 1:N) {
  v_pin <- v_pin + dc*(n1+n2)*(kappa*p+(1-kappa)*cval(i))*(G(cval(i))^(n1+n2-1))*dG(cval(i))
}

v_mtia = 0
for (i in 1:N) {
  v_mtia <- v_mtia + dc*p*n1*(kappa + (1-kappa)*cval(i))*(G(cval(i))^(n1-1))*(G(Inverse_Fn(kappa + (1-kappa)*cval(i), beta))^n2)*dG(cval(i))
  v_mtia <- v_mtia + dc*p*n2*(kappa + (1-kappa)*cval(i))*(G(max((beta[i]-kappa)/(1-kappa),0))^n1)*(G(cval(i))^(n2-1))*dG(cval(i))
  
  v_mtia <- v_mtia + dc*(1-p)*n1*(1-kappa)*cval(i)*(G(cval(i))^(n1-1))*(G(Inverse_Fn((1-kappa)*cval(i), beta))^n2)*dG(cval(i))
  v_mtia <- v_mtia + dc*(1-p)*n2*(1-kappa)*cval(i)*(G(min(beta[i]/(1-kappa),1))^n1)*(G(cval(i))^(n2-1))*dG(cval(i))
}

v_ct = v_pin

print("Compute common--value conversion rates:")
print(paste("v_pin =", v_pin, "v_mtia =", v_mtia, "v_ct =", v_ct))

################################################################################
# Compute independent--value publisher revenues

w_pin = 0
for (i in 1:N) {
  w_pin <- w_pin + dc*(n1+n2)*(n1+n2-1)*cval(i)*(1-Gtilde(cval(i)))*(Gtilde(cval(i))^(n1+n2-2))*dGtilde(cval(i))
}

w_mtia = 0
for (i in 1:N) {
  if (n1 > 1) {
    w_mtia <- w_mtia + dc*n1*(n1-1)*cval(i)*(1-Gtilde(cval(i)))*(Gtilde(cval(i))^(n1-2))*(G(min(max((cval(i)-kappa*p)/(1-kappa),0),1))^n2)*dGtilde(cval(i))
  }
  w_mtia <- w_mtia + dc*n1*n2*cval(i)*(1-G(min(max((cval(i)-kappa*p)/(1-kappa),0),1)))*(Gtilde(cval(i))^(n1-1))*(G(min(max((cval(i)-kappa*p)/(1-kappa),0),1))^(n2-1))*dGtilde(cval(i))
  w_mtia <- w_mtia + dc*n1*n2*(kappa*p + (1-kappa)*cval(i))*(1-Gtilde(kappa*p+(1-kappa)*cval(i)))*(Gtilde(kappa*p+(1-kappa)*cval(i))^(n1-1))*(G(cval(i))^(n2-1))*dG(cval(i))
  if (n2 > 1) {
    w_mtia <- w_mtia + dc*n2*(n2-1)*(kappa*p + (1-kappa)*cval(i))*(1-G(cval(i)))*(Gtilde(kappa*p+(1-kappa)*cval(i))^n1)*(G(cval(i))^(n2-2))*dG(cval(i))
  }
}

w_ct = 0
for (i in 1:N) {
  w_ct <- w_ct + dc*(n1+n2)*(n1+n2-1)*(kappa*p+(1-kappa)*cval(i))*(1-G(cval(i)))*(G(cval(i))^(n1+n2-2))*dG(cval(i))
}

print("Compute independent--value publisher revenues:")
print(paste("w_pin =", w_pin, "w_mtia =", w_mtia, "w_ct =", w_ct))

################################################################################
# Compute independent--value conversion rates

v_pin = 0
# for (i in 1:N) {
#   v_pin <- v_pin + dc*(n1+n2)*cval(i)*(Gtilde(cval(i))^(n1+n2-1))*dGtilde(cval(i))
# }
for (i in 1:N) {
  v_pin <- v_pin + dc*(1-p)*(n1+n2)*(1-kappa)*cval(i)*(((1-p)*G(cval(i)))^(n1+n2-1))*dG(cval(i))
  v_pin <- v_pin + dc*p*(n1+n2)*(kappa+(1-kappa)*cval(i))*((p*G(cval(i))+1-p)^(n1+n2-1))*dG(cval(i))
}

v_mtia = 0
for (i in 1:N) {
  v_mtia <- v_mtia + dc*n1*cval(i)*(Gtilde(cval(i))^(n1-1))*(G(min(max((cval(i)-kappa*p)/(1-kappa),0),1))^n2)*dGtilde(cval(i))
  v_mtia <- v_mtia + dc*n2*(kappa*p + (1-kappa)*cval(i))*(Gtilde(kappa*p+(1-kappa)*cval(i))^n1)*(G(cval(i))^(n2-1))*dG(cval(i))
}

v_ct = 0
for (i in 1:N) {
  v_ct <- v_ct + dc*(n1+n2)*(kappa*p+(1-kappa)*cval(i))*(G(cval(i))^(n1+n2-1))*dG(cval(i))
}

print("Compute independent--value conversion rates:")
print(paste("v_pin =", v_pin, "v_mtia =", v_mtia, "v_ct =", v_ct))

################################################################################
