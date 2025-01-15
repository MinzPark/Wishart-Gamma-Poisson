#install.packages("LaplacesDemon")


#############################################
# generate data
#############################################


# S: scale param pxp
# nu: degree freedom


set.seed(1)
N <- 100000
d <- 3


#-------------------------------------------
# 1) S with d=4
#-------------------------------------------

delta <- 4
S <- matrix(c( 0.25, 0.15,  0.15, 
                0.15,  0.25, 0.15, 
                0.15, 0.15, 0.25), nrow = 3, byrow = TRUE)



#lambs <- rep(c(0.3,0.6,0.9),N)
lambs <- rep(c(0.5), N*d)
Lambs <- matrix(lambs, nrow = N, ncol = d, byrow = TRUE)
dim(Lambs)




Us <- c()
for(i in 1:N){
  U <- LaplacesDemon::rwishart( delta, S)  
  #U <- rwishart(delta, p = nrow(S), S, SqrtSigma = diag(p))
  Us <- c(Us, diag(U))
}


Us.mat <- matrix(Us, nrow = N, byrow = TRUE)


cor(Us.mat)


Ns <- rpois(N * d, lambs * Us)
table(Ns)
N.vec = Ns
Ns <- matrix(N.vec, nrow = N, byrow = TRUE)


#############################################
# optim off- diag elt
#############################################
Sig_hat <- matrix(rep(NA, d*d), nrow = d, byrow = TRUE)
Sig_hat




# -------------------------------------------
# 1) optim diag
# -------------------------------------------


nb_ll <- function(lambs, N.vec, param){
  delta <- param[1]
  p = delta/ (2 * lambs + delta)
  r = delta/2
  
  return( -sum( dnbinom( N.vec, size = r, prob = p, log = TRUE)))
}




delta.hat <- optim(par = c(0), fn = nb_ll, N.vec = N.vec, lambs = lambs, method = 'Brent', lower = 0, upper = 10)$par
diag(Sig_hat) <- 1/delta.hat


# -------------------------------------------
# 1) optim diag
# -------------------------------------------


cal_num <- function(N, delta){
  delta2 <- 2/delta
  ret = 0
  if(N == 0){
    ret = 0
  }else{
    for(n in 1:N){
      ret = ret + log(1 +(n -1) * delta2)
    }
  }
  return(ret)
}


delta_f <- function( params, Lambs, Ns){
  delta <- params[1]
  Num <- apply(Ns, c(1:2), function(N) cal_num(N, delta))
  Den <- (delta/2 + Ns ) * log(1 + 2/delta * Lambs)
  l = -sum( (Num - Den) )
  
  
  return(l)
}


(delta.hat <- optim(par = c(0), fn = delta_f, Ns = Ns, Lambs = Lambs,
                    method = 'Brent', lower = 0, upper = 100)$par)
diag(Sig_hat) <- 1/delta.hat




# -------------------------------------------
# 2) optim off- diag elt
# -------------------------------------------
prop1 <- function(a,b,c, k1, k2, delta, lamb1, lamb2){
  ret = 0
  
  
  if(k1 * k2 == 0){
    
    ret = (exp(lgamma(delta/2 + k1 + k2)) * (a+c*lamb2)^(k1) * (b+c*lamb1)^(k2)) /
      (exp(lgamma(delta/2)) * (1 + a*lamb1 + b*lamb2 + c*lamb1*lamb2)^(delta/2 +k1 +k2))
    
  }else {
    ret = exp(lgamma(k1 + 1) + lgamma(k2 +1) - lgamma(delta/2))
    sum_t = 0
    for(k in 0:(min(k1,k2))){
      tmp.num =  ( (-c)^k * exp(lgamma(delta/2 + k1 +k2 -k)) * (a+c*lamb2)^(k1-k) * (b + c*lamb1)^(k2-k))
      tmp.den = exp(lgamma(k2 - k + 1) + lgamma(k1 - k +1) + lgamma(k+1)) *(1 + a*lamb1 + b*lamb2 + c*lamb1*lamb2)^(delta/2 + k1 + k2 -k)
      sum_t = sum_t + tmp.num/tmp.den
    }
    ret = ret * sum_t
  }
  return(ret)
}




## sigma j1j2


off_diag_elt <- function(param, j1, j2, Ns, Lamb, delta){
  
  l = c()

  
  sig_j1_j2 <- param[1]
  
  a <- 2 * 1/delta
  b <- 2 * 1/delta
  c <- (4/(delta^2) - 4*(sig_j1_j2)^2)
  
  for(i in 1:nrow(Ns)){
    k1 <- Ns[i,j1]
    k2 <- Ns[i,j2]
    lamb1 <- Lambs[i,j1]
    lamb2 <- Lambs[i,j2]
    
    pr1 <- prop1(a,b,c, k1, k2, delta, lamb1, lamb2)
    t.pr1 = (lamb1^(k1) * lamb2^(k2)  *  pr1)
    l <- c(l, t.pr1)
  }
  return(-sum(log(l)))
}




# pairwise estimation index
combination <-c()
for(i in 1:(d-1)){
  for(j in i:d){
    if(i!=j){
      combination <- rbind(combination,c(i,j))
    }
  }
}


# optimize off diag elt for Sigma
start.time <- Sys.time()
for(idx in 1:nrow(combination)){
  j1 = combination[idx, 1]
  j2 = combination[idx, 2]
  sigma.hat <- optim(par = c(0), fn = off_diag_elt,
                     j1 = j1, j2 = j2, Ns = Ns, Lamb = Lambs, delta = delta.hat,
                     method = 'Brent', lower = -1/delta.hat, upper = 1/delta.hat)$par
  Sig_hat[j1,j2] <- sigma.hat
  Sig_hat[j2,j1] <- sigma.hat
}
end.time <- Sys.time()
(exe.time <- end.time - start.time)




# -------------------------------------------
# 3) select the best Sigma estimation among sgn(Sigma)
# -------------------------------------------


Sig_ll <- function(M, Lambs, Ns, delta.hat, Sig_hat){
  if(min(eigen(Sig_hat)$value < 0) > 0){
    return("Some eigen value is neg")
  }else{  
    l_total = 0
    # Sigma : -179408.2  Sgn * Sigma -358873.8
    start.time <- Sys.time()
    
    Us <- c()
    d <- dim(Sig_hat)[1]
    for(j in 1:M){
      #U <- diag(rwishart(delta.hat, Sig_hat))
      U <- diag(dlm::rwishart(delta.hat, p = nrow(Sig_hat),Sig_hat, SqrtSigma = diag(p)))
      Us <- c(Us, U)
    }
    Us.mat <- matrix(Us, nrow = M, byrow = TRUE)
    
    for(k in 1:N){
      
      lamb <- Lambs[k,]
      N.p <- Ns[k,]
      ll = exp(-apply(sweep(Us.mat, 2, lamb, `*`),1,sum)) * apply(sweep(Us.mat, 2, N.p, `^`),1, function(x) prod(x))
      
      l_total = l_total + log(mean(ll))
    }
    
    end.time <- Sys.time()
    execution.time <- end.time - start.time
    print(execution.time)
    
    return(l_total)}
}






# find non D-similar Sign Sigma matrices ( # = 2^{(d-2)(d-1)/2} )


Sign.scen <- function(d, Sigma){
  
  sign.mat <- expand.grid(rep(list(c(-1, 1)), (d-2)*(d-1)/2))
  scen.l <- list()
  
  for(i in 1:nrow(sign.mat)){
    
    idx <- 1
    b.mat <- matrix(rep(1, d*d), nrow=d)
    sign.vec <- unlist(sign.mat[i,])
    
    for(j in (d-2):1){
      
      b.mat[(d-j)    , (d-j+1):d] <- sign.vec[idx:(idx+j-1)]
      b.mat[(d-j+1):d, (d-j)    ] <- sign.vec[idx:(idx+j-1)]
      
      idx <- idx+j
    }
    
    scen.l[[i]] <- b.mat * Sigma
  }
  return(scen.l)
}








####################




#m.set <- Sign.scen(d, S)
rep = 10
tmp <- c()

for(scn in 1:rep){
  m.set <- Sign.scen(d, Sig_hat)
  
  M<-1000
  sig.l.val <- c()
  i <- 1
  
  for(mat in m.set){
    cat(i, 'th case\n')
    sig.ll <- Sig_ll(M, Lambs, Ns, delta.hat, mat)
    sig.l.val <- c(sig.l.val, sig.ll)
    i = i +1 
    cat("\n")
  }
  tmp <- rbind(tmp, sig.l.val)
}


result <- cbind(seq(1,rep),tmp,tmp[,1] - tmp[,2])
colnames(result) <- c('iter',"Matrix1", "Matrix2", 'difference')

# save the result for Sigma scenario and likelihood ftn of each Sigma
# save sigma scenario
sigma_scenario <- "Sigma Scenario:\n"
for (i in seq_along(m.set)) {  
  sigma_scenario <- paste0(
    sigma_scenario, 
    "Matrix ", i, ":\n", 
    paste(capture.output(print(m.set[[i]])), collapse = "\n"), 
    "\n\n"  
  )
}
writeLines(sigma_scenario, paste0(filename,"/Sim_sigma_scenario.txt"))

# save the result of Simulation
filename = file.path(path,'wishart/result',gsub("-", "", Sys.Date()))
ifelse(!file.exists(filename), dir.create(filename, recursive = TRUE) ,0)

write.csv(result, paste0(filename,'/Sim_result.csv'), row.names = FALSE)



