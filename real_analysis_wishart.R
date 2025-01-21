
# install.packages("CASdatasets", repos = "http://dutangc.free.fr/pub/RRepos/", type="source")
# #or 
# install.packages("CASdatasets", repos = "http://dutangc.perso.math.cnrs.fr/RRepository/", type="source")
# #or 
# install.packages("CASdatasets", repos = "https://cas.uqam.ca/pub/", type="source")
# library(CASdatasets)
library(mvtnorm)

filepath <- "/Users/baengminji/Desktop"
#############################################
# import real data set
#############################################

#-------------------------------------------
# 1) S with d=3
#-------------------------------------------
path <- '/Users/baengminji/Desktop/wishart/data'


Ns <- as.matrix(read.csv(paste0(path,'/N_mat.csv')))
N.vec <- c(t(Ns))

N <- dim(Ns)[1]
d <- dim(Ns)[2]

Lambs <- as.matrix(read.csv(paste0(path,'/lamb_mat.csv')))
lambs <- c(t(Lambs))

dim(Lambs)



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



delta.result <- optim(par = c(0), fn = nb_ll, N.vec = N.vec, lambs = lambs, method = 'Brent', lower = 0, upper = 10)
delta.hat <- optim(par = c(0), fn = nb_ll, N.vec = N.vec, lambs = lambs, method = 'Brent', lower = 0, upper = 10)$par
diag(Sig_hat) <- 1/delta.hat



# -------------------------------------------
# 0) plot optim delta with nb-likelihood ftn
# -------------------------------------------

x.vec <- seq(0.0001,4,5/10000)
y.vec <-c()
for(i in x.vec){
  l <-  nb_ll(lambs, N.vec, i)
  y.vec <- c(y.vec,l)
}

plot(x.vec, y.vec, cex = 0.5)
points(as.numeric(delta.result$par),as.numeric(delta.result$value), col = "red", pch = 19)


filename = file.path(filepath,'wishart/result',gsub("-", "", Sys.Date()))
ifelse(!file.exists(filename), dir.create(filename, recursive = TRUE) ,0)
png(paste0(filename,"/delta.png"), width = 800, height = 600) 

par(mar = c(6, 4, 3, 2))  
plot(x.vec[which(x.vec > 0.3)], -y.vec[which(x.vec > 0.3)], 
     cex = 2,type = 'l' ,
     xlab = "",ylab = "",
     main =bquote("log-likelihood of " ~ f(N)),
     cex.lab = 2,  
     cex.axis = 1.5, 
     cex.main = 2)

mtext(expression(delta), side = 1, line = 4, cex = 2.5)

points(as.numeric(delta.result$par),-as.numeric(delta.result$value), col = "black", pch = 19,cex = 2)

dev.off()
# select delta = 1 by using nb likelihood val plot of diag elt

diag(Sig_hat) <- 1

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
  if(sum(eigen(Sig_hat)$value < 0) > 0){
    return("Some eigen value is neg")
  }else{  
    l_total = 0
    # Sigma : -179408.2  Sgn * Sigma -358873.8
    start.time <- Sys.time()
    
    Us <- c()
    d <- dim(Sig_hat)[1]
    for(j in 1:M){
      # U <- diag(dlm::rwishart(delta.hat, p = nrow(Sig_hat),Sig_hat, SqrtSigma = diag(p)))
      Z <- rmvnorm(3, mean = rep(0,d), sigma = Sig_hat)
      U <- Z * Z
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


result <- cbind(seq(1,rep),tmp,tmp[,2] - tmp[,1])
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
writeLines(sigma_scenario, paste0(filename,"/Real_sigma_scenario.txt"))



write.csv(result, paste0(filename,'/Real_result.csv'), row.names = FALSE)


