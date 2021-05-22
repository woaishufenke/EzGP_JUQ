###Example 4.1 code for the EZGP model
###This code uses a parallel computing for nr replications;
###this code use Low-storage BFGS optimization implemented by R package nloptr
###Each replication may take up to several minutes varied by computers; 
###To further reduce computing time, you can set "maxeval" = 50 in the following codes

# clear menory
rm(list=ls())

###R pacakge required
library('nloptr')
library(parallel)

#set number of quanititave factors
p = 3 
#set number of qualitative factors
q = 3
#set vector containing numbers of levels in qualitative factors
m=c(3,3,3) 
#set number of training data
n = 81
#set number of testing data
n2 = 1215
#set nugget if needed
tau = 0
#number of replications, (we use nr = 100 in the paper on a 26 core cluster)
#nr = 100
nr=4 #for a simple illustration

###read all training data for 100 replications
tradata = as.matrix(read.table("eg4_1train.txt", h =T))
# dim(tradata)
# head(tradata)

###read all testing data for 100 replications
testdata = as.matrix(read.table("eg4_1test.txt", h =T))
# head(testdata)
# dim(testdata)


### EzGP model function used in parallel computing for relications. 
parafunc <- function(ii, p, q, m, n, n2, tau) ##iith replication.
{
  ##total number of parameters in the model
  npar = 1 + q + p + p*sum(m)
  
  ## a help function
  psum <- function(x1,x2, par2)
  {
    return(sum(-par2*(x1-x2)^2))
  }
  
  # calculating covariance between two inputs w1 and w2
  covx <- function(w1,w2, parv){
    #variance parameter sigma^2
    par1 = parv[1:(q+1)]
    #correlation parameter in G0
    par2 = parv[(q+2):(q+1+p)]
    #correlation parameter in G1 to Gq
    par3 = parv[(q+2+p): npar]
    x1 = w1[1:p]
    z1 = w1[(p+1):(p+q)]
    x2 = w2[1:p]
    z2 = w2[(p+1):(p+q)]
    res1 = par1[1]*exp(psum(x1,x2,par2))
    for (i in 1:q){
      if(z1[i] != z2[i]){
        res1 = res1+0
      }
      else{
        l = z1[i]
        res1 = res1 + par1[i+1]*exp(psum(x1,x2, par3[(9*(i-1)+3*(l-1)+1) : (9*(i-1)+3*(l-1)+3)]))
      }
    }
    return(res1)
  }
  
  ### a modified version of covx where w.12 = (w1, w2)
  covx.m <- function(w.12, parv){
    return( covx(w1 = w.12[1:(p+q)], w2 = w.12[(p+q+1):(2*p+2*q)], parv) )
  }
  
  #expand grid to avoid for looping
  rcoord <- cbind(
    rep(seq_len(n - 1), times = rev(seq_len(n - 1))), 
    unlist(lapply(
      X = rev(seq_len(n - 1)), 
      FUN = function(nn, nm) seq_len(nn) + nm - nn, nm = n)))
  
  ### logliklihood function and the analytical gradients
  eval_f_list <- function(parv) {
    #check if all parameters are all positive
    if (min(parv) < 0)
    {
      return(list( "objective" = 10000000000,
            "gradient" = rep(NA, npar)))
    }
    #covariance matrix
    covm = matrix(0,n,n)
    # first compute the vector of elements in covariance matrix
    Rtemp <- apply(cbind(dat[rcoord[, 1], ], dat[rcoord[, 2], ]), 1, FUN = covx.m, parv=parv)
    covm[rcoord] <- Rtemp
    covm <- covm + t(covm)
    diag(covm) <- sum(parv[1:(q+1)]) + tau
    #response vector
    y = as.matrix(y)
    Tm = try(chol(covm), silent=TRUE)
    #round(t(Tm)%*%Tm,2) == round(covm,2)
    if ('try-error' %in% class(Tm)) {
      return(list( "objective" = 10000000000,
                   "gradient" = rep(NA, npar)))
    }
    m1 = as.matrix(c(rep(1,n)))
    invT = backsolve(Tm, diag(dim(Tm)[1]))
    invc = invT%*%t(invT)
    #liklihood function value
    MLE_result = 2*sum(log(diag(Tm))) + (t(y) %*% invc %*% y) - 1/sum(invc)*(t(m1) %*% invc %*% y)^2
    
    #calcaulte analytical gradients
    # trend estimator
    mu = as.numeric(1/sum(invc)*(t(m1) %*% invc %*% y))
    
    ### derivative for variance parameter sigma^2_0
    grad_var0 <- function()
    {
      ### derative function of sigma_0 for wi, wj 
      gradf_var <- function(w12)
      {
        x1 = w12[1:p]
        x2 = w12[(p+q+1):(2*p+q)]
        #correlation parameter in G0
        par2 = parv[(q+2):(q+1+p)]
        gx = exp(psum(x1,x2,par2))
        return(as.numeric(gx))
      }
      der_m = matrix(0,n,n)
      Rtemp_var1 <- apply(cbind(dat[rcoord[, 1], ], dat[rcoord[, 2], ]), 1, FUN = gradf_var)
      der_m[rcoord] <- Rtemp_var1
      der_m <- der_m + t(der_m)
      diag(der_m) = 1
      result = sum(diag(invc %*% der_m)) - t(y-mu) %*% invc %*% der_m %*% invc %*% (y-mu)
      return(result)
    }
    #grad_var0()
    
    ### derivative for variance parameter sigma^2_h h =1, ..., q
    grad_var <- function(h)
    {
      ### derative function of sigma_h for wi, wj 
      gradf_var <- function(w12,h)
      {
        x1 = w12[1:p]
        z1 = w12[(p+1):(p+q)]
        x2 = w12[(p+q+1):(p+q+p)]
        z2 = w12[(p+q+p+1):(2*p+2*q)]
        if(z1[h] != z2[h]){
          return(0)
        }
        else{
          l = as.numeric(z1[h])
          par3 = parv[(q+2+p): npar]
          gx = exp(psum(x1,x2, par3[(sum(m)*(h-1)+m[h]*(l-1)+1) : (sum(m)*(h-1)+m[h]*l)]))
          return(as.numeric(gx))
        }
      }
      der_m = matrix(0,n,n)
      Rtemp_var1 <- apply(cbind(dat[rcoord[, 1], ], dat[rcoord[, 2], ]), 1, FUN = gradf_var, h=h)
      der_m[rcoord] <- Rtemp_var1
      der_m <- der_m + t(der_m)
      diag(der_m) = 1
      result = sum(diag(invc %*% der_m)) - t(y-mu) %*% invc %*% der_m %*% invc %*% (y-mu)
      return(result)
    }
    
    ### derivative for correlation parameter theta_(0)_s for s=1,...,p
    grad_cor0 <- function(s)
    {
      ### derative function of theta_0_s for wi, wj 
      gradf_cor0 <- function(w12,s)
      {
        x1 = w12[1:p]
        #z1 = w12[(p+1):(p+q)]
        x2 = w12[(p+q+1):(p+q+p)]
        #z2 = w12[(p+q+p+1):(2*p+2*q)]
        #correlation parameter in G0
        par2 = parv[(q+2):(q+1+p)]
        gx = -parv[1] * (x1[s] - x2[s])^2 * exp(psum(x1,x2,par2))
        return(as.numeric(gx))
      }
      der_m = matrix(0,n,n)
      Rtemp_cor <- apply(cbind(dat[rcoord[, 1], ], dat[rcoord[, 2], ]), 1, FUN = gradf_cor0, s=s)
      der_m[rcoord] <- Rtemp_cor
      der_m <- der_m + t(der_m)
      diag(der_m) = 0
      result = sum(diag(invc %*% der_m)) - t(y-mu) %*% invc %*% der_m %*% invc %*% (y-mu)
      return(result)
    }
    
    ### derivative for correlation parameter theta_(h)_l_s for h=1,...,q, l=1,...,m_h, s=1,...,p
    grad_cor <- function(h,l,s)
    {
      ### derative function of theta_0_s for wi, wj 
      gradf_corhs <- function(w12,h,l,s)
      {
        x1 = w12[1:p]
        z1 = w12[(p+1):(p+q)]
        x2 = w12[(p+q+1):(p+q+p)]
        z2 = w12[(p+q+p+1):(2*p+2*q)]
        if((z1[h] != l) | (z2[h] != l)){
          return(0)
        }
        else{
          #variance parameter sigma^2
          par1 = parv[1:(q+1)]
          #correlation parameter in G1 to Gq
          par3 = parv[(q+2+p): npar]
          gx = -par1[h+1] * (x1[s] - x2[s])^2 * exp(psum(x1,x2, par3[(sum(m)*(h-1)+m[h]*(l-1)+1) : (sum(m)*(h-1)+m[h]*l)]))
          return(as.numeric(gx))
        }
      }
      der_m = matrix(0,n,n)
      Rtemp_cor <- apply(cbind(dat[rcoord[, 1], ], dat[rcoord[, 2], ]), 1, FUN = gradf_corhs, h=h, l=l, s=s)
      der_m[rcoord] <- Rtemp_cor
      der_m <- der_m + t(der_m)
      diag(der_m) = 0
      result = sum(diag(invc %*% der_m)) - t(y-mu) %*% invc %*% der_m %*% invc %*% (y-mu)
      return(result)
    }
    return( list( "objective" = MLE_result,
                  "gradient" = c(grad_var0(), grad_var(1), grad_var(2), grad_var(3),
                                 grad_cor0(1), grad_cor0(2), grad_cor0(3),
                                 grad_cor(1,1,1), grad_cor(1,1,2), grad_cor(1,1,3),
                                 grad_cor(1,2,1), grad_cor(1,2,2), grad_cor(1,2,3),
                                 grad_cor(1,3,1), grad_cor(1,3,2), grad_cor(1,3,3),
                                 grad_cor(2,1,1), grad_cor(2,1,2), grad_cor(2,1,3),
                                 grad_cor(2,2,1), grad_cor(2,2,2), grad_cor(2,2,3),
                                 grad_cor(2,3,1), grad_cor(2,3,2), grad_cor(2,3,3),
                                 grad_cor(3,1,1), grad_cor(3,1,2), grad_cor(3,1,3),
                                 grad_cor(3,2,1), grad_cor(3,2,2), grad_cor(3,2,3),
                                 grad_cor(3,3,1), grad_cor(3,3,2), grad_cor(3,3,3))                    ))
  }
  
  # calculate rmse
  rmse <- function(parv,xn,yn){
    #covariance matrix
    covm = matrix(0,n,n)
    # first compute the vector of elements in covariance matrix
    Rtemp <- apply(cbind(dat[rcoord[, 1], ], dat[rcoord[, 2], ]), 1, FUN = covx.m, parv=parv)
    covm[rcoord] <- Rtemp
    covm <- covm + t(covm)
    diag(covm) <- sum(parv[1:(q+1)]) + tau
    #response vector
    y = as.matrix(y)
    Tm = try(chol(covm), silent=TRUE)
    #round(t(Tm)%*%Tm,2) == round(covm,2)
    if ('try-error' %in% class(Tm)) {
      return(NULL)
    }
    m1 = as.matrix(c(rep(1,n)))
    invT = backsolve(Tm, diag(dim(Tm)[1]))
    invc = invT%*%t(invT)
    mu = as.numeric(1/sum(invc)*(t(m1) %*% invc %*% y))
    
    # prediction function
    prey <- function(wn){
      covv = matrix(0,n)
      for(i in 1:n){
        if (sum(round(wn,5)!=round(dat[i,],5)) > 0)
        {
          covv[i] = covx(wn,dat[i,],parv)
        } else {
          covv[i] = sum(parv[1:(q+1)]) + tau
        }
      }
      gamma = as.matrix(covv)
      result = mu  + t(gamma) %*% invc %*% (y - mu * m1) 
      return(result)
    }
    dif = c()
    nn = length(yn)
    for ( i in 1:nn){
      dif[i] = prey(xn[i,]) - yn[i]
    }
    result = sqrt(sum(dif^2)/nn)
    return(result)
  }
  
  # training data design part
  dat = tradata[c((n*(ii-1)+1):(n*ii)), c(1:(p+q))]
  ### training response vector
  y = tradata[c((n*(ii-1)+1):(n*ii)), (p+q+1)]
  # testing data design part
  ndat = testdata[c((n2*(ii-1)+1):(n2*ii)), c(1:(p+q))]
  ### testing response vector
  ny = testdata[c((n2*(ii-1)+1):(n2*ii)), (p+q+1)]
  
  #lower and upper bounds of control
  lb <- rep(0.1, npar)
  ub <- c(rep(100, (q+1)), rep(10, (npar-q-1)))
  #initial settings
  x0 = as.vector((lb+ub)/2)
  
  opts <- list("algorithm" = "NLOPT_LD_LBFGS",
                "xtol_rel" = 1.0e-5,
                "maxeval" = 100
                )
  
  res <- nloptr( x0=x0,
                 eval_f=eval_f_list,
                 lb=lb,
                 ub=ub,
                 opts=opts)
  ##mle of parameters
  solpar = as.vector(res$solution)
  ##gradients
  #grad = as.vector(eval_f_list(solpar)$gradient)
  ### calculating MSE
  rmseresult = rmse(solpar, ndat, ny)
  ### calculating NSE
  nse = 1 - n2*rmseresult^2/var(ny)/(n2-1)
  return(c(rmseresult,nse))
  #return(c(rmseresult,nse, solpar, grad))
}

# Calculate the number of cores
no_cores <- detectCores()
cl <- makeCluster(no_cores)
#cl <- makeCluster(26)

clusterExport(cl, "tradata") 
clusterExport(cl, "testdata") 
clusterEvalQ(cl, library('nloptr'))

# Initiate cluster
result = parLapply(cl, 1:nr, parafunc, p=p, q=q, m=m, n=n, n2=n2, tau=tau)
stopCluster(cl)

##number of parameters
#npar = 1 + q + p + p*sum(m)

rmser = c()
nser = c()
#solpar = matrix(0,nr,npar)
#grad = matrix(0,nr,npar)
for(i in 1:nr)
{
  rmser[i] = result[[i]][1]
  nser[i] = result[[i]][2]
  #solpar[i,] = result[[i]][3:(npar+2)]
  #grad[i,] = result[[i]][(npar+3):(npar+2+npar)]
}

###output results
write.table(as.vector(rmser), file = "eg4_1_ezgp_rmse.txt", sep="\t")
write.table(as.vector(nser), file = "eg4_1_ezgp_nse.txt", sep="\t")
