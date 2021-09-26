

print("Load Packages needed")
library(mvtnorm)
library(rstiefel)
library(matrixcalc)

PRKM  <- function(data, k=4, q=2, max.iter = 100, likelihood = F, show = F, init_A = "orthogonal",  init_eps = 0.1, ortho = F,seed = 1, print_partition = F,limit = 0.01, init_partition = c(), early_stop = F){
    #print("#############################################")
    #print("############### Let's Start #################")
    #print("#############################################")
    # set the seed number 
    likeli_list = list()
    n = nrow(data)
    j = ncol(data)
    # number of groups 
    
    group = 1:k
    iter = 0
    
    #get initial 
    ## A : we may use PCA initializer here 
    if (init_A =="pca"){
        pca = prcomp(data)
        A = pca$rotation[,1:q]
    }
    else if (init_A == "orthogonal"){
        set.seed(seed)
        A= t(rmf.matrix(matrix(rep(0,j*j), ncol = j)))[,1:q]
    }
    
    ## mu : just pick among t(A)X randomly
    #XA = data %*% A # this is same as XA = t(t(A) %*% t(data))
    if (length(init_partition) == 0){
        set.seed(seed)
        random = sample(1:n, k)
    }
    else {
        random = init_partition
    }
    upper_mu = data[random,]
    
    mu = matrix(rep(0, q*k), nrow = k)
    for (g in group){
        mu[g,] = t(A) %*% upper_mu[g,]
    }
    
    ## sigma : set 1 uniformly since usually the data is standardized and we have no information
    sigma = matrix(1)
    
    ## eps : set 0.1 as initial subspace transformation error
    eps = matrix(init_eps)
    
    # j : the proportion of each cluster (set it uniformly)
    p = matrix(rep(1/k, k), ncol = 1)
    like = 0
    
    # iteration steps 
    for (iteration in 1:max.iter){
        
        ###############################
        # E-step : fill missing values 
        ###############################
        ##1. get posterior of Z r_k(x_i)
        # calculate posterior 
        pf = matrix(rep(0, k*n), ncol = k)
        for(g in group){
            pf[,g] = dmvnorm(data, mean = upper_mu[g,], sigma = sigma[1] * A %*% t(A) + diag(eps[1],j,j)) * p[g,]
        }
        
        # prevent posterior becoming too small. it will cause some errors 
        pf = pf + 0.1e-300
        # calculate marginal of pf (nominator)
        sum_pf  = apply(pf,1,sum)
        # calculate posterior probability 
        r = pf/sum_pf
        
        ##2. Fill in the data in lower dimension; Y
        e_array = array(dim = c(n,q,k))
        M = sigma[1] * A%*%t(A) + diag(eps[1],j,j)
        for (g in group){
            x_sub_Amu = sweep(data ,2, upper_mu[g,])
            e_array[,,g] = t(mu[g,] + sigma[1] * t(A) %*% solve(M + diag(0.1e-30,j,j)) %*% t(x_sub_Amu))
        }
        
        V_array = array(dim = c(q,q,n,k))
        for (g in group){
            for (i in 1:n){
                V_array[,,i,g] = diag(sigma[1],q,q) - sigma[1]^2 * t(A) %*% solve(M + diag(0.1e-30,j,j))  %*% A + e_array[i,,g] %*% t(e_array[i,,g])
            }
        }
        
        ###############################
        # M-step : Estimate Parameters 
        ###############################
        n_k = as.matrix(apply(r, 2,sum))
        
        # update the proportion of partition 
        p = as.matrix(apply(r, 2,sum))/n
        # update the mu
        
        mu = matrix(rep(0, q*k), nrow = k)
        for (g in group){
            mu[g,] = apply(sweep(e_array[,,g], 1, r[,g], "*"), 2, sum)/n_k[g]
        }
        
        ## update the sigma 
        #sigma = matrix(0)
        #for (g in group){
        #    E_tr_V = apply(V_array[,,,g],3,tr)
        #    sigma[1] =  sigma[1] + (r[,g] %*% E_tr_V - 2 * sum(sweep(e_array[,,g],1,r[,g],"*") %*% mu[g,]) + n_k[g] * t(mu[g,]) %*% mu[g,]) 
        #}
        #sigma[1] = sigma[1] / (q * n)
        
        ## update the A 
        #1st part 
        first = matrix(rep(0, q^2), ncol = q)
        for (g in group){
            first = first + apply(sweep(V_array[,,,g], 3, r[,g], "*" ), 1:2, sum)
        }
        
        second = matrix(rep(0,j*q), ncol = q)
        for (g in group){
            second = second + t(apply(data,2, function(x){return(x*r[,g])})) %*% e_array[,,g]
        }
        
        A = second %*% solve(first+diag(0.1e-30, q,q))
        #A = varimax(A)$loading
        
        #t(A) %*% A orthogonality check 
        if (ortho == T){
            print("t(A)%*%A")
            print(t(A) %*% A)
        }
        
        ## Update eps
        eps = matrix(0) 
        for (g in group){
            eps[1] = eps[1] +  (tr(apply(data,2, function(x){return(x*r[,g])}) %*% t(data)) - 2 * tr(apply(data,2, function(x){return(x*r[,g])}) %*% A %*% t(e_array[,,g]) )
                                + r[,g] %*% apply(V_array[,,,g],3, function(x){return(tr(x %*% t(A) %*% A))}))/(j*n)
        }
        
        ###############################
        #Calculate the log-likelihood##
        ###############################
        
        new_like = 0
        data = as.matrix(data)
        temp = matrix(rep(0, n*k), ncol = k)
        
        # update upper_mu
        upper_mu = t(A %*% t(mu))
        
        for (g in group){
            set.seed(seed)
            temp[,g] = dmvnorm(data, mean = upper_mu[g,], sigma = sigma[1] * A %*% t(A) + diag(eps[1],j,j)) * p[g,]
        }
        
        new_like = sum(log(apply(temp, 1, sum)))
        if (likelihood == T){
            print(paste0(iteration,"-th iteration"))
            print(new_like)
        }
        
        # Terminate if the likelihood do not improve
        if (abs(like - new_like)<limit){
            break
        }
        
        
        # update the likelihood
        like = new_like
        
        
        if (show == T){
            XA = data %*% A
            plot(XA[,1],XA[,2],)
        }
        
        likeli_list = cbind(likeli_list,new_like)
        
        
        if (print_partition == T){
            set.seed(seed)
            Partition = apply(r, 1, function(x){return(which.max(rank(x, ties.method = "random")))})  
            print(Partition)
        }
        iter = iter + 1
        if (early_stop == T){
            if (length(unique(prkm$Partition)) != 4){
                break
            }
        }
        
        
        
        
    }
    Partition = apply(r, 1, which.max)
    return(list("mu" = mu, "sigma" = sigma,"A" = A, "p" = p,"epsilon" = eps, "Partition"=Partition,"Posterior" = r, "Likelihood" =new_like, "Likelihood_list"=likeli_list, "iter"= iter, "e" = e_array, "upper_mu" = upper_mu ))
}






