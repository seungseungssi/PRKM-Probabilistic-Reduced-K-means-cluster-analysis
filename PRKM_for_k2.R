
######3############
###############

print("Load Packages needed")
library(mvtnorm)
library(rstiefel)
library(matrixcalc)



PRKM_slow  <- function(data, k=4, q=2, max.iter = 1000, likelihood = T, show = F, init_A = "random", init_sigma = 1, init_eps = 0.1,seed = 1, print_partition = F,limit = 0.001, init_partition = c()){
    #print("#############################################")
    #print("############### Let's Start #################")
    #print("#############################################")
    # set the seed number 
    #addTaskCallback(function(...) {set.seed(seed);TRUE})
    likeli_list = list()
    n = nrow(data)
    j = ncol(data)
    group = 1:k # number of groups 
    
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
    else if (init_A == "const"){
        set.seed(seed)
        A = matrix(rep(1, j*q), ncol = q)
    }
    else if (init_A == "random"){
        set.seed(seed)
        A = matrix(sample( c(1:10),j*q, replace=T), ncol = q)/10
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
    for (g in 1:k){
        mu[g,] = t(A) %*% upper_mu[g,]
    }
    
    ## sigma : set 1 uniformly since usually the data is standardized and we have no information
    sigma = matrix(rep(init_sigma,k), ncol = 1)
    
    ## eps : set 0.1 as initial subspace transformation error
    eps = matrix(init_eps)
    
    # j : the proportion of each cluster (set it uniformly)
    p = matrix(rep(1/k, k), ncol = 1)
    like = 0
    iter = 0
    
    # iteration steps 
    for (iteration in 1:max.iter){
        
        ###############################
        # E-step : fill missing values 
        ###############################
        ##1. get posterior of Z r_k(x_i)
        # calculate posterior 
        pf = matrix(rep(0, k*n), ncol = k)
        for(g in group){
            pf[,g] = dmvnorm(data, mean = upper_mu[g,], sigma = sigma[g] * A %*% t(A) + diag(eps[1],j,j)) * p[g,]
        }
        
        # prevent posterior becoming too small. it will cause some errors 
        pf = pf + 0.1e-300
        # calculate marginal of pf (nominator)
        sum_pf  = apply(pf,1,sum)
        # calculate posterior probability 
        r = pf/sum_pf
        
        ##2. Fill in the data in lower dimension; Y
        e_array = array(dim = c(n,q,k))
        for (g in group){
            M = sigma[g] * A%*%t(A) + diag(eps[1],j,j)
            x_sub_Amu = sweep(data ,2, upper_mu[g,])
            e_array[,,g] = t(mu[g,] + sigma[g] * t(A) %*% solve(M  + diag(0.1e-30,j,j)) %*% t(x_sub_Amu))
        }
        
        V_array = array(dim = c(q,q,n,k))
        for (g in group){
            for (i in 1:n){
                M = sigma[g] * A %*% t(A) + diag(eps[1],j,j)
                V_array[,,i,g] = diag(sigma[g],q,q) - sigma[g]^2 * t(A) %*% solve(M + diag(0.1e-30,j,j))  %*% A + e_array[i,,g] %*% t(e_array[i,,g])
            }
        }
        
        ###############################
        # M-step : Estimate Parameters 
        ###############################
        
        # update the proportion of partition 
        p = as.matrix(apply(r, 2,sum))/n
        
        # update the mu
        mu = matrix(rep(0, q*k), nrow = k)
        for (i in 1:n){
            for (g in 1:k){
                mu[g,] = mu[g,] + r[i,g]*e_array[i,,g]
            }
        }
        
        for (g in group){
            mu[g,] = mu[g,]/ as.matrix(apply(r, 2,sum))[g,]
        }
        
        
        
       
        ## update the A 
        #1st part 
        first = matrix(rep(0, q^2), ncol = q)
        for (i in 1:n){
            for (g in group){
                first = first + r[i,g] * V_array[,,i,g]
            }
        }
        
        second = matrix(rep(0,j*q), ncol = q)
        for (i in 1:n){
            for (g in group){
                second = second + r[i,g] * data[i,] %*% t(e_array[i,,g])
            }
        }
        
        A = second %*% solve(first+diag(0.1e-30, q,q))
        
      
        ## Update eps
        
        eps = matrix(0)
        for (g in 1:k){
            for (i in 1:n){
                eps = eps + r[i,g] * (t(data[i,]) %*% data[i,] - 2 * t(data[i,]) %*% A %*% e_array[i,,g] + tr(V_array[,,i,g] %*% t(A) %*% A))/(j*n)
            }
        }
        
        ###############################
        #Calculate the log-likelihood##
        ###############################
        
        new_like = 0
        data = as.matrix(data)
        temp = matrix(rep(0, n*k), ncol = k)
        
        #update upper_mu
        upper_mu = t(A %*% t(mu))
        
        
        for (g in 1:k){
            temp[,g] = dmvnorm(data, mean = A %*% mu[g,], sigma = sigma[g] * A %*% t(A) + diag(eps[1],j,j)) * p[g,]
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
        
        # calculate mu in upper dimention 
        if (print_partition == T){
            print(apply(r, 1, which.max))
        }
        likeli_list = cbind(likeli_list,new_like)
        iter = iter + 1
    }
    # removing the seed number 

    Partition = apply(r, 1, which.max)
    return(list("mu" = mu, "sigma" = sigma,"A" = A, "p" = p,"epsilon" = eps, "Partition"=Partition,"Posterior" = r, "Likelihood" =new_like, "Likelihood_list"=likeli_list, "e" = e_array, "iter" = iter))
}

