library(pracma)

make_simul_data <- function(A_sim=factor_loading, c= 3, n = 25, sigma = 1, PSR = 0.15, noise = 6, seed = 1){
    #controls the overlapping
    #c=3; n = 25; sigma = 1; PSR = 0.15;A_sim=factor_loading
    
    F_sim= t(matrix(c(c,c,-c,-c,c,-c,-c,c), byrow = F, nrow = 2)) #C X Q centroid matrix 
    A_null_sim  = nullspace(t(A_sim)) # J X (J-Q) = 8  X (8-2)
    
    
    #making the U matix 
    U_sim = matrix(rep(0, n*4*4), ncol = 4) # I X C
    U_sim[1:n,1] = 1;U_sim[(n+1):(2*n),2] = 1;U_sim[(2*n+1):(3*n),3] = 1;U_sim[(3*n+1):(4*n),4] = 1;
    
    #subspace residual
    set.seed(seed)
    E_sim =matrix(rnorm(n*4*2, mean = 0, sd = sqrt(sigma)), ncol = 2)
    
    #complement residual 
    sigma_ortho = sigma/PSR-sigma
    set.seed(seed+1)
    E_ortho_sim = matrix(rnorm(n*4*noise, mean = 0, sd = sqrt(sigma_ortho)), ncol = noise)
    
    X_sim = U_sim%*%F_sim%*%t(A_sim) + E_sim%*%t(A_sim) + E_ortho_sim%*%t(A_null_sim)
    return(list("X" = X_sim, "U" = U_sim))
}
