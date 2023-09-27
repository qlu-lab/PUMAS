######### Implement Y~PRS Elastic Net Model with Sumstats ##########

# Combine and sort SNP weights for all PRS methods
sort_sumstats <- function(xty.snp){
    tmp.files <- list.files(paste0(weight_path))
    tmp.models <- grep(x=tmp.files,pattern=paste0("^ite1\\.(.*)\\.txt"),value=T)
    tmp.params <- gsub(x=tmp.models,pattern=paste0("ite1(.*).txt"),replacement="\\1")
    
    weight.mat <- mclapply(1:k, function(ite) {
        weight.out <- c()
        for (tmp.param in tmp.params){
            weight.out.tmp <- rep(0,nrow(xty.snp))
            snp.tmp <- xty.snp$A1
            prs.weight <- as.data.frame(fread(paste0(weight_path,"/ite",ite,tmp.param,".txt"),header=F))
            weight.out.tmp[match(prs.weight$V1,xty.snp$SNP)] <- prs.weight$V3
            snp.tmp[match(prs.weight$V1,xty.snp$SNP)] <- prs.weight$V2
            weight.out.tmp <- ifelse(snp.tmp==xty.snp$A1,weight.out.tmp,-weight.out.tmp)
            weight.out <- cbind(weight.out,weight.out.tmp)
        }
        return(weight.out)
    }, mc.cores=4)

    return(weight.mat)
}


# Calculate linear regression weights
ss_lm <- function(prsty, prs.cov){
    beta.ridge <- solve(prs.cov) %*% prsty
    return(as.numeric(beta.ridge))
}


# Calculate ridge regression weights
ss_ridge <- function(prsty, prs.cov, lambda){
    beta.ridge <- solve(prs.cov + lambda*diag(nrow(prs.cov))) %*% prsty
    return(as.numeric(beta.ridge))
}


# Calculate EN regression weights
EN_ss_weight <- function(prsty, prs.cov, param, thr){

    # iterate through lambda and alpha
    EN_weight <- c()
    prev_alpha <- 0
    for (i in 1:nrow(param)) {
        alpha <- as.numeric(param[i,2])
        if (alpha == 0){
            EN_weight_tmp <- ss_ridge(prsty=prsty, prs.cov=prs.cov, lambda=as.numeric(param[i,1]))
        } else {
            # warm start for lambda
            if (alpha == prev_alpha) {
              beta <- EN_weight_tmp
            } else {
              beta <- rep(0,nrow(prs.cov))
            }
            EN_weight_tmp <- wrap_coord_des(prsty, prs.cov, nrow(prs.cov), alpha, as.numeric(param[i,1]), beta, thr)
            #EN_weight_tmp <- wrap_coord_des(prsty, prs.cov, alpha, as.numeric(param[i,1]), beta)
            prev_alpha <- alpha
        }
        #cat("\nnrow(prs.cov): ", nrow(prs.cov), "\n")
        #cat("\nalpha: ", alpha, "\n")
        #cat("\nEN_weight_tmp: ", EN_weight_tmp, "\n")
        
        EN_weight <- rbind(EN_weight,EN_weight_tmp)
    }
    
    return(EN_weight)
}


# Calculate EN PRS r2
EN_ss_r2 <- function(prsty, prs.cov, prs.weight){
    
    cov.Y_Y.hat <- sum(prs.weight * prsty)
    var.Y.hat <- matrix(prs.weight,nrow=1) %*% prs.cov %*% matrix(prs.weight,ncol=1)
    sum.R2 <- cov.Y_Y.hat^2/(as.numeric(var.Y.hat))
    
    # return values
    return(sum.R2)
}


# Calculate single PRS R2
ss_r2 <- function(prsty, prs.cov){
    
    cov.Y_Y.hat <- prsty
    var.Y.hat <- prs.cov
    sum.R2 <- cov.Y_Y.hat^2/(as.numeric(var.Y.hat))
    
    # return values
    return(sum.R2)
}
