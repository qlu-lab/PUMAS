######### Helpers for Y~PRS Elastic Net and Super Learning Model with Sumstats ##########

# Combine and sort SNP weights for all PRS methods
sort_sumstats <- function(prs_methods,xty.snp,iterations,full=FALSE) {
    weight.mat <- mclapply(1:iterations, function(ite) {
        weight.out <- c()
        for (prs_method in prs_methods) {
            if (full) {
                prs.method.weights <- as.data.frame(fread(paste0(full_weight_path,trait_name,".",prs_method,".ite1.txt"),header=F))
            } else {
                prs.method.weights <- as.data.frame(fread(paste0(weight_path,trait_name,".",prs_method,".ite",ite,".txt"),header=F))
            }
            for (index in 3:ncol(prs.method.weights)) {
                weight.out.tmp <- rep(0,nrow(xty.snp))
                snp.tmp <- xty.snp$A1
                match_indices <- match(prs.method.weights$V1, xty.snp$SNP)
                non_na_indices <- which(!is.na(match_indices))
                weight.out.tmp[match_indices[non_na_indices]] <- prs.method.weights[non_na_indices,index]
                snp.tmp[match_indices[non_na_indices]] <- prs.method.weights[non_na_indices,2]
                weight.out.tmp <- ifelse(snp.tmp==xty.snp$A1,weight.out.tmp,-weight.out.tmp)
                weight.out <- cbind(weight.out,weight.out.tmp)
                colnames(weight.out)[ncol(weight.out)] <- paste0("model", index - 2, ".", prs_method)
            }
        }
        return(weight.out)
    }, mc.cores=cores)
    
    return(weight.mat)
}

# Calculate ridge regression weights
ss_ridge <- function(prsty, prs.cov, lambda) {
    beta.ridge <- tryCatch({
        solve(prs.cov + lambda*diag(nrow(prs.cov))) %*% prsty
    }, error = function(e) {
        message("Error for lambda = ", lambda, ": ", e$message)
        matrix(0, nrow = nrow(prsty), ncol = ncol(prsty))
    })

    return(as.numeric(beta.ridge))
}

# Calculate EN regression weights
EN_ss_weight <- function(prsty, prs.cov, param, thr) {

    # iterate through lambda and alpha
    EN_weight <- c()
    prev_alpha <- 0
    EN_weight_tmp <- ss_ridge(prsty=prsty, prs.cov=prs.cov, lambda=as.numeric(param[1,1]))
    for (i in 1:nrow(param)) {
        alpha <- as.numeric(param[i,2])
        # warm start for lambda
        if (alpha == prev_alpha) {
            beta <- EN_weight_tmp
        } else {
            beta <- rep(0,nrow(prs.cov))
        }
        EN_weight_tmp <- wrap_coord_des(prsty, prs.cov, nrow(prs.cov), alpha, as.numeric(param[i,1]), beta, thr)
        prev_alpha <- alpha
        
        EN_weight <- rbind(EN_weight,EN_weight_tmp)
    }
    
    return(EN_weight)
}

# Calculate EN PRS r2
EN_ss_r2 <- function(prsty, prs.cov, prs.weight) {
    
    cov.Y_Y.hat <- sum(prs.weight * prsty)
    var.Y.hat <- matrix(prs.weight,nrow=1) %*% prs.cov %*% matrix(prs.weight,ncol=1)
    sum.R2 <- cov.Y_Y.hat^2/(as.numeric(var.Y.hat))
    
    # return values
    return(sum.R2)
}

# Calculate single PRS R2
ss_r2 <- function(prsty, prs.cov) {
    
    cov.Y_Y.hat <- prsty
    var.Y.hat <- prs.cov
    sum.R2 <- cov.Y_Y.hat^2/(as.numeric(var.Y.hat))
    
    # return values
    return(sum.R2)
}

# Calculate SNP weights from ensemble learning prs weights
ensemble_snp_weights <- function(ensemble.results,weight.col,prs.weight,xty.snp,ite) {

    ensemble.weight <- c()
    for (i in 1:ite){
        ensemble.weight <- cbind(ensemble.weight, t(ensemble.results[[i]][[weight.col]]))
    }
    avg.ensemble.weight <- rowMeans(ensemble.weight)
    ensemble.snp.weight <- prs.weight %*% avg.ensemble.weight
    ensemble.snp.weight <- cbind(xty.snp$SNP, xty.snp$A1, ensemble.snp.weight)
    colnames(ensemble.snp.weight) <- c("SNP", "A1", "Ensemble_Weight")
    
    return(ensemble.snp.weight)
}

# Calculate SNP weights from ensemble learning prs weights
prs_model_weights <- function(ensemble.results,weight.col,ite) {

    ensemble.weight <- c()
    for (i in 1:ite){
        ensemble.weight <- cbind(ensemble.weight, t(ensemble.results[[i]][[weight.col]]))
    }
    avg.ensemble.weight <- rowMeans(ensemble.weight)

    prs.model.weights <- data.frame(model = colnames(ensemble.results[[1]][[weight.col]]), weight = avg.ensemble.weight)
    
    return(prs.model.weights)
}

# create the start time for a function and print start time
start_time <- function(name) {
    begin.time = Sys.time()
    cat("\n", name, " begins at ", format(begin.time, format = "%F %R %Z"), ":", sep = "", "\n\n")

    return(begin.time)
}

# get the end time for a function and print how long the function took
end_time <- function(name, begin) {
    end.time <- Sys.time()
    total.time <- difftime(time1=end.time,time2=begin,units="sec")
    hrs <- floor(floor(total.time)/3600)
    mins <- floor(floor(floor(total.time) - hrs * 3600)/60)
    secs <- total.time - hrs*3600 - mins*60
    cat("\n", name, " ended at ", format(end.time, format = "%F %R %Z"), ". \n  Total time: ", hrs, " hours, ", mins, " minutes and ", secs, " seconds.\n\n", sep = "")
}
