######### Helpers for Y~PRS Elastic Net and Super Learning Model with Sumstats ##########

# Combine and sort SNP weights for all PRS methods
sort_sumstats <- function(prs_methods,xty.snp,iterations,full=FALSE){
    weight.mat <- mclapply(1:iterations, function(ite) {
        weight.out <- c()
        for (prs_method in prs_methods) {
            if (full) {
                prs.method.weights <- as.data.frame(fread(paste0(full_weight_path,trait_name,".",prs_method,".txt"),header=F))
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
            }
        }
        return(weight.out)
    }, mc.cores=cores)
    return(weight.mat)
}


# calculate r2 for single PRS
single_ss_main <- function(X.ref, snp.weight, ite){
    # read statistics
    xty.dat <- fread(paste0(xty_path,trait_name,".xty.omnibus.ite",ite,".txt"))
    stat.dat <- fread(paste0(stats_path, trait_name, ".omnibus.forEVAL.txt"))
    set.seed(ite)
    ran.ord <- sample(1:1000,size=1000)
    ref.train <- ran.ord[1:round(nrow(X.ref)/3)]
    ref.tune <- ran.ord[(round(nrow(X.ref)/3)+1):round(nrow(X.ref)*2/3)]
    ref.final <- ran.ord[(round(nrow(X.ref)*2/3)+1):nrow(X.ref)]
    
    # get snp weights info
    snp.weight.sd <- apply(snp.weight,2,sd)
    snp.weight <- scale(snp.weight,center=F)
    
    # get standardized PRStY
    Y.hat <- X.ref[ref.train,] %*% snp.weight
    #Y.hat <- X.ref %*% snp.weight
    sd.yhat <- as.numeric(apply(Y.hat,2,sd))
    Y.hat.std <- scale(Y.hat,center=F,scale=sd.yhat)
    cov.yhat <- cov(Y.hat.std)
    prsty.train <- t(snp.weight) %*% xty.dat$test / (sqrt(stat.dat$var.Y) * stat.dat$N.t * sd.yhat)

    model.r2.train <- c()

    for (i in 1:ncol(snp.weight)){
        model.r2.tmp <- ss_r2(prsty=prsty.train[i,], prs.cov=cov.yhat[i,i])
        model.r2.train <- c(model.r2.train, model.r2.tmp)
    }
       
    # calculate testing r2
    Y.hat <- X.ref[ref.final,] %*% snp.weight
    sd.yhat <- as.numeric(apply(Y.hat,2,sd))
    Y.hat.std <- scale(Y.hat,center=F,scale=sd.yhat)
    cov.yhat <- cov(Y.hat.std)
    prsty.test <- t(snp.weight) %*% xty.dat$validation_test / (sqrt(stat.dat$var.Y) * stat.dat$N.vt * sd.yhat)

    model.r2.test <- c()

    for (i in 1:ncol(snp.weight)){
        model.r2.tmp <- ss_r2(prsty=prsty.test[i,], prs.cov=cov.yhat[i,i])
        model.r2.test <- c(model.r2.test, model.r2.tmp)
    }
    
    #write.table(data.frame(r2=model.r2.test),paste0(output_path,"/ite",ite,".single.r2.txt"),col.names = F, row.names = F,sep = "\t",quote = F)
    return(data.frame(r2=as.numeric(model.r2.train),beta.sd=snp.weight.sd))
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
            prev_alpha <- alpha
        }
        
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


# Calculate SNP weights from ensemble learning prs weights
ensemble_snp_weights <- function(ensemble.results,weight.col,prs.weight,xty.snp,ite){

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
