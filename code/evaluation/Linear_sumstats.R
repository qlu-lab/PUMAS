## calculate sumstats-based R2 for a single PRS method
get_sumR2 <- function(XtY.t, weight, var.Y, N.t, X.ref){
    cat("\nXtY.t: ", dim(XtY.t))
    cat("\nweight: ", dim(weight))
    cat("\nN.t: ", dim(N.t))
  weight <- t((t(weight) - colMeans(weight)) / apply(weight, 2, sd)) # idk why we cannot use colSds (from matrixStats)
  cat("\nweight after transform : ", dim(weight))
  cov.Y_Y.Hat <- colSums(weight*(XtY.t/N.t))
  Y.Hat <- X.ref %*% weight
  var.Y.Hat <- apply(Y.Hat, 2, var)
  sum.R2 <- (unlist(cov.Y_Y.Hat))^2/(var.Y*unlist(var.Y.Hat))
  sum.R2[is.na(sum.R2)] <-0
  return(sum.R2)
}

## calculate sumstats-based omnibus R2
get_omnibus_sumR2 <- function(XtY.vt, weights, prs.weights, var.Y, N.vt, X.ref){

    weights.std <- apply(weights,2,function(s){return((s-mean(s))/sd(s))})
    Y.Hats <- X.ref %*% weights.std

    # calculate omnibus R2
    cov.Y_Y.Hat.vt <- t(prs.weights) %*% t(weights.std) %*% XtY.vt / N.vt
    var.Y.Hat <- var(as.numeric(Y.Hats %*% prs.weights))
    sum.R2 <- cov.Y_Y.Hat.vt^2/(var.Y*unlist(var.Y.Hat))
    
    # return values
    return(list(sum.R2=sum.R2))
}

get_omnibus_weights <- function(XtY.vtr, weights, N.vtr, X.ref){

    weights.std <- apply(weights,2,function(s){return((s-mean(s))/sd(s))})
    Y.Hats <- X.ref %*% weights.std
    
    # calculate PRS weights
    cov.Y_Y.Hat.vtr <- t(weights.std) %*% XtY.vtr
    Sigma.Y.hats <- N.vtr*cov(Y.Hats)
    prs.weights <- as.numeric(solve(Sigma.Y.hats) %*% cov.Y_Y.Hat.vtr)
    prs.weights[prs.weights<0] <- 0

    # return values
    return(list(prs.weights=prs.weights))
}

## main function for single method R2
single_prs <- function(prs_method){
    snp.w <- as.data.frame(fread(paste0(weight_path,trait_name,".",prs_method,".ite1.txt"),header=T))
    suffix <- colnames(snp.w)[-c(1:4)]
    xty <- as.data.frame(fread(paste0(xty_path,trait_name,".xty.omnibus.ite1.txt"),header=T))
    stats.pumas <- as.data.frame(fread(paste0(stats_path,trait_name,".omnibus.forEVAL.txt"),header=T))
    
    pumas.cor2 <- c()
    for (j in 1:k) {
        xty <- fread(paste0(xty_path,trait_name,".xty.omnibus.ite",j,".txt"),header=T)
        snp.w <- as.data.frame(fread(paste0(weight_path,trait_name,".",prs_method,".ite",j,".txt"),header=T))
        snp.weight <- snp.w[, suffix]
        snp.weight[snp.weight==Inf|snp.weight==-Inf] <- 0
        snp.weight[is.na(snp.weight)] <- 0
        pumas.cor2.tmp <- get_sumR2(XtY.t=xty$test, weight=snp.weight, var.Y=stats.pumas$var.Y, N.t=stats.pumas$N.t, X.ref=ref.geno)
        pumas.cor2 <- rbind(pumas.cor2,pumas.cor2.tmp)
    }
    colnames(pumas.cor2) <- suffix
    # store R2
    write.table(pumas.cor2,paste0(output_path,trait_name,".",prs_method,".tuning.txt"),col.names = T,row.names=F,quote=F,sep="\t")
}


single_prs_test <- function(prs_method){
    snp.w <- as.data.frame(fread(paste0(weight_path,trait_name,".",prs_method,".ite1.txt"),header=T))
    suffix <- colnames(snp.w)[-c(1:4)]
    xty <- as.data.frame(fread(paste0(xty_path,trait_name,".xty.omnibus.ite1.txt"),header=T))
    stats.pumas <- as.data.frame(fread(paste0(stats_path,trait_name,".omnibus.forEVAL.txt"),header=T))
    
    pumas.cor2 <- c()
    for (j in 1:k) {
        xty <- fread(paste0(xty_path,trait_name,".xty.omnibus.ite",j,".txt"),header=T)
        snp.w <- as.data.frame(fread(paste0(weight_path,trait_name,".",prs_method,".ite",j,".txt"),header=T))
        snp.weight <- snp.w[, suffix]
        snp.weight[snp.weight==Inf|snp.weight==-Inf] <- 0
        snp.weight[is.na(snp.weight)] <- 0
        pumas.cor2.tmp <- get_sumR2(XtY.t=xty$validation_test, weight=snp.weight, var.Y=stats.pumas$var.Y, N.t=stats.pumas$N.vt, X.ref=ref.geno)
        pumas.cor2 <- rbind(pumas.cor2,pumas.cor2.tmp)
    }
    colnames(pumas.cor2) <- suffix
    # store R2
    write.table(pumas.cor2,paste0(output_path,trait_name,".",prs_method,".testing.txt"),col.names = T,row.names=F,quote=F,sep="\t")
}


## main function for omnibus R2
omnibus_prs <- function(prs_methods){
    # first find the best tuning parameter for each PRS method
    best_param <- c()
    for (method in prs_methods){
        method_r2 <- as.matrix(fread(paste0(output_path,trait_name,".",method,".tuning.txt"),header=T))
        best_param <- c(best_param,  colnames(method_r2)[which.max(apply(method_r2,2,function(s){mean(s,na.rm=T)}))])
    }
    
    # coordinate genotype and stats
    xty <- as.data.frame(fread(paste0(xty_path,trait_name,".xty.omnibus.ite1.txt"),header=T))
    stats.pumas <- as.data.frame(fread(paste0(stats_path,trait_name,".omnibus.forEVAL.txt"),header=T))
    
    # calculate PRS weights
    pumas.weights <- c()
    for (j in 1:k) {
        xty <- fread(paste0(xty_path,trait_name,".xty.omnibus.ite",j,".txt"),header=T)
        
        # make a combined m*l matrix of m SNPs' weights for l methods
        snp.w <- as.data.frame(fread(paste0(weight_path,trait_name,".",prs_methods[1],".ite",j,".txt"),header=T,select=best_param[1]))
        for (i in 2:length(prs_methods)){
            snp.w <- cbind(snp.w, as.data.frame(fread(paste0(weight_path,trait_name,".",prs_methods[i],".ite",j,".txt"),header=T,select=best_param[i])))
        }
        snp.w[snp.w==Inf|snp.w==-Inf] <- 0
        snp.w[is.na(snp.w)] <- 0
        pumas.tmp <- get_omnibus_weights(XtY.vtr=xty$validation_train, weights=snp.w, N.vtr=stats.pumas$N.vtr, X.ref=ref.geno)
        pumas.weights <- rbind(pumas.weights,pumas.tmp$prs.weights)
    }
    pumas.weights.avg <- colMeans(pumas.weights)
    
    
    # calculate omnibus PRS r2
    pumas.cor2 <- c()
    for (j in 1:k){
        xty <- fread(paste0(xty_path,trait_name,".xty.omnibus.ite",j,".txt"),header=T)
        
        # make a combined m*l matrix of m SNPs' weights for l methods
        snp.w <- as.data.frame(fread(paste0(weight_path,trait_name,".",prs_methods[1],".ite",j,".txt"),header=T,select=best_param[1]))
        for (i in 2:length(prs_methods)){
            snp.w <- cbind(snp.w, as.data.frame(fread(paste0(weight_path,trait_name,".",prs_methods[i],".ite",j,".txt"),header=T,select=best_param[i])))
        }
        snp.w[snp.w==Inf|snp.w==-Inf] <- 0
        snp.w[is.na(snp.w)] <- 0
        pumas.tmp <- get_omnibus_sumR2(XtY.vt=xty$validation_test, weights=snp.w, prs.weights=pumas.weights.avg, var.Y=stats.pumas$var.Y, N.vt=stats.pumas$N.vt, X.ref=ref.geno)
        
        pumas.cor2 <- c(pumas.cor2,unlist(pumas.tmp$sum.R2))
    }
    
    pumas.cor2 <- as.data.frame(pumas.cor2)
    colnames(pumas.cor2) <- "Omnibus_PRS"
    colnames(pumas.weights) <- prs_methods
    
    # get r2 from single PRS method on validation_test
    
    # store results
    write.table(pumas.weights,paste0(output_path,trait_name,".omnibus.weights.txt"),col.names = T,row.names=F,quote=F,sep="\t")
    write.table(pumas.cor2,paste0(output_path,trait_name,".omnibus.r2.txt"),col.names = T,row.names=F,quote=F,sep="\t")
}


linear_ensemble_learning <- function() {
    # get single PRS method r2 separately, model tuning
    for (single_method in prs_method) {
        single_prs(single_method)
    }

    # get single PRS method r2 separately, testing
    for (single_method in prs_method) {
        single_prs_test(single_method)
    }

    # get omnibus r2
    omnibus_prs(prs_method)
} 
