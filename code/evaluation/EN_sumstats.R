######### Implement Y~PRS Elastic Net Model with Sumstats ##########
# Obtain weights for ensemble PRS
EN_ss_main <- function(X.ref, param, weight, exclude, thr){
    stat.dat <- fread(paste0(stats_path, trait_name, ".omnibus.forEVAL.txt"))
    
    en_results <- mclapply(1:k, function(ite) {
        # read statistics
        xty.dat <- fread(paste0(xty_path,trait_name,".xty.omnibus.ite",ite,".txt"))
        snp.weight <- apply(weight[[ite]][,-exclude],2,scale,center=F)
        
        # partition reference data to subset
        set.seed(ite)
        ran.ord <- sample(1:1000,size=1000)
        ref.train <- ran.ord[1:round(nrow(X.ref)/3)]
        ref.tune <- ran.ord[(round(nrow(X.ref)/3)+1):round(nrow(X.ref)*2/3)]
        ref.final <- ran.ord[(round(nrow(X.ref)*2/3)+1):nrow(X.ref)]
        
        # get standardized PRStY
        Y.hat <- X.ref[ref.train,] %*% snp.weight
        #Y.hat <- X.ref %*% snp.weight
        sd.yhat <- as.numeric(apply(Y.hat,2,sd))
        Y.hat.std <- scale(Y.hat)
        cov.yhat <- cov(Y.hat.std)
        prsty.train <- t(snp.weight) %*% xty.dat$test / (sqrt(stat.dat$var.Y) * stat.dat$N.t * sd.yhat)
        

        # calculate EN regression weights
        model.weights <- EN_ss_weight(prsty=prsty.train, prs.cov=cov.yhat, param=param, thr=thr)
        lm.weights <- ss_lm(prsty=prsty.train, prs.cov=cov.yhat)
        
        # calculate EN r2 for fine-tuning lambda and alpha
        model.r2.tmp <- c()
        Y.hat <- X.ref[ref.tune,] %*% snp.weight
        sd.yhat <- as.numeric(apply(Y.hat,2,sd))
        Y.hat.std <- scale(Y.hat)
        cov.yhat <- cov(Y.hat.std)
        prsty.vali <- t(snp.weight) %*% xty.dat$validation_train / (sqrt(stat.dat$var.Y) * stat.dat$N.vtr * sd.yhat)
        
        for (i in 1:nrow(model.weights)){
            model.r2.tmp2 <- EN_ss_r2(prsty=prsty.vali, prs.cov=cov.yhat, prs.weight=as.numeric(model.weights[i,]))
            model.r2.tmp <- c(model.r2.tmp, model.r2.tmp2)
        }
        best_param <- matrix(param[which.max(model.r2.tmp),],1,2)
        
        # ensemble weights calculation (final)
        Y.hat <- X.ref[c(ref.train,ref.tune),] %*% snp.weight
        sd.yhat <- as.numeric(apply(Y.hat,2,sd))
        Y.hat.std <- scale(Y.hat)
        cov.yhat <- cov(Y.hat.std)
        prsty.total <- t(snp.weight) %*% (xty.dat$validation_train + xty.dat$test) / (sqrt(stat.dat$var.Y) * (stat.dat$N.vtr + stat.dat$N.t) * sd.yhat)
        ensemble.weights.tmp <- EN_ss_weight(prsty=prsty.total, prs.cov=cov.yhat, param=best_param, thr=thr)
        
        # ensemble R2 calculation
        Y.hat <- X.ref[ref.final,] %*% snp.weight
        sd.yhat <- as.numeric(apply(Y.hat,2,sd))
        Y.hat.std <- scale(Y.hat)
        cov.yhat <- cov(Y.hat.std)
        prsty.test <- t(snp.weight) %*% xty.dat$validation_test / (sqrt(stat.dat$var.Y) * stat.dat$N.vt * sd.yhat)
        
        ensemble.model.r2.best <- EN_ss_r2(prsty=prsty.test, prs.cov=cov.yhat, prs.weight=as.numeric(ensemble.weights.tmp))
        
        # lm R2 calculation
        lm.model.r2 <- EN_ss_r2(prsty=prsty.test, prs.cov=cov.yhat, prs.weight=as.numeric(lm.weights))

        return(list(en.weights=ensemble.weights.tmp,lm.weights=t(lm.weights),en.r2=ensemble.model.r2.best,lm.r2=lm.model.r2))
    }, mc.cores=cores)

    if (!is.null(full.snp.weight)) {
        mclapply(1:2, function(i) {
            if (i == 1) {
                en.weights <- ensemble_snp_weights(ensemble.results=en_results,weight.col="en.weights",prs.weight=full.snp.weight,xty.snp=xty.tmp,ite=k)
                write.table(en.weights,paste0(output_path,trait_name,".en.weights.txt"),col.names = T, row.names = F,sep = "\t",quote = F)
            } else {
                lm.weights <- ensemble_snp_weights(ensemble.results=en_results,weight.col="lm.weights",prs.weight=full.snp.weight,xty.snp=xty.tmp,ite=k)
                write.table(lm.weights,paste0(output_path,trait_name,".lm.weights.txt"),col.names = T, row.names = F,sep = "\t",quote = F)
            }
        }, mc.cores=cores)
    }
    en.r2 <- c()
    for (i in 1:k) {
        en.r2 <- rbind(en.r2, en_results[[i]]$en.r2)
    }
    write.table(en.r2,paste0(output_path,trait_name,".en.r2.txt"),col.names = F, row.names = F,sep = "\n",quote = F)
    lm.r2 <- c()
    for (i in 1:k) {
        lm.r2 <- rbind(lm.r2, en_results[[i]]$lm.r2)
    }
    write.table(lm.r2,paste0(output_path,trait_name,".lm.r2.txt"),col.names = F, row.names = F,sep = "\n",quote = F)
}


############# main function ##############
# k, ref_path, trait_name, prs_method, xty_path, stats_path, weight_path, output_path
EN_ensemble_learning <- function() {
    # initialize tuning parameters
    lambda <- 10^seq(-4,0,0.2)
    alpha <- seq(0,1,by=0.25)
    param <- data.frame(lambda=rep(rev(lambda),n=length(alpha)),alpha=rep(alpha,each=length(lambda)))    
    # calculate SS-ensemble weights
    EN_ss_main(X.ref=ref.geno, param=param, weight=snp.weights.mat, exclude=prs_exclude, thr=threshold)
}