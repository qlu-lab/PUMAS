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
        write.table(data.frame(alpha=best_param[,2],lambda=best_param[,1],r2=ensemble.model.r2.best),paste0(output_path,"/ite",ite,".en.r2.txt"),col.names = F, row.names = F,sep = "\t",quote = F)
        write.table(ensemble.weights.tmp,paste0(output_path,"/ite",ite,".en.weights.txt"),col.names = F, row.names = F,sep = "\n",quote = F)
        
        # lm R2 calculation
        lm.model.r2 <- EN_ss_r2(prsty=prsty.test, prs.cov=cov.yhat, prs.weight=as.numeric(lm.weights))
        write.table(lm.model.r2,paste0(output_path,"/ite",ite,".lm.r2.txt"),col.names = F, row.names = F,sep = "\n",quote = F)
        write.table(lm.weights,paste0(output_path,"/ite",ite,".lm.weights.txt"),col.names = F, row.names = F,sep = "\n",quote = F)

        return(data.frame(alpha=best_param[,2],lambda=best_param[,1],r2=ensemble.model.r2.best))
    }, mc.cores=cores)
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