######### Implement Y~PRS Elastic Net Super Learning Model with Sumstats ##########
EN_sl_main <- function(X.ref, param, weight, exclude, thr){

    stat.dat <- fread(paste0(stats_path, trait_name, ".omnibus.forEVAL.txt"))
    
    sl_results <- mclapply(1:k, function(ite) {
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
        
        # calculate EN regression weights
        model.weights <- list()
        for (ite2 in 1:k){
            xty.dat.inner <- fread(paste0(stats_path,trait_name,".innerloop.ite",ite,".xty.ite",ite2,".txt"))
            prsty.train <- t(snp.weight) %*% xty.dat.inner$test / (sqrt(stat.dat$var.Y) * 0.5 * stat.dat$N.t * sd.yhat)

            ensemble.weights.tmp <- EN_ss_weight(prsty=prsty.train, prs.cov=cov.yhat, param=param, thr=thr)
            model.weights[[ite2]] <- ensemble.weights.tmp
        }

        # calculate EN r2 for fine-tuning lambda and alpha
        model.r2 <- rep(0,nrow(param))
        Y.hat <- X.ref[ref.tune,] %*% snp.weight
        sd.yhat <- as.numeric(apply(Y.hat,2,sd))
        Y.hat.std <- scale(Y.hat)
        cov.yhat <- cov(Y.hat.std)
        for (ite2 in 1:k){
            model.r2.tmp <- c()
            xty.dat.inner <- fread(paste0(stats_path,trait_name,".innerloop.ite",ite,".xty.ite",ite2,".txt"))
            prsty.vali <- t(snp.weight) %*% (xty.dat$test - xty.dat.inner$test) / (sqrt(stat.dat$var.Y) * 0.5 * stat.dat$N.t * sd.yhat)
            
            for (i in 1:nrow(model.weights[[ite2]])){
                model.r2.tmp2 <- EN_ss_r2(prsty=prsty.vali, prs.cov=cov.yhat, prs.weight=as.numeric(model.weights[[ite2]][i,]))
                model.r2.tmp <- c(model.r2.tmp, model.r2.tmp2)
            }
            model.r2 <- model.r2 + model.r2.tmp/k
        }
        best.ridge <- matrix(param[1:21,][which.max(model.r2[1:21]),],1,2)
        best.lasso <- matrix(param[85:105,][which.max(model.r2[85:105]),],1,2)
        
        # ensemble weights calculation (final)
        Y.hat <- X.ref[c(ref.train,ref.tune),] %*% snp.weight
        sd.yhat <- as.numeric(apply(Y.hat,2,sd))
        Y.hat.std <- scale(Y.hat)
        cov.yhat <- cov(Y.hat.std)
        prsty.EN <- t(snp.weight) %*% xty.dat$test / (sqrt(stat.dat$var.Y) * stat.dat$N.t * sd.yhat)
        
        weights.ridge <- EN_ss_weight(prsty=prsty.EN, prs.cov=cov.yhat, param=best.ridge, thr=thr)
        weights.lasso <- EN_ss_weight(prsty=prsty.EN, prs.cov=cov.yhat, param=best.lasso, thr=thr)
        weights.lm <- ss_lm(prsty=prsty.EN, prs.cov=cov.yhat)
        
        # calculate sl weights
        sl.input.weights <- rbind(weights.ridge,weights.lasso,weights.lm)
        Y.hat <- X.ref[ref.final,] %*% snp.weight %*% t(sl.input.weights)
        sd.yhat <- as.numeric(apply(Y.hat,2,sd))
        Y.hat.std <- scale(Y.hat)
        cov.yhat <- cov(Y.hat.std)
        prsty.sl <- sl.input.weights %*% t(snp.weight) %*% xty.dat$validation_train / (sqrt(stat.dat$var.Y) * stat.dat$N.vtr * sd.yhat)
        
        sl.weights <- ss_lm(prsty=prsty.sl, prs.cov=cov.yhat)
        final.weights <- t(sl.weights) %*% sl.input.weights

        # ensemble R2 calculation
        Y.hat <- X.ref[ref.final,] %*% snp.weight
        sd.yhat <- as.numeric(apply(Y.hat,2,sd))
        Y.hat.std <- scale(Y.hat)
        cov.yhat <- cov(Y.hat.std)
        prsty.test <- t(snp.weight) %*% xty.dat$validation_test / (sqrt(stat.dat$var.Y) * stat.dat$N.vt * sd.yhat)
        
        sl.model.r2 <- EN_ss_r2(prsty=prsty.test, prs.cov=cov.yhat, prs.weight=as.numeric(final.weights))
        
        return(list(weights=final.weights,r2=sl.model.r2))
    }, mc.cores=cores)

    if (!is.null(full.snp.weight)) {
        sl.weights <- ensemble_snp_weights(ensemble.results=sl_results,weight.col="weights",prs.weight=full.snp.weight,xty.snp=xty.tmp,ite=k)
        write.table(sl.weights,paste0(output_path,trait_name,".sl.weights.txt"),col.names = T, row.names = F,sep = "\t",quote = F)
    }
    sl.r2 <- c()
    for (i in 1:k) {
        sl.r2 <- rbind(sl.r2, sl_results[[i]]$r2)
    }
    write.table(sl.r2,paste0(output_path,trait_name,".sl.r2.txt"),col.names = F, row.names = F,sep = "\n",quote = F)
}

############# main function ##############
# k, ref_path, trait_name, prs_method, xty_path, stats_path, weight_path, output_path
EN_super_learning <- function() {
    # initialize tuning parameters
    lambda <- 10^seq(-4,0,0.2)
    alpha <- seq(0,1,by=0.25)
    param <- data.frame(lambda=rep(rev(lambda),n=length(alpha)),alpha=rep(alpha,each=length(lambda)))

    # calculate SS-ensemble weights
    EN_sl_main(X.ref=ref.geno, param=param, weight=snp.weights.mat, exclude=prs_exclude, thr=threshold)
}