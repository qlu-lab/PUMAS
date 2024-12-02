######### Implement Y~PRS Elastic Net Model with Sumstats ##########
# Obtain weights for ensemble PRS
EN_ss_en <- function(X.ref, param, weight, thr){
    stat.dat <- fread(paste0(stats_path, trait_name, ".omnibus.forEVAL.txt"))
    
    en_results <- mclapply(1:k, function(ite) {
        # read statistics
        xty.dat <- fread(paste0(xty_path,trait_name,".xty.omnibus.ite",ite,".txt"))
        snp.weight <- apply(weight[[ite]],2,scale,center=F)
        
        # partition reference data to subset
        set.seed(ite)
        ran.ord <- sample(1:nrow(X.ref),size=nrow(X.ref))
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
        colnames(ensemble.weights.tmp) = colnames(snp.weight)

        # ensemble R2 calculation
        Y.hat <- X.ref[ref.final,] %*% snp.weight
        sd.yhat <- as.numeric(apply(Y.hat,2,sd))
        Y.hat.std <- scale(Y.hat)
        cov.yhat <- cov(Y.hat.std)
        prsty.test <- t(snp.weight) %*% xty.dat$validation_test / (sqrt(stat.dat$var.Y) * stat.dat$N.vt * sd.yhat)
        
        ensemble.model.r2.best <- EN_ss_r2(prsty=prsty.test, prs.cov=cov.yhat, prs.weight=as.numeric(ensemble.weights.tmp))

        # fwrite(as.data.frame(ensemble.weights.tmp),paste0(output_path,trait_name,".ite", ite, ".prs.weights.txt.gz"),col.names = T, row.names = F,sep = "\t",quote = F)
        return(list(en.weights=ensemble.weights.tmp, en.r2=ensemble.model.r2.best))
    }, mc.cores=cores)

    en.r2 <- c()
    for (i in 1:k) {
        en.r2 <- rbind(en.r2, en_results[[i]]$en.r2)
    }
    if (all(!is.na(en.r2))) {
        write.table(en.r2,paste0(output_path,trait_name,".en.r2.txt"),col.names = F, row.names = F,sep = "\n",quote = F)
    }

    if (!is.null(full.snp.weight)) {
        en.weights <- ensemble_snp_weights(ensemble.results=en_results,weight.col="en.weights",prs.weight=full.snp.weight,xty.snp=xty.tmp,ite=k)
        return(en.weights)
    } else {
        cat("\nNo full snp weights provided, outputting average PRS weight vector instead of SNP-level weights.")
        en.prs.model.weights <- prs_model_weights(ensemble.results=en_results,weight.col="en.weights",ite=k)
        write.table(en.prs.model.weights,paste0(output_path,trait_name,".en.prs.model.weights.txt"),col.names = T, row.names = F,sep = "\t",quote = F)
    }
}

EN_ss_SL <- function(X.ref, param, weight, thr){
    stat.dat <- fread(paste0(stats_path, trait_name, ".omnibus.forEVAL.txt"))

    sl_results <- mclapply(1:k, function(ite) {
        # read statistics
        xty.dat <- fread(paste0(xty_path,trait_name,".xty.omnibus.ite",ite,".txt"))
        snp.weight <- apply(weight[[ite]],2,scale,center=F)
        
        # partition reference data to subset
        set.seed(ite)
        ran.ord <- sample(1:nrow(X.ref),size=nrow(X.ref))
        ref.train <- ran.ord[1:round(nrow(X.ref)/3)]
        ref.tune <- ran.ord[(round(nrow(X.ref)/3)+1):round(nrow(X.ref)*2/3)]
        ref.final <- ran.ord[(round(nrow(X.ref)*2/3)+1):nrow(X.ref)]
        
        # get standardized PRStY
        Y.hat <- X.ref[ref.train,] %*% snp.weight
        sd.yhat <- as.numeric(apply(Y.hat,2,sd))
        Y.hat.std <- scale(Y.hat)
        cov.yhat <- cov(Y.hat.std)
        prsty.train <- t(snp.weight) %*% (xty.dat$test) / (sqrt(stat.dat$var.Y) * (stat.dat$N.t) * sd.yhat)
        
        # calculate EN regression weights
        model.weights <- list()
        for (ite2 in 1:k){
            xty.dat.inner <- fread(paste0(xty_path, trait_name, ".innerloop.ite",ite,".xty.ite",ite2,".txt"))
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
            xty.dat.inner <- fread(paste0(xty_path, trait_name, ".innerloop.ite",ite,".xty.ite",ite2,".txt"))
            prsty.vali <- t(snp.weight) %*% (xty.dat$test - xty.dat.inner$test) / (sqrt(stat.dat$var.Y) * 0.5 * stat.dat$N.t * sd.yhat)
            
            for (i in 1:nrow(model.weights[[ite2]])){
                model.r2.tmp2 <- EN_ss_r2(prsty=prsty.vali, prs.cov=cov.yhat, prs.weight=as.numeric(model.weights[[ite2]][i,]))
                model.r2.tmp <- c(model.r2.tmp, model.r2.tmp2)
            }
            model.r2 <- model.r2 + model.r2.tmp/k
        }
        best_param <- matrix(param[22:84,][which.max(model.r2.tmp[22:84]),],1,2)
        best.ridge <- matrix(param[1:21,][which.max(model.r2[1:21]),],1,2)
        best.lasso <- matrix(param[85:105,][which.max(model.r2[85:105]),],1,2)
        
        # ensemble weights calculation (train+tune)
        Y.hat <- X.ref[c(ref.train,ref.tune),] %*% snp.weight
        sd.yhat <- as.numeric(apply(Y.hat,2,sd))
        Y.hat.std <- scale(Y.hat)
        cov.yhat <- cov(Y.hat.std)
        prsty.total <- t(snp.weight) %*% xty.dat$test / (sqrt(stat.dat$var.Y) * stat.dat$N.t * sd.yhat)
        weights.en <- EN_ss_weight(prsty=prsty.total, prs.cov=cov.yhat, param=best_param, thr=thr)
        weights.ridge <- EN_ss_weight(prsty=prsty.total, prs.cov=cov.yhat, param=best.ridge, thr=thr)
        weights.lasso <- EN_ss_weight(prsty=prsty.total, prs.cov=cov.yhat, param=best.lasso, thr=thr)

        # calculate sl weights
        sl.input.weights <- rbind(weights.ridge,weights.lasso,weights.en)
        Y.hat <- X.ref[ref.final,] %*% snp.weight %*% t(sl.input.weights)
        sd.yhat <- as.numeric(apply(Y.hat,2,sd))
        Y.hat.std <- scale(Y.hat)
        cov.yhat <- cov(Y.hat.std)
        prsty.sl <- sl.input.weights %*% t(snp.weight) %*% xty.dat$validation_train / (sqrt(stat.dat$var.Y) * stat.dat$N.vtr * sd.yhat)
        
        super.weights <- as.numeric(EN_ss_weight(prsty=prsty.sl, prs.cov=cov.yhat, param=matrix(c(0,0),nrow=1), thr=thr))
        final.weights <- t(super.weights) %*% sl.input.weights
        colnames(final.weights) = colnames(snp.weight)

        # ensemble R2 calculation
        Y.hat <- X.ref[ref.final,] %*% snp.weight
        sd.yhat <- as.numeric(apply(Y.hat,2,sd))
        Y.hat.std <- scale(Y.hat)
        cov.yhat <- cov(Y.hat.std)
        prsty.test <- t(snp.weight) %*% xty.dat$validation_test / (sqrt(stat.dat$var.Y) * stat.dat$N.vt * sd.yhat)
        
        superlearn.model.r2.best <- EN_ss_r2(prsty=prsty.test, prs.cov=cov.yhat, prs.weight=as.numeric(final.weights))
        
        return(list(sl.weights=final.weights, sl.r2=superlearn.model.r2.best))
    }, mc.cores=cores)

    sl.r2 <- c()
    for (i in 1:k) {
        sl.r2 <- rbind(sl.r2, sl_results[[i]]$sl.r2)
    }
    if (all(!is.na(sl.r2))) {
        write.table(sl.r2,paste0(output_path,trait_name,".sl.r2.txt"),col.names = F, row.names = F,sep = "\n",quote = F)
    }

    if (!is.null(full.snp.weight)) {
        suplearn.weights <- ensemble_snp_weights(ensemble.results=sl_results,weight.col="sl.weights",prs.weight=full.snp.weight,xty.snp=xty.tmp,ite=k)
        return(suplearn.weights)
    } else {
        cat("\nNo full snp weights provided, outputting average PRS weight vector instead of SNP-level weights.")
        sl.prs.model.weights <- prs_model_weights(ensemble.results=sl_results,weight.col="sl.weights",ite=k)
        write.table(sl.prs.model.weights,paste0(output_path,trait_name,".sl.prs.model.weights.txt"),col.names = T, row.names = F,sep = "\t",quote = F)
    }
}

############# main function ##############
# k, ref_path, trait_name, prs_method, xty_path, stats_path, weight_path, output_path
EN_en_learning <- function() {
    # initialize tuning parameters
    lambda <- 10^seq(-4,0,0.2)
    alpha <- seq(0,1,by=0.25)
    param <- data.frame(lambda=rep(rev(lambda),n=length(alpha)),alpha=rep(alpha,each=length(lambda)))
    # calculate SS-ensemble weights
    EN_ss_en(X.ref=ref.geno, param=param, weight=snp.weights.mat, thr=1e-3)
}

EN_super_learning <- function() {
    # initialize tuning parameters
    lambda <- 10^seq(-4,0,0.2)
    alpha <- seq(0,1,by=0.25)
    param <- data.frame(lambda=rep(rev(lambda),n=length(alpha)),alpha=rep(alpha,each=length(lambda)))
    # calculate SS-ensemble weights
    EN_ss_SL(X.ref=ref.geno, param=param, weight=snp.weights.mat, thr=1e-3)
}