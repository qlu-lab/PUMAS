######### Implement Y~PRS Elastic Net Model with Sumstats ##########
list.files()
source("./code/evaluation/PUMACUBS_helpers.R")
sourceCpp("./code/evaluation/CoordDescent.cpp")

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
        
        # lm R2 calculation
        lm.model.r2 <- EN_ss_r2(prsty=prsty.test, prs.cov=cov.yhat, prs.weight=as.numeric(lm.weights))
        write.table(lm.model.r2,paste0(output_path,"/ite",ite,".lm.r2.txt"),col.names = F, row.names = F,sep = "\n",quote = F)

        return(data.frame(alpha=best_param[,2],lambda=best_param[,1],r2=ensemble.model.r2.best))
    }, mc.cores=4)
    print(en_results)
}


# calculate r2 for single PRS
single_ss_main <- function(X.ref, param, snp.weight, ite){
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
    
    write.table(data.frame(r2=model.r2.test),paste0(output_path,"/ite",ite,".single.r2.txt"),col.names = F, row.names = F,sep = "\t",quote = F)
    return(data.frame(r2=as.numeric(model.r2.train),beta.sd=snp.weight.sd))
}


############# main function ##############
# k, ref_path, trait_name, prs_method, xty_path, stats_path, weight_path, output_path
EN_ensemble_learning <- function() {

    # initialize tuning parameters
    lambda <- 10^seq(-4,0,0.2)
    alpha <- seq(0,1,by=0.25)
    param <- data.frame(lambda=rep(rev(lambda),n=length(alpha)),alpha=rep(alpha,each=length(lambda)))

    # read reference genotype
    cat("\n\nReading reference genotypes \n")
    begin.time = Sys.time()
    cat("\nReading reference genotypes ... begins at ", format(begin.time, format = "%F %R %Z"), ":", sep = "", "\n\n")
    ref.geno <- BEDMatrix(ref_path)
    rs.geno <- fread(paste0(ref_path,".bim"),header=F)$V2
    xty.tmp <- fread(paste0(xty_path,trait_name,".xty.omnibus.ite1.txt"),header=T)
    ref.geno <- as.matrix(ref.geno[,match(xty.tmp$SNP,rs.geno)])
    # report computing times for reading reference genotypes
    end.time <- Sys.time()
    total.time <- difftime(time1=end.time,time2=begin.time,units="sec")
    hrs <- floor(floor(total.time)/3600)
    mins <- floor(floor(floor(total.time) - hrs * 3600)/60)
    secs <- total.time - hrs*3600 - mins*60
    cat("\nReading reference genotypes ... ended at ", format(end.time, format = "%F %R %Z"), ". \n  Total time: ", hrs, " hours, ", mins, " minutes and ", secs, " seconds.\n\n", sep = "")

    # sort PRS weights and make weight matrix for pumas
    cat("\n\nSorting PRS weights ...\n")
    begin.time = Sys.time()
    cat("\nSorting PRS weights ... begins at ", format(begin.time, format = "%F %R %Z"), ":", sep = "", "\n\n")
    snp.weights.mat <- sort_sumstats(xty.snp=xty.tmp)
    # report computing times for sorting PRS weights
    end.time <- Sys.time()
    total.time <- difftime(time1=end.time,time2=begin.time,units="sec")
    hrs <- floor(floor(total.time)/3600)
    mins <- floor(floor(floor(total.time) - hrs * 3600)/60)
    secs <- total.time - hrs*3600 - mins*60
    cat("\nSorting PRS weights ended at ", format(end.time, format = "%F %R %Z"), ". \n  Total time: ", hrs, " hours, ", mins, " minutes and ", secs, " seconds.\n\n", sep = "")

    # calculate single PRS r2 and remove some PRS methods
    cat("\n\nCalculate single PRS and remove PRS methods ... \n")
    begin.time = Sys.time()
    cat("\nsingle_ss_main ... begins at ", format(begin.time, format = "%F %R %Z"), ":", sep = "", "\n\n")
    prs_exclude_list <- mclapply(1:k, function(ite) {
        tmp <- single_ss_main(X.ref=ref.geno, param=param, snp.weight=snp.weights.mat[[ite]], ite=ite)
        tmp.ratio <- (tmp[,2])/(tmp[,1])
        return(which(tmp.ratio >= (median(tmp.ratio) + 1.5*IQR(tmp.ratio))))
    }, mc.cores=4)
    prs_exclude <- unique(unlist(prs_exclude_list))
    # report computing times for single PRS calcs
    end.time <- Sys.time()
    total.time <- difftime(time1=end.time,time2=begin.time,units="sec")
    hrs <- floor(floor(total.time)/3600)
    mins <- floor(floor(floor(total.time) - hrs * 3600)/60)
    secs <- total.time - hrs*3600 - mins*60
    cat("\n  single_ss_main ended at ", format(end.time, format = "%F %R %Z"), ". \n  Total time: ", hrs, " hours, ", mins, " minutes and ", secs, " seconds.\n\n", sep = "")

    write.table(prs_exclude,paste0(output_path,"/prs_toExclude.txt"),col.names = F, row.names = F,sep = "\n",quote = F)

    # calculate SS-ensemble weights
    cat("\n\nEnsemble learning to calculate weights ... \n")
    begin.time = Sys.time()
    cat("\n Ensemble learning ... begins at ", format(begin.time, format = "%F %R %Z"), ":", sep = "", "\n\n")
    EN_ss_main(X.ref=ref.geno, param=param, weight=snp.weights.mat, exclude=prs_exclude, thr=1e-3)
    # report computing times for Ensemble learning
    end.time <- Sys.time()
    total.time <- difftime(time1=end.time,time2=begin.time,units="sec")
    hrs <- floor(floor(total.time)/3600)
    mins <- floor(floor(floor(total.time) - hrs * 3600)/60)
    secs <- total.time - hrs*3600 - mins*60
    cat("\n  Ensemble learning ended at ", format(end.time, format = "%F %R %Z"), ". \n  Total time: ", hrs, " hours, ", mins, " minutes and ", secs, " seconds.\n\n", sep = "")
}