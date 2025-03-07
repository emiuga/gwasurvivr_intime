survFitInt_intime <- function(SNP,
                       cox.params,
                       cov.interaction,
                       print.covs) 
    {
    ## creating model matrix
    X <- cbind(INTER.TERM=cox.params$pheno.file[,cov.interaction]*SNP,
               SNP,
               cox.params$pheno.file)
    
    ## remove NA samples in genotype data
    X <- X[!is.na(SNP),]
    Y <- cox.params$Y[!is.na(SNP)]
    ROWNAMES <- cox.params$ROWNAMES[!is.na(SNP)]
    
    ## run fit with pre-defined parameters including INIT
    ## emiuga: Use aggre.fit for left-truncated data (type!="right") and method!='exact'
    type <- attr(Y, "type")
    if(type=="right") {
       stop("Data is not in left-truncated format")
    } else if(cox.params$METHOD=="breslow" || cox.params$METHOD =="efron") {
    fit <- agreg.fit(X,
                     Y,
                     cox.params$STRATA,
                     cox.params$OFFSET,
                     cox.params$INIT, 
                     cox.params$CONTROL,
                     cox.params$WEIGHTS,
                     cox.params$METHOD, 
                     ROWNAMES)
    } else {
       stop("Method is not 'breslow' or 'efron'.")
    }
    
    ## extract statistics
    if(print.covs=="only") {
        coef <- fit$coefficients[1]
        serr <- sqrt(diag(fit$var)[1])
        n.sample <- nrow(X)
        n.event <- sum(!grepl("[+]", as.character(Y)))
        res <- cbind(coef=coef, serr=serr, n.sample=n.sample, n.event=n.event)
        return(res)
        
    } else if(print.covs=="some"){
        coef <- fit$coefficients
        serr <- sqrt(diag(fit$var))
        res <- cbind(coef, serr)
        drop <- colnames(cox.params$pheno.file)[!colnames(cox.params$pheno.file)
                                                %in% cov.interaction]
        res <- res[!rownames(res) %in% drop,]
        res.names <- dimnames(res)
        n.sample <- nrow(X)
        n.event <- sum(!grepl("[+]", as.character(Y)))
        res <- c(res, n.sample, n.event)
        names(res) <- c(paste(toupper(res.names[[2]][1]),
                              res.names[[1]],
                              sep="_"),
                        paste(toupper(res.names[[2]][2]),
                              res.names[[1]], 
                              sep="_"),
                        "n.sample",
                        "n.event")
        return(res)
        
    } else if(print.covs=="all"){
        coef <- fit$coefficients
        serr <- sqrt(diag(fit$var))
        res <- cbind(coef, serr)
        res.names <- dimnames(res)
        n.sample <- nrow(X)
        n.event <- sum(!grepl("[+]", as.character(Y)))
        res <- c(res, n.sample, n.event)
        names(res) <- c(paste(toupper(res.names[[2]][1]),
                              res.names[[1]],
                              sep="_"),
                        paste(toupper(res.names[[2]][2]),
                              res.names[[1]], 
                              sep="_"),
                        "n.sample",
                        "n.event")
        return(res)
    }
}
