## emiuga (2023-10-23): 
## added 'cox.params$' to 'covariates' (line 19) because error: "object 'covariates' not found"
## added 'getSnpSpikeCoxOut_intime' to call interval time Surv. functions.

getGenotypesCoxOut_intime <- function(inter.term, genotypes, cl, cox.params,
                               print.covs) {
  
  if(is.null(inter.term)){
    if(is.matrix(genotypes)){
      cox.out <- t(parApply(cl=cl,
                            X=genotypes, 
                            MARGIN=1, 
                            FUN=survFit_intime, 
                            cox.params=cox.params,
                            print.covs=print.covs))
    } else {
      cox.out <- survFit_intime(genotypes,
                         cox.params=cox.params,
                         print.covs=print.covs) 
    }
  } else if(inter.term %in% cox.params$covariates) {
    if(is.matrix(genotypes)){
      cox.out <- t(parApply(cl=cl,
                            X=genotypes,
                            MARGIN=1,
                            FUN=survFitInt_intime, 
                            cox.params=cox.params, 
                            cov.interaction=inter.term,
                            print.covs=print.covs))
    } else {
      cox.out <- survFitInt_intime(genotypes,
                            cox.params=cox.params,
                            cov.interaction=inter.term, 
                            print.covs=print.covs)
    }
  }
  
  return(cox.out)
  
}

getSnpSpikeCoxOut_intime <- function(inter.term, snp.spike, cox.params, print.covs){
  
  if(is.null(inter.term)){
    cox.out <- t(apply(snp.spike, 1, survFit_intime,
                       cox.params=cox.params,
                       print.covs=print.covs))
  } else {
    cox.out <- t(apply(snp.spike,
                       1,
                       survFitInt_intime,
                       cox.params=cox.params,
                       cov.interaction=inter.term, 
                       print.covs=print.covs) )
  }
  
}
