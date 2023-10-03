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
    
    return(cox.out)
    
  } else if(inter.term %in% covariates) {
  	stop("Interaction terms have not been implemented for left-truncated data!")
  }
