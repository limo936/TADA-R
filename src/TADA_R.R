#########################################################################################
# This is a rewrite of the original TADA software (He et al 2015). It focus on the calculation of TADA test statistics, Pvalues, and Bayes Qvalues.
# Many of the previous extraneous functions are no longer supported in this software but can be found in the original version on 
# the old software page. 
# Version 2:
#     - Does combine the functionality of the TADA and TADAdenovo into one function.
#     - Allows for different mutation rates for different variants (ie denovo and missense)
#     - No longer supports gene specific denovo-only statistics, denovo-only (T or F) now applies to all 
#                genes for a mutation-variant.
# 
#
# To use this software:
#   1) set up the data in the required format
#   2) source("TADA_v2.R")
#   3) TADA(.....) to obtain the test statistics
#   4) Bayesian.FDR(....) to obtain the Qvalue for the test statistics
#   when pvalues are desired
#   5) TADAnull(....) to obtain the test statistics for the data set under the null hypothesis
#   6) bayesFactor.pvalue(....) to obtain the permuations based p-values
#
# Bert Klei
# Computational Genetics
# WPIC-UPMC
# October 2015
#########################################################################################

TADA <- function(tada.counts, sample.counts, mu, hyperpar, denovo.only, mu.frac=1, pi.gene=1) {
  # Genome-wide application of TADA for K classes of variants
  # tada.counts: list of K data frames in which each dataframe consists of vectors for counts for denovo, case and control mutation counts
  # sample.counts: list of K data frames in which each dataframe consists of a vector with three entries of total sample counts, 
  #                 one for denovo (# trios), cases (# cases + # trios) and controls (# controls + # trios)
  # mu: data frame with K vectors of mutation rates. One for each mutation category.
  # mu.frac: data frame with the fraction to use for each mutation category K
  # hyperpar: list of K data frames. Each data frame consists of entries for gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0
  # denovo.only: data frame with K Boolean variables indicating whether only denovo counts should be used (T) or whether both denovo and case-control counts be used (F).
  # pi.gene: data frame with K vectors of estimated fractions of causal variants, one for each class of variants. These fractions will be used to set gene-specific RR (case-control)
  # Output: dataframe with BF for each of the K classes of variants as well as BF.total. One entry for each of the genes.
  
  ### make sure every list and dataframe has the same elements
  print("CALCULATION OF TADA TEST STATISTICS")
  print("checking the input for consistent variable names"); flush.console()
  mutation.types=names(tada.counts)
  n.mutation=length(mutation.types)
  ### make sure mu.frac and pi.gene are data.frames
  if(!is.data.frame(mu.frac)){
    mu.frac=data.frame(matrix(mu.frac,1,n.mutation))
    names(mu.frac)=mutation.types
  }
  if(!is.data.frame(pi.gene)){
    pi.gene=data.frame(matrix(pi.gene,nrow(tada.counts[[1]]),n.mutation))
    colnames(pi.gene)=mutation.types
  }
  
  if(sum(mutation.types %in% names(mu)) != n.mutation | sum(mutation.types %in% names(mu.frac)) != n.mutation | 
     sum(mutation.types %in% names(hyperpar)) != n.mutation | sum(mutation.types %in% names(pi.gene)) != n.mutation |
     sum(mutation.types %in% names(denovo.only)) != n.mutation | sum(mutation.types %in% names(sample.counts)) != n.mutation){
    return("miss-match in names for the different variables")
  }
  names.trans.categories=names(sample.counts)
  names.N=c("dn","ca","cn","ca2","cn2","ar")
  for(mutation in mutation.types){
    if(sum(names.N %in% names(tada.counts[[mutation]])) != 6){
      return(paste("columns of ",mutation," do not match the required 'dn' 'ca' 'cn' 'ca2' 'cn2' 'ar'",sep=""))
    }
  }
  
  ### find the number of genes and the number of different kinds of mutations
  n.gene=nrow(mu)  ### was m
  n.mutation=length(mutation.types) ### was K
  
  ###
  BF=NULL
  for(mutation in mutation.types){
    print(paste("working on :: ",mutation)); flush.console()
    print(system.time(
      BF.mut<-unlist(lapply(1:nrow(tada.counts[[mutation]]),calculateBF,counts=tada.counts[[mutation]],n=sample.counts[[mutation]],
                            mu=mu[,mutation],mu.frac=mu.frac[,mutation],hyperpar=hyperpar[[mutation]],denovo.only=denovo.only[,mutation],
                            pi.gene=pi.gene[,mutation]))
    )); flush.console()
    BF=cbind(BF,BF.mut)
  }
  colnames(BF)=mutation.types
  rownames(BF)=rownames(tada.counts[[mutation]])
  
  ### calculate the overall BF
  BF.total=exp(rowSums(log(BF)))
  return (list(BF=BF, BF.total=BF.total))
}

TADAnull <- function(tada.counts, sample.counts, mu, hyperpar, denovo.only, mu.frac=1, nrep=100, dn.max=20, ca.max = 200, cn.max = 200, ar.max=5, max.gap=50){
  # Genome-wide application of TADA for K classes of variants
  # This function determines the distribution of the null hypothesis test statistics which in turn can be used to determine approximate p-values
  # tada.counts: list of K data frames in which each dataframe consists of vectors for counts for denovo, case and control mutation counts
  # sample.counts: list of K data frames in which each dataframe consists of a vector with three entries of total sample counts, 
  #                 one for denovo (# trios), cases (# cases + # trios) and controls (# controls + # trios)
  # mu: data frame with K vectors of mutation rates. One for each mutation category.
  # mu.frac: data frame with the fraction to use for each mutation category K
  # hyperpar: list of K data frames. Each data frame consists of entries for gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0
  # denovo.only: data frame with K Boolean variables of indicating whether only denovo counts should be used (T) or whether both denovo and case-control 
  # counts be used (F).
  # nrep: number of repetitions to use recommended to be at least 100. For smaller numbers of genes nrep should be increased.
  # dn.max: number of denovo events for which the BF is pre-computed and stored in a table. This speeds up the simulation process. The function will use dn.max or
  #         the maximum number of denovo events for a gene, whichever is smaller.
  # ca.max and cn.max: as dn.max this is used to pre-compute a table of common case-control count events to speed up processing. The larger these numbers, the longer the
  #                     the longer the pre computation step takes. For very large values of ca.max and cn.max an intergration error might occur.
  # max.gap: this is used internally to control the size of the pre-compution of the case-control BF matrix. It represents the gap between two genes in 
  #           case-control count when ordered from smallest to largest count. This essentially identifies outlying values that hardly ever happen and do not have to be
  #           pre-computed.
  # Output: dataframe with BFnull for each of the K classes of variants as well as BFnull.total. One entry for each of the genes times nrep.
  
  print("CALCULATION OF TADA TEST STATISTICS UNDER THE NULL HYPOTHESIS")
  ### determine the variants being analyzed
  mutations=names(hyperpar)
  ### make sure mu.frac is a data.frame
  if(!is.data.frame(mu.frac)){
    mu.frac=data.frame(matrix(mu.frac,1,length(mutations)))
    names(mu.frac)=mutations
  }
  
  # Pre-compute the bayes-factors for the denovo data
  table.BF.dn=list()
  for(mutation in mutations){
    print(paste("working on creating DN table for :: ",mutation)); flush.console()
    x=0:dn.max
    param=hyperpar[[mutation]]
    n=sample.counts[[mutation]]
    
    print(system.time(
      BF<-do.call(cbind.data.frame,lapply(1:length(x),table.BF.dn.wrapper,x=x,n.dn=n$dn,mu=mu[,mutation]*mu.frac[,mutation],
                                          gamma.dn=param$gamma.mean.dn,beta.dn=param$beta.dn))
    )); flush.console()
    colnames(BF)=paste("X",x,sep="")
    table.BF.dn[[mutation]]=BF
    
  }
  
  table.BF.cc=list()
  if(!denovo.only[mutation]){
    # Pre-compute the bayes-factors of the case-control data
    for (mutation in mutations) {
      print(paste("working on creating CC table for :: ",mutation)); flush.console()
      tada.counts[[mutation]]$Ncc=rowSums(tada.counts[[mutation]][,c("ca","cn")])
      Ncc=sort(tada.counts[[mutation]]$Ncc)
      Ncc.gaps=Ncc[-1]-Ncc[-length(Ncc)]
      i.gap=which(Ncc.gaps > max.gap)[1]
      n.ca=min(ca.max,Ncc[i.gap],na.rm=T)
      n.cn=min(cn.max,Ncc[i.gap],na.rm=T)
      x=expand.grid(ca=0:n.ca,cn=0:n.cn)
      param=hyperpar[[mutation]]
      n=sample.counts[[mutation]]
      dim(system.time(
        BF<-unlist(lapply(1:nrow(x),table.BF.cc.wrapper,x=x,n.cc=n[,c("ca","cn")],gamma.cc=param$gamma.mean.CC,beta.cc=param$beta.CC,
                          rho1=param$rho1,nu1=param$nu1,rho0=param$rho0,nu0=param$nu0))
      )); flush.console()
      table.BF.cc[[mutation]]=matrix(NA,max(x$ca)+1,max(x$cn)+1)
      table.BF.cc[[mutation]][cbind(x$ca+1,x$cn+1)]=BF
    }
  }
  
  table.BF.ar=list()
  if(!denovo.only[mutation]){
    # Pre-compute the bayes-factors of the case-control data
    for (mutation in mutations) {
      print(paste("working on creating ar table for :: ",mutation)); flush.console()
      x=0:ar.max
      param=hyperpar[[mutation]]
      n=sample.counts[[mutation]]
      dim(system.time(
        BF<-do.call(cbind.data.frame,lapply(1:length(x),table.BF.ar.wrapper,x=x,n.ar=n$ar,gamma.ar=param$gamma.mean.ar,beta.ar=param$beta.ar,
                                            rho1=param$rho1,nu1=param$nu1,rho0=param$rho0,nu0=param$nu0))
      )); flush.console()
      colnames(BF)=paste("X",x,sep="")
      table.BF.ar[[mutation]]=BF
    }
  }
  
  ### determine BF under the null distribution through permuations.s
  BF=NULL
  for(mutation in mutations){
    print(paste("working on creating null data for :: ",mutation)); flush.console()
    print(system.time(
      
      BF<-cbind(BF,unlist(lapply(1:nrow(tada.counts[[mutation]]),permute.gene,mu.rate=mu[,mutation]*mu.frac[,mutation], counts=tada.counts[[mutation]], 
                                 n=n, nrep=nrep, param=hyperpar[[mutation]], denovo.only=denovo.only[mutation],
                                 table.cc=table.BF.cc[[mutation]], table.dn=table.BF.dn[[mutation]], table.ar=table.BF.ar[[mutation]])))
    )); flush.console()
  }
  
  # calculate BF total
  colnames(BF)=mutations
  rownames(BF)=1:nrow(BF)
  BFtotal=exp(rowSums(log(BF)))
  
  return (list(BFnull=BF, BFnull.total=BFtotal))
  
}  

calculateBF <- function(i.gene,counts,n,mu,mu.frac,hyperpar,denovo.only,pi.gene){
  # wrapper so that lapply can be used to determine the BF for a gene and a particular mutation variant
  # i.gene: gene of interest
  # counts: counts for a particular variant, dataframe with vectors for dn, ca, cn
  # n: total samples counts, dataframe with entries for dn, ca, cn
  # mu: vector with mutation rates for the variant of interest for each gene
  # mu.frac: fraction to multiple mu with for the variant of interest
  # hyperpar: dataframe with entries for gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0 for the 
  #           variant of interst
  # denovo.only: Booelean vector indicating whether only denovo contribution (denovo.only = T) or a combination of 
  #              denovo and case-control contributions is to be used (denovo.only=F)
  # pi.gene: vector with K vectors of estimated fractions of causal variants, one for each class of variants. 
  #          These fractions will be used to set gene-specific RR (case-control)
  
  # if(i.gene%%100 == 0 | i.gene == nrow(counts)){
  #   pb <- txtProgressBar(min = 0, max = nrow(counts), style = 3)
  #   setTxtProgressBar(pb, i.gene)
  # }
  
  ### set the hyper parameters for this gene
  hyperpar.gene <- hyperpar
  RR.product <- hyperpar.gene$gamma.mean.CC*hyperpar.gene$beta.CC
  hyperpar.gene$gamma.mean.CC <- hyperpar.gene$gamma.mean.CC*pi.gene[i.gene] + (1-pi.gene[i.gene])
  hyperpar.gene$beta.CC <- RR.product/hyperpar.gene$gamma.mean.CC
  
  #determine the BAYES factor
  BF <- bayes.factor(x=counts[i.gene,], n=n, mu=mu[i.gene]*mu.frac, param=hyperpar.gene, denovo.only=denovo.only)
  
  return(BF=BF)
  
}

bayes.factor <- function(x,n,mu,param,denovo.only){
  # Bayes factor of the gene combining de novo and case-control
  # x: a list of (dn, ca, cn), counts in de novo, cases and controls
  # n: a list of (dn, ca, cn), sample sizes
  # param: (gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0)
  # denovo.only: Boolean vector indicating whether only denovo contribution (denovo.only = T) or a combination of 
  #              denovo and case-control contributions is to be used (denovo.only=F)
  # Prior distribution of RR in de novo: gamma.dist.dn ~ Gamma(gamma.mean.dn*beta.dn, beta.dn)
  # Prior distribution of RR in C/C data: gamma.dist.cc ~ Gamma(gamma.mean.CC*beta.CC, beta.CC)
  # Prior distribution of q|H1: Gamma(rho1, nu1)
  # Prior distribution of q|H0: Gamma(rho0, nu0)
  
  #contribution of denovo variants in families
  BF.dn=bayes.factor.dn(x.dn=x$dn,n.dn=n$dn,mu=mu,gamma.dn=param$gamma.mean.dn,beta.dn=param$beta.dn)
  if(denovo.only == F){
    #contribution of variants in cases and controls
    BF.cc=bayes.factor.cc(x.cc=x[,c("ca","cn")],n.cc=n[,c("ca","cn")],gamma.cc=param$gamma.mean.CC,beta.cc=param$beta.CC,
                          rho1=param$rho1,nu1=param$nu1,rho0=param$rho0,nu0=param$nu0)
    BF.ar=bayes.factor.ar(x.ar=x$ar,n.ar=n$ar,gamma.ar=param$gamma.mean.ar,beta.ar=param$beta.ar,
                          rho1=param$rho1,nu1=param$nu1,rho0=param$rho0,nu0=param$nu0)
  }else{
    BF.cc=1
  }
  # combine the pieces of information
  BF <- BF.dn * BF.cc * BF.ar
  return(BF=BF)
  
}

bayes.factor.dn <- function(x.dn,n.dn,mu,gamma.dn,beta.dn){
  # Bayes factor of de novo counts of a gene 
  # x.dn: the de novo count
  # n.dn: the sample size (number of families)
  # mu: the mutation rate (of this type of mutational events)
  # Prior distribution of RR: gamma ~ Gamma(gamma.dn*beta.dn, beta.dn)
  marg.lik0 <- dpois(x.dn, 2*n.dn*mu)
  marg.lik1 <- dnbinom(x.dn, gamma.dn*beta.dn, beta.dn/(beta.dn+2*n.dn*mu))
  BF <- marg.lik1/marg.lik0
  if(marg.lik0==0)
    BF=1
  return (BF=BF)
}

bayes.factor.cc <- function(x.cc,n.cc,gamma.cc,beta.cc,rho1,nu1,rho0,nu0){
  # Bayes factor of the case-control data
  # BF.cn and BF.ca: contribution from control and case data, respectively
  # Input: the count data x.cc, the sample size n.cc and the parameters gamma.cc, beta.cc, rho1 and nu1
  # Prior distribution of RR: gamma ~ Gamma(gamma.cc*beta.cc, beta.cc)
  # Prior distribution of q|H1: Gamma(rho1, nu1)
  # Prior distribution of q|H0: Gamma(rho0, nu0)
  marglik0.cc <- evidence.null.cc(x.cc, n.cc, rho0, nu0)
  marglik1.cc <- evidence.alt.cc(x.cc, n.cc, gamma.cc, beta.cc, rho1, nu1)
  BF.cn <- marglik1.cc$cn / marglik0.cc$cn
  BF.ca <- marglik1.cc$ca / marglik0.cc$ca
  if(marglik0.cc$cn==0)
    BF.cn=1
  if(marglik0.cc$ca==0)
    BF.ca=1
  BF <- BF.cn * BF.ca
  return(BF=BF)
}


evidence.null.cc <- function(x.cc, n.cc, rho0, nu0) {
  # model evidence of case-control data: P(x_1,x_0|H_0) 
  # Input: the count data x.cc, the sample size n.cc and the rho0 and nu0
  # Prior distribution of q|H0: Gamma(rho0, nu0)
  marglik0.ctrl.log <- log(dnbinom(x.cc$cn, rho0, nu0/(nu0+n.cc$cn)))
  marglik0.case.log <- log(dnbinom(x.cc$ca, rho0+x.cc$cn, (nu0+n.cc$cn)/(nu0+n.cc$cn+n.cc$ca)))
  marglik0.log <- marglik0.ctrl.log + marglik0.case.log
  
  return (list(cn=exp(marglik0.ctrl.log), ca=exp(marglik0.case.log), total=exp(marglik0.log)))
}

evidence.alt.cc <- function(x.cc, n.cc, gamma.cc, beta.cc, rho1, nu1, q.lower=1e-8, q.upper=0.1, debug=FALSE) {
  # model evidence of case-control data: P(x_1,x_0|H_1) 
  # Input: the count data x.cc, the sample size n.cc and the parameters gamma.cc, beta.cc, rho1 and nu1
  # Prior distribution of RR: gamma ~ Gamma(gamma.cc*beta.cc, beta.cc)
  # Prior distribution of q|H1: Gamma(rho1, nu1)
  integrand <- function(u) {
    q <- exp(u)
    return (dnbinom(x.cc$ca, gamma.cc*beta.cc, beta.cc/(beta.cc+n.cc$ca*q)) * dgamma(q, rho1+x.cc$cn, nu1+n.cc$cn) * exp(u))
  }
  
  marglik1.ctrl <- dnbinom(x.cc$cn, rho1, nu1/(nu1+n.cc$cn))
  marglik1.case <- integrate(integrand, lower=log(q.lower), upper=log(q.upper))$value
  
  marglik1 <- marglik1.ctrl * marglik1.case
  
  return (list(cn=marglik1.ctrl, ca=marglik1.case, total=marglik1))
}


bayes.factor.ar <- function(x.ar,n.ar,gamma.ar,beta.ar,rho1,nu1,rho0,nu0,q.lower=1e-8, q.upper=0.1){
  integrand.lik0 <- function(u) {
    q <- exp(u)
    return (dpois(x.ar, n.ar*q^2/4) * dgamma(q, rho0, nu0) * exp(u))
  }
  
  integrand.lik1 <- function(u) {
    q <- exp(u)
    return (dnbinom(x.ar, gamma.ar*beta.ar, beta.ar/(beta.ar+n.ar*q^2/4)) * dgamma(q, rho1, nu1) * exp(u))
  }
  
  marg.lik0 <- integrate(integrand.lik0, lower=log(q.lower), upper=log(q.upper))$value
  marg.lik1 <- integrate(integrand.lik1, lower=log(q.lower), upper=log(q.upper))$value
  BF <- marg.lik1/marg.lik0
  if(marg.lik0==0)
    BF=1
  return(BF=BF)
}



table.BF.dn.wrapper <- function(i.x,x,n.dn,mu,gamma.dn,beta.dn){
  # a wrapper function used for generating the denovo table, see bayes.factor.dn for a description of the variables
  BF=bayes.factor.dn(x[i.x],n.dn=n.dn,mu=mu,gamma.dn=gamma.dn,beta.dn=beta.dn)
  return(BF)
}


table.BF.cc.wrapper <- function(i.x,x,n.cc,gamma.cc,beta.cc,rho1,nu1,rho0,nu0){
  # a wrapper function used for generating the case-control table, see bayes.factor.cc for a description of the variables
  # if((i.x%%100 == 0 | i.x == nrow(x)) & nrow(x) > 1000){
  #   pb <- txtProgressBar(min = 0, max = nrow(x), style = 3)
  #   setTxtProgressBar(pb, i.x)
  # }
  BF=bayes.factor.cc(x[i.x,],n.cc=n.cc,gamma.cc=gamma.cc,beta.cc=beta.cc,rho1=rho1,nu1=nu1,rho0=rho0,nu0=nu0)
  return(BF)
}

table.BF.ar.wrapper <- function(i.x,x,n.ar,gamma.ar,beta.ar,rho1,nu1,rho0,nu0){
  # a wrapper function used for generating the case-control table, see bayes.factor.cc for a description of the variables
  # if((i.x%%100 == 0 | i.x == nrow(x)) & nrow(x) > 1000){
  #   pb <- txtProgressBar(min = 0, max = nrow(x), style = 3)
  #   setTxtProgressBar(pb, i.x)
  # }
  BF=bayes.factor.ar(x[i.x],n.ar=n.ar,gamma.ar=gamma.ar,beta.ar=beta.ar,rho1=rho1,nu1=nu1,rho0=rho0,nu0=nu0)
  return(BF)
}

permute.gene <- function(i.gene,mu.rate, counts, n, nrep, param, denovo.only, table.cc, table.dn, table.ar) {
  # Compute permutation BFs of one gene
  # mu.rate: the mutation rate of a gene for the variant of interest
  # counts: dn, ca, and co counts for the variant of interest to be permuted, this also has a column for Ncc=ca+cn
  # n: sample size, for values for de novo, case, control, and case+control
  # nrep: number of permutations
  # param: set of hyper parameters for the variant of interest. 
  # table.cc: table of precomputed BFs for case control events of size max.ca by max.cn for the variant of interest
  # table.dn: table of precomputed BFs for denovo events of size number of genes by max.dn for the variant of interest
  # Output: vector of nrep BF generated under the null hypothesis
  
  if(i.gene%%100 == 0 | i.gene == nrow(counts)){
    pb <- txtProgressBar(min = 0, max = nrow(counts), style = 3)
    setTxtProgressBar(pb, i.gene)
  }
  
  # generate permutation data for denovo events
  sample.dn = rpois(nrep, 2*n$dn*mu.rate[i.gene])
  # look up the BF value in the table.
  BFdn=table.dn[i.gene,sample.dn+1]
  
  
  if(!denovo.only){
    # when both denovo and case-control BF are needed
    #generate permuation data for case-control events
    max.ca=nrow(table.cc)
    max.cn=ncol(table.cc)
    sample.ca = rhyper(nrep, counts$Ncc[i.gene], n$ca+n$cn-counts$Ncc[i.gene], n$ca)
    sample.cn = counts$Ncc[i.gene] - sample.ca
    
    ## find the generated counts that are outside of the pre-computed table
    i.na=which(sample.ca + 1 > max.ca | sample.cn + 1 > max.cn)
    if(length(i.na > 0 )){
      # calculate their BF on a case by case basis
      BF.na=unlist(lapply(i.na,table.BF.cc.wrapper,x=data.frame(ca=sample.ca,cn=sample.cn),n.cc=n[,c("ca","cn")],gamma.cc=param$gamma.mean.CC,beta.cc=param$beta.CC,
                          rho1=param$rho1,nu1=param$nu1,rho0=param$rho0,nu0=param$nu0))
    }
    
    ## set the counts outside the range to missing
    sample.ca[ sample.ca+1 > max.ca]=NA
    sample.cn[ sample.cn+1 > max.cn]=NA
    
    ## gather the BF values that can be taken from the pre-computed table
    BFcc=table.cc[cbind(sample.ca+1,sample.cn+1)]
    
    ## replace the missing values with the pre-computed ones
    i.na=which(is.na(BFcc))
    if(length(i.na) > 0){
      BFcc[i.na]=BF.na
    }
  }else{
    ## if denovo only needed then set BFdn to 1
    BFcc=1
  }
  
  if(!denovo.only){
    # generate permutation data for denovo events\
    sample.rate = rgamma(nrep, param$rho0, param$nu0)
    sample.ar = rpois(nrep, n$ar*sample.rate^2/4)
    # look up the BF value in the table.
    BFar=table.ar[i.gene,sample.ar+1]
  }else{
    ## if denovo only needed then set BFdn to 1
    BFar=1
  }
  
  
  #determine the total BF from the two components
  BF=BFcc*BFdn*BFar
  #BF=BFcc*BFdn
  
  return(BF=BF)
}



Bayesian.FDR <- function(BF, pi0) {
  # Bayesian FDR control (PMID:19822692, Section2.3)
  # BF: a svector of BFs 
  # pi0: the prior probability that the null model is true
  # Return: the q-value of each BF, and the number of findings with q below alpha.
  
  # order the BF in decreasing order, need to retain order to get results back in proper order 
  i.order=order(BF, decreasing = T)
  BF=BF[i.order]
  # convert BFs to PPA (posterior probability of alternative model)
  pi <- 1-pi0
  q <- pi*BF/(1-pi+pi*BF) # PPA
  q0 <- 1 - q # posterior probability of null model
  
  # the FDR at each PPA cutoff
  FDR=cumsum(q0)/(1:length(BF))
  
  # reorder to the original order
  FDR[i.order]=FDR
  
  return (FDR=FDR)
}

bayesFactor.pvalue <- function(BF,BF.null){
  ## determines the pvalue for the BF using permutations under the null hypothesis BF.null
  # BF : vector with bayes factors  based on the data
  # BF.null : vector with bayes factors based on permuted data 
  
  BF.null <- sort(BF.null, decreasing=TRUE)
  pval <- findInterval(-BF, -BF.null)/length(BF.null)
  pval[pval==0] <- 0.5/length(BF.null)
  
  return(pval=pval)
}

denovo.MOM <- function(k, N, mu, C, beta, d=2, S=100, max.kvec = NULL) {
  # Estimating relative risk and the number of multiple hits from de novo data
  # Input:  k - number of disease genes
  #         N - sample size
  #         mu - mutation rate for all genes
  #         C - oberved number of denovo events
  #         beta - parameter of the prior distribution of gamma
  #         d - number of events to use (1 is 1 or more, 2 is 2 or more)
  #         S - number of samples to generate per gene
  #         max.kvec - used to generate a time line.
  # Output: gamma.mean  - the average relative risk, 
  #         M           - the expected number of multi-hit genes
  if(!is.null(max.kvec)){
    if(k%%100 == 0 | k == max.kvec){
      pb <- txtProgressBar(min = 0, max = max.kvec, style = 3)
      setTxtProgressBar(pb, k)
    }
  }
  
  m <- length(mu) # number of genes
  
  # enrichment of de novo events
  nu <- C / (2 * N * sum(mu))
  
  # MOM estimator of gamma.mean
  gamma.mean <- (nu-1)*m/k +1
  
  # expected M (choose d = 2)
  rs <- count.multihit(N, mu, k/m, gamma.mean, beta, d=d, S=S)
  M <- sum(rs$M1) + sum(rs$M0)
  
  return (list(gamma.mean=gamma.mean, M=M))
}


count.multihit <- function(N, mu, pi, gamma.mean, beta, d, S) {
  # Estimate the number of multihit genes in a genome. 
  # Input:  N - sample size
  #         mu - mutation rate for all genes
  #         pi - ratio of number of risk genes and total number of genes
  #         gamma.mean  - the average relative risk
  #         beta - parameter of the prior distribution of gamma
  #         d - number of events to use (1 is 1 or more, 2 is 2 or more)
  #         S - number of samples to generate per gene
  # Output: M0 - number of mutliple hit genes for the non-risk genes 
  #         M1 - number of multiple hit genes for risk genes
  m <- length(mu)
  
  # M1: the number of causal genes having d or more de novo mutations
  p.alt=NULL
  for(j in 1:length(d)){
    p.alt=cbind(p.alt,unlist(lapply(mu, multihit.prob, N, gamma.mean, beta, d=d[j], S=S)))
  }
  M1 <- m * pi * colMeans(as.matrix(p.alt))
  
  # M0: the number of non-causal genes having d or more de novo mutations
  p.null=NULL
  for (j in 1:length(d)) {
    p.null <-cbind(p.null,(1 - ppois(d[j], 2*N*mu)))
  }
  M0 <- m * (1-pi) * colMeans(as.matrix(p.null))
  
  result <- data.frame(d=d, M0=M0, M1=M1)
  return (result)
}

multihit.prob <- function(mu, N, gamma.mean, beta, d, S) {
  # Prob. of having d or more de novo mutations under H1
  # Use simulation, but could also use analytic form 
  # Input:  mu  - mutation rate for a gene
  #         N   - sample size
  #         gamma.mean  - the average relative risk
  #         beta - parameter of the prior distribution of gamma
  #         d - number of events to use (1 is 1 or more, 2 is 2 or more)
  #         S - number of samples to generate per gene
  # Output: p - avarage probabiltiy of having d or more de novo mutations
  gamma <- rgamma(S, gamma.mean*beta, rate=beta)
  p <- 1 - ppois(d, 2*N*mu*gamma)
  return (p=mean(p))
}

