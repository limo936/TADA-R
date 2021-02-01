source("../src/TADA_R.R")

### read mutation data
tada.file="test_data.txt"
tada.data=read.table(tada.file,header=T)

prior_dmis=read.table("../data/posterior_dmis.txt", header = T)
prior_lof=read.table("../data/posterior_lof.txt", header = T)
prior_lof <- prior_lof[match(tada.data$gene.id, prior_lof$Gene),]
prior_dmis <- prior_dmis[match(tada.data$gene.id, prior_dmis$Gene),]

  
### specify the number of families and the number of cases and control samples included in the analysis
n.family = 3805
n = data.frame(dn=n.family, ca=n.family, cn=n.family, ca2=n.family, cn2=n.family, ar=n.family)
sample.counts <- list(cls1=n, cls2=n) #cls1=lof, cls2=dmis
  
# create the mutational data used by TADA
cls1.counts=data.frame(dn=tada.data$dn.lof, ca=tada.data$trans.lof, cn=tada.data$ntrans.lof, ca2=0, cn2=0, ar=tada.data$ar.lof)
rownames(cls1.counts)=tada.data$gene.id
cls2.counts=data.frame(dn=tada.data$dn.dmis, ca=tada.data$trans.dmis, cn=tada.data$ntrans.dmis, ca2=0, cn2=0, ar=tada.data$ar.dmis)
rownames(cls2.counts)=tada.data$gene.id
tada.counts=list(cls1=cls1.counts,cls2=cls2.counts)
  
### set up mutation rates
mu=data.frame(cls1=tada.data$mut.lof,cls2=tada.data$mut.dmis)
  
### set up denovo only TRUE/FALSE, here we do not want to restrict ourselves to de novo only analyses
denovo.only=data.frame(cls1=FALSE,cls2=FALSE)
  
### set up parameters
cls1= data.frame(gamma.mean.dn=34.0,beta.dn=1,gamma.mean.CC=9.5,beta.CC=4.0,gamma.mean.ar=16.0,beta.ar=4.0,rho1=prior_lof$Parameter1,nu1=prior_lof$Parameter2,rho0=prior_lof$Parameter1,nu0=prior_lof$Parameter2)
cls2= data.frame(gamma.mean.dn=20.0,beta.dn=1,gamma.mean.CC=4,beta.CC=1000,gamma.mean.ar=16.0,beta.ar=1000,rho1=prior_dmis$Parameter1,nu1=prior_dmis$Parameter2,rho0=prior_dmis$Parameter1,nu0=prior_dmis$Parameter2)
hyperpar=list(cls1=cls1,cls2=cls2)
  
# running TADA
re.TADA <- do.call(cbind.data.frame, TADA(tada.counts=tada.counts, sample.counts=sample.counts, mu=mu, hyperpar=hyperpar, denovo.only=denovo.only))

# Bayesian FDR control
re.TADA$qval=Bayesian.FDR(re.TADA$BF.total, pi0 = 0.95)
  
# run permutation to get the null distributions to use for calculating p-values for TADA
re.TADA.null=do.call(cbind.data.frame, TADAnull(tada.counts=tada.counts, sample.counts=sample.counts, mu=mu, hyperpar=hyperpar, denovo.only=denovo.only, nrep=2e4))
re.TADA$pval=bayesFactor.pvalue(re.TADA$BF.total,re.TADA.null$BFnull.total)

#print result  
print(re.TADA)

