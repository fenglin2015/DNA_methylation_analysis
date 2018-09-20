###################
methy_eigen=read.table("methy_set_QC_1.qassoc", as.is = T, header = T) #EigenGWAS
#full_eigen=read.table("../gwas/als.pc1.qassoc", as.is = T, header = T) #
#mg=merge(x=methy_eigen, y=full_eigen, by.x = "SNP", by.y = "SNP", sort = F)
#plot(-log10(mg$P.x), -log10(mg$P.y), cex=0.5, pch=16, bty="l")

load("../QCed_methylation_data/norm.beta.Robj")
y=array(NA, dim = ncol(norm.beta))
yPC=array(NA, dim = ncol(norm.beta))
for(i in 1:length(y)) {
  idx=which(evec[,1] == colnames(norm.beta)[i])
  if (length(idx)>0) {
    if(evec[idx,3]>0) {
      y[i] = 1
    } else {
      y[i] = 0
    }
    yPC[i] = evec[idx,3]
  }
}

bReg=matrix(0, nrow(norm.beta), 4)
bRegPC=matrix(0, nrow(norm.beta), 4)
library(RcppArmadillo)
for(i in 1:nrow(norm.beta)) {
  idx=which(!is.na(y))
#  mod=lm(y[idx]~norm.beta[i,idx])
  mod=fastLmPure(cbind(1,norm.beta[i,idx]),y[idx])
  bReg[i,1]=mod$coefficients[2,1]
  bReg[i,2]=mod$stderr[2,1]
  bReg[i,3]=bReg[i,1]/bReg[i,2]
  bReg[i,4]=2*pt(abs(bReg[i,3]), df=mod$df.residual, lower.tail = F)

  modPC=fastLmPure(cbind(1, norm.beta[i,idx]), yPC[idx])
  bRegPC[i,1]=modPC$coefficients[2,1]
  bRegPC[i,2]=modPC$stderr[2,1]
  bRegPC[i,3]=bRegPC[i,1]/bRegPC[i,2]
  bRegPC[i,4]=2*pt(abs(bRegPC[i,3]), df=modPC$df.residual, lower.tail = F)
}
hist(bReg[,4])
print(median(bReg[,4]))
mGC=qchisq(median(bReg[,4]), 1, lower.tail = F)/qchisq(0.5, 1, lower.tail = F)
print(mGC)