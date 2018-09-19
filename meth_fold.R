rm(list=ls())
setwd("/Users/mac/desktop/Methylation_data_N_S")
load("meth_data/norm.beta.Robj")
read.table("EigenGWAS_meth/meth_qc.eigenvec",header = T) -> pca
library(plyr)
dimnames(norm.beta)[[2]] -> FID
meth_sample <- cbind(FID,FID)
dimnames(norm.beta)[[1]] -> meth_name
write.table(meth_sample,"meth_sample.txt",col.names = F,row.names = F,quote = F,sep="\t")
join(as.data.frame(FID),pca) -> PHE_pca
write.table(PHE_pca,"PHE_pca.txt",row.names = F,col.names = T,quote=F,sep="\t")
sample(FID,689,replace = F) ->FID_random
write.table(FID_random,"FID_random.txt",row.names = F,col.names = T,sep="\t",quote = F)
#########################################################################################

list(c(1:137), c(138:275), c(276:413), c(414:551), c(552:689)) -> fold_order

for(temp_m in 1:5)
  {
   join(as.data.frame(FID),pca) -> PHE_pca
   for(i in fold_order[[temp_m]])
      {
        PHE_pca$PC1[which(PHE_pca$FID==FID_random[i])] <- NA
      }

   P  <-  rep(0,dim(norm.beta)[1])
   b0 <-  rep(0,dim(norm.beta)[1])
   b1 <-  rep(0,dim(norm.beta)[1])
  for(i in 1:dim(norm.beta)[1])
     {
       summary(lm(PHE_pca$PC1~norm.beta[i,])) -> beta
       P[i]    <- beta$coefficients[2,4]
       b0[i]   <- beta$coefficients[1,1]
       b1[i]   <- beta$coefficients[2,1]
     }
  paste("P_fold_",temp_m,".txt",sep = "") ->file_out
  write.table(cbind(b0,b1,P),file_out,row.names = F,col.names = F,sep="\t",quote = F)

####################################################################################

  Threshold <- c(1,0.5,0.1,1e-2,1e-3,1e-4)
  for(temp_i in 1:length(Threshold))
      {
        y_predit <- matrix(0,length(which(P < Threshold[temp_i])),length(which(is.na(PHE_pca$PC1)==TRUE)))
        temp_j = 1
        for(i in which(P<Threshold[temp_i]))
            {
              temp=1
 
              for(j in c(which(is.na(PHE_pca$PC1)==TRUE)))
                  {
                    y_predit[temp_j,temp] = b1[i] * norm.beta[i,j]
                   temp=temp+1
                  }
             temp_j=temp_j+1
               }
#y_predit <- b1[which(P <- Threshold[temp_i] )] * norm.beta[which(P <- Threshold[temp_i]), which(is.na(PHE_pca$PC1)==TRUE)]
paste("y_predict_fold",temp_m,"_",Threshold[temp_i],".txt",sep = "") ->file_out
write.table(y_predit,file_out,col.names = F,row.names = F,quote = F,sep="\t")
apply(y_predit,2,mean) -> pri_phe
FID_predict<-data.frame(PHE_pca[which(is.na(PHE_pca$PC1)==TRUE),]$FID,pri_phe)
dimnames(FID_predict)[[2]][1] <- "FID"
join(pca,FID_predict,by="FID") -> pca_predict
na.omit(pca_predict) -> pca_predict
cor.test(pca_predict$PC1,pca_predict$pri_phe,use="pairwise.complete.obs")$p.value->corr
cbind(Threshold[temp_i],corr) -> corr
cbind(temp_m,corr)            -> corr
write.table(corr,"corr_real_predict_fold.txt",append = T,quote = F,col.names = F,row.names = F,sep="\t")
      }

}

#########################################################################################