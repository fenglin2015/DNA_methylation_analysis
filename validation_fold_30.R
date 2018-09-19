 rm(list=ls())
 setwd("E:/linfeng//Methylation_data_N_S")
 load("meth_data/norm.beta.Robj")
 read.table("EigenGWAS_meth/meth_qc.eigenvec",header = T) -> pca
 library(plyr) 
 library(RcppArmadillo)
 library(ggplot2)
 colnames(norm.beta)  -> FID
 rownames(norm.beta)  -> meth_name
 join(as.data.frame(FID),pca) -> PHE_pca


#########################################################################################
 for(num in 1:30)
+ {
+   sample(FID,689,replace = F)  -> FID_random
+   
+   join(as.data.frame(FID),pca) -> PHE_pca
+   for(i in 1:137)
+   {
+     PHE_pca$PC1[which(PHE_pca$FID==FID_random[i])] <- NA
+   }
+   
+   P   <-  rep(0,nrow(norm.beta))
+   std <-  rep(0,nrow(norm.beta))
+   b1  <-  rep(0,nrow(norm.beta))
+   idx <- which(!is.na(PHE_pca$PC1))
+   for(i in 1:nrow(norm.beta))
+   {
+     mod <- fastLmPure(cbind(1,norm.beta[i,idx]),PHE_pca$PC1[idx])
+     
+     b1[i]  <- mod$coefficients[2,1]
+     std[i] <- mod$stderr[2,1]
+     temp   <- b1[i]/std[i]
+     P[i]   <- 2*pt(abs(temp), df=mod$df.residual, lower.tail = F)
+     
+   }
+   paste("P_fold_",num,".txt",sep = "") -> file_out
+   write.table(cbind(b1,P),file_out,row.names = F,col.names = F,sep="\t",quote = F)
+   
+   ####################################################################################
+   
+   Threshold <- c(1,0.5,0.1,1e-2,1e-3,1e-4)
+   for(temp_i in 1:length(Threshold))
+   {
+     y_predit <- matrix(0,length(which(P < Threshold[temp_i])),
+                        length(which(is.na(PHE_pca$PC1)==TRUE)))
+     y_predit <- b1[which(P<Threshold[temp_i] )] * norm.beta[which(P < Threshold[temp_i]), which(is.na(PHE_pca$PC1)==TRUE)]
+     paste("y_predict_fold",num,"_",Threshold[temp_i],".txt",sep = "") -> file_out
+     write.table(y_predit,file_out,col.names = F,row.names = F,quote = F,sep="\t")
+     apply(y_predit,2,mean) -> pri_phe
+     FID_predict<-data.frame(PHE_pca[which(is.na(PHE_pca$PC1)==TRUE),]$FID,pri_phe)
+     dimnames(FID_predict)[[2]][1] <- "FID"
+     join(pca,FID_predict,by="FID") -> pca_predict
+     na.omit(pca_predict) -> pca_predict
+     cor(pca_predict$PC1,pca_predict$pri_phe) -> corr
+     cbind(Threshold[temp_i],corr) -> corr
+     cbind(num,corr)               -> corr
+     write.table(corr,"corr_real_predict_fold.txt",append = T,quote = F,col.names = F,
+                 row.names = F,sep="\t")
+   }
+   
+ }
axisSetting <- theme(  # remove grid line
+   panel.background = element_rect(fill = "transparent",colour = NA), 
+   panel.border = element_blank(),
+   axis.line.y = element_line(colour = "black"),
+   axis.line.x = element_line(colour = "black"),
+   panel.grid.major = element_blank(),
+   panel.grid.minor = element_blank(),
+   axis.title.x = element_text(size=12,color='black',face='bold'),
+   axis.title.y = element_text(size=12,color='black',face='bold'),
+   axis.text.x = element_text(size=10,color='black',face='bold'),
+   axis.text.y = element_text(size=10,color='black',face='bold'),
+   #axis.ticks.x = element_blank(),
+   #legend.position ="top",
+   legend.position ="none",
+   #legend.position = c(1.0, 0.5),
+   legend.title=element_blank(),
+   #legend.text=element_text(size=12,color='black',face='bold'),
+   #legend.background=element_rect(size=1),
+   #legend.key.size=unit(0.2,'cm'),
+   #legend.key.width=unit(0.6,'cm'),
+ )
 tiff(filename = "validation_30times.tiff",
+      width = 2500, height = 2500, units = "px",res = 300,compression="lzw")
 
 plot <- ggplot(corr_par, aes(x=factor(Sig_level), y=corr_mean)) + geom_bar(position="dodge",width=0.5,stat="identity") +
+         geom_errorbar(aes(ymax=corr_mean + corr_sd,ymin=corr_mean - corr_sd),position = position_dodge(0.9), width = 0.15)
 
 plot <- plot + axisSetting
 
 plot <- plot + xlab("Sig_level") + ylab("Correlation")
 
 
 print (plot)

 dev.off()