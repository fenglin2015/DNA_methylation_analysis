rm(list=ls())
setwd("E://linfeng//Methylation_data_N_S")
library(plyr)
library(gridExtra)
library(ggplot2)
load("meth_data/norm.beta.Robj")
read.table("Meth_anno.txt",header = TRUE) -> meth_anno
read.table("meth_predict_age_513.txt",header=TRUE) -> meth_age
#read.table("meth_predict_age.txt",header = T) -> meth_age
read.table("age_info.txt",header = T) ->age_info
read.table("EigenGWAS_meth/meth_qc.eigenvec",header = T) -> loc_pca
join(meth_age,meth_anno) -> meth_predict_age_adjust
read.table("loc_info.txt",header = T,fileEncoding = "UTF-8") ->loc_info
read.table("sex_info.txt",header=T) -> sex_info

for(i in 1:nrow(loc_pca))
{
  if(as.numeric(loc_pca$PC1[i])>0.0){
    loc_pca$PC_N_S[i]="S"
    
  }else{
    
    loc_pca$PC_N_S[i]="N"
  }
  
}
loc_pca <- data.frame(loc_pca$FID,loc_pca$PC_N_S)
colnames(loc_pca) <- c("FID","PC_N_S")
dimnames(norm.beta)[[2]] -> FID
dimnames(norm.beta)[[1]] -> meth_name
as.data.frame(norm.beta) -> norm.beta
norm.beta$meth_name <- meth_name
join(meth_predict_age_adjust[,1:2],norm.beta) -> norm_beta
write.table(norm_beta,"norm_beta.txt",col.names = T,row.names = T,quote = F,sep="\t")
for(i in 3:dim(norm_beta)[[2]])
{
  norm_beta[,i]=norm_beta[,i]*norm_beta[,2]
  
  
}
apply(norm_beta[,c(3:dim(norm_beta)[[2]])],2,sum,na.rm=T) -> age_predict
age_predict + norm_beta[1,2] -> age_predict
as.data.frame(age_predict)->age_predict
cbind(FID,age_predict) ->age_predict
dimnames(age_predict)[[2]] <- list(c("FID"),c("meth_age"))
join(age_predict,age_info) -> age
join(age,loc_pca) -> age
join(age,sex_info) -> age
age$age_den <- age$meth_age/age$age
#############################################################################################
axisSetting <- theme(  # remove grid line
  panel.background = element_rect(fill = "transparent",colour = NA), 
  panel.border = element_blank(),
  axis.line.y = element_line(colour = "black"),
  axis.line.x = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size=12,color='black',face='bold'),
  axis.title.y = element_text(size=12,color='black',face='bold'),
  axis.text.x = element_text(size=8,color='black',face='bold'),
  axis.text.y = element_text(size=10,color='black',face='bold'),
  #axis.ticks.x = element_blank(),
  #legend.position ="top",
  legend.position ="none",
  #legend.position = c(1.0, 0.5),
  legend.title=element_blank(),
  plot.title = element_text(hjust = 0.5,color='black',face='bold')
  #legend.text=element_text(size=12,color='black',face='bold'),
  #legend.background=element_rect(size=1),
  #legend.key.size=unit(0.2,'cm'),
  #legend.key.width=unit(0.6,'cm'),
)
####################################################################################################
age[-which(is.na(age$age)==TRUE),] -> age
age[-which(is.na(age$PC_N_S)==TRUE),] -> age_pc
age_pc_N <- age_pc[which(age_pc$PC_N_S=="N"),]
cor(age_pc_N$meth_age,age_pc_N$age)^2 -> cor_N
age_pc_S <- age_pc[which(age_pc$PC_N_S=="S"),]
cor(age_pc_S$meth_age,age_pc_S$age)^2 -> cor_S
age[-which(is.na(age$sex)==TRUE),]    -> age_sex
age_sex_M <- age_sex[which(age_sex$sex=="Male"),]
cor(age_sex_M$meth_age,age_sex_M$age)^2 -> cor_M
age_sex_F <- age_sex[which(age_sex$sex=="Female"),]
cor(age_sex_F$meth_age,age_sex_F$age)^2 -> cor_F
tiff(filename = "predict_age_513.tiff",
     width = 3000, height = 3000, units = "px",res = 400,compression="lzw")
###############
plot_N <- ggplot(age_pc_N,aes(age,meth_age))+geom_point(shape=20,size=2,col="Blue") 
plot_N <- plot_N + axisSetting
plot_N <- plot_N + geom_abline(slope=1,intercept=0)
plot_N <- plot_N + xlab("Age") + scale_x_continuous(limits = c(20,90)) + scale_y_continuous(limits = c(10,90)) +
  ylab("DNAm age") + ggtitle("N")
plot_N <- plot_N +annotate("text", x = 60, y = 25,
                           label = "paste(italic(R) ^ 2, \" = 0.74\")", parse = TRUE,size=5)
################
plot_S <- ggplot(age_pc_S,aes(age,meth_age))+geom_point(shape=20,size=2,col="red") 
plot_S <- plot_S  + axisSetting
plot_S <- plot_S  + geom_abline(slope=1,intercept=0)
plot_S <- plot_S  + xlab("Age") + scale_x_continuous(limits = c(20,90)) + scale_y_continuous(limits = c(10,90))+
  ylab("DNAm_age") + ggtitle("S")
plot_S <- plot_S  + annotate("text", x = 60, y = 25,
                             label = "paste(italic(R) ^ 2, \" = 0.83\")", parse = TRUE,size=5)
#################
plot_den_loc <- ggplot(age_pc,aes(x=age_den, color=factor(PC_N_S))) + geom_line(stat = "density")
plot_den_loc <- plot_den_loc + axisSetting  
plot_den_loc <- plot_den_loc + scale_color_manual(values = c("blue","red"))
plot_den_loc <- plot_den_loc + ggtitle("Density") + xlab("DNAm_age/Age ")
#################
plot_M  <-  ggplot(age_sex_M,aes(age,meth_age))+geom_point(shape=20,size=2,col="blue") 
plot_M  <-  plot_M  + axisSetting
plot_M  <-  plot_M  + xlab("Age")+scale_x_continuous(limits = c(20,90))+scale_y_continuous(limits = c(10,90))+
  ylab("DNAm_age") +  ggtitle("Male")
plot_M  <-  plot_M  + annotate("text", x = 60, y =25,
                               label = "paste(italic(R) ^ 2, \" = 0.76\")", parse = TRUE,size=5)
plot_M  <-  plot_M  + geom_abline(slope=1,intercept=0)
########################
plot_F  <-  ggplot(age_sex_F,aes(age,meth_age))+geom_point(shape=20,size=2,col="red") 
plot_F  <-  plot_F  + axisSetting
plot_F  <-  plot_F  + xlab("Age")+scale_x_continuous(limits = c(20,90))+scale_y_continuous(limits = c(10,90))+
  ylab("DNAm_age") +  ggtitle("Female")
plot_F  <-  plot_F  + annotate("text", x = 60, y = 25,
                               label = "paste(italic(R) ^ 2, \" = 0.76\")", parse = TRUE,size=5)
plot_F  <-  plot_F  + geom_abline(slope=1,intercept=0)
###################
plot_den_sex <- ggplot(age_sex,aes(x=age_den, color=factor(sex))) + geom_line(stat = "density")
plot_den_sex <- plot_den_sex + axisSetting  
plot_den_sex <- plot_den_sex + scale_color_manual(values = c("red","blue"))
plot_den_sex <- plot_den_sex +  ggtitle("Density") + xlab("DNAm_age/Age ")
######################
grid.arrange(plot_N,plot_S,plot_den_loc,plot_M,plot_F,plot_den_sex,nrow=2,ncol=3)
dev.off()
####################################################################################################