rm(list=ls())

library(ggplot2)

setwd("E:/linfeng/Methylation_data_N_S/comm_1KG_ALS/maf0.03")
#########################################
#read.table("ALS.PC1.qassoc",header = T)      -> PC1

#read.table("CHB_CHS_recode12.txt",header=F) -> CHB_CHS

#read.table("CHB_CHS.eigenvec",header = T)    -> CHB_CHS_PC
######################################### 
read.table("CHB_ChS.PC1.qassoc",header = T)   -> PC1

read.table("ALS_redcode12.txt",header=F)      -> CHB_CHS

read.table("ALS.eigenvec",header = T)         -> CHB_CHS_PC

PC1$logp <- -log10(PC1$P)

cor <- rep(0,6)

y_predict <- matrix(NA,nrow = 6,ncol = ncol(CHB_CHS)-2)

for(i in 1:6)
{
 
  
  y_predict[i,] <- t(as.matrix(PC1$BETA[which(PC1$logp > i)])) %*%
               as.matrix(CHB_CHS[which(PC1$logp > i),c(3:4087)])
  
  #y_predict <- y_predict/length(which(PC1$logp > i))
  
  cor(y_predict[i,],CHB_CHS_PC$PC1) -> cor[i]
  
}

#write.table(y_predict,"y_predict.txt",col.names = F,row.names = F,quote = F,sep="\t")

Thresold <- c(1,2,3,4,5,6)

t(rbind(Thresold,cor)) -> corr_par

colnames(corr_par) <- c("sig_level","corr")

as.data.frame(corr_par) -> corr_par

##################################################################
axisSetting <- theme(  # remove grid line
  panel.background = element_rect(fill = "transparent",colour = NA), 
  panel.border = element_blank(),
  axis.line.y = element_line(colour = "black"),
  #axis.line.x = element_line(colour = "black"),
  axis.line.x = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size=12,color='black',face='bold'),
  axis.title.y = element_text(size=12,color='black',face='bold'),
  axis.text.x = element_text(size=12,color='black',face='bold'),
  axis.text.y = element_text(size=12,color='black',face='bold'),
  #axis.ticks.x = element_blank(),
  #legend.position ="top",
  legend.position ="none",
  #legend.position = c(1.0, 0.5),
  legend.title=element_blank(),
  #legend.text=element_text(size=12,color='black',face='bold'),
  #legend.background=element_rect(size=1),
  #legend.key.size=unit(0.2,'cm'),
  #legend.key.width=unit(0.6,'cm'),
)
##############################################################
tiff(filename = "Train_CHB_CHS_cor.tiff",
          width = 2500, height = 2500, units = "px",res = 300,compression="lzw")

plot <- ggplot(corr_par, aes(x=factor(sig_level), y=corr)) + 
        geom_bar(position="dodge",width=0.5,stat="identity") 

plot <- plot + axisSetting

plot <- plot + xlab(expression(-log[10](italic(p)))) + ylab("Correlation") +
         scale_y_continuous(expand=c(0,0))


print (plot)

dev.off()

