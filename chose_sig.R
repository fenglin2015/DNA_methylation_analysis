rm(list=ls())
library(plyr)
setwd("E:\\linfeng\\QCed_methylation_data\\")
load("norm.beta.Robj")
dimnames(norm.beta)[[2]]->indi
as.data.frame(indi)     ->indi
dimnames(norm.beta)[[1]]->rowname
read.table("Meth_pca_maf0.03.eigenvec",header = T)  -> pca
read.table("age.txt",header=T) ->age
read.table("chr14_snp_sig.ped",header = F) ->sex
sex[,-c(2:4,6:dim(sex)[[2]])] ->sex
dimnames(sex)[[2]]<-c("indi","sex")
for(i in 1:length(pca$indi))
{
  if(as.numeric(pca$pca1[i])>0.0){
    pca$phe[i]=0
    
  }else{
    
    pca$phe[i]=1
  }
  
}
join(indi,pca) ->indi_adjust
join(indi_adjust,age) ->indi_adjust_age
join(indi_adjust_age,sex) ->indi_adjust_age_sex
write.table(indi_adjust_age_sex,"indi_pca_N_S_age_sex.phe",row.names =F,col.names = T,quote = F,sep="\t")

P<-rep(0,dim(norm.beta)[1])
for(i in 1:dim(norm.beta)[1])
{
  glm(indi_adjust_age_sex$phe~norm.beta[i,])->temp
  summary(temp)$coefficients[2,4] ->P[i]
  
  
}
row.names(norm.beta)->meth_name
as.data.frame(meth_name) -> meth_name
cbind(meth_name,P) ->Man_meth


read.table("meth_anno.txt",header = T)->meth_anno

join(meth_anno,Man_meth) ->Man_draw
Man_draw$logp            <- -log10(Man_draw$P)



write.table(Man_draw,"indi_case_control_glm.txt",row.names =F,col.names = T,quote = F,sep="\t")
tiff(filename = "Manhadun_pca_glm_delsnp.tiff",
     width = 2500, height = 2000, units = "px",res = 300,compression="lzw")
if(!require("ggplot2")) {
  print("Please install the ggplot plackage")
  install.packages('ggplot2', dep=T)
}

require(ggplot2)

d <- read.table("indi_case_control_delsnp_glm.txt",header=T,stringsAsFactors = F)

d <- na.omit(d)
d$CHR[which(d$CHR=="X")]<-23
d$CHR[which(d$CHR=="Y")]<-23
d$CHR <- as.numeric(d$CHR)
d$BP  <- as.numeric(d$BP)
d$logp<- as.numeric(d$logp)
d <- subset(d,as.numeric(d$CHR) <23)
maxy <- 30

sig <- -log10(0.05/length(d$P))
gc<-qchisq(median(d$P),1, lower.tail = F)/0.455

qchisq(d$P,1,lower.tail = F)/gc->d$P_adjust

pchisq(d$P_adjust,1,lower.tail = F)->d$P_adjust


d$logp <- -log10(d$P_adjust)
#d$logp <- -log10(d$P)
d$x <- NA
ticks <- NULL
lastbase <- 0

order <- order(d$CHR,d$BP)
d <- d[order,]
d[d$logp > maxy,'logp'] <- maxy

numchrs <- length(unique(d$CHR))

for( i in unique(d$CHR) ) {
  if( i==1) {
    d[d$CHR==i,]$x <- d[d$CHR==i,]$BP
  } else {
    lastbase <- lastbase+tail(subset(d,CHR==i-1)$BP,1)
    d[d$CHR==i,]$x <- d[d$CHR==i,]$BP+lastbase
  }
  ticks <- c(ticks,d[d$CHR==i,]$x[floor(length(d[d$CHR==i,]$x)/2)+1])
}
ticklim <- c(min(d$x),max(d$x))

if(numchrs > 1) {
  cols <- rep(c("red","blue"),max(d$CHR))
} else {
  cols <- "gray"
}

if (maxy=="max") maxy=ceiling(max(d$logp)) else maxy=maxy
if (maxy<8) maxy=8
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
  axis.text.x = element_text(size=8,color='black',face='bold'),
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

plot<- ggplot(d,aes(x=x,y=logp,colour=factor(CHR)))+geom_point(size=1)

plot <- plot+ylab(expression(-log[10](italic(p))) )
plot <- plot + scale_x_continuous(expand=c(0.02,0),name="Chromosome", limits=ticklim,
                                  breaks=ticks,
                                  labels=(unique(d$CHR)))
plot <- plot + scale_y_continuous(expand=c(0,0),limits = c(0,max(d$logp)+0.2))
plot <- plot + scale_colour_manual(values=cols)
#plot <- plot + ggtitle("Manhattan Plot")
plot <- plot + axisSetting
plot <- plot + geom_abline(intercept=sig,slope=0, col="black")

print(plot)
dev.off()
subset(Man_draw,logp>sig) ->chose_P
write.table(chose_P,"chose.txt",row.names = F,col.names = T,quote = F,sep="\t")

