rm(list=ls())
setwd("E:\\linfeng\\Methylation_data_N_S\\QCed_GWAS_plink")
read.table("als.pc1.qassoc",header = T) -> d
xlab <- paste("Theoretical", expression(italic(symbol(c))^2),"value")
ylab <- paste("Observed",expression(italic(symbol(c))^2),"value")
tiff(filename = "qqplot_P_adjust_eigenvalue.tiff",
     width = 2500, height = 2500, units = "px",res = 300,compression="lzw")
if(!require("ggplot2")) {
  print("Please install the ggplot plackage")
  install.packages('ggplot2', dep=T)
}

require(ggplot2)
axisSetting <- theme(  # remove grid line
  panel.background = element_rect(fill = "transparent",colour = NA), 
  panel.border = element_blank(),
  axis.line.y = element_line(colour = "black"),
  axis.line.x =  element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size=12,color='black',face='bold'),
  axis.title.y = element_text(size=12,color='black',face='bold'),
  axis.text.x = element_text(size=12,color='black',face='bold'),
  axis.text.y = element_text(size=12,color='black',face='bold'),
  #axis.ticks.x = element_blank(),
  #legend.position ="top",
  #legend.position ="none",
   legend.position = c(0.2, 0.5),
  legend.title=element_blank(),
  #legend.text=element_text(size=12,color='black',face='bold'),
  #legend.background=element_rect(size=1),
  #legend.key.size=unit(0.2,'cm'),
  #legend.key.width=unit(0.6,'cm'),
)


d <- na.omit(d)

d <- subset(d,as.numeric(d$CHR) <23)
maxy <- 30

sig <- -log10(0.05/length(d$P))
gc  <- qchisq(median(d$P),1, lower.tail = F)/0.455

qchisq(d$P,1,lower.tail = F)/gc -> d$P_adjust_GC

d$logp_GC <- d$P_adjust_GC

d$logp_eig <- qchisq(d$P,1,lower.tail = F)/6.19


d$logp     <- qchisq(d$P,1,lower.tail = F)

rep("P_adjust_GC",length(d$logp_GC)) -> P1

rep("P_adjust_eigenvalue",length(d$logp_eig)) -> P2

rep("P_noadjust",length(d$logp_eig)) -> P3


rbind(as.matrix(sort(d$logp_GC)),as.matrix(sort(d$logp_eig)),as.matrix(sort(d$logp))) -> P

rbind(as.matrix(P1),as.matrix(P2),as.matrix(P3)) -> P4

qchisq((1:length(d$logp_GC)/length(d$logp_GC)),1) -> chisq_vale

rbind(as.matrix(chisq_vale),as.matrix(chisq_vale),as.matrix(chisq_vale)) -> chisq



data.frame(P,P4,chisq) -> P_adjust


plot <- ggplot(P_adjust, aes(y = P,x = chisq, colour = factor(P4))) + geom_point()
       

plot <- plot + axisSetting

plot <- plot + xlab(quote(paste("Theoretical ", italic(symbol(c))^2," value"))) +
               ylab(quote(paste("Observed ",italic(symbol(c))^2," value")))

plot <- plot + geom_abline(slope = 1,intercept = 0)
print(plot)

dev.off()

