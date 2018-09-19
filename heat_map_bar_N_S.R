rm(list=ls())
setwd("E:\\linfeng\\map")
legend_title<- expression(log[10](italic(p)))
legend_title<- "-log10(p-value)"
#############################################################################
tiff(filename = "map.tiff",
     width = 2500, height = 2500, units = "px",res = 300,compression="lzw")
read.table("E:\\linfeng\\pca_meth_case_control.txt",header = T,fileEncoding = "UTF-8")->loc_pca
read.table("E:\\linfeng\\long_lat.txt",header=T,fileEncoding = "UTF-8")->long_lat
as.character(long_lat$lg_la)->long_lat$lg_la
as.factor(long_lat$NAME)->temp
#loc_pca$pca1<--loc_pca$pca1
fun<-function(x){sample(x,1,replace = T)}
for(i in 1:length(loc_pca$indi))
{
  tapply(long_lat$lg_la,temp,fun)->longlat
  for(j in 1:length(levels(temp)))
   {
    
    if(loc_pca$loc[i]==levels(temp)[j])
    {
      loc_pca$long[i]<-as.numeric(strsplit(longlat[j],"_")[[1]][1])
      loc_pca$lat[i] <-as.numeric(strsplit(longlat[j],"_")[[1]][2])
      
    }
  
  }
  
}

for(i in 1:length(loc_pca$indi))
{
  if(as.numeric(loc_pca$pca1[i])>0.0){
    loc_pca$pca3_[i]=1
    
  }else{
    
    loc_pca$pca3_[i]=0
  }
  
}
library(RColorBrewer)
library(gridExtra)
library(rgeos)
library(ggplot2)
library(maps)
library(maptools)
library(mapdata)
library(plyr)
cols<-brewer.pal(8,"Dark2")[1:8]
newpalette<-colorRampPalette(brewer.pal(9,"Blues"))(10)
china_map<-readShapePoly("bou2_4p.shp")
province_city<-read.csv("chinaprovincecity.csv",header=T,as.is=T)
data1<- china_map@data      
data2<- data.frame(id=row.names(data1),data1) 
china_map1 <- fortify(china_map) 
china_data <- join(china_map1,data2, type = "full") 

subset(loc_pca,pca3_==0)->loc_pca_0
subset(loc_pca,pca3_==1)->loc_pca_1
plot1<-ggplot(china_data,aes(long,lat))+
  geom_polygon(aes(group=group),fill="white",colour="grey60")+
  geom_point(data=loc_pca_0,aes(x=long,y=lat,colour="red"))+
  coord_map("polyconic") + labs(x="PCA<0")+
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )
plot2<-ggplot(china_data,aes(long,lat))+
  geom_polygon(aes(group=group),fill="white",colour="grey60")+
  geom_point(data=loc_pca_1,aes(x=long,y=lat,colour="red"))+
  coord_map("polyconic") + labs(x="PCA>0")+
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )
grid.arrange(plot1,plot2,nrow=1,ncol=2)

dev.off()

###########################################################################

tiff(filename = "map1.tiff",
     width = 2500, height = 2500, units = "px",res = 300,compression="lzw")

 ggplot(china_data,aes(long,lat))+
  geom_polygon(aes(group=group),fill="white",colour="grey60")+
  geom_point(data=loc_pca,aes(x=long,y=lat,colour=factor(pca3_)))+
  coord_map("polyconic") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )
dev.off()

##############################################################################

tiff(filename = "bar_plot.tiff",
     width = 2500, height = 2500, units = "px",res = 300,compression="lzw")

loc  <-as.factor(loc_pca$loc)
tapply(loc_pca$pca3_,loc,count)->loc_count
loc  <-rep(0,length(loc_count))
pca_0<-rep(0,length(loc_count))
pca_1<-rep(0,length(loc_count))
pca_count<-data.frame(loc,pca_0,pca_1)

for(i in 1:length(loc_count))
{
  pca_count$loc[i]<-names(loc_count)[i]
  if(loc_count[[i]]$x[1]==0)
  {
    pca_count$pca_0[i] <-loc_count[[i]]$freq[1]
  }
  if(loc_count[[i]]$x[1]==1)
  {
    pca_count$pca_1[i] <-loc_count[[i]]$freq[1]
  }
  if(length(loc_count[[i]]$freq)==2)
  {
    pca_count$pca_1[i] <-loc_count[[i]]$freq[2]
  }
    
  
}
rbind(as.matrix(pca_count$loc),as.matrix(pca_count$loc))    ->loc
rbind(as.matrix(pca_count$pca_0),as.matrix(pca_count$pca_1))->num
rbind(as.matrix(rep("S",length(pca_count$loc))),as.matrix(rep("N",length(pca_count$loc))))->temp
pca_count_adjust<-data.frame(loc,temp,num)
loc_pca[order(loc_pca$loc),] ->loc_pca_sort
order                        <- seq(0,(length(levels(loc_pca_sort$loc)))*4,4)
ticks<-seq(0.5,(length(levels(loc_pca_sort$loc)))*4,4)

x<-rep(0,length(loc_pca_sort$pca1))
count_i<-rep(0,length(levels(loc_pca_sort$loc)))
for(i in 1:length(levels(loc_pca_sort$loc)))
{
  length(loc_pca_sort[loc_pca_sort$loc==levels(loc_pca_sort$loc)[i],]$pca1)->count_i[i]
  temp<-rep(0,count_i[i])
  sample(seq(order[i],order[i+1]-3,0.001),count_i[i])->temp 
  x[head(which(x==0),1):(head(which(x==0),1)+count_i[i]-1)]<-temp
  
}

 
loc_pca_sort$x <-x

as.data.frame(loc_pca_sort)->loc_pca_sort

subset(loc_pca,pca3_==0)->loc_pca_0
subset(loc_pca,pca3_==1)->loc_pca_1
as.factor(loc_pca_0$loc) ->factor_temp
tapply(loc_pca_0$pca1,factor_temp,mean)->pca0_sum
as.data.frame(pca0_sum)->pca0_sum
as.factor(loc_pca_1$loc) ->factor_temp
tapply(loc_pca_1$pca1,factor_temp,mean)->pca1_sum
as.data.frame(pca1_sum)->pca1_sum
loc  <-rep(0,length(loc_count))
loc_pca_mean<-data.frame(loc)
loc_pca_mean$loc<-pca_count$loc
loc_pca_mean$pca0_mean<-pca0_sum$pca0_sum
loc_pca_mean$pca1_mean<-pca1_sum$pca1_sum
loc_pca_mean$pca1_mean[is.na(loc_pca_mean$pca1_mean)]<-NA
loc_pca_mean$pca0_mean[is.na(loc_pca_mean$pca0_mean)]<-NA
rbind(as.matrix(loc_pca_mean$loc),as.matrix(loc_pca_mean$loc))    ->loc
rbind(as.matrix(loc_pca_mean$pca0_mean),as.matrix(loc_pca_mean$pca1_mean))->num
rbind(as.matrix(rep("S",length(loc_pca_mean$loc))),as.matrix(rep("N",length(loc_pca_mean$loc))))->temp
rbind(as.matrix(ticks),as.matrix(ticks)) ->ticks_com
loc_pca_mean_adjust<-data.frame(loc,temp,num,ticks_com)
join(loc_pca_mean_adjust,loc_pca_sort)->loc_pca_mean_adjust
axisSetting <- theme(  # remove grid line
  panel.background = element_rect(fill = "transparent",colour = NA), 
  panel.border = element_blank(),
  axis.line.y = element_line(colour = "black"),
  axis.line.x = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size=12,color='black',face='bold'),
  axis.title.y = element_text(size=12,color='black',face='bold'),
  axis.text.x = element_text(angle=90,size=10,color='black', hjust=1, vjust=.5),
  axis.text.y=element_text(size=10,color='black'),
  #axis.ticks.x = element_blank(),
  legend.position =c(0.95,0.8),
  legend.title = element_blank()
)
en_name<-read.csv("Loc_Engname.csv",header=T,as.is=T,fileEncoding = "UTF-8")
pca_count_adjust<-join(pca_count_adjust,en_name)
#pca_count_adjust<-within(pca_count_adjust,{city<-factor(city,levels=c("Anhui","Beijing","Fujian","Guangdong","Guangxi","Guizhou","Hebei",""))})
plot <-ggplot(pca_count_adjust,aes(x=city,y=num,fill=temp))+geom_bar(position="dodge",width=0.5,stat="identity")
plot <-plot+axisSetting
plot <-plot+scale_y_continuous(expand=c(0,0),limits = c(0,max(pca_count_adjust$num)))
plot <-plot+scale_fill_manual(values=c(cols[1],cols[4]))
plot <-plot+xlab("Province")+ ylab(" ")
plot
loc_pca_mean_adjust<-join(loc_pca_mean_adjust,en_name)
loc_pca_sort<-join(loc_pca_sort,en_name)
loc_pca_mean_adjust[order(loc_pca_mean_adjust$city),]->loc_pca_mean_adjust
plot1 <-ggplot(loc_pca_mean_adjust,aes(x=x,y=pca1,color=factor(loc_pca_mean_adjust$pca3_)))+geom_point(size=0.5)
plot1 <-plot1+scale_x_continuous(expand = c(0.01,0),
  limits = c(order[1],order[length(order)]-3),
   breaks=ticks,
  labels=(unique(loc_pca_mean_adjust$city)))
plot1 <-plot1+axisSetting
plot1 <-plot1+scale_colour_manual(values=c(cols[4],cols[1],cols[1],cols[4]))
plot1 <-plot1+xlab("Province")+ ylab(" ")
plot1 <-plot1+theme(legend.position = "none")
plot1 <-plot1+geom_point(aes(x=ticks_com,y=num,color=factor(loc_pca_mean_adjust$temp)),size=4,shape=3)
plot1 <-plot1+axisSetting
plot1 <-plot1+xlab("Province")+ ylab(" ")
plot1 <-plot1+theme(legend.position = "none")
plot1 <-plot1+geom_hline(aes(yintercept=0))
grid.arrange(plot,plot1,nrow=2,ncol=1)

dev.off()

################################################################################################
tiff(filename = "heat_map_bar.tiff",
     width = 2500, height = 2500, units = "px",res = 300,compression="lzw")
subset(loc_pca,pca3_==0) -> pca0_t
subset(loc_pca,pca3_==1) -> pca1_t
NAME<-rep(0,length(levels(pca0_t$loc)))
pca_t<-data.frame(NAME)
pca_t$NAME<-levels(pca0_t$loc)
#temp2<-rep(0,length(levels(pca0_t$loc)))
for(i in 1:length(levels(pca0_t$loc)))
{
  if(length(pca0_t$pca1[which(pca0_t$loc==levels(pca0_t$loc)[i])])>1&length(pca1_t$pca1[which(pca1_t$loc==levels(pca1_t$loc)[i])])>1)
{   t.test(-pca0_t$pca1[which(pca0_t$loc==levels(pca0_t$loc)[i])],
       pca1_t$pca1[which(pca1_t$loc==levels(pca1_t$loc)[i])])[3][1]->temp1
    t.test(-pca0_t$pca1[which(pca0_t$loc==levels(pca0_t$loc)[i])],
           pca1_t$pca1[which(pca1_t$loc==levels(pca1_t$loc)[i])])[1][1]->temp2
          pca_t$t_pvalue[i]<--log10(temp1$p.value)*sign(temp2$statistic)
    
  }
  else if(length(pca0_t$pca1[which(pca0_t$loc==levels(pca0_t$loc)[i])])<2&length(pca1_t$pca1[which(pca1_t$loc==levels(pca1_t$loc)[i])])>1){
    pca_t$t_pvalue[i]<--12
  }else if(length(pca0_t$pca1[which(pca0_t$loc==levels(pca0_t$loc)[i])])>1&length(pca1_t$pca1[which(pca1_t$loc==levels(pca1_t$loc)[i])])<2){
    pca_t$t_pvalue[i]<-12
  } else {
    pca_t$t_pvalue[i]<-NA
  }
  
}
na.omit(pca_t)->pca_t
pca_t[pca_t$t_pvalue>5,]$t_pvalue<-5
pca_t[pca_t$t_pvalue<(-5),]$t_pvalue<--5
pca_t$province<-pca_t$NAME
join(province_city,pca_t)->province_city_adjust
province_city_adjust[-which(is.na(province_city_adjust$t_pvalue)),]->province_city_adjust
rename(pca_t,c(t_pvalue=legend_title))->pca_t
china_map_data <- join(china_data, pca_t,type="full")
as.data.frame(china_map_data)->china_map_data
ggplot(china_map_data,aes(long,lat))+
  geom_polygon(aes(group=group,fill=`-log10(p-value)`),colour="grey",size=0.01)+
  scale_fill_gradient(low="green",high="red")+
  coord_map("polyconic")+
  geom_text(aes(x=jd,y=wd,label=city),data=province_city_adjust,colour="black",size=2)+
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.border = element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    axis.title=element_blank(),
    legend.position=c(0.9,0.5)
  )
dev.off()

##################################################################################################
tiff(filename = "heat_map.tiff",
     width = 2500, height = 2500, units = "px",res = 300,compression="lzw")


sum(pca_count$pca_0)        ->pca_0
sum(pca_count$pca_1)        ->pca_1
pca_0/(pca_1+pca_0)         ->freq_0
pca_count$freq      <-pca_count$pca_0*(pca_count$pca_0+pca_count$pca_1)/(pca_count$pca_0+pca_count$pca_1)
Fun<-function(x){(x-freq_0*pca_0)/sqrt(freq_0*(1-freq_0))}
apply(as.matrix(pca_count$freq),2,Fun)->pca_count$pca_t
pca_count$pca_t<-pca_count$pca_t/sqrt(pca_count$pca_0+pca_count$pca_1)
pt(pca_count$pca_t,as.matrix(pca_count$pca_0+pca_count$pca_1)-1)->pca_count$p_value
pca_count$NAME<-pca_count$loc
pca_count$p_value[is.nan(pca_count$p_value)]<-NA
pca_count$p_value<--log10(pca_count$p_value)
pca_count$p_value<-pca_count$p_value*sign((pca_count$pca_0-pca_count$pca_1)+0.1)
rename(pca_count,c(p_value=legend_title))->pca_count
 
china_map_data <- join(china_data, pca_count, type="full")
pca_count$province<-pca_count$NAME
join(province_city,pca_count)->province_city_adjust
province_city_adjust[-which(is.na(province_city_adjust$`-log10(p-value)`)),]->province_city_adjust

ggplot(china_map_data,aes(long,lat))+
  geom_polygon(aes(group=group,fill=`-log10(p-value)`),colour="grey",size=0.01)+
  scale_fill_gradient(low="green",high="red")+
  coord_map("polyconic")+
  geom_text(aes(x=jd,y=wd,label=city),data=province_city_adjust,colour="black",size=2)+
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.border = element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    axis.title=element_blank(),
    legend.position=c(0.9,0.5)
  )
dev.off()

###################################################################################
