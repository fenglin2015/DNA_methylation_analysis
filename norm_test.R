normal_test<- function(input.data,alpha=0.05,picplot=TRUE){
  if(picplot==TRUE){#画图形
    dev.new()#新建窗口画图
    #par(mfrow=c(2,1))
    
    #Q-Q图法
   # qqnorm(input.data,main="qq图")
   # qqline(input.data)
    
    #概率密度曲线比较法
    hist(input.data,freq=F,main="直方图和密度估计曲线")
    #如果画出的图缺少尖端部分则使用下面这句代码
    #hist(input.data,freq=F,main="直方图和密度估计曲线",ylim = c(0,0.5))#使用合适的值来避免红蓝线缺少尖端部分，这里根据已经跑出来的图像我得出0.5
    lines(density(input.data),col="blue") #密度估计曲线
    x<-seq(min(input.data),max(input.data),0.0001)
    #使用seq(),若取0.0000001太密集跑大一点的数据就容易死机，建议0.0001
    lines(x,dnorm(x,mean(input.data),sd(input.data)),col="red") 
    #正态分布曲线，思想是根据求每个x应该对应的标准正态y值，然后将x与求出的y放在一起做出所求数据如果按照正态分布应该是怎样的，并于实际密度曲线（蓝线）对比 
  }#sd标准差 mean平均值
  
  #夏皮罗-威尔克（Shapiro-Wilk）检验法【数据不能过大，范围为3~5000，假如有一个300*300的矩阵那么这个方法运行函数时作废】
  shapiro_result<- shapiro.test(input.data)
  if(shapiro_result$p.value>alpha){
    print(paste("success:服从正态分布,p.value=",shapiro_result$p.value,">",alpha))    
  }else{
    print(paste("error:不服从正态分布,p.value=",shapiro_result$p.value,"<=",alpha))
 }
  shapiro_result
}