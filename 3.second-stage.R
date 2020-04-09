
################################################################################
# This is the code with sample data for the analysis in:
#
#   "Mortality burden attributable to ambient temperature and temperature variability: 
#    results from 364 Chinese locations"
#
#   Author:
#     Tao Liu, Chunliang Zhou, Yize Xiao, Biao Huang, Yanjun Xu, Xiaojun Xu, 
#     Lifeng Lin, JianXiong Hu, Jianpeng Xiao, Weilin Zeng, Xing Li, Siqi Chen, 
#     Lingchuan Guo, Yonghui Zhang, Cunrui Huang, Yaodong Du, Min Yu, Maigeng Zhou, 
#     Wenjun Ma
#
#   Affiliation:
#     Guangdong Provincial Institute of Public Health
#     Guangdong Provincial Center for Disease Control and Prevention
#   
#   Available at:
################################################################################

################################################################################
# 3.second-stage analysis
################################################################################

###(1)Analysis of TM ####
mvtm <- mvmeta(coeftm,vcovtm,method = "reml")

##extact coef and vcov after meta analysis
mvcoeftm <- coef(mvtm)
mvvcovtm <- vcov(mvtm)

##get the cumulative curve
#Generate temperature series
perctm <- quantile(data$tm,c(0:1000)/1000,na.rm=T)

#Set parameters and temperature range
tmknots2 <-  quantile(data$tm,c(10,50,90)/100,na.rm=TRUE)
bound <- range(data$tm,na.rm = T)

#genarate one-basis
argvar1 <- list(fun=varfun,degree=vardegree,
                knots=tmknots2,Boundary.knots=bound)
cb2 <- crossbasis(perctm,lag=lag,argvar=argvar1,
                  arglag=arglag)
bvar <- do.call("onebasis",c(list(x=perctm),attr(cb2,"argvar")))

#make predictions based on the results of meta analysis
  #Note:The curve is inconsistent with the results 
  #     in the article because of limited locations
tm.cp <- crosspred(bvar,coef=mvcoeftm,vcov=mvvcovtm,at=perctm)

#Get MMT and MMP
MMT <- round(as.numeric(names(which.min(tm.cp$allfit))),2)
MMT

MMP <- names(which(round(quantile(data$tm,c(0:1000)/1000,na.rm = T),2)==MMT))
MMP
#Set MMT and re-predict
tm.cpn <- crosspred(bvar,coef=mvcoeftm,vcov=mvvcovtm,at=data$tm,cen = MMT)
plot(tm.cpn,"overall",ylab="YLL Rate(per 10^5 populations)",
     ylim=c(-2,18),col="red",lwd=4,xlab="Temperature",cex.lab=1.3)

#Get the attributable YLL rate of the corresponding TM
tmyllr <- data.frame(tm=tm.cpn$predvar,tmyllr=tm.cpn$allfit,
                    tmyllrL=tm.cpn$alllow,tmyllrH=tm.cpn$allhigh)
tmyllr$tm <- round(tmyllr$tm,1)
tmyllr <- tmyllr[!duplicated(tmyllr$tm),]
data$tm <- round(data$tm,1)
data <- left_join(data,tmyllr,by="tm")
#show result
 #Note: It is a pity that this example has no statistical significance
apply(data[c("tmyllr","tmyllrL","tmyllrH")], 2, sum)




###(2)Analysis of TV ####
mvtv <- mvmeta(coeftv,vcovtv,method = "reml")

##extact coef and vcov after meta analysis
mvcoeftv <- coef(mvtv)
mvvcovtv <- vcov(mvtv)

##get the cumulative curve
#Generate temperature variability series
perctv <- quantile(data$tv,c(0:1000)/1000,na.rm=T)

#Set parameters and temperature range
tvknots2 <-  quantile(data$tv,c(10,50,90)/100,na.rm=TRUE)
boundv <- range(data$tv,na.rm = T)

#genarate one-basis
argvar1v <- list(fun=varfun,degree=vardegree,
                knots=tvknots2,Boundary.knots=boundv)
cb2v <- crossbasis(perctv,lag=lag,argvar=argvar1v,
                  arglag=arglag)
bvarv <- do.call("onebasis",c(list(x=perctv),attr(cb2v,"argvar")))

#make predictions based on the results of meta analysis
      #Note:The curve is inconsistent with the results 
      #     in the article because of limited locations
tv.cp <- crosspred(bvarv,coef=mvcoeftv,vcov=mvvcovtv,at=perctv)

#Get MMTv and MMPv
MMTv <- round(as.numeric(names(which.min(tv.cp$allfit))),2)
MMTv

MMPv <- names(which(round(quantile(data$tv,c(0:1000)/1000,na.rm = T),2)==MMTv))
MMPv

#Set MMT and re-predict
tv.cpn <- crosspred(bvarv,coef=mvcoeftv,vcov=mvvcovtv,at=data$tv,cen = MMTv)
plot(tv.cpn,"overall",ylab="YLL Rate(per 10^5 populations)",
     ylim=c(-8,15),col="red",lwd=4,xlab="Temperature variability",cex.lab=1.3)

#Get the attributable YLL rate of the corresponding TV
tvyllr <- data.frame(tv=tv.cpn$predvar,tvyllr=tv.cpn$allfit,
                     tvyllrL=tv.cpn$alllow,tvyllrH=tv.cpn$allhigh)
tvyllr$tv <- round(tvyllr$tv,1)
tvyllr <- tvyllr[!duplicated(tvyllr$tv),]
data$tv <- round(data$tv,1)
data <- left_join(data,tvyllr,by="tv")
#show result
 #Note: It is a pity that this example has no statistical significance
apply(data[c("tvyllr","tvyllrL","tvyllrH")], 2, sum)
