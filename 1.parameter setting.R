
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
#     https://github.com/gztt2002/YLL-of-Tm-and-Tv
################################################################################

################################################################################
# 1.parameters setting
################################################################################

##parameters in crossbasis
#exposure function
varfun = "bs"
vardegree = 2
varper <- c(10,50,90) 
cenper <- 75

#lag function
lag <- 21
arglag <- list(fun="ns",knots=logknots(21,df=4,int=T))


##parameters in model
#degree of seasonality
dft <- 7

#define percentiles of TM
per <- t(sapply(datalist,function(x) 
  quantile(x$tm,c(2.5,10,25,50,75,90,97.5)/100,na.rm=T)))
per

#define percentiles of TV
perv <- t(sapply(datalist,function(x) 
  quantile(x$tv,c(2.5,10,25,50,75,90,97.5)/100,na.rm=T)))
perv


##model formula
formula <- yllr ~ cb+av+dow+PM10+ns(t,df=dft*length(unique(year)))+ ns(rh,3)


