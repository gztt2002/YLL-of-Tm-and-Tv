
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
# 2.first-stage analysis
################################################################################

##create the objects to store the results
#coefficients and vocv for each location
coeftm <- coeftv <- matrix(NA,length(datalist),length(varper)+vardegree,
                           dimnames=list(location,paste("b",seq(5),sep="")))

vcovtm <- vcovtv <- vector("list",length(datalist))
names(vcovtm) <- names(vcovtv) <- location



##run the loop
tic("Time to complete") #Used to record the running time

for(i in seq(length(datalist))) {
  
  # print
  print(paste(i,location[i],sep = "-"))
  
  # extract single location data
  sub <- datalist[[i]]
  sub$tm <- round(sub$tm,1)
  sub$tv <- round(sub$tv,1)
  
  # define cross-basis of TM
  argvar <- list(fun=varfun,degree=vardegree,
                 knots=quantile(sub$tm,varper/100,na.rm=T))
  cb <- crossbasis(sub$tm,lag=lag,argvar=argvar,
                     arglag=arglag)
  #summary(cb)
  
  # define cross-basis of TV
  argvarv <- list(fun=varfun,degree=vardegree,
                  knots=quantile(sub$tv,varper/100,na.rm=T))
  av <- crossbasis(sub$tv,lag=lag,argvar=argvarv,
                     arglag=arglag)
  #summary(av)

  # run the model
  model <- glm(formula,family=gaussian,sub,na.action="na.exclude")
  
  # reduction to overall cumulative
    #note: The setting of cen does not affect coef and vcov
  redtm <- crossreduce(cb,model,cen=quantile(sub$tm,cenper/100,na.rm=T))
  coeftm[i,] <- coef(redtm)
  vcovtm[[i]] <- vcov(redtm)
  
  redtv <- crossreduce(av,model,cen=quantile(sub$tv,cenper/100,na.rm=T))
  coeftv[i,] <- coef(redtv)
  vcovtv[[i]] <- vcov(redtv)
}

toc()

