
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
# 0.data preparation
################################################################################


##load packages
library(dlnm); library(mvmeta) 
library(splines);library(reshape2)
library(lubridate);library(tictoc)
library(dplyr)

##load the dataset
  #note: We randomly selected 20 locations 
  #from the national database as examples
data <- read.csv("./dataset/Sample dataset.csv",stringsAsFactors = F)
data$date <- as.Date(data$date)


##generate variable of temperature variability(TV)
#Define the lag function for dataframe
lagdf <- function(x, k) {
  c(rep(NA, k), x)[1 : length(x)] 
}

#caculate TV in each location
location <- unique(data$site)
tvdata <- c()
for (i in 1:length(location)) {
  print(paste("Calculating TV in",location[i]))
  sub <- subset(data,site==location[i])
  sub <- sub[order(sub$date),]
  sub$lagtmin <- lagdf(sub$tmin,1)
  sub$lagtmax <- lagdf(sub$tmax,1)
  sub$tv <- apply(sub[c("tmin","tmax","lagtmin","lagtmax")],1,sd,na.rm=T)
  tvdata <- c(tvdata,sub$tv)
}
data$tv <- tvdata


##genarate variable of YLL rate
data$yllr <- data$yll/data$pop*100000

##arrange the data as a list of dataset
datalist <- lapply(location, function(x){
  data[data$site==x,]
})

