
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
# 4.Attributable YLL per death and components analysis
################################################################################

###(1)Attributable YLL per death ####
#Calculate the mortality rate per 100,000 people
data$mor <- data$death/data$pop*100000

#Calculate the loss of life for each death caused by TM
apply(data[c("tmyllr","tmyllrL","tmyllrH")], 2, sum)/sum(data$mor)

#Calculate the loss of life for each death caused by TV
apply(data[c("tvyllr","tvyllrL","tvyllrH")], 2, sum)/sum(data$mor)


###(2)components analysis ####
#Divided into 4 components
data$indtm <- ifelse(data$tm <= quantile(data$tm,0.025),1,
                     ifelse(data$tm > quantile(data$tm,0.025) & data$tm <= MMT,2,
                            ifelse(data$tm > MMT & data$tm <= quantile(data$tm,0.975),3,4)))
#1=extremcold,2=moderatecold,3=moderateheat,4=extremheat

data$indtv <- ifelse(data$tv <= MMTv,1,2)
#1=lowtv,2=hightv


#attributable YLL rate of different composition TM and TV
Atmc <- aggregate(data[c("tmyllr","tmyllrL","tmyllrH")],mor)
Atmc
Atvc <- aggregate(data[c("tvyllr","tvyllrL","tvyllrH")],by=list(comtm=data$indtv),sum)
Atvc


#attributable YLL per death of different composition TM and TV
Atmc[c("tmyllr","tmyllrL","tmyllrH")]/aggregate(data["mor"],by=list(comtm=data$indtm),sum)$mor
Atvc[c("tvyllr","tvyllrL","tvyllrH")]/aggregate(data["mor"],by=list(comtm=data$indtv),sum)$mor


