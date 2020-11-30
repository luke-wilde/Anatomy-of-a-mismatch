#Script for godwit chick growth curves

##First, make sure you can calculate the age of the chick at each time point. To do this, you need to subtract their hatch date from their 'banding' date, which is the date on which they were captured each time

setwd("C:/Users/14064/Dropbox/Chapter 2/Data/Input Files/")
x <- c("nlme", "AICcmodavg", "ggplot2", "dplyr", "ggthemes")
lapply(x, require, character.only = TRUE)

#load the data
chick<-read.csv(file="HUGO_growth.csv") 
head(chick)

max(chick$Mass)

#only include chicks that are at least 1 day old
chick <- chick[chick$Include == 1,]
chick$Year <- as.numeric(chick$Year)
chick <- chick[chick$Year != 2014,]

#make sure that year is factor
chick$Year <- as.factor(chick$Year)

head(chick)
unique(chick$Code)

levels(chick$Year)

#view the data, visual inspection
plot(Mass ~ day, data=chick)

#can use 
fit <- nls(Mass ~ SSlogis(day, Asym, xmid, scal), data = chick); summary(fit)

curve(predict(fit, newdata = data.frame(day=x)), add=TRUE)


chick$hdate<-as.numeric(as.Date(as.character(chick$Hatching.Date), format="%m/%d/%y")) #make hatch date a number for subtraction
chick$bdate<-as.numeric(as.Date(as.character(chick$Banding.Date), format="%m/%d/%y")) #make measuring date a number
chick$day<-(chick$bdate-chick$hdate) #subtract the two for the age in days


##Next, make sure that each chick has a unique identifier.

chick<-chick[!is.na(chick$Mass),]
chick$fband<-as.factor(chick$Code)

options(scipen = 999)

head(chick)

##Then, this is the logistic and gompertz equations we created for Beluga, where we fixed the asymptotic mass to the average mass of all adults (i.e., not differentiating between males and females):
##Then you need to group the data:
chickg<-groupedData(Mass~day|Species, chick, order.groups=F) #nlme required grouped Data for whatever reason

# use an lm to determine how many coeffecients you will need - with these data you only need 8 since there are 6 K*Year combinations (since 2014 was excluded)
test <- length(coef(lm(Mass ~ day + Year, data = chick))); test

#body mass curves

#Build all candidate models: logistic
# 
#   la0.b0.b1<-nlme(Mass~249/(1+exp(-b0*(age-b1))), data=chickg, fixed=list(b0+b1~1), random=b0+b1~1|Code,   na.action=na.omit, start=c(b0=0.12, b1=14.04), control=list(maxIter=500))
#   
# 
# 
# #now do with gompertz equation
#   ga0.b0.b1<-nlme(Mass~249*exp(-exp(-b0*(age-b1))), data=chickg, fixed=list(b1~1, b0~Year), random=b0+b1~1|Code,   na.action=na.omit, start=c(b0=0.17,0,0,0,0,0, b1=11.465), control=list(maxIter=500))
#  
#   
# #from Loonstra et al. 2017 - We included chickID as a random effect, to account for pseudoreplication (Pinheiro & Bates 2000). However, models including a random effect for all three different growth parameters (y???, k and Ti) did not converge. A closer examination of the correlation between the estimated random effects revealed that they were highly correlated and that the model was overfitted (negative variances of the random effects). Exploratory analyses showed that the convergence problems were solved when individuals were only allowed to vary randomly for asymptotic size (y???). We therefore decided to only include a random effect for the asymptotic growth parameter (y???). We then tested, using both growth functions (Gompertz and logistic), for an effect of sex on y???, k and Ti.  
  
# logistic function
logisfix<-deriv(~a0/(1+exp(-K*(day-i))), c("a0","i","K"), function(a0,day,i,K) {})


# test candidate models
logfix1<-nlme(Mass~logisfix(a0,day,i,K), chickg, start=c(249, 10.7, .12), fixed=c(a0~1,i~1,K~1), random=a0~1|Code, verbose=T, control=list(maxIter=100))

logfix2<-nlme(Mass~logisfix(a0,day,i,K), chickg, start=c(249, 10.7,10.7,10.7,10.7,10.7,10.7, .12), fixed=c(a0~1,i~Year,K~1), random=a0~1|Code, verbose=T, control=list(maxIter=100))

logfix3<-nlme(Mass~logisfix(a0,day,i,K), chickg, start=c(249, 10.7, .12,.12,.12,.12,.12,.12), fixed=c(a0~1,i~1,K~Year), random=a0~1|Code, verbose=T, control=list(maxIter=100))

logfix4<-nlme(Mass~logisfix(a0,day,i,K), chickg, start=c(249, 10.7,10.7,10.7,10.7,10.7,10.7, .12,.12,.12,.12,.12,.12), fixed=c(a0~1,i~Year,K~Year), random=a0~1|Code, verbose=T, control=list(maxIter=100))

logfix5<-nlme(Mass~logisfix(a0,day,i,K), chickg, start=c(249,249,249,249,249,249, 10.7,.12), fixed=c(a0~Year,i~1,K~1), random=a0~1|Code, verbose=T, control=list(maxIter=100)); summary(logfix5); AIC(logfix5)

logfix6<-nlme(Mass~logisfix(a0,day,i,K), chickg, start=c(249,249,249,249,249,249, 10.7,10.7,10.7,10.7,10.7,10.7, .12), fixed=c(a0~Year,i~Year,K~1), random=a0~1|Code, verbose=T, control=list(maxIter=100))

logfix7<-nlme(Mass~logisfix(a0,day,i,K), chickg, start=c(249,249,249,249,249,249, 10.7, .12,.12,.12,.12,.12,.12), fixed=c(a0~Year,i~1,K~Year), random=a0~1|Code, verbose=T, control=list(maxIter=100))

sapply(list(logfix1,logfix2,logfix3,logfix4, logfix5, logfix6, logfix7), AIC)

summary(logfix4)

Year = c(2009,2010,2011,2015,2016,2019)
a0 <- c(165.32896,165.32896,165.32896,165.32896,165.32896,165.32896)
i <- c(11.46413, 12.23624,13.16011,11.23496,13.64994,9.06131)
k <- c(0.17215, 0.17578,0.1672,0.16872,0.13724,0.23168)

curves <- data.frame("Year" = Year, "a" = a0, "Ti" = i, "k" = k)
curves
# logistic with fixed assymptote

logfix.assym <- deriv(~249/(1+exp(-K*(day-i))), c("i","K"), function(day,i,K) {})
#
#
logfix.a1<-nlme(Mass~logfix.assym(day,i,K), chickg, start=c(10.7, .12), fixed=c(i~1,K~1), random=i~1|Code, verbose=T, control=list(maxIter=100)); summary(logfix.a1); AIC(logfix.a1) 

logfix.a2<-nlme(Mass~logfix.assym(day,i,K), chickg, start=c(10.7, .12), fixed=c(i~1,K~1), random=i+K~1|Code, verbose=T, control=list(maxIter=100)); summary(logfix.a2); AIC(logfix.a2)

logfix.a3<-nlme(Mass~logfix.assym(day,i,K), chickg, start=c(10.7, .12), fixed=c(i~1,K~1), random=K~1|Code, verbose=T, control=list(maxIter=100)); summary(logfix.a3); AIC(logfix.a3)

logfix.a4<-nlme(Mass~logfix.assym(day,i,K), chickg, start=c(10.7,10.7,10.7,10.7,10.7,10.7, .12), fixed=c(i~Year,K~1), random=i~1|Code, verbose=T, control=list(maxIter=100)); summary(logfix.a4); AIC(logfix.a4)

logfix.a5<-nlme(Mass~logfix.assym(day,i,K), chickg, start=c(10.7,10.7,10.7,10.7,10.7,10.7, .12), fixed=c(i~Year,K~1), random=K~1|Code, verbose=T, control=list(maxIter=100)); summary(logfix.a5); AIC(logfix.a5)

logfix.a6<-nlme(Mass~logfix.assym(day,i,K), chickg, start=c(10.7,10.7,10.7,10.7,10.7,10.7, .12), fixed=c(i~Year,K~1), random=i+K~1|Code, verbose=T, control=list(maxIter=100)); summary(logfix.a6); AIC(logfix.a6)

logfix.a7<-nlme(Mass~logfix.assym(day,i,K), chickg, start=c(10.7,.12,.12,.12,.12,.12,.12), fixed=c(i~1,K~Year), random=K~1|Code, verbose=T, control=list(maxIter=100)); summary(logfix.a7); AIC(logfix.a7)

logfix.a8<-nlme(Mass~logfix.assym(day,i,K), chickg, start=c(10.7,.12,.12,.12,.12,.12,.12), fixed=c(i~1,K~Year), random=i~1|Code, verbose=T, control=list(maxIter=100)); summary(logfix.a8); AIC(logfix.a8)

logfix.a9<-nlme(Mass~logfix.assym(day,i,K), chickg, start=c(10.7,.12,.12,.12,.12,.12,.12), fixed=c(i~1,K~Year), random=i+K~1|Code, verbose=T, control=list(maxIter=100)); summary(logfix.a9); AIC(logfix.a9)

logfix.a10<-nlme(Mass~logfix.assym(day,i,K), chickg, start=c(10.7,10.7,10.7,10.7,10.7,10.7,.12,.12,.12,.12,.12,.12), fixed=c(i~Year,K~Year), random=i~1|Code, verbose=T, control=list(maxIter=100)); summary(logfix.a10); AIC(logfix.a10)

logfix.a11<-nlme(Mass~logfix.assym(day,i,K), chickg, start=c(10.7,10.7,10.7,10.7,10.7,10.7,.12,.12,.12,.12,.12,.12), fixed=c(i~Year,K~Year), random=K~1|Code, verbose=T, control=list(maxIter=100)); summary(logfix.a11); AIC(logfix.a11)

logfix.a12<-nlme(Mass~logfix.assym(day,i,K), chickg, start=c(10.7,10.7,10.7,10.7,10.7,10.7,.12,.12,.12,.12,.12,.12), fixed=c(i~Year,K~Year), random=i+K~1|Code, verbose=T, control=list(maxIter=100)); summary(logfix.a12); AIC(logfix.a12)

list <- sapply(list(logfix.a1,logfix.a2,logfix.a3,logfix.a4,logfix.a5,logfix.a7,logfix.a8,logfix.a10,logfix.a11), AICc)

list - AICc(logfix.a2)



summary(logfix4)
summary(logfix.a2)



# write gompertz function
gompfix<-deriv(~249*exp(-exp(-K*(day-i))), c("i","K"), function(day,i,K) {})

#now construct gompertz candidate models
gompfix1<-nlme(Mass~ gompfix(day,i,K), chickg, start=c(10.7, .12), fixed=c(i~1,K~1), random=i~1|Code, verbose=T, control=list(maxIter=100)); summary(gompfix1); AIC(gompfix1)

gompfix5<-nlme(Mass~ gompfix(day,i,K), chickg, start=c(10.7, .12), fixed=c(i~1,K~1), random=K~1|Code, verbose=T, control=list(maxIter=100)); summary(gompfix5); AIC(gompfix5)

gompfix2<-nlme(Mass~ gompfix(day,i,K), chickg, start=c(10.7,10.7,10.7,10.7,10.7,10.7, .12), fixed=c(i~Year,K~1), random=i~1|Code, verbose=T, control=list(maxIter=100)); summary(gompfix2); AIC(gompfix2)

gompfix3<-nlme(Mass~ gompfix(day,i,K), chickg, start=c(10.7,10.7,10.7,10.7,10.7,10.7, .12), fixed=c(i~Year,K~1), random=K~1|Code, verbose=T, control=list(maxIter=100)); summary(gompfix3); AIC(gompfix3)

gompfix4<-nlme(Mass~ gompfix(day,i,K), chickg, start=c(10.7, .12,.12,.12,.12,.12,.12), fixed=c(i~1,K~Year), random=K~1|Code, verbose=T, control=list(maxIter=100)); summary(gompfix4); AIC(gompfix4)

gompfix6<-nlme(Mass~ gompfix(day,i,K), chickg, start=c(10.7,10.7,10.7,10.7,10.7,10.7, .12,.12,.12,.12,.12,.12), fixed=c(i~Year,K~Year), random=i~1|Code, verbose=T, control=list(maxIter=100)); summary(gompfix6); AIC(gompfix6)

summary(logfix4)

#visualize

curves <- read.csv("C:/Users/14064/Dropbox/Chapter 2/Data/Output/growth_curve_estimates_v1.csv", fileEncoding="UTF-8-BOM")
head(curves)

head(chick)

curves$Year <- as.factor(curves$Year)
chick$Year <- as.factor(chick$Year)

p <- ggplot(data = chick, mapping = aes(x = x)) 

fun.1 <- function(x) (newdata$a/(1+exp((-1*k)*(x-Ti))))

p + stat_function(fun = fun.1, geom="line")

cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00")
  
ggplot() + geom_line(data = curves, aes(day, exp.mass, color = year)) + theme_classic() + geom_point(data = chick, aes(day, Mass, color = Year)) + xlab("Age") + ylab("Mass") + theme(axis.text.x = element_text(colour = "grey30", size = 12),  axis.text.y = element_text(colour = "grey30", size = 12), axis.title.x = element_text(colour = "grey30", size = 16), axis.title.y = element_text(colour = "grey30", size = 16)) + scale_y_continuous(expand = c(0,0))+ coord_cartesian(ylim=c(0,200), xlim=c(0,28))+ scale_x_continuous(expand = c(0,0)) + scale_fill_manual(values= cbp2)
  


#
new <- right_join(curves, chick, by = c("Year","day"))
new

new$BCI <- (new$Mass - new$exp.mass)

write.csv(new, "C:/Users/14064/Dropbox/Chapter 2/Data/Output/BCI_estimates.csv")
