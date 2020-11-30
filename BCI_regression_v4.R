# script to run linear regression on chick BCI values

#load packages
x <- c("lme4", "MuMIn", "arm", "ggplot2", "mgcv::gamlss", "mgcv", "ggpubr", "visreg", "plyr", "gridExtra")
lapply(x, require, character.only = TRUE)

#load the dataset
data <- read.csv("C:/Users/14064/Dropbox/Chapter 2/Data/Data for Analyses/BCI Regression/BCI_regression_v2.csv")
head(data)


data <- na.omit(data)

data$YEAR <- as.factor(data$YEAR)

data$YEAR

data$BROOD <- as.factor(data$BROOD)
data$CODE <- as.factor(data$CODE)


vif.gam <- function(object){
  
  obj.sum <- mgcv::summary.gam(object)
  
  s2 <- object$sig2 # estimate of standard deviation of residuals
  X <- object$model # data used to fit the model
  n <- nrow(X) # how many observations were used in fitting?
  v <- -1 # omit the intercept term, it can't inflate variance
  varbeta <- obj.sum$p.table[v,2]^2 # variance in estimates
  selected_col <- row.names(obj.sum$p.table)[v]
  selected_col <- gsub("TRUE", "", selected_col)
  varXj <- apply(X=X[, selected_col],MARGIN=2, var) # variance of all the explanatory variables
  VIF <- varbeta/(s2/(n-1)*1/varXj) # the variance inflation factor, obtained by rearranging
  # var(beta_j) = s^2/(n-1) * 1/var(X_j) * VIF_j
  
  VIF.df <- tibble::tibble(variable=names(VIF),
                           vif=VIF)
  
  return(VIF.df)
}


sd(data$Biomass.7); sd(data$Hatch.Date)

31

#need to first rescale all predictors
data$zAge <- rescale(data$Age)
data$zBCI.fixed <- rescale(data$BCI.fixed)
data$zHatch <- rescale(data$Hatch.Date)
data$zPreySize <- rescale(data$PreySize)
data$zBiomass <- rescale(data$Biomass)
data$zmatch <- rescale(data$match)

data$zPreySize.7 <- rescale(data$PreySize.7)
data$zBiomass.7 <- rescale(data$Biomass.7)

data$zPreySize.3 <- rescale(data$PreySize.3)
data$zBiomass.3 <- rescale(data$Biomass.3)

data$zPreySize.1 <- rescale(data$PreySize.1)
data$zBiomass.1 <- rescale(data$Biomass.1)


#test models with and without random effects
mod <- mgcv::gam(BCI.fixed ~ zHatch +  zPreySize.7 + zBiomass.7, family = "gaussian", data = data, na.action = na.fail)
mod.ryear <- mgcv::gam(BCI.fixed ~ zHatch +  zPreySize.7 + zBiomass.7 + s(YEAR, bs = "re"), family = "gaussian", data = data, na.action = na.fail)
mod.ryear.rnest <- mgcv::gam(BCI.fixed ~ zHatch +  zPreySize.7 + zBiomass.7 + s(YEAR, bs = "re") + s(BROOD, bs = "re"), family = "gaussian", data = data, na.action = na.fail)
mod.rnest <- mgcv::gam(BCI.fixed ~ zHatch +  zPreySize.7 + zBiomass.7 + s(BROOD, bs = "re"), family = "gaussian", data = data, na.action = na.fail)
mod.ryear.rindiv <- mgcv::gam(BCI.fixed ~ zHatch +  zPreySize.7 + zBiomass.7 + s(YEAR, bs = "re") + s(CODE, bs = "re"), family = "gaussian", data = data, na.action = na.fail)

mod.ryear.rindiv.rnest <- mgcv::gam(BCI.fixed ~ zHatch +  zPreySize.7 + zBiomass.7 + s(YEAR, bs = "re") + s(CODE/BROOD, bs = "re"), family = "gaussian", data = data, na.action = na.fail)

AICc(mod.rnest)

sapply(list(mod, mod.ryear, mod.ryear.rnest, mod.ryear.rindiv, mod.rnest, mod.ryear.rindiv.rnest), AICc)

runif(7, 80.6, 105)

#so now will go on running the others with nested and year as additive

mod.day <- lm(BCI.fixed ~ zPreySize + zBiomass, data)
mod.day.1 <- lmer(BCI.fixed ~ zPreySize.1 + zBiomass.1 +  (1|YEAR), data)
mod.day.3 <- lm(BCI.fixed ~ zPreySize.3 + zBiomass.3, data)
mod.day.7 <- lm(BCI.fixed ~ zBiomass.7 + zPreySize.7, data)

sapply(list(mod.day, mod.day.1, mod.day.3, mod.day.7), AIC)



####################
#Begin model tests
##################
global.0 <- mgcv::gam(BCI.fixed ~ zHatch + PreySize.7 + Biomass.7 + s(zAge), family = "gaussian", data = data, na.action = na.fail); summary(global.0)

dd<-dredge(global.0, trace = 2)
data.frame(subset(dd, delta<300))
model.avg(dd, subset = delta < 4)
avgmod.95p <- model.avg(dd, delta<4)
summary(avgmod.95p);confint(avgmod.95p)
summary(get.models(dd, 1)[[1]])

global.1 <- mgcv::gam(BCI.fixed ~ zBiomass.7 + zHatch + s(zAge) + s(Year, bs = "re") + s(Brood, bs = "re"), family = "gaussian", data, na.action = na.fail); summary(global.1); AICc(global.1); confint(global.1)

#the top model
options(scipen = 999)

top.mod <- mgcv::gam(BCI.fixed ~ Hatch.Date + Biomass.7 + s(Age), family = "gaussian", data = data, na.action = na.fail); summary(top.mod); AICc(top.mod); confint(top.mod)

min(data$Biomass.7);max(data$Biomass.7)

match <- visreg(global.0, "zmatch", gg=T, alpha = 0.1, line=list(col="black", size = 2), fill=list(fill="grey", alpha=0.3), points=list(size=2, pch=1)) + theme_classic() + geom_hline(aes(yintercept = 1), lty = 2, alpha = 0.6); match

match_group1 <- factor(data$match_group, levels=c("> -10", "-10 to -7", "-7 to -1", "0", "1 to 7", "> 10"))


mod <- lm(BCI.fixed ~ match_group, data); summary(mod)

eff<- as.data.frame(summary(mod)$coefficients[, 1])
se <- as.data.frame(summary(mod)$coefficients[, 2])
as.data.frame(rbind(eff[, 1], se[, 1]))

ggplot() + geom_smooth(aes(match, BCI.fixed), data, method = "gam") + geom_point(aes(match, BCI.fixed, group = YEAR, color = YEAR), data) + coord_cartesian(ylim=c(0.25,1.75), xlim = c(-17,23)) + theme(axis.text.x = element_text(colour = "grey30", size = 10),  axis.text.y = element_text(colour = "grey30", size = 10), axis.title.x = element_text(colour = "grey30", size = 12), axis.title.y = element_text(colour = "grey30", size = 12)) + scale_x_continuous(expand = c(0,0), limits = c(-17,23), breaks = seq(-17, 23, by = 5)) + scale_y_continuous(expand = c(0,0), limits = c(0.25,1.75), breaks = seq(0.25,1.75, by = .25)) + theme_classic() + geom_hline(aes(yintercept = 1), lty = 2, alpha = 0.6) + geom_vline(aes(xintercept = 1), lty = 1, size = 1.5, alpha = 0.3, color = "grey") + labs(x = "Match (days)", y = "Body Condition Index")



#ggplot() + geom_boxplot(aes(match_group1, BCI.fixed), data = data) + theme_classic() + geom_hline(aes(yintercept = 1), lty = 2, alpha = 0.6) + labs(x = "Match (days)", y = "Body Condition Index")


biomass <- visreg(top.mod, "Biomass.7", gg=T, alpha = 0.1, line=list(col="black", size = 2), fill=list(fill="grey", alpha=0.3), points=list(size=2, pch=1)) + theme_classic() + geom_hline(aes(yintercept = 1), lty = 2, alpha = 0.6) + labs(x = "Invertebrate Biomass (mg)", y = "Body Condition Index", title = "a") + coord_cartesian(ylim=c(0.25,1.75), xlim = c(22, 138)) + theme(axis.text.x = element_text(colour = "grey30", size = 10),  axis.text.y = element_text(colour = "grey30", size = 10), axis.title.x = element_text(colour = "grey30", size = 12), axis.title.y = element_text(colour = "grey30", size = 12)) + scale_x_continuous(expand = c(0,0), limits = c(24, 136), breaks = seq(24, 136, by = 14)) + scale_y_continuous(expand = c(0,0), limits = c(0.25,1.75), breaks = seq(0.25,1.75, by = .25)) ; biomass



plot <- ggplot() + geom_smooth(formula = y ~ x + I(x^2), aes(Age, BCI.fixed), method = "lm", data = data, color = "black", fill = "grey") + labs(x = "Invertebrate Biomass (mg)", y = "Body Condition Index") + coord_cartesian(ylim=c(0.25,1.75), xlim = c(22, 138)) + theme(axis.text.x = element_text(colour = "grey30", size = 10),  axis.text.y = element_text(colour = "grey30", size = 10), axis.title.x = element_text(colour = "grey30", size = 12), axis.title.y = element_text(colour = "grey30", size = 12)) + scale_x_continuous(expand = c(0,0), limits = c(24, 136), breaks = seq(24, 136, by = 14)) + scale_y_continuous(expand = c(0,0), limits = c(0.25,1.75), breaks = seq(0.25,1.75, by = .25)) + theme_classic() + geom_hline(aes(yintercept = 1), lty = 2, alpha = 0.6)

plot



min(data$Hatch.Date);max(data$Hatch.Date)
hatch <- visreg(global.0, "PreySize.7", gg=T, alpha = 0.1, line=list(col="white", size = 0, alpha=0.0000001), fill=list(fill="white", alpha=0.0000001), points=list(size=2, pch=1)) + theme_classic() + geom_hline(aes(yintercept = 1), lty = 2, alpha = 0.6) + labs(x = "Median Invert. Body Size (mg)", y = "", title = "b") + coord_cartesian(ylim=c(0.25,1.75), xlim = c(1.9,4.75)) + theme(axis.text.x = element_text(colour = "grey30", size = 10),  axis.text.y = element_text(colour = "grey30", size = 10), axis.title.x = element_text(colour = "grey30", size = 12), axis.title.y = element_text(colour = "grey30", size = 12)) + scale_x_continuous(expand = c(0,0), limits = c(1.9,4.8), breaks = seq(1.9,4.7, by = .4)) + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) ; hatch

#+ theme(axis.line.y = element_blank())

jpeg(filename = "C:/Users/14064/Dropbox/Chapter 2/Final Materials/Figures/Fig2_v5.jpeg", res = 600, width = 6, height = 4, units = 'in')
ggarrange(biomass,hatch)
dev.off()
