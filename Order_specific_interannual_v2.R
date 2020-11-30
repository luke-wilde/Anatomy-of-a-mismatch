data <- read.csv("C:/Users/14064/Dropbox/Chapter 2/Data/Data for Analyses/Invertebrate Biomass Data/daily_avgbiomass_size_final.v1.csv"); head(data)
levels(as.factor(data$Year))

data1 <- subset(data, Order != "COLLEMBOLA")
data1 <- subset(data1, Year != 2017)
head(data1)

require(lme4); require(ggplot2); library(dplyr); library(arm)

x <- c("gamlss", "mgcv", "ggplot2", "MuMIn", "arm", "sfsmisc", "ggpubr"); lapply(x, require, character.only = T)

data1$z.avg.Biomass <- rescale(data1$avg.Biomass)
data1$z.sum.Biomass <- rescale(data1$sum.Biomass)

mod.z.size <- lmer(z.avg.Biomass ~ Year + (1+Year|SD), data1); summary(mod.z.size); confint(mod.z.size)
mod.z.abund <- lm(z.sum.Biomass ~ Year, data1); summary(mod.z.abund)
mod.z.advance <- lmer(z.sum.Biomass ~ Year + (1+Year|SD), data1); summary(mod.z.advance); confint(mod.z.advance)

ggplot() + geom_point(aes(Year, sum.Biomass), data1) + geom_smooth(aes(Year, sum.Biomass), method="lm", data1)

mod.size<- lmer(avg.Biomass ~ Year + (1+Year|SD), data1); summary(mod.size); confint(mod.size)
mod.abun <- lm(sum.Biomass ~ Year, data1); summary(mod.abun); confint(mod.abun)
mod.adv <- lmer(sum.Biomass ~ Year + (1+Year|SD), data1); summary(mod.adv); confint(mod.adv)

#mod.p.abund <- lmer(perc.avg.Biomass ~ Year*Order + (1|SD), data); summary(mod.p.abund)
#mod.p.size <- lmer(perc.sum.Biomass ~ Year*Order + (1|SD), data); summary(mod.p.size)

x <- data1 %>%
  group_by(Year, SD) %>%
  summarise(avg.BIOMASS = mean(sum.Biomass))
counted <- as.data.frame(x)
counted

early <- mean(counted[1:4,2]); late <- mean(counted[5:8,2])
late/early

ggplot() + geom_smooth(aes(Year, z.avg.Biomass, group = stage, color = stage), method = "lm", se = F, data1, alpha = 0.3) + geom_smooth(aes(Year, z.avg.Biomass), method = "lm", se = T, alpha = 0.15, color = "black", data1) + theme_classic()
ggplot() + geom_smooth(aes(Year, z.sum.Biomass, group = stage, color = stage), method = "lm", se = F, data1, alpha = 0.3) + geom_smooth(aes(Year, z.sum.Biomass), method = "lm", se = T, alpha = 0.15, color = "black", data1) + theme_classic()


effects <- read.csv("C:/Users/14064/Dropbox/Chapter 2/Data/Output/resource_change_effects.csv")
head(effects)
levels(effects$Test)

sum <- as.data.frame(effects[1:7,])
avg <- as.data.frame(effects[8:14,])
peak <- as.data.frame(effects[15:21,])
obs.peak <- as.data.frame(effects[22:28,])

sum$Order <- factor(sum$Order, levels=rev(unique(sum$Order)))
avg$Order <- factor(avg$Order, levels=rev(unique(avg$Order)))
peak$Order <- factor(peak$Order, levels=rev(unique(peak$Order)))
obs.peak$Order <- factor(obs.peak$Order, levels=rev(unique(obs.peak$Order)))


abundance <- ggplot()+ geom_pointrange(aes(x = Estimate, y = factor(Order), xmin = Lower, xmax = Upper), shape = 1, size = .3, data = sum) + theme_classic() + geom_vline(xintercept = 0, lty = 1, color = "grey", alpha = .5, lwd = .8) + coord_cartesian(xlim = c(-0.2,0.2)) + scale_x_continuous(expand = c(0,0), limits = c(-0.2,0.2), breaks = seq(-0.2,0.2, by = 0.1)) + labs(x = "Standardized Estimate", y = "", title = "Invertebrate Biomass") + theme(axis.text.x = element_text(colour = "grey30", size = 8),  axis.text.y = element_text(colour = "grey30", size = 8), axis.title.x = element_text(colour = "grey30", size = 8), axis.title.y = element_text(colour = "grey30", size = 8), title = element_text(colour = "grey30", size = 7)) + theme(legend.position = "none")





avg.biomass <- ggplot()+ geom_pointrange(aes(x = Estimate, y = factor(Order), xmin = Lower, xmax = Upper), shape = 1, size = .3, data = avg) + theme_classic() + geom_vline(xintercept = 0, lty = 1, color = "grey", alpha = .5, lwd = .8) + coord_cartesian(xlim = c(-0.2,0.2)) + scale_x_continuous(expand = c(0,0), limits = c(-0.2,0.2), breaks = seq(-0.2,0.2, by = 0.1)) + labs(x = "", y = "", title = "Invert. Body Size") + theme(axis.text.x = element_text(colour = "grey30", size = 8),  axis.text.y = element_text(colour = "grey30", size = 8), axis.title.x = element_text(colour = "grey30", size = 8), axis.title.y = element_text(colour = "grey30", size = 8), title = element_text(colour = "grey30", size = 7)) + theme(legend.position = "none")


#ggarrange(abundance, avg.biomass, ncol = 2, common.legend = T)  

#dodge <- position_dodge(width = 0.65)

peakday <- ggplot()+ geom_pointrange(aes(x = Estimate, y = Order, xmin = Lower, xmax = Upper), shape = 1, size = .3, data = peak) + theme_classic() + geom_vline(xintercept = 0, lty = 1, color = "grey", alpha = .5, lwd = .8) + coord_cartesian(xlim = c(-0.2,0.2)) + scale_x_continuous(expand = c(0,0), limits = c(-0.5,0.5), breaks = seq(-0.2,0.2, by = 0.1)) + labs(x = "", y = "", title = "Seasonal Peak") + theme(axis.text.x = element_text(colour = "grey30", size = 8),  axis.text.y = element_text(colour = "grey30", size = 8), axis.title.x = element_text(colour = "grey30", size = 8), axis.title.y = element_text(colour = "grey30", size = 8), title = element_text(colour = "grey30", size = 7)) + theme(legend.position = "none")

ggarrange(peakday, abundance, avg.biomass, ncol = 3, common.legend = T)  
 
# new.data <- read.csv("C:/Users/14064/Dropbox/Chapter 2/Data/Data for Analyses/Invert Assemblage Analysis/peak_interannual_change.csv")
# old.data <- read.csv("C:/Users/14064/Dropbox/Chapter 2/Data/Data for Analyses/Invert Assemblage Analysis/shape_trends.csv")
# mod.peak.change <- lm(peak ~ (year + I(year^2)) *order, data = new.data); summary(mod.peak.change)
# old.data$zPeak <- rescale(old.data$Peak)
# mod.peak.shape.change <- lm(zPeak ~ (Year + I(Year^2))*Order, data = old.data); summary(mod.peak.shape.change); confint(mod.peak.shape.change)


require(ggpubr)

jpeg(filename = "C:/Users/14064/Dropbox/Chapter 2/Final Materials/Figures/Figure_interannual_changes_v4.jpeg", res = 600, width = 6, height = 4, units = 'in')
ggarrange(peakday, abundance, avg.biomass, ncol = 3, common.legend = T)  
dev.off()

