# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(ggplot2)
library(doBy)
library(gridExtra)

#################################################

envDF = read.table('01-EncFac.txt', header = T, sep = '\t', stringsAsFactors = F)
envDF = na.omit(envDF)

vfaDF = read.table('01-VFAs.txt', header = T, sep = '\t', stringsAsFactors = F)
vfaDF = na.omit(vfaDF)

#################################################

DF.mean = summaryBy(.~EnvFac + EnvFacLong + TimeGroup + GP, data = envDF, FUN = mean)
DF.mean$lngroup = paste0(DF.mean$EnvFac, '_',DF.mean$GP)

DF.mean$EnvFacLong = factor(DF.mean$EnvFacLong, 
                        levels = c("Sulfate concentration (mM)",
                                   "Sulfite concentration (mM)",
                                   "Total Fe concentration (mg/L)",
                                   "TOC concentration (mg/L)",
                                   "TIC concentration (mg/L)",
                                   "Total phosphorus (mM)","pH"))

DF.sd = summaryBy(.~EnvFac + EnvFacLong + TimeGroup + GP, data = envDF, FUN = sd)
DF.mean = cbind(DF.mean, SD = round(DF.sd$values.sd, digits = 4))

p1 = ggplot(DF.mean, aes(x = TimeGroup, y = values.mean, group = lngroup, color = GP)) + 
  geom_line(linetype = "solid",  size = 0.3, alpha = 0.3) +
  geom_point(shape = 16,  size = 1, alpha = 0.75) + 
  geom_smooth(method = "lm", se = F, size = 0.5, alpha = 0.75) + 
  geom_errorbar(aes(ymin = values.mean - SD, ymax = values.mean + SD), 
                alpha = 0.5, size = 0.2, width = .2) +
  scale_color_manual(values = c('#E69F00', '#0072B2')) + 
  facet_wrap(.~EnvFacLong, ncol = 4, scales = 'free_y') + 
  labs(x = "", y = "") + theme_bw()+
  theme(panel.grid = element_blank(), 
        strip.text = element_text(size = 8),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 6, hjust = 0.5, vjust = 0.5),
        axis.text.x = element_text(size = 8, color = 'black'),
        legend.position = "none")

# ggsave("03-EnvFac.pdf", p1, units = "in", height = 4, width = 8)

write(c(), '03-Env_lmSmooth.txt')
DF.mean$Time = as.integer(substr(DF.mean$TimeGroup, 2, 3))
for (gp in unique(DF.mean$GP)) {
  DF.mean.f = DF.mean[DF.mean$GP == gp, ]
  for (fac in unique(DF.mean.f$EnvFac)) {
    DF.mean.lm = DF.mean.f[DF.mean.f$EnvFac == fac, ]
    res = summary(lm(DF.mean.lm$values.mean ~ DF.mean.lm$Time))
    write('----------------------', '03-Env_lmSmooth.txt', append = T)
    write(paste0('Group: ', gp), '03-Env_lmSmooth.txt', append = T)
    write(paste0('Factor: ', fac), '03-Env_lmSmooth.txt', append = T)
    write(paste0('r.squared: ', res$r.squared), '03-Env_lmSmooth.txt', append = T)
    write(paste0('p_value: ', res$coefficients[,'Pr(>|t|)'][2]), '03-Env_lmSmooth.txt', append = T)
  }
}

#################################################

DF.mean = summaryBy(.~VFAs + VFAsLong + TimeGroup + GP, data = vfaDF, FUN = mean)
DF.mean$lngroup = paste0(DF.mean$VFAs, '_',DF.mean$GP)

DF.mean$VFAsLong = factor(DF.mean$VFAsLong, 
                          levels = c(unique(DF.mean$VFAsLong)[1],
                                     unique(DF.mean$VFAsLong)[6],
                                     unique(DF.mean$VFAsLong)[2],
                                     unique(DF.mean$VFAsLong)[4],
                                     unique(DF.mean$VFAsLong)[7],
                                     unique(DF.mean$VFAsLong)[5],
                                     unique(DF.mean$VFAsLong)[3]))

DF.sd = summaryBy(.~VFAs + VFAsLong + TimeGroup + GP, data = vfaDF, FUN = sd)
DF.mean = cbind(DF.mean, SD = round(DF.sd$values.sd, digits = 4))

p2 = ggplot(DF.mean, aes(x = TimeGroup, y = values.mean, group = lngroup, color = GP)) + 
  geom_line(linetype = "solid",  size = 0.3, alpha = 0.3) +
  geom_point(shape = 16,  size = 1, alpha = 0.75) + 
  geom_smooth(method = "lm", se = F, size = 0.5, alpha = 0.75) + 
  geom_errorbar(aes(ymin = values.mean - SD, ymax = values.mean + SD), 
                alpha = 0.5, size = 0.2, width = .2) +
  scale_color_manual(values = c('#E69F00', '#0072B2')) + 
  facet_wrap(.~VFAsLong, ncol = 4, scales = 'free_y') + 
  labs(x = "", y = "") + theme_bw()+
  theme(panel.grid = element_blank(), 
        strip.text = element_text(size = 8),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 6, hjust = 0.5, vjust = 0.5),
        axis.text.x = element_text(size = 8, color = 'black'),
        legend.position = "none")

# ggsave("03-VfaFac.pdf", p2, units = "in", height = 4, width = 8)

write(c(), '03-Vfa_lmSmooth.txt')
DF.mean$Time = as.integer(substr(DF.mean$TimeGroup, 2, 3))
for (gp in unique(DF.mean$GP)) {
  DF.mean.f = DF.mean[DF.mean$GP == gp, ]
  for (fac in unique(DF.mean.f$VFAs)) {
    DF.mean.lm = DF.mean.f[DF.mean.f$VFAs == fac, ]
    res = summary(lm(DF.mean.lm$values.mean ~ DF.mean.lm$Time))
    write('----------------------', '03-Vfa_lmSmooth.txt', append = T)
    write(paste0('Group: ', gp), '03-Vfa_lmSmooth.txt', append = T)
    write(paste0('Factor: ', fac), '03-Vfa_lmSmooth.txt', append = T)
    write(paste0('r.squared: ', res$r.squared), '03-Vfa_lmSmooth.txt', append = T)
    write(paste0('p_value: ', res$coefficients[,'Pr(>|t|)'][2]), '03-Vfa_lmSmooth.txt', append = T)
  }
}

#################################################

ggsave("04-MergeEnvVFAs.pdf", grid.arrange(p1, p2, nrow = 2), 
       units = "in", height = 7.5, width = 8)
