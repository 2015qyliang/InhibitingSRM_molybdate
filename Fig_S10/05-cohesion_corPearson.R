# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(tidyverse)
library(reshape2)
library(ggplot2)
library(grid)
library(doBy)

#####################################################

fnscoh = list.files(pattern = '03-cohesion_')

values.cohesion = Group = Type = c()
for (fn in fnscoh) {
  gp = substr(fn, 13, 16)
  dfCoh = readRDS(fn)
  pos.cohes = dfCoh$`Positive Connectedness`
  neg.cohes = dfCoh$`Negative Connectedness`
  values.cohesion = append(values.cohesion, c(pos.cohes, neg.cohes))
  Type = append(Type, c(rep('Positive', length(pos.cohes)), 
                        rep('Negative', length(neg.cohes))))
  Group = append(Group, rep(gp, length(c(pos.cohes, neg.cohes))))
}

DFcohesion = data.frame(Group, Type, values.cohesion)
DFcohesion$GP = substr(DFcohesion$Group, 4, 4)
DFcohesion$Time = as.integer(substr(DFcohesion$Group, 2, 3))

DFcohesion.m = DFcohesion[DFcohesion$Group == 'D00C', ]
DFcohesion.m$Group = 'D00M'
DFcohesion.m$GP = 'M'
DFcohesion = rbind(DFcohesion, DFcohesion.m)

DFcohesion = DFcohesion[, c('Group', 'values.cohesion', 'Type')]

sumDF = summaryBy(values.cohesion~Group+Type, data = DFcohesion, FUN = mean)

sumDF = sumDF[!(sumDF$Group == 'D00M'), ]

sumDF.r = dcast(data = sumDF, Type ~ Group, value.var = 'values.cohesion.mean')

sumDF.r = sumDF.r[, 3:ncol(sumDF.r)]

#####################################################

rawdf = read.table('04-complexity+stability.txt', header = T, sep = '\t', 
                   stringsAsFactors = F, row.names = 1)

rawdf = rbind(rawdf, sumDF.r)
rownames(rawdf)[12:13] = c('negative.c', 'positive.c')

rawdf.c = t(rawdf[, str_detect(colnames(rawdf), 'C')])
rawdf.m = t(rawdf[, str_detect(colnames(rawdf), 'M')])

cc = cor(rawdf.c)
cm = cor(rawdf.m)

cc = as.data.frame(cc[1:11, 12:13])
cm = as.data.frame(cm[1:11, 12:13])

cc$complex = rownames(cc)
cm$complex = rownames(cm)

cc.r = melt(cc, variable.name = "Stability", measure.name = colnames(cc)[1:2],  
            value.name ="cor")
cc.r$cor = round(cc.r$cor, 2)
cm.r = melt(cm, variable.name = "Stability", measure.name = colnames(cm)[1:2],  
            value.name ="cor")
cm.r$cor = round(cm.r$cor, 2)

cc.r$complex = factor(cc.r$complex, levels = rev(rownames(cc)))
cm.r$complex = factor(cm.r$complex, levels = rev(rownames(cm)))

cc.r$text = cc.r$cor
cc.r$text[abs(cc.r$text) < 0.7] = NA
pc = ggplot(cc.r, aes(Stability, complex, fill = cor))+
  geom_raster(color = "white")+
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",  
                       midpoint = 0, limit = c(-1,1), space = "Lab") +
  geom_text(aes(Stability, complex, label = text), color = "black", size = 2) +
  theme( panel.background = element_blank(),
         panel.grid = element_blank(), 
         axis.title = element_blank(), 
         axis.ticks = element_blank(), 
         axis.text.x = element_text(size = 8, angle = 30, hjust = 1, vjust = 1, color = 'black'),
         axis.text.y = element_text(size = 8, hjust = 1, vjust = 0.5, color = 'black'),
         legend.key.width = unit(.07, 'in'),
         legend.key.height = unit(.2, 'in'),
         legend.position = "right",
         legend.background = element_blank(),
         legend.key = element_blank())

cm.r$text = cm.r$cor
cm.r$text[abs(cm.r$text) < 0.7] = NA
pm = ggplot(cm.r, aes(Stability, complex, fill = cor))+
  geom_raster(color = "white")+
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",  
                       midpoint = 0, limit = c(-1,1), space = "Lab") +
  geom_text(aes(Stability, complex, label = text), color = "black", size = 2) +
  theme( panel.background = element_blank(),
         panel.grid = element_blank(), 
         axis.title = element_blank(), 
         axis.ticks = element_blank(), 
         axis.text.x = element_text(size = 8, angle = 30, hjust = 1, vjust = 1, color = 'black'),
         axis.text.y = element_text(size = 8, hjust = 1, vjust = 0.5, color = 'black'),
         legend.key.width = unit(.07, 'in'),
         legend.key.height = unit(.2, 'in'),
         legend.position = "right",
         legend.background = element_blank(),
         legend.key = element_blank())

# ggsave('06-cor_complex_stability.pdf', grid.arrange(pc, pm, nrow = 1), 
#        units = 'in', width = 5, height = 2.5)

##################################################################################
##################################################################################
##################################################################################
##################################################################################

fnscoh = list.files(pattern = '03-cohesion_')

values.cohesion = Group = Type = c()
for (fn in fnscoh) {
  gp = substr(fn, 13, 16)
  dfCoh = readRDS(fn)
  pos.cohes = dfCoh$`Positive Connectedness`
  neg.cohes = dfCoh$`Negative Connectedness`
  values.cohesion = append(values.cohesion, c(pos.cohes, neg.cohes))
  Type = append(Type, c(rep('Positive', length(pos.cohes)), 
                        rep('Negative', length(neg.cohes))))
  Group = append(Group, rep(gp, length(c(pos.cohes, neg.cohes))))
}

DFcohesion = data.frame(Group, Type, values.cohesion)
DFcohesion$GP = substr(DFcohesion$Group, 4, 4)
DFcohesion$Time = as.integer(substr(DFcohesion$Group, 2, 3))

DFcohesion.m = DFcohesion[DFcohesion$Group == 'D00C', ]
DFcohesion.m$Group = 'D00M'
DFcohesion.m$GP = 'M'
DFcohesion = rbind(DFcohesion, DFcohesion.m)

head(DFcohesion)

##########################################################

p = ggplot(DFcohesion, aes(x = Time, y = values.cohesion, group = GP, color = GP)) + 
  geom_point(shape = 16,  size = 1, alpha = 0.2) + 
  geom_smooth(method = "lm", se = F, size = 0.5, alpha = 0.75) + 
  scale_color_manual(values = c('#E69F00', '#0072B2')) + 
  facet_wrap(.~Type, ncol = 2, scales = 'free_y') + 
  labs(x = "", y = "") + theme_bw()+
  theme(panel.grid = element_blank(),
        # panel.grid = element_line(size = 0.1, color = 'grey90'), 
        strip.text = element_text(size = 10),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 8, hjust = 0.5, vjust = 0.5, angle = 90),
        axis.text.x = element_text(size = 10, color = 'black'),
        legend.position = "none")

# ggsave("06-plot_cohesion.pdf", p, units = "in", height = 2.5, width = 4)

# set the layout of result
pdf(width = 8, height = 2.5, file = '06-cohesion+cor.pdf' )
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 11)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}
print(p, vp = vplayout(1,1:5))
print(pc, vp = vplayout(1,6:8))
print(pm, vp = vplayout(1,9:11))
dev.off()


##########################################################

write(c(), '07-cohesion_lmSmooth.txt')
for (gp in unique(DFcohesion$GP)) {
  DFcohesion.f = DFcohesion[DFcohesion$GP == gp, ]
  for (fac in unique(DFcohesion.f$Type)) {
    DFcohesion.lm = DFcohesion.f[DFcohesion.f$Type == fac, ]
    res = summary(lm(DFcohesion.lm$values.cohesion ~ DFcohesion.lm$Time))
    write('----------------------', '07-cohesion_lmSmooth.txt', append = T)
    write(paste0('Group: ', gp), '07-cohesion_lmSmooth.txt', append = T)
    write(paste0('Factor: ', fac), '07-cohesion_lmSmooth.txt', append = T)
    write(paste0('slope: ', res$coefficients[,'Estimate'][1]), '07-cohesion_lmSmooth.txt', append = T)
    write(paste0('r.squared: ', res$r.squared), '07-cohesion_lmSmooth.txt', append = T)
    write(paste0('p_value: ', res$coefficients[,'Pr(>|t|)'][2]), '07-cohesion_lmSmooth.txt', append = T)
  }
}


