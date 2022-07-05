# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(vegan)
library(tidyverse)
library(gridExtra)
library(microeco)
library(ape)
library(ggplot2)
library(ggsci)

# #############################################################################

sample_info = read.table('02-samplesGroup.txt', header = T, 
                         sep = '\t', stringsAsFactors = F)
rownames(sample_info) = sample_info$SampleID
otu_table = read.table('03-OTUtable.txt', header = T, sep = '\t',
                       row.names = 1, stringsAsFactors = F)
taxonomy_table = read.table('03-comTax.txt',  header = T, sep = '\t', 
                            row.names = 1, stringsAsFactors = F)
phylo_tree = read.tree('02-FastTree.txt')
dataset = microtable$new(sample_table = sample_info, 
                         otu_table = otu_table, 
                         tax_table = taxonomy_table, 
                         phylo_tree = phylo_tree)

dataset$cal_abund()

# #############################################################################
diyplotBar = function(dataset, tax, numRank) {
  t1 = trans_abund$new(dataset = dataset, taxrank = tax,
                       ntaxa = numRank, groupmean = "Group")
  p = t1$plot_bar(use_colors = colorRampPalette(pal_npg()(10))(numRank),
                  others_color = "grey70",
                  xtext_keep = FALSE, legend_text_italic = T) + 
    theme(panel.background = element_rect(fill="white",color="black"),
          panel.grid = element_blank(), 
          axis.title.y.left =  element_text(size = 10),
          axis.text.x = element_text(size = 8, color = 'black', angle = 45, hjust = 1, vjust = 1), 
          axis.text.y = element_text(size = 8, color = 'black', angle = 90), 
          legend.key = element_blank(), 
          legend.position = 'bottom',
          legend.title = element_blank(),
          legend.text = element_text(size = 4)) + 
    guides(fill = guide_legend( ncol = 3,
                                keywidth = unit(1, 'mm'),
                                keyheight = unit(1, 'mm'),
                                label.hjust = 0,
                                label.vjust = 0.5,
                                label.position = 'right',
                                title.position = 'top', 
                                direction = 'horizonal'))
  return(p)
}

# #############################################################################

newtax = dataset$taxa_abund
for (i in 2:length(newtax)) {
  df1 = newtax[[i]]
  newRN = c()
  for (rn in rownames(df1)) {
    rn = strsplit(rn, split =  '|', fixed = T)[[1]][i]
    newRN = append(newRN, rn)
  }
  rownames(newtax[[i]]) = newRN
}
dataset$taxa_abund = newtax

tax = 'Family'
dpb = diyplotBar(dataset, tax, 20)

# #############################################################################
# #############################################################################
# #############################################################################
# #############################################################################
# #############################################################################

dclassDF = read.table('02-samplesGroup.txt', header = T, 
                      sep = '\t', stringsAsFactors = F)

otuSampleDF = read.table('01-Res_aFRI.txt', 
                         header = T,sep = '\t',row.names = 1)
otuSampleDF = t(otuSampleDF[, 1:180])

nmds = metaMDS(vegdist(otuSampleDF, method = 'jaccard', binary = F), 
               k = 2, trymax = 100, wascores = TRUE)
info.stress = round(nmds$stress, digits = 2)

data.scores = as.data.frame(scores(nmds))
dt.scores = data.frame(data.scores, dclass = dclassDF$Group)

newdf = data.frame() 
newdf = rbind(newdf, dt.scores)

# write.table(newdf, '02-NMDS.txt', append = F, quote = F,
#             sep = '\t', row.names = F, col.names = T)
# newdf = read.table('02-NMDS.txt', header = T, sep = '\t')

c.newdf = newdf[which(newdf$dclass %in% c('D00C', 'D05C', 'D12C', 'D21C', 'D30C')), ]
c.newdf$dclass = factor(c.newdf$dclass,
                        labels = c('D00C', 'D05C', 'D12C', 'D21C', 'D30C'),
                        levels = c('D00C', 'D05C', 'D12C', 'D21C', 'D30C'))
m.newdf = newdf[which(newdf$dclass %in% c('D05M', 'D12M', 'D21M', 'D30M')), ]
m.newdf$dclass = factor(m.newdf$dclass,
                        labels = c('D05M', 'D12M', 'D21M', 'D30M'),
                        levels =  c('D05M', 'D12M', 'D21M', 'D30M'))

nmds.p = ggplot(newdf, aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = c.newdf, aes(x = NMDS1, y = NMDS2, color = dclass),
             size = 1.5, alpha = 0.7, shape = 16) + # point
  geom_point(data = m.newdf, aes(x = NMDS1, y = NMDS2, color = dclass),
             size = 1.5, alpha = 0.7, shape = 17) + # triangle
  labs(caption = paste0('Stress: ', paste(info.stress, collapse = ', '))) +
  scale_colour_manual(values = c(pal_npg()(5), pal_npg()(5)[2:5]),
                      breaks = c('D00C', 'D05C', 'D12C', 'D21C', 'D30C',
                                 'D05M', 'D12M', 'D21M', 'D30M')) +
  geom_hline(yintercept=0, colour="grey60", linetype="dashed")+
  geom_vline(xintercept=0, colour="grey60", linetype="dashed")+
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(),
        plot.caption =  element_text(size = 8), 
        axis.title = element_text(size = 10),
        axis.text.y = element_text(size = 8, hjust = 0.5, vjust = 0.5, angle = 90, color = 'black'),
        axis.text.x = element_text(size = 8, color = 'black'),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 4),
        legend.key = element_blank())  +
  guides(color = guide_legend(nrow = 2,
                              keywidth = unit(1, 'mm'),
                              keyheight = unit(1, 'mm'),
                              label.position = 'bottom',
                              title.position = 'top',
                              direction = 'horizonal'))

# #############################################################################

ggsave('05-NetCommunity_CompsFRI.pdf', grid.arrange(dpb, nmds.p, nrow = 1), 
       units = 'in', width = 7, height = 4)

