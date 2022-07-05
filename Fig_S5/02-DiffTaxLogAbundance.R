# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(RColorBrewer)
library(tidyverse)
library(reshape2)
library(microeco)
library(ape)
library(ggplot2)
library(grid)
library(doBy)

##############################################################

diy.aovSum = function(abundDF.all, Taxons, comp.gps) {
  p.adj = c()
  diff.g = c()
  G1mean = G2mean = c()
  for (taxon in Taxons) {
    abundDF = abundDF.all[abundDF.all$Taxonomy %in% taxon, ]
    abundDF = abundDF[abundDF$Group %in% comp.gps, ]
    abundDF$Group = factor(abundDF$Group, levels = comp.gps)
    result = t.test(Abundance ~ Group, data = abundDF)
    p.adj = append(p.adj, result$p.value)
    diff.g = append(diff.g, result$statistic)
    ag = abundDF[which(abundDF$Group == comp.gps[1]), 'Abundance' ]
    G1mean = append(G1mean, as.numeric(summary(ag)[4]))
    ag = abundDF[which(abundDF$Group == comp.gps[2]), 'Abundance' ]
    G2mean = append(G2mean, as.numeric(summary(ag)[4]))
  }
  newdf =data.frame(TaxLevel = rep(tax, length(Taxons)), Taxs = Taxons, 
                    P.adj = p.adj, Diff = diff.g, G1mean, G2mean )
  colnames(newdf)[5:6] = paste0(c('C', 'M'), 'mean')
  newdf$Groups = rep(paste(comp.gps, collapse = '_vs_'), nrow(newdf))
  return(newdf)
}

##############################################################

sample_info = read.table('01-samplesGroup.txt', header = T, 
                         sep = '\t', stringsAsFactors = F)
rownames(sample_info) = sample_info$SampleID
otu_table = read.table('01-OTUtable.txt', header = T, sep = '\t',
                       row.names = 1, stringsAsFactors = F)
taxonomy_table = read.table('01-comTax.txt',  header = T, sep = '\t', 
                            row.names = 1, stringsAsFactors = F)
phylo_tree = read.tree('01-FastTree.txt')
dataset = microtable$new(sample_table = sample_info, 
                         otu_table = otu_table, 
                         tax_table = taxonomy_table, 
                         phylo_tree = phylo_tree)
dataset$cal_abund()

##############################################################

fgs = rownames(dataset$taxa_abund$Genus)
Fname = Gname = c()
for (i in 1:length(fgs)) {
  fm = strsplit(fgs[i], split =  '|', fixed = T)[[1]][5]
  gm = strsplit(fgs[i], split =  '|', fixed = T)[[1]][6]
  Fname = append(Fname, fm)
  Gname = append(Gname, gm)
}

Gname = gsub('g_unclassified', 'g_u', Gname)
fgDF = data.frame(Fname, Gname)

##############################################################

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

##############################################################

tax = 'Genus'
cutoff.abund = 0.5

##############################################################

numRank = nrow(dataset$taxa_abund[[tax]])
t1 = trans_abund$new(dataset = dataset, taxrank = tax, ntaxa = numRank)
abundDF.all = t1$abund_data
Taxons = t1$use_taxanames

gps.list = list(T05 = c('T05C', 'T05M'), T12 = c('T12C', 'T12M'),
                T21 = c('T21C', 'T21M'), T30 = c('T30C', 'T30M'))
newdf = data.frame()
for (gps in gps.list) {
  newdf = rbind(newdf, diy.aovSum(abundDF.all = abundDF.all, 
                                  Taxons = Taxons, comp.gps = gps))
}

newdf$Log2Change = round(log2((newdf$Cmean+0.0001)/(newdf$Mmean+0.0001)), digits = 3)
newdf = newdf[!(is.infinite(newdf$Log2Change)), ]
newdf = newdf[!(newdf$Log2Change == 0), ]

write.table(newdf, '03-GenusDiffTaxsAmongGroups.txt', 
            append = F, quote = F, sep = '\t', 
            row.names = F, col.names = T)

##############################################################

DiffTaxsDF = read.table('03-GenusdiffTaxsAmongGroups.txt', header = T, 
                        sep = '\t', stringsAsFactors = F)

cgtm.df = DiffTaxsDF[DiffTaxsDF$Diff > 0 & DiffTaxsDF$P.adj <= 0.05, ]
cgtm.df = cgtm.df[which(cgtm.df$Cmean >= cutoff.abund), ]
cgtm.df = cgtm.df[!(str_detect(cgtm.df$Taxs, '_norank_') | 
                      str_detect(cgtm.df$Taxs, 'g_Clostridium_sensu') |
                      str_detect(cgtm.df$Taxs, 'g_livecontrol') |
                      str_detect(cgtm.df$Taxs, 'g_RBG') | 
                      str_detect(cgtm.df$Taxs, 'g_EBM') | 
                      str_detect(cgtm.df$Taxs, 'g_unclassified_c_Clostridia') | 
                      str_detect(cgtm.df$Taxs, 'g_unclassified_c_Bacteroidia') | 
                      str_detect(cgtm.df$Taxs, 'g_unclassified_c_Alphaproteobacteria') | 
                      str_detect(cgtm.df$Taxs, 'g_unclassified_o_Peptostreptococcales') | 
                      str_detect(cgtm.df$Taxs, 'g_unclassified_o_Sphingobacteriales')), ]

clem.df = DiffTaxsDF[DiffTaxsDF$Diff < 0 & DiffTaxsDF$P.adj <= 0.05, ]
clem.df = clem.df[which(clem.df$Mmean >= cutoff.abund), ]
clem.df = clem.df[!(str_detect(clem.df$Taxs, '_norank_') | 
                      str_detect(clem.df$Taxs, 'g_Clostridium_sensu') |
                      str_detect(clem.df$Taxs, 'g_livecontrol') |
                      str_detect(clem.df$Taxs, 'g_RBG') | 
                      str_detect(clem.df$Taxs, 'g_EBM') | 
                      str_detect(clem.df$Taxs, 'g_unclassified_c_Clostridia') | 
                      str_detect(clem.df$Taxs, 'g_unclassified_c_Bacteroidia') | 
                      str_detect(clem.df$Taxs, 'g_unclassified_c_Alphaproteobacteria') | 
                      str_detect(clem.df$Taxs, 'g_unclassified_o_Peptostreptococcales') | 
                      str_detect(clem.df$Taxs, 'g_unclassified_o_Sphingobacteriales')), ]


###############################################################

newdf = rbind(cgtm.df, clem.df)
newdf = rbind(cgtm.df, clem.df)
newdf$Taxs = gsub('g_unclassified', 'g_u', newdf$Taxs)

fgDF.f = fgDF[fgDF$Gname %in% newdf$Taxs, ]
taxlev = fgDF.f[order(fgDF.f$Fname, decreasing = T), 'Gname']
newdf$Taxs = factor(newdf$Taxs, levels = taxlev)

filter.taxs =  sort(unique(newdf$Taxs))

newdf$Groups = factor(newdf$Groups, levels = sort(unique(newdf$Groups)))
newdf$signif = rep(' ', nrow(newdf))
newdf$signif[which(newdf$P.adj <= 0.001)] = '***'
newdf$signif[which(newdf$P.adj > 0.001 & newdf$P.adj <= 0.01)] = '**'
newdf$signif[which(newdf$P.adj > 0.01 & newdf$P.adj <= 0.05)] = '*'


p1 = ggplot(newdf, aes(x = Groups, y = Taxs, fill = Log2Change)) + 
  geom_raster() + geom_text(label = newdf$signif, size = 3, nudge_y = -0.25) +
  scale_fill_gradient2(low = "#D55E00", mid = "grey80", 
                       high = "#0072B2", midpoint = 0) +
  theme_bw()+
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_text(size = 8, color = 'black'), 
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key =  element_blank(),
        legend.key.height = unit(1, 'mm'),
        legend.text = element_text(size = 6))

###############################################################

DiffTaxsDF = read.table('03-GenusdiffTaxsAmongGroups.txt', header = T,
                        sep = '\t', stringsAsFactors = F)

cgtm.df = DiffTaxsDF[DiffTaxsDF$Diff > 0 & DiffTaxsDF$P.adj <= 0.05, ]
cgtm.df = cgtm.df[which(cgtm.df$Cmean >= cutoff.abund), ]
cgtm.df = cgtm.df[!(str_detect(cgtm.df$Taxs, '_norank_') | 
                      str_detect(cgtm.df$Taxs, 'g_Clostridium_sensu') |
                      str_detect(cgtm.df$Taxs, 'g_livecontrol') |
                      str_detect(cgtm.df$Taxs, 'g_RBG') | 
                      str_detect(cgtm.df$Taxs, 'g_EBM')  | 
                      str_detect(cgtm.df$Taxs, 'g_unclassified_c_Clostridia') | 
                      str_detect(cgtm.df$Taxs, 'g_unclassified_c_Bacteroidia') | 
                      str_detect(cgtm.df$Taxs, 'g_unclassified_c_Alphaproteobacteria') | 
                      str_detect(cgtm.df$Taxs, 'g_unclassified_o_Peptostreptococcales') | 
                      str_detect(cgtm.df$Taxs, 'g_unclassified_o_Sphingobacteriales')), ]

clem.df = DiffTaxsDF[DiffTaxsDF$Diff < 0 & DiffTaxsDF$P.adj <= 0.05, ]
clem.df = clem.df[which(clem.df$Mmean >= cutoff.abund), ]
clem.df = clem.df[!(str_detect(clem.df$Taxs, '_norank_') | 
                      str_detect(clem.df$Taxs, 'g_Clostridium_sensu') |
                      str_detect(clem.df$Taxs, 'g_livecontrol') |
                      str_detect(clem.df$Taxs, 'g_RBG') | 
                      str_detect(clem.df$Taxs, 'g_EBM') | 
                      str_detect(clem.df$Taxs, 'g_unclassified_c_Clostridia') | 
                      str_detect(clem.df$Taxs, 'g_unclassified_c_Bacteroidia') | 
                      str_detect(clem.df$Taxs, 'g_unclassified_c_Alphaproteobacteria') | 
                      str_detect(clem.df$Taxs, 'g_unclassified_o_Peptostreptococcales') | 
                      str_detect(clem.df$Taxs, 'g_unclassified_o_Sphingobacteriales')), ]

newdf = rbind(cgtm.df, clem.df)

DFmelt = c()
for (gp in unique(newdf$Groups)) {
  gpc = unlist(strsplit(gp, split = '_vs_'))[1]
  gpm = unlist(strsplit(gp, split = '_vs_'))[2]
  gDF = newdf[which(newdf$Groups == gp), ]
  gDF = melt(gDF, measure.vars = c('Cmean', 'Mmean'),
             variable.name = 'Grp', value.name = 'Abundance')
  gDF$Grp = as.character(gDF$Grp)
  gDF$Grp[which(gDF$Grp == 'Cmean')] = gpc
  gDF$Grp[which(gDF$Grp == 'Mmean')] = gpm
  DFmelt = rbind(DFmelt, gDF)
}

DFmelt$Groups = factor(DFmelt$Groups, levels = sort(unique(DFmelt$Groups)))
DFmelt$Grp = factor(DFmelt$Grp, levels = sort(unique(DFmelt$Grp)))

DFmelt = DFmelt[, c('Grp', 'Taxs', 'Abundance')]

##############################

meanDF = summaryBy(Abundance~Taxonomy+Group, data = abundDF.all, FUN = mean)
D00.df = meanDF[meanDF$Taxonomy %in% DFmelt$Taxs, c('Group', 'Taxonomy', 'Abundance.mean')]
D00.df = D00.df[D00.df$Group == 'T00C', ]
colnames(D00.df) = colnames(DFmelt)
DFmelt = rbind(D00.df, DFmelt)
DFmelt$Abundance[DFmelt$Abundance < cutoff.abund] = NA

DFmelt$Taxs = gsub('g_unclassified', 'g_u', DFmelt$Taxs)
DFmelt = DFmelt[DFmelt$Taxs %in% filter.taxs, ]

fgDF.f = fgDF[fgDF$Gname %in% DFmelt$Taxs, ]
taxlev = fgDF.f[order(fgDF.f$Fname, decreasing = T), 'Gname']
DFmelt$Taxs = factor(DFmelt$Taxs, levels = taxlev)

##############################

DFmelt.abund = DFmelt
DFmelt.abund$Grp = gsub('T', 'D', DFmelt.abund$Grp)
abundaceDF = dcast(DFmelt.abund, Taxs~ Grp, value.var = 'Abundance')
rownames(fgDF.f) = as.character(fgDF.f$Gname)
write.table(data.frame(Family = fgDF.f[as.character(abundaceDF$Taxs), 'Fname'], abundaceDF),
            na = '<0.5', '04-GenusAbundaceDF.txt', append = F, quote = F,
            sep = '\t', row.names = F, col.names = T)

##############################

outDF = fgDF.f[order(fgDF.f$Fname, decreasing = F), ]
write(paste0(outDF$Gname, '_', outDF$Fname), '05-Genus_familyGenus.txt',
      append = F, sep = '\n')

##############################

p2 = ggplot(DFmelt, aes(x = Grp, y = Taxs, fill = Abundance)) + 
  geom_raster() + 
  scale_fill_gradient2(low = brewer.pal(9,'Spectral')[9], 
                       mid = brewer.pal(9,'Spectral')[5], 
                       high = brewer.pal(9,'Spectral')[1], 
                       midpoint = 12, na.value = NA) +
  scale_y_discrete(position = "right") +
  theme_bw()+
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.text.y.right = element_text(size = 8, color = 'black'), 
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key =  element_blank(),
        legend.key.height = unit(1, 'mm'),
        legend.text = element_text(size = 6))


# set the layout of result
pdf(width = 8, height = 5.5, file = '06-GenusDiffTaxsLog2Abundance.pdf' )
grid.newpage()  
pushViewport(viewport(layout = grid.layout(1, 10))) 
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p1, vp = vplayout(1,(1:4)+1))
print(p2, vp = vplayout(1,(5:8)+1))
dev.off()

