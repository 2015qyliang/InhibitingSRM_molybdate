# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

##################################################################

library(tidyverse)
library(doBy)
library(reshape2)


##################################################################
# filter Bins by diff-taxs
binTaxDF = read.table('PD.Bin_TopTaxon.csv', 
                      header = T, sep = ',', stringsAsFactors = F)
binsDF.f = binTaxDF[, c("Bin", "TopTaxonID", "TopTaxon.Phylum", 'BinRA')]
colnames(binsDF.f) = c("Bin", "OTUID", 'Phylum', 'BinRA')
rownames(binsDF.f) = binsDF.f$OTUID

confirm.phylum = c('p_Myxococcota','p_Desulfobacterota','p_Bdellovibrionota',
                   'p_Bacteroidota','p_Fusobacteriota','p_Proteobacteria',
                   'p_Acidobacteriota','p_Firmicutes','p_Actinobacteriota',
                   'p_Spirochaetota','p_Chloroflexi','p_Synergistota',
                   'p_Planctomycetota','p_Verrucomicrobiota','p_Cyanobacteria')

binsDF.f$Phylum[!(binsDF.f$Phylum %in% confirm.phylum)] = 'Other'
binsDF.f$Phylum = gsub(binsDF.f$Phylum, pattern = 'p_', replacement = '')

##################################################################
# filter contribution of Bins -- HoS & DL
binContriDF = read.table('PD.BinContributeToProcess_EachGroup.csv', 
                         header = T, sep = ',', stringsAsFactors = F)
gps = unique(binContriDF$Group[1:45])
binContriDF.f = binContriDF[binContriDF$Group %in% gps, 3:ncol(binContriDF)]
colnames(binContriDF.f)[3:ncol(binContriDF.f)] = paste0('Bin', 1:(ncol(binContriDF.f)-2))
binContriDF.f = binContriDF.f[, c(1:2, which(colnames(binContriDF.f) %in% binsDF.f$Bin))]

binContriDF.f$GP = substr(binContriDF.f$Group, 4, 4)
binContriDF.f = binContriDF.f[, 2:ncol(binContriDF.f)]

binContriDF.f.melt = melt(binContriDF.f, variable.name = 'Bins',
                          measure.vars = colnames(binContriDF.f)[2:(ncol(binContriDF.f)-1)])
binContriDF.long = summaryBy(.~GP+Bins+Process, data = binContriDF.f.melt, FUN = mean)

bindf = binContriDF.long[binContriDF.long$Process %in% c('HoS', 'DL'), ]
bindf$Process = paste0(bindf$Process, '.', bindf$GP)
bindf = bindf[, 2:4]
bindf = dcast(bindf, Bins~Process, value.var = 'value.mean')
bindf$HoS = 100*(bindf$HoS.M - bindf$HoS.C)
bindf$DL = 100*(bindf$DL.M - bindf$DL.C)

write.table(data.frame(Bins = binsDF.f$Bin, 
                       Phylum = binsDF.f$Phylum, 
                       abundace = 100*binsDF.f$BinRA, 
                       Hos = bindf$HoS, DL = bindf$DL), 
            '04-M-simpleBar.txt', append = F, 
            row.names = F, col.names = T, sep = '\t', quote = F)

##################################################################

