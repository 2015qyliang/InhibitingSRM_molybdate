# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

##################################################################

library(ggsci)
library(tidyverse)
library(doBy)
library(reshape2)

##################################################################
# filter Bins by diff-taxs
binTaxDF = read.table('PD.Bin_TopTaxon.csv', 
                      header = T, sep = ',', stringsAsFactors = F)
binsDF.f = binTaxDF[, c("Bin", "TopTaxonID", "TopTaxon.Family", "TopTaxon.Phylum")]
colnames(binsDF.f) = c("Bin", "OTUID", "Family", "Phylum")

#######################################

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

#######################################

as.pre = c("HeS", "HoS", "HD", "DL", "DR")
as.pre.cols = pal_npg(alpha = 0.75)(5)

bindf.m = binContriDF.long[binContriDF.long$GP == 'M', ]
bindf.m = dcast(bindf.m[, 2:4], Bins~Process, value.var = 'value.mean' )

#######################################

bindf.r = bindf.m[, 2:6]
rownames(bindf.r) = bindf.m$Bins
rsums = rowSums(bindf.r)
bindf.r$HeS = bindf.r$HeS/rsums
bindf.r$HoS = bindf.r$HoS/rsums
bindf.r$DL = bindf.r$DL/rsums
bindf.r$HD = bindf.r$HD/rsums
bindf.r$DR = bindf.r$DR/rsums
bindf.r = bindf.r[, c("HeS", "HoS", "HD", "DL", "DR")]

write('DATASET_MULTIBAR', '02-M-multipleBar.txt', append = F)
write('SEPARATOR	TAB', '02-M-multipleBar.txt', append = T)
write('DATASET_LABEL	Assembly', '02-M-multipleBar.txt', append = T)
write('COLOR	#ff0000', '02-M-multipleBar.txt', append = T)
write('LEGEND_TITLE	legend', '02-M-multipleBar.txt', append = T)
write('LEGEND_SHAPES	1	1	1	1	1', '02-M-multipleBar.txt', append = T)
write(paste0('FIELD_COLORS	', paste(as.pre.cols, collapse = '\t')), 
      '02-M-multipleBar.txt', append = T)
write(paste0('FIELD_LABELS	', paste(as.pre, collapse = '\t')), 
      '02-M-multipleBar.txt', append = T)
write(paste0('LEGEND_COLORS	', paste(as.pre.cols, collapse = '\t')), 
      '02-M-multipleBar.txt', append = T)
write(paste0('LEGEND_LABELS	', paste(as.pre, collapse = '\t')), 
      '02-M-multipleBar.txt', append = T)
write('DATA', '02-M-multipleBar.txt', append = T)
write.table(bindf.r, '02-M-multipleBar.txt', append = T, 
            row.names = T, col.names = F, sep = '\t', quote = F)

##################################################################

