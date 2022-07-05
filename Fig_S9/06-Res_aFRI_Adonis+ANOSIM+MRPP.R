# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(vegan)
library(tidyverse)

#############################################################

compares = list(c('D00C', 'D05C'), c('D00C', 'D12C'), c('D00C', 'D21C'), 
                c('D00C', 'D30C'), c('D00C', 'D05M'), c('D00C', 'D12M'),
                c('D00C', 'D21M'), c('D00C', 'D30M'), c('D05C', 'D05M'),
                c('D12C', 'D12M'), c('D21C', 'D21M'), c('D30C', 'D30M'))

###############################
otuSampleDF = read.table('01-Res_aFRI.txt', 
                         header = T,sep = '\t',row.names = 1)
otuSampleDF = t(otuSampleDF[, 1:180])

otuSampleDF = as.matrix(otuSampleDF)
# otuSampleDF[otuSampleDF != 0] = 1

adonis.F = adonis.p = anosim.r = anosim.p = mrpp.delta = mrpp.p = gpPairs = c()
for (i in 1:length(compares)) {
  choseSamples = str_detect(rownames(otuSampleDF), compares[[i]][1]) | 
    str_detect(rownames(otuSampleDF), compares[[i]][2])
  choseDF = otuSampleDF[choseSamples, ]
  choseDF = choseDF[, colSums(choseDF) != 0]
  res.adonis = adonis(choseDF ~ substr(rownames(choseDF), 1, 4),
                      permutations = 999, method = 'jaccard')
  res.anosim = anosim(choseDF, substr(rownames(choseDF), 1, 4), 
                      permutations = 999, distance = 'jaccard')
  res.mrpp = mrpp(choseDF, substr(rownames(choseDF), 1, 4), 
                      permutations = 999, distance = 'jaccard')
  adonis.F = append(adonis.F, res.adonis$aov.tab$F.Model[1])
  adonis.p = append(adonis.p, res.adonis$aov.tab$`Pr(>F)`[1])
  anosim.r = append(anosim.r, res.anosim$statistic)
  anosim.p = append(anosim.p, res.anosim$signif)
  mrpp.delta = append(mrpp.delta, res.mrpp$delta)
  mrpp.p = append(mrpp.p, res.mrpp$Pvalue)
  gpPairs = append(gpPairs, paste(compares[[i]], collapse = '_vs_'))
}

write.table(data.frame(gpPairs, adonis.F, adonis.p, anosim.r, anosim.p, mrpp.delta, mrpp.p), 
            '07-CompositionTest_aFRI_pair.txt', append = F, quote = F, sep = '\t', 
            row.names = F, col.names = T)

###############################

