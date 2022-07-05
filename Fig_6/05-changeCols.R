
library(tidyverse)

for (gp in c("M")) {
  dirgp = paste0(gp, '_simpleBar')
  abfn = readLines(file.path(dirgp, 'iTOL_simplebar-abundace.txt'))
  hosfn = readLines(file.path(dirgp, 'iTOL_simplebar-HoS.txt'))
  dlfn = readLines(file.path(dirgp, 'iTOL_simplebar-DL.txt'))
  
  ln = str_detect(abfn, '#000000')
  abfn[ln] = gsub(abfn[ln], pattern = '#000000', replacement = '#ADACA5')
  write(abfn, file.path(dirgp, 'iTOL_simplebar-abundace.txt'))
  
  ln = str_detect(hosfn, '#000000')
  hosfn[ln] = gsub(hosfn[ln], pattern = '#000000', replacement = '#4DBBD5BF')
  write(hosfn, file.path(dirgp, 'iTOL_simplebar-HoS.txt'))
  
  ln = str_detect(dlfn, '#000000')
  dlfn[ln] = gsub(dlfn[ln], pattern = '#000000', replacement = '#3C5488BF')
  write(dlfn, file.path(dirgp, 'iTOL_simplebar-DL.txt'))
}


# abundance - #000000 #ADACA5
# HoS - #000000 #4DBBD5BF
# DL - #000000 #3C5488BF

