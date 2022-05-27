
library(INLA)
library(rgdal)

if(!file.exists('Ohio_data/Data.txt')){
    if(!file.exists('Ohio_data.zip'))
        download.file(
            'https://sites.google.com/a/r-inla.org/stbook/Ohio_data.zip?attredirects=0&d=1',
            'Ohio_data.zip')
    unzip('Ohio_data.zip')
}

ohio <- readOGR('Ohio_data/', 'tl_2010_39_county00')
head(ohio@data)

odata <- read.delim('Ohio_data/Data.txt', skip=16)

names(odata)[which(names(odata)=='n')] <- 'N'

head(odata)

ograph <- inla.read.graph('Ohio_data/Ohio.graph')
apropos('graph')

save(list=c('ograph', 'odata', 'ohio'),
     file='ohio.RData', compress='xz')
