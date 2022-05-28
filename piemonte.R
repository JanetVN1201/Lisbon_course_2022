### fit the model in
## https://link.springer.com/article/10.1007/s10182-012-0196-3

library(INLA)

inla.setOption(
    inla.mode='experimental',
    num.threads='8:-1',
    smtp='pardiso', 
    pardiso.license='~/.pardiso.lic')

ctri <- list(
    int.strategy='ccd',
    parallel.linesearch=TRUE)

### get the dataset
u0 <- paste0(
    'http://inla.r-inla-download.org/',
    'r-inla.org/case-studies/Cameletti2012/')
coofl <- 'coordinates.csv'
datafl <- 'Piemonte_data_byday.csv'
bordersfl <- 'Piemonte_borders.csv'

if(!file.exists(bordersfl)) 
    download.file(paste0(u0, bordersfl), bordersfl)

dim(pborders <- read.csv(bordersfl))

if(!file.exists(coofl)) 
    download.file(paste0(u0, coofl), coofl)
dim(locs <- read.csv(coofl))

if(!file.exists(datafl)) 
    download.file(paste0(u0, datafl), datafl)
dim(pdata <- read.csv(datafl))

head(pdata)

### prepare and select time 
range(pdata$Date <- as.Date(pdata$Date, '%d/%m/%y'))
pdata$time <- as.integer(difftime(
    pdata$Date, min(pdata$Date), units='days'))+1

table(pdata$time)
if(!any(ls()=='nt'))
    nt <- 182 ### max number of time points to be used
pdata <- pdata[pdata$time<=nt,]

sd(pdata$PM10, na.rm=TRUE)

### as in Cameleti et. al. 2012
smesh <- inla.mesh.2d(
        cbind(locs[,2], locs[,3]),
    loc.domain=pborders, 
    max.edge=c(50, 300), 
    offset=c(10, 140), 
    cutoff=5, 
    min.angle=c(26, 21))

if(FALSE) {## alternative mesh

    bnd <- inla.nonconvex.hull(
        rbind(cbind(locs[,2], locs[,3]),
              cbind(pborders[,1], pborders[,2])), 
        convex=10, concave=50, resolution=100)
    
    smesh <- inla.mesh.2d(
        boundary=bnd, 
        max.edge=c(20, 100), 
        offset=c(50, 150), 
        cutoff=10, 
        min.angle=c(26, 21))

}

smesh$n

par(mar=c(3,3,1,0))
plot(smesh, asp=1)
points(locs[,2:3], pch=8, col='red')
lines(pborders, lwd=2, col='blue')

### define the stationary spatial model
psigma <- c(20, 0.1) ## P(sigma > 20) = 0.1
prs <- c(100, 0.5) ## P(range < 100) = 0.5
spde <- inla.spde2.pcmatern(
    mesh=smesh, alpha=2,
    prior.sigma=psigma, 
    prior.range=prs)
smesh$n

### define the projector matrix
A <- inla.spde.make.A(
    smesh, as.matrix(pdata[c('UTMX', 'UTMY')]),
    group=pdata$time, n.group=nt)
dim(A)
stopifnot(sum(A)==nrow(pdata))

### prepare the data including the covariates
xnames <- c('A', 'WS', 'TEMP', 'HMIX', 'PREC', 'EMI')
xmean <- colMeans(pdata[, xnames])
xsd <- sapply(pdata[xnames], sd)
xx.scaled <- scale(pdata[xnames], xmean, xsd)


### spacetime index set
idx.st <- inla.spde.make.index(
    name='s', n.spde=spde$n.spde, n.group=nt)

sapply(idx.st, range)

dsstack <- inla.stack(
    tag='e',
    data=list(y=log(pdata$PM10)), ### not agree
    effects=list(
        data.frame(b0=1, xx.scaled),
        idx.st), 
    A=list(1, A))

### define the separable vector AR1 model
prho <- list(
    theta=list(
        prior='pc.cor1',
        param=c(0,0.9)) ## P(cor>0)=0.9
)

### the linear predictor formula
formulae <- update(
    y~0+b0+f(s, model=spde, group=s.group,
             control.group=list(model='ar1', hyper=prho)),
    paste('.~.+', paste(xnames, collapse='+')))

### likelihood precision prior
lkprec <- list(
    prec=list(initial=-4, fixed=FALSE,
              prior='pcprec',
              param=c(1, 0.01)) ## P(sigma>1)=0.01
)

### fit the separable model
result <- inla(formulae,
           data=inla.stack.data(dsstack),
           control.predictor=list(
               A=inla.stack.A(dsstack)),
           control.fixed=list(prec=c(b0=1)),
           control.inla=ctri, 
           control.family=list(hyper=lkprec))

result$cpu

result$summary.fix

##             mean          sd  0.025quant    0.5quant   0.975quant mode kld
## b0    3.08777713 0.400555350  2.30270307  3.08777713  3.872851194   NA   0
## A    -0.16770320 0.046827214 -0.25948285 -0.16770320 -0.075923544   NA   0
## WS   -0.06000185 0.008389926 -0.07644580 -0.06000185 -0.043557898   NA   0
## TEMP -0.11966910 0.034778672 -0.18783404 -0.11966910 -0.051504156   NA   0
## HMIX -0.02405629 0.013023105 -0.04958111 -0.02405629  0.001468527   NA   0
## PREC -0.05353958 0.008549584 -0.07029646 -0.05353958 -0.036782704   NA   0
## EMI   0.03785135 0.014613458  0.00920950  0.03785135  0.066493204   NA   0

result$summary.hy[,1:2]

##                                                mean           sd
## Precision for the Gaussian observations  30.7281502  1.313238533
## Range for s                             278.6025572 17.974170176
## Stdev for s                               1.1745614  0.110575884
## GroupRho for s                            0.9624957  0.007290129

bb <- t(apply(pborders, 2, range))
bb

dxy <- apply(bb, 1, diff)
dxy

gridp <- inla.mesh.projector(
    smesh,
    xlim=bb[1,], ylim=bb[2,], 
    dims=round(2*dxy))

### make maps from here...
st.m <- matrix(result$summary.random$s$mean, smesh$n)
dim(st.m)

### identify pixels outside the borders map
spborders <- SpatialPolygons(
    list(Polygons(list(Polygon(pborders)), '0')))
id.out <- which(is.na(over(
    SpatialPoints(gridp$lattice$loc),
    spborders)))

### project it
z <- inla.mesh.project(gridp, st.m[,1])

### set pixels outside as NA
z[id.out] <- NA

### visualize the first time on the map
library(fields)
plot(pborders, asp=1)
image.plot(gridp$x, gridp$y,
           z, add=TRUE)

### make a unique breaks (overall) for colors
bk <- seq(min(st.m)-1e-5, max(st.m)+1e-5, length=51)
colorsk <- tim.colors(length(bk)-1)

### visualize some days on the map
par(mfrow=c(5,8), mar=c(0,0,0,0))
for(k in round(seq(1, nt, length=40))) {
    z <- inla.mesh.project(gridp, st.m[,k])
    z[id.out] <- NA
    plot(spborders)
    image.plot(gridp$x, gridp$y, z, add=TRUE,
               breaks=bk, col=colorsk)
    legend('topleft', '', bty='n',
           title=paste(sort(unique(pdata$Date))[k]))
}

