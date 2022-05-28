library(INLA)

inla.setOption(
    inla.mode='experimental',
    smtp='pardiso',
   pardiso.license='~/.pardiso.lic')

## inla.pardiso()

ctri <- list(
    strategy='gaussian', ### run faster for non-gaussian likelihood
    int.strategy='ccd',
    control.vb = list(
        enable = TRUE)
)

ctrc <- list(
    dic = TRUE, waic = TRUE, cpo = TRUE)

source('get_ohio_data.R')
ls()
head(odata,2)

summary(odata)

odata$gender <- odata$gender - 1
odata$race <- odata$race - 1

formula0 <- y ~ gender + race +
    f(year, model='iid') +
    f(county, model='besag', graph=ograph)

result0 <- inla(
    formula = formula0,
    family = 'poisson',
    data = odata,
    control.inla = ctri,
    control.compute = ctrc)

result0$cpu

formula0b <- y ~ gender + race +
    f(year, model='iid') +
    f(year2, model='rw2') + 
    f(county, model='bym2', graph=ograph)

odata$year2 <- odata$year

result0b <- inla(
    formula = formula0b,
    family = 'poisson',
    data = odata,
    control.inla = ctri,
    control.compute = ctrc)

result0b$cpu

c(result0$dic$dic,
  result0b$dic$dic)

### set up an 'iid' spacetime term
### interaction type I (Knorr-Held, 2000)
formula1 <- update(
    formula0b,
    .~.+f(st, model='iid'))

nt <- length(unique(odata$year))
nt

nc <- length(unique(odata$county))
nc

nt * nc
dim(odata)

odata$st <- odata$county +
    (odata$year-1)*nc
length(unique(odata$st))

result1 <- inla(
    formula = formula1,
    family = 'poisson',
    data = odata,
    control.inla = ctri,
    control.compute = ctrc)

c(-sum(log(result0b$cpo$cpo)),
  -sum(log(result1$cpo$cpo)))

### Martinez-Beneito (2008) model
formula1b <- update(
    formula0b,
    .~.+f(county2, model='besag',
          graph=ograph, 
          group=year3,
          control.group=list(model='ar1')))

odata$county2 <- odata$county
odata$year3 <- odata$year

result1b <- inla(
    formula = formula1b,
    family = 'poisson',
    data = odata,
    verbose=TRUE,
    control.inla = ctri,
    control.compute = ctrc)

c(result0$dic$dic,
  result0b$dic$dic,
  result1$dic$dic,
  result1b$dic$dic)

c(-sum(log(result0b$cpo$cpo)),
  -sum(log(result1$cpo$cpo)),
  -sum(log(result1b$cpo$cpo)))


### interaction type IV (Knorr-Held, 2000)
## help(inla.knmodels)
tgraph <- sparseMatrix(
    i=c(2:nt, 1:(nt-1)),
    j=c(1:(nt-1), 2:nt), x=1)
tgraph

formula0c <- y ~ gender + race +
    f(year, model='bym2', graph=tgraph) +
    f(county, model='bym2', graph=ograph)

result4 <- inla.knmodels(
    formula0c,
    control.st=list(
        time=year, space=county,
        spacetime=st,
        graph=ograph, type='4'),
    family = 'poisson',
    data = odata,
    verbose=TRUE,
    control.inla = ctri,
    control.compute = ctrc)
    

c(-sum(log(result0b$cpo$cpo)),
  -sum(log(result1$cpo$cpo)),
  -sum(log(result1b$cpo$cpo)),
  -sum(log(result4$cpo$cpo)))

plot(ohio)

names(result4$summary.random)
st4 <- result4$summary.random$st
dim(st4)
head(st4)

ohio@data <-
    data.frame(
        ohio@data,
        st=matrix(st4$mean, nrow=nc))

names(ohio)            

spplot(ohio, 'st.1')

spplot(ohio,
       c('st.1', 'st.7', 'st.13',
         'st.21'))
