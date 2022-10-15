################################################################################
source("Dataset_Preparation_New.R")
################################################################################
### Smooth the curves SET BY SET (Preferable for testing the "PRACTICE EFFECT")
### Consider the object "subject.dfs.imp". For each subject create a nested list
### whose elements are the dataframes for each set

rangeval <- c(1,set.len)
basis.set <- create.bspline.basis(rangeval = rangeval, norder = 4,
                                  nbasis = 8)

# Define a Penalty object (penalization of the curvature, i.e., 2nd derivative)
pen.obj <- fdPar(fdobj = basis.set, Lfdobj = 2, lambda = 10^(1))
# pen.obj <- fdPar(fdobj = basis.set)

# Same argval for each subject
abscissa <- 1:set.len

# Smooth the curves: group together the curves for all the subjects in a set
variable.dfs.FD.byset <- vector("list", n.vars)
variable.dfs.FD.byset <- lapply(variable.dfs.byset, function(var.df)
{
  lapply(var.df, function(set.df)
  {
    smooth.basis(argvals = abscissa, 
                 y = as.matrix(set.df), 
                 fdParobj = pen.obj)$fd
  })
})

x11()
par(mfrow = c(4,2))
for(i in 1:8) 
{
  plot(variable.dfs.FD.byset$AbsDirErr[[i]][1,])
  points(1:set.len, as.matrix(variable.dfs.byset$AbsDirErr[[i]][,1]), 
         pch = 20, col = "blue")
}

# Join the sets in a plot
smooth.values1 <- lapply(variable.dfs.FD.byset$AbsDirErr, function(x)
  {eval.fd(1:set.len, x[1,])})
aa<-do.call("c", smooth.values1) %>% as.numeric()
plot(1:784, aa, type="l", lwd = 2)
abline(v = set.delim, col = "red")
points()

loglam = seq(4,9,0.25)
nlam = length(loglam)
dfsave = rep(NA,nlam)
gcvsave = rep(NA,nlam)
for (ilam in 1:nlam) 
{
  cat(paste("log10 lambda =",loglam[ilam], "\n"))
  lambda = 10^loglam[ilam]
  fdParobj = fdPar(daybasis, harmaccelLfd, lambda)
  smoothlist = smooth.basis(day.5, logprecav,
                            fdParobj)
  dfsave[ilam] = smoothlist$df
  gcvsave[ilam] = sum(smoothlist$gcv)
}

################################################################################
### SMOOTHING OF THE ENTIRE TIME SERIES
# Due to numeric problems with linear regression, we decided to smooth the data
# in order to just consider their global trend. 
# So we do smoothing with strong penalization.

# Define the basis system
rangeval <- c(1,data.len)
abscissa <- 1:data.len
# # Select n.knots internal knots for each set
# n.knots <- 5
# knotsvec <- lapply(set.delim[-length(set.delim)], function(x)
# {
#   knots <- seq(1, set.len, length.out = n.knots + 2)
#   # Drop external knots
#   knots <- knots[-c(1,n.knots+2)]
#   # Add the set position
#   knots <- knots + x
#   return(knots)
# })
# knotsvec <- unlist(knotsvec)
# knotsvec <- sort(c(knotsvec, rep(set.delim[-c(1, length(set.delim))], 3)))
# basis.set <- create.bspline.basis(rangeval = rangeval, norder = 4,
#                                   breaks = knotsvec)
# pen.obj <- fdPar(fdobj = basis.set)
basis.set <- create.bspline.basis(rangeval = rangeval, norder = 4,
                                  nbasis = 100)
pen.obj <- fdPar(fdobj = basis.set, Lfdobj = 2, lambda = 10^(6))

# To solve numerical issues with fRegress function we decided to fit the model
# on the exponentially transformed data
variable.dfs.exp <- lapply(variable.dfs, exp)

variable.dfs.FD <- lapply(variable.dfs.exp, function(x)
{
  smooth.basis(argvals = abscissa, 
               y = as.matrix(x), 
               fdParobj = pen.obj)$fd
               # fdParobj = basis.set)$fd
})
# n.b. This is the list of FD covariates for the model

for(i in 1:N)
{
  plot(variable.dfs.FD$BetaMod[i,])
  abline(v = set.delim, col = "red", lty = 2, lwd = 1.5)
  # points(1:data.len, as.matrix(variable.dfs$BetaMod[,i]), pch = 16, col = "blue")
  Sys.sleep(3)
}
