########################
### BUILDING THE FUNCTIONAL LINEAR MODEL

# N.B. For the moment we stick with a CONCURRENT flm. 
# (For theoretical analysis, see 10.2 book)

# N.B. We consider BetaMod as the response variable
# So we need to drop ERS

# List of Covariates
cov.fd.list <- variable.dfs.FD[-5]
# cov.fd.list <- cov.fd.list$AbsArea
# cov.fd.list <- list(Intercept = rep(1, N), AbsDirErr = cov.fd.list)
cov.fd.list <- append(cov.fd.list, NA, 0)
cov.fd.list[[1]] <- rep(1, N)
names(cov.fd.list)[1] <- "Intercept"
names(cov.fd.list)
n.covars <- 5 # intercept included

# Response variable
BetaMod.fd <- variable.dfs.FD$BetaMod
# ERS.fd <- variable.dfs.FD$ERS

# List of Beta parameters for each variable
beta.fdPar <- fdPar(basis.set)
beta.list <- list(beta.fdPar, beta.fdPar, beta.fdPar, beta.fdPar, beta.fdPar)
# beta.list <- list(beta.fdPar, beta.fdPar)
# glimpse(beta.list)

# Fit the model
model.list <- fRegress(y = BetaMod.fd, 
                       xfdlist = cov.fd.list, 
                       betalist = beta.list)

# Get the fitted values
BetaMod.fd.fitted <- model.list$yhatfdobj # (y_hat)_i, i = 1,...,23
plot(BetaMod.fd.fitted)

# Get the estimated beta
beta.est.list <- model.list$betaestlist
# glimpse(beta.est.list)

par(mfrow = c(2,3))
for(i in 1:5)
{
  plot(beta.est.list[[i]])
}

# Compute some model diagnostics metrics
# 1) Residual variance-covariance matrix estimate
BetaMod.mat <- eval.fd(evalarg = abscissa, fdobj = BetaMod.fd) # smoothed values
BetaMod.mean <- rowMeans(BetaMod.mat)
BetaMod.fitted.mat <- eval.fd(evalarg = abscissa, fdobj = BetaMod.fd.fitted) # fitted values
# res.mat <- variable.dfs$BetaMod - BetaMod.fitted.mat # residuals wrt original values
res.mat <- BetaMod.mat - BetaMod.fitted.mat # MODEL RESIDUALS
SigmaE <- cov(t(res.mat))

# Compute SSE and R.sq
res.mat0 <- BetaMod.mat - BetaMod.mean %*% matrix(1, 1, N) # SMOOTHING RESIDUALS
SSE0 <- apply((res.mat0)^2, 1, sum)
SSE1 <- apply((res.mat)^2, 1, sum)
Rsqr <- (SSE0 - SSE1)/SSE0
plot(Rsqr, type="l")
abline(h=0, lty=2)

# Rsqr<0 significa che "la variabilità dei dati fittati con il modello lineare
# è maggiore rispetto alla variabilità originale dei dati".

BetaMod.y2cMap <- smooth.basis(argvals = abscissa, 
                               y = as.matrix(variable.dfs$BetaMod), 
                               fdParobj = basis.set)
BetaMod.y2cMap <- BetaMod.y2cMap$y2cMap
model2 <- fRegress.stderr(model.list, y2cMap = BetaMod.y2cMap, SigmaE = SigmaE)
beta.stderr.list <- model2$betastderrlist

# Build confidence intervals manually to do better plots
beta.se.eval.list <- lapply(beta.stderr.list, function(x)
{
  x <- eval.fd(evalarg = abscissa, fdobj = x) %>% as.numeric()
})
beta.eval.list <- lapply(beta.est.list, function(x)
{
  x <- eval.fd(evalarg = abscissa, fdobj = x$fd) %>% as.numeric() 
})

# Compute the confidence intervals
beta.CI.list <- mapply(x = beta.eval.list, y = beta.se.eval.list, 
  function(x,y)
  {
    tibble(value = x, 
           lower = x - qnorm(0.975)*y, 
           upper = x + qnorm(0.975)*y)
  }, SIMPLIFY = F) 
glimpse(beta.CI.list)

# Build a dataframe for ggplot2
beta.CI.vect <- bind_rows(beta.CI.list)
cov.names <- c("Intercept", 
               "Absolute Area", 
               "Absolute Directional Error", 
               "Reaction Time", 
               "Peak Velocity")
beta.gg.df <- tibble(beta.CI.vect, 
                     abscissas = rep(1:data.len, n.covars),
                     coeff = factor(rep(cov.names, rep(data.len, n.covars)), 
                                    levels=cov.names))
beta.gg.df
# Create ggplot object
gggg <- ggplot(beta.gg.df, aes(x = abscissas, y = value)) + 
  geom_hline(yintercept = 0, col = "red", size = 0.5, linetype = "dashed") +
  # geom_line(aes(y = lower)) +
  # geom_line(aes(y = upper)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = coeff), alpha=0.4) +
  geom_line(aes(y = value)) +
  facet_wrap(facets = vars(coeff), scales = "free") + 
  labs(title = "Confidence Intervals (95%) for the Model Coefficients",
       subtitle = "Complete model",
       x = NULL, y=NULL) +
  theme(legend.position = "bottom")+
  theme(plot.title=element_text(size=24, face="bold")) +
  theme(plot.subtitle=element_text(size=18, face ="italic")) +
  theme(plot.title=element_text(hjust=0)) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=15)) +
  theme(strip.text=element_text(size=12,face ="bold")) +
  theme(legend.key.size = unit(1,"cm")) +
  scale_fill_discrete(labels=paste("β",0:(n.covars-1),sep="")) +
  guides(colour = guide_legend(override.aes = list(size=5)))
gggg


################################################################################
### MODELS WITH ONLY ONE COVARIATE

### 1) Only AbsArea
# List of Covariates
cov.fd.list <- variable.dfs.FD[-5]
cov.fd.list <- cov.fd.list$PVel
cov.fd.list <- list(Intercept = rep(1, N), PVel = cov.fd.list)
# cov.fd.list <- append(cov.fd.list, NA, 0)
# cov.fd.list[[1]] <- rep(1, N)
# names(cov.fd.list)[1] <- "Intercept"
names(cov.fd.list)
n.covars <- 2 # intercept included

# List of Beta parameters for each variable
beta.list <- list(beta.fdPar, beta.fdPar)

# Fit the model
model.list <- fRegress(y = BetaMod.fd, 
                       xfdlist = cov.fd.list, 
                       betalist = beta.list)

# Get the fitted values
BetaMod.fd.fitted <- model.list$yhatfdobj # (y_hat)_i, i = 1,...,23
plot(BetaMod.fd.fitted)

# Get the estimated beta
beta.est.list <- model.list$betaestlist
# glimpse(beta.est.list)

par(mfrow = c(1,2))
for(i in 1:2)
{
  plot(beta.est.list[[i]])
}

# Compute some model diagnostics metrics
# 1) Residual variance-covariance matrix estimate
BetaMod.mat <- eval.fd(evalarg = abscissa, fdobj = BetaMod.fd) # smoothed values
BetaMod.mean <- rowMeans(BetaMod.mat)
BetaMod.fitted.mat <- eval.fd(evalarg = abscissa, fdobj = BetaMod.fd.fitted) # fitted values
# res.mat <- variable.dfs$BetaMod - BetaMod.fitted.mat # residuals wrt original values
res.mat <- BetaMod.mat - BetaMod.fitted.mat # MODEL RESIDUALS
SigmaE <- cov(t(res.mat))

# Compute SSE and R.sq
res.mat0 <- BetaMod.mat - BetaMod.mean %*% matrix(1, 1, N) # SMOOTHING RESIDUALS
SSE0 <- apply((res.mat0)^2, 1, sum)
SSE1 <- apply((res.mat)^2, 1, sum)
Rsqr <- (SSE0 - SSE1)/SSE0
plot(Rsqr, type="l")
abline(h=0, lty=2)

# Comppute the std err
BetaMod.y2cMap <- smooth.basis(argvals = abscissa, 
                               y = as.matrix(variable.dfs$BetaMod), 
                               fdParobj = basis.set)
BetaMod.y2cMap <- BetaMod.y2cMap$y2cMap
model2 <- fRegress.stderr(model.list, y2cMap = BetaMod.y2cMap, SigmaE = SigmaE)
beta.stderr.list <- model2$betastderrlist

# Build confidence intervals manually to do better plots
beta.se.eval.list <- lapply(beta.stderr.list, function(x)
{
  x <- eval.fd(evalarg = abscissa, fdobj = x) %>% as.numeric()
})
beta.eval.list <- lapply(beta.est.list, function(x)
{
  x <- eval.fd(evalarg = abscissa, fdobj = x$fd) %>% as.numeric() 
})

# Compute the confidence intervals
beta.CI.list <- mapply(x = beta.eval.list, y = beta.se.eval.list, 
                       function(x,y)
                       {
                         tibble(value = x, 
                                lower = x - qnorm(0.975)*y, 
                                upper = x + qnorm(0.975)*y)
                       }, SIMPLIFY = F) 
glimpse(beta.CI.list)

# Build a dataframe for ggplot2
beta.CI.vect <- bind_rows(beta.CI.list)
cov.names <- c("Intercept", 
               "Peak Velocity")
beta.gg.df <- tibble(beta.CI.vect, 
                     abscissas = rep(1:data.len, n.covars),
                     coeff = factor(rep(cov.names, rep(data.len, n.covars)), 
                                    levels=cov.names))
beta.gg.df
# Create ggplot object
gggg <- ggplot(beta.gg.df, aes(x = abscissas, y = value)) + 
  geom_hline(yintercept = 0, col = "red", size = 0.5, linetype = "dashed") +
  # geom_line(aes(y = lower)) +
  # geom_line(aes(y = upper)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = coeff), alpha=0.4) +
  geom_line(aes(y = value)) +
  facet_wrap(facets = vars(coeff), scales = "free") + 
  labs(title = "Confidence Intervals (95%) for the Model Coefficients",
       subtitle = "One Covariate - Peak Velocity",
       x = NULL, y=NULL) +
  theme(legend.position = "bottom")+
  theme(plot.title=element_text(size=24, face="bold")) +
  theme(plot.subtitle=element_text(size=18, face ="italic")) +
  theme(plot.title=element_text(hjust=0)) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=15)) +
  theme(strip.text=element_text(size=12,face ="bold")) +
  theme(legend.key.size = unit(1,"cm")) +
  scale_fill_manual(labels=paste("β",0:(n.covars-1),sep=""),
    values = c(hue_pal()(5)[1], hue_pal()(5)[5])) +
  guides(colour = guide_legend(override.aes = list(size=5)))
gggg




