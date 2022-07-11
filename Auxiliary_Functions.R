library(dplyr)
library(openxlsx)
library(fda)
library(readxl)
library(gtools)
library(tibble)
library(tictoc)
library(imputeTS)
library(fdatest)
library(progress)
library(pbapply)
library(parallel)
library(roahd)

split.sets <- function(data, set.length)
{
  n.sets <- nrow(data) / set.length
  n <- ncol(data) # number of subjects
  
  # Transform the data in the "SET STRUCTURE"
  # Factor identifying the groups (in the sense of the observations)
  group.fct <- factor(rep(1:n.sets, rep(set.length, n.sets)))
  
  # Separate the dataframe in sub-dfs, one for each set
  data.split <- split(data, group.fct) 
  names(data.split) <- paste("Set.", 1:n.sets, sep = "")
  
  return(data.split)
}

perm.func.anova <- function(data, niter = 10000, n.by.group, l)
{
  n.groups <- length(n.by.group)
  n <- ncol(data)
  
  # Transform the data in the structure described above
  # Factor identifying the groups (in the sense of the observations)
  group.fct <- factor(rep(1:n.groups, rep(l, n.groups)))
  # Separate the dataframe in sub-dfs, one for each set (this will simplify permutations)
  data.split <- split(data, group.fct) 
  data.trsf <- lapply(data.split, t) # each row corresponds to a function
  data.trsf <- lapply(data.trsf, as.data.frame)
  data.trsf <- lapply(data.trsf, 
                      function(x) {names(x) <- paste("v",1:l, sep=""); return(x)})
  subj.names <- row.names(data.trsf[[1]])
  data.all <- bind_rows(data.trsf)
  
  # Compute the test statistic for the original sample  
  set.fct <- rep(1:n.groups, rep(n, n.groups)) %>% as.factor()
  fit <- manova(as.matrix(data.all) ~ set.fct)
  T0 <- - summary.manova(fit, test="Wilks")$stats[1,2]
  
  T.stat <- numeric(niter)
  pb <- progress_bar$new(total = niter)
  pb$tick(0)
  for(k in 1:niter)
  {
    # Permutation of data
    # We have to switch the units between the groups
    # e.g. units#1 in group 1 -> unit#1 in group 5
    #      unit#2 in group 1 -> unit#1 in group 7 (and so on...)
    perm.idxs <- replicate(n = n, expr = sample(1:n.groups))
    perm.data <- vector("list", n.groups)
    # Initialize dfs inside perm.data object
    perm.data <- lapply(perm.data, function(x) 
    { as.data.frame(matrix(0, nrow = n, ncol = l), row.names = r.names) })
    for(i in 1:n)
    {
      for(j in 1:n.groups)
      {
        curr.idx <- perm.idxs[j,i] 
        perm.data[[curr.idx]][i,] <- data.trsf[[j]][i,]
      }
    } 
    perm.data.all <- bind_rows(perm.data)
    
    # Compute the Test Stat for the permuted values
    fit.perm <- manova(as.matrix(perm.data.all) ~ set.fct)
    T.stat[k] <- - summary.manova(fit.perm, test="Wilks")$stats[1,2]
    
    pb$tick()
  }
  
  return(list(p.val = sum(T.stat >= T0)/niter,
              T.stats = T.stat,
              T0 = T0))
}

IWT1 <- function(data,mu=0,B=1000,dx=NULL,recycle=TRUE){
  pval.correct <- function(pval.matrix){
    matrice_pval_2_2x <- cbind(pval.matrix,pval.matrix)
    p <- dim(pval.matrix)[2]
    matrice_pval_2_2x <- matrice_pval_2_2x[,(2*p):1]
    corrected.pval.matrix <- matrix(nrow=p,ncol=p)
    corrected.pval.matrix[p,] <- pval.matrix[p,p:1]
    for(var in 1:p){
      pval_var <- matrice_pval_2_2x[p,var]
      inizio <- var
      fine <- var #inizio fisso, fine aumenta salendo nelle righe
      for(riga in (p-1):1){
        fine <- fine + 1
        pval_cono <- matrice_pval_2_2x[riga,inizio:fine]
        pval_var <- max(pval_var,pval_cono,na.rm=TRUE)
        corrected.pval.matrix[riga,var] <- pval_var
      }
    }
    corrected.pval.matrix <- corrected.pval.matrix[,p:1]
    return(corrected.pval.matrix)
  }
  
  # data preprocessing
  if(is.fd(data)){ # data is a functional data object
    rangeval <- data$basis$rangeval
    if(is.null(dx)){
      dx <- (rangeval[2]-rangeval[1])*0.01
    }
    abscissa <- seq(rangeval[1],rangeval[2],by=dx)
    coeff <- t(eval.fd(fdobj=data,evalarg=abscissa))
  }else if(is.matrix(data)){
    coeff <- data
  }else{
    stop("First argument must be either a functional data object or a matrix.")
  }
  
  if (is.fd(mu)){ # mu is a functional data
    rangeval.mu <- mu$basis$rangeval
    if(sum(rangeval.mu == rangeval)!=2){
      stop("rangeval of mu must be the same as rangeval of data.")
    }
    if(is.null(dx)){
      dx <- (rangeval.mu[2]-rangeval.mu[1])*0.01
    }
    abscissa <- seq(rangeval.mu[1],rangeval.mu[2],by=dx)
    mu.eval <- t(eval.fd(fdobj=mu,evalarg=abscissa))
  }else if(is.vector(mu)){
    mu.eval <- mu
  }else{
    stop("Second argument must be either a functional data object or a numeric vector.")
  }
  
  n <- dim(coeff)[1]
  p <- dim(coeff)[2]
  data.eval <- coeff <- coeff - matrix(data=mu.eval,nrow=n,ncol=p,byrow=TRUE)
  
  #univariate permutations
  print('Point-wise tests')
  T0 <- abs(colMeans(coeff))^2  #sample mean
  T_coeff <- matrix(ncol=p,nrow=B)
  for (perm in 1:B){
    signs <- rbinom(n,1,0.5)*2 - 1
    coeff_perm <- coeff*signs
    T_coeff[perm,] <- abs(colMeans(coeff_perm))^2
  }
  pval <- numeric(p)
  for(i in 1:p){
    pval[i] <- sum(T_coeff[,i]>=T0[i])/B
  }
  
  #combination
  print('Interval-wise tests')
  
  #asymmetric combination matrix:
  matrice_pval_asymm <- matrix(nrow=p,ncol=p)
  matrice_pval_asymm[p,] <- pval[1:p]
  T0_2x <- c(T0,T0)
  T_coeff_2x <- cbind(T_coeff,T_coeff)
  
  maxrow <- 1
  # con parametro scale
  #maxrow <- p-scale+1
  
  if(recycle==TRUE){
    for(i in (p-1):maxrow){ # rows
      for(j in 1:p){ # columns
        inf <- j
        sup <- (p-i)+j
        T0_temp <- sum(T0_2x[inf:sup])
        T_temp <- rowSums(T_coeff_2x[,inf:sup])
        pval_temp <- sum(T_temp>=T0_temp)/B
        matrice_pval_asymm[i,j] <- pval_temp
      }
      print(paste('creating the p-value matrix: end of row ',as.character(p-i+1),' out of ',as.character(p),sep=''))
    }
  }else{ # without recycling
    for(i in (p-1):maxrow){ # rows
      for(j in 1:i){ # columns
        inf <- j
        sup <- (p-i)+j
        T0_temp <- sum(T0_2x[inf:sup])
        T_temp <- rowSums(T_coeff_2x[,inf:sup])
        pval_temp <- sum(T_temp>=T0_temp)/B
        matrice_pval_asymm[i,j] <- pval_temp
      }
      print(paste('creating the p-value matrix: end of row ',as.character(p-i+1),' out of ',as.character(p),sep=''))
    }
  }
  
  corrected.pval.matrix <- pval.correct(matrice_pval_asymm)
  corrected.pval <- corrected.pval.matrix[1,]
  
  print('Interval-Wise Testing completed')
  IWT.result <- list(
    test = '1pop', mu = mu.eval,
    adjusted_pval = corrected.pval,
    unadjusted_pval = pval,
    pval_matrix = matrice_pval_asymm,
    data.eval=data.eval)
  class(IWT.result) = 'IWT1'
  return(IWT.result)
}

IWT2 <- function(data1,data2,mu=0,B=1000,paired=FALSE,dx=NULL,recycle=TRUE,alternative="two.sided"){
  pval.correct <- function(pval.matrix){
    matrice_pval_2_2x <- cbind(pval.matrix,pval.matrix)
    p <- dim(pval.matrix)[2]
    matrice_pval_2_2x <- matrice_pval_2_2x[,(2*p):1]
    corrected.pval <- numeric(p)
    corrected.pval.matrix <- matrix(nrow=p,ncol=p)
    corrected.pval.matrix[p,] <- pval.matrix[p,p:1]
    for(var in 1:p){
      pval_var <- matrice_pval_2_2x[p,var]
      inizio <- var
      fine <- var #inizio fisso, fine aumenta salendo nelle righe
      for(riga in (p-1):1){
        fine <- fine + 1
        pval_cono <- matrice_pval_2_2x[riga,inizio:fine]
        pval_var <- max(pval_var,pval_cono,na.rm=TRUE)
        corrected.pval.matrix[riga,var] <- pval_var
      }
      corrected.pval[var] <- pval_var
    }
    corrected.pval <- corrected.pval[p:1]
    corrected.pval.matrix <- corrected.pval.matrix[,p:1]
    return(corrected.pval.matrix)
  }
  
  possible_alternatives <- c("two.sided", "less", "greater")
  if(!(alternative %in% possible_alternatives)){
    stop(paste0('Possible alternatives are ',paste0(possible_alternatives,collapse=', ')))
  }
  
  # data preprocessing
  if(is.fd(data1)){ # data1 is a functional data object
    rangeval1 <- data1$basis$rangeval
    rangeval2 <- data2$basis$rangeval
    if(is.null(dx)){
      dx <- (rangeval1[2]-rangeval1[1])*0.01
    }
    if(sum(rangeval1 == rangeval2)!=2){
      stop("rangeval of data1 and data2 must coincide.")
    }
    abscissa <- seq(rangeval1[1],rangeval1[2],by=dx)
    coeff1 <- t(eval.fd(fdobj=data1,evalarg=abscissa))
    coeff2 <- t(eval.fd(fdobj=data2,evalarg=abscissa))
    
  }else if(is.matrix(data1)){
    coeff1 <- data1
    coeff2 <- data2
  }else{
    stop("First argument must be either a functional data object or a matrix.")
  }
  
  if (is.fd(mu)){ # mu is a functional data
    rangeval.mu <- mu$basis$rangeval
    if(sum(rangeval.mu == rangeval1)!=2){
      stop("rangeval of mu must be the same as rangeval of data.")
    }
    if(is.null(dx)){
      dx <- (rangeval.mu[2]-rangeval.mu[1])*0.01
    }
    abscissa <- seq(rangeval.mu[1],rangeval.mu[2],by=dx)
    mu.eval <- t(eval.fd(fdobj=mu,evalarg=abscissa))
  }else if(is.vector(mu)){
    mu.eval <- mu
  }else{
    stop("Second argument must be either a functional data object or a numeric vector.")
  }
  
  n1 <- dim(coeff1)[1]
  n2 <- dim(coeff2)[1]
  p <- dim(coeff1)[2]
  n <- n1+n2
  etichetta_ord <- c(rep(1,n1),rep(2,n2))
  coeff1 <- coeff1 - matrix(data=mu.eval,nrow=n1,ncol=p)
  
  #print('First step: basis expansion')
  #splines coefficients:
  eval <- coeff <- rbind(coeff1,coeff2)
  
  data.eval <- eval
  data.eval[1:n1,] <- data.eval[1:n1,] + matrix(data=mu.eval,nrow=n1,ncol=p)
  
  print('Point-wise tests')
  #univariate permutations
  meandiff <- colMeans(coeff[1:n1,,drop=FALSE],na.rm=TRUE) - colMeans(coeff[(n1+1):n,,drop=FALSE],na.rm=TRUE)
  sign.diff <- sign(meandiff)
  sign.diff[which(sign.diff==-1)] <- 0
  T0 <- switch(alternative,
               two.sided =  (meandiff)^2,
               greater   =  (meandiff*sign.diff)^2,
               less      =  (meandiff*(sign.diff-1))^2)
  
  T_coeff <- matrix(ncol=p,nrow=B)
  for (perm in 1:B){
    if(paired){
      if.perm <- rbinom(n1,1,0.5)
      coeff_perm <- coeff
      for(couple in 1:n1){
        if(if.perm[couple]==1){
          coeff_perm[c(couple,n1+couple),] <- coeff[c(n1+couple,couple),]
        }
      }
    }else{
      permutazioni <- sample(n)
      coeff_perm <- coeff[permutazioni,]
    }
    meandiff <- colMeans(coeff_perm[1:n1,,drop=FALSE],na.rm=TRUE) - colMeans(coeff_perm[(n1+1):n,,drop=FALSE],na.rm=TRUE)
    sign.diff <- sign(meandiff)
    sign.diff[which(sign.diff==-1)] <- 0
    T_coeff[perm,] <- switch(alternative,
                             two.sided =  (meandiff)^2,
                             greater   =  (meandiff*sign.diff)^2,
                             less      =  (meandiff*(sign.diff-1))^2)
  }
  pval <- numeric(p)
  for(i in 1:p){
    pval[i] <- sum(T_coeff[,i]>=T0[i])/B
  }
  
  #combination
  print('Interval-wise tests')
  
  #asymmetric combination matrix:
  matrice_pval_asymm <- matrix(nrow=p,ncol=p)
  matrice_pval_asymm[p,] <- pval[1:p]
  T0_2x <- c(T0,T0)
  T_coeff_2x <- cbind(T_coeff,T_coeff)
  
  maxrow <- 1
  # con parametro scale
  #maxrow <- p-scale+1
  
  if(recycle==TRUE){
    for(i in (p-1):maxrow){ # rows
      for(j in 1:p){ # columns
        inf <- j
        sup <- (p-i)+j
        T0_temp <- sum(T0_2x[inf:sup])
        T_temp <- rowSums(T_coeff_2x[,inf:sup])
        pval_temp <- sum(T_temp>=T0_temp)/B
        matrice_pval_asymm[i,j] <- pval_temp
      }
      print(paste('creating the p-value matrix: end of row ',as.character(p-i+1),' out of ',as.character(p),sep=''))
    }
  }else{ # without recycling
    for(i in (p-1):maxrow){ # rows
      for(j in 1:i){ # columns
        inf <- j
        sup <- (p-i)+j
        T0_temp <- sum(T0_2x[inf:sup])
        T_temp <- rowSums(T_coeff_2x[,inf:sup])
        pval_temp <- sum(T_temp>=T0_temp)/B
        matrice_pval_asymm[i,j] <- pval_temp
      }
      print(paste('creating the p-value matrix: end of row ',as.character(p-i+1),' out of ',as.character(p),sep=''))
    }
  }
  
  corrected.pval.matrix <- pval.correct(matrice_pval_asymm)
  corrected.pval <- corrected.pval.matrix[1,]
  
  print('Interval-Wise Testing completed')
  IWT.result <- list(
    test = '2pop', mu = mu.eval,
    adjusted_pval = corrected.pval,
    unadjusted_pval = pval,
    pval_matrix = matrice_pval_asymm,
    data.eval=data.eval,
    ord_labels = etichetta_ord)
  class(IWT.result) = 'IWT2'
  return(IWT.result)
}

plot.IWT1 <- function(x, xrange = c(0,1),
                      alpha1 = 0.05, alpha2 = 0.01,
                      ylab = 'Functional Data', main = NULL, 
                      lwd = 1, col = 1, 
                      ylim = NULL,type='l', ...) {
  if (class(x) != "IWT1") stop("First argument is not a IWT1 object.")
  if (alpha1 < alpha2) {
    temp <- alpha1
    alpha1 <- alpha2
    alpha2 <- temp
  }
  object <- x
  
  devAskNewPage(ask = TRUE)
  p <- length(object$unadjusted_pval)
  n <- dim((object$data.eval))[1]
  xmin <- xrange[1]
  xmax <- xrange[2]
  abscissa_pval <- seq(xmin, xmax, len = p)
  main_data <- paste(main, ': Functional Data')
  main_data <- sub("^ : +", "", main_data)
  n_coeff <- dim(object$data.eval)[2]
  data_eval <- object$data.eval
  if (is.null(ylim)) ylim <- range(data_eval)
  matplot(abscissa_pval, t(data_eval), type='l', main = main_data, 
          ylab = ylab, col = col, lwd = lwd, ylim = ylim, ...)
  difference1 <- which(object$adjusted_pval < alpha1)
  if (length(difference1) > 0) {
    for (j in 1:length(difference1)) {
      min_rect <- abscissa_pval[difference1[j]] - (abscissa_pval[2] - abscissa_pval[1]) / 2
      max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
      rect(min_rect, par("usr")[3], max_rect, par("usr")[4], col = "gray90", 
           density = -2, border = NA)
    }
    rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], 
         col = NULL, border = "black")
  }
  difference2 <- which(object$adjusted_pval < alpha2)
  if (length(difference2) > 0) {
    for (j in 1:length(difference2)) {
      min_rect <- abscissa_pval[difference2[j]] - (abscissa_pval[2] - abscissa_pval[1]) / 2
      max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
      rect(min_rect, par("usr")[3], max_rect, par("usr")[4], col = "gray80", 
           density = -2, border = NA)
    }
    rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], 
         col = NULL, border = "black")
  }
  matplot(abscissa_pval, t(data_eval), type = 'l', main = main_data, 
          ylab = ylab, col = col, lwd = lwd, add = TRUE, ...)
  if (length(object$mu) == 1) { # mu is a constant function
    mu_eval <- rep(object$mu, p)
  } else { # mu is a functional data with no constant coefficients
    mu <- object$mu
    mu_eval <- mu
  }
  abscissa_mu <- abscissa_pval
  lines(abscissa_mu, mu_eval, col = 'gray', lwd = 2)
  # Plot adjusted p-values
  main_p <- paste(main,': Adjusted p-values')
  main_p <- sub("^ : +", "", main_p)
  plot(abscissa_pval, object$adjusted_pval, ylim = c(0, 1), 
       main = main_p, ylab = 'p-value', type=type,lwd=lwd, ...)
  difference1 <- which(object$adjusted_pval < alpha1)
  if (length(difference1) > 0) {
    for (j in 1:length(difference1)) {
      min_rect <- abscissa_pval[difference1[j]] - (abscissa_pval[2] - abscissa_pval[1]) / 2
      max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
      rect(min_rect, par("usr")[3], max_rect, par("usr")[4], col = "gray90", 
           density = -2, border = NA)
    }
    rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], 
         col = NULL, border = "black")
  }
  difference2 <- which(object$adjusted_pval < alpha2)
  if (length(difference2) > 0) {
    for (j in 1:length(difference2)) {
      min_rect <- abscissa_pval[difference2[j]] - (abscissa_pval[2] - abscissa_pval[1]) / 2
      max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
      rect(min_rect, par("usr")[3], max_rect, par("usr")[4], col = "gray80", 
           density = -2, border = NA)
    }
    rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], 
         col = NULL, border = "black")
  }
  for (j in 0:10) {
    abline(h = j/10, col = 'lightgray', lty = "dotted")
  }
  points(abscissa_pval, object$adjusted_pval,type=type,lwd=lwd)
  
  
  devAskNewPage(ask = FALSE)
}

# Since for some patients we have consecutive missing values, in order to avoid
# major fitting mistakes we place the knots adaptively.
# SCHEME: sequentially group available taskNo values in group of k. Put a knot 
# at the average of each k values.
# Exception: if the size of the last group is less than k, it is merged 
# the previous group.
# Notation: - x is the vector on which the knots are taken;
#           - k is the size of the groups
#           - h is the number of point put at set delimiters
adapt.knotsvec <- function(x, k, h = 2)
{
  # Split the vector in groups of size k
  groups <- split(x, ceiling(seq_along(x)/k))
  n <- length(groups)
  
  # Check the size of the last group
  if(length(groups[[n]]) < k)
  {
    groups[[n-1]] <- c(groups[[n-1]], groups[[n]])
    groups <- groups[-n]
  }
  
  # Compute the mean inside each group to get the knots location
  knots <- lapply(groups, mean) %>% unlist()
  
  # Add additional knots in correspondence to set delimiters avoid smooth junction
  knots <- sort(c(knots, rep(set.delim[-c(1, length(set.delim))], 2)))
  
  return(knots)
}

plot.IWT2 <- function(x, xrange = c(0,1),
                      alpha1 = 0.05, alpha2 = 0.01,
                      ylab = 'Functional Data', main = NULL, 
                      lwd = 0.5, col=c(1,2), 
                      ylim = NULL, type='l', ...) {
  if (class(x) != "IWT2") stop("First argument is not a IWT2 object.")
  if (alpha1 < alpha2) {
    temp <- alpha1
    alpha1 <- alpha2
    alpha2 <- temp
  }
  object <- x
  n <- dim(t(object$data.eval))[1]
  
  colors <- numeric(n)
  id_pop1 <- unique(object$ord_labels)[1]
  id_pop2 <- unique(object$ord_labels)[2]
  colors[which(object$ord_labels == id_pop1)] <- col[1]
  colors[which(object$ord_labels == id_pop2)] <- col[2]
  
  devAskNewPage(ask = TRUE) 
  
  p <- length(object$unadjusted_pval)
  xmin <- xrange[1]
  xmax <- xrange[2]
  abscissa_pval = seq(xmin, xmax, len = p)
  main_data <- paste(main, ': Functional Data')
  main_data <- sub("^ : +", "", main_data)
  n_coeff <- dim(object$data.eval)[2]
  data_eval <- object$data.eval
  if (is.null(ylim)) ylim <- range(data_eval,na.rm=TRUE)
  matplot(abscissa_pval, t(data_eval), type = 'l', main = main_data, 
          ylab = ylab, col = colors, lwd = lwd, ylim = ylim, ...)
  mean1 = colMeans(object$data.eval[which(object$ord_labels==id_pop1),],na.rm=TRUE)
  mean2 = colMeans(object$data.eval[which(object$ord_labels==id_pop2),],na.rm=TRUE)
  
  difference1 <- which(object$adjusted_pval < alpha1)
  if (length(difference1) > 0) {
    for (j in 1:length(difference1)) {
      min_rect <- abscissa_pval[difference1[j]] - (abscissa_pval[2] - abscissa_pval[1])/2
      max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
      rect(min_rect, par("usr")[3], max_rect, par("usr")[4], col = "gray90", 
           density = -2, border = NA)
    }
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
         col = NULL, border = "black")
  }
  difference2 <- which(object$adjusted_pval < alpha2)
  if (length(difference2) > 0) {
    for (j in 1:length(difference2)) {
      min_rect <- abscissa_pval[difference2[j]] - (abscissa_pval[2] - abscissa_pval[1])/2
      max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
      rect(min_rect, par("usr")[3], max_rect, par("usr")[4], col = "gray80", 
           density = -2, border = NA)
    }
    rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], 
         col = NULL, border = "black")
  }
  matplot(abscissa_pval, t(data_eval), type = 'l', main = main_data,
          ylab = ylab, col = colors, lwd = lwd, add = TRUE, ...)
  #matlines(abscissa_pval,cbind(mean1,mean2),col=col,lwd=2,lty=1)
  
  #  adjusted p-values
  main_p <- paste(main,': Adjusted p-values')
  main_p <- sub("^ : +", "", main_p)
  plot(abscissa_pval, object$adjusted_pval, ylim = c(0, 1),
       main = main_p, ylab = 'p-value', type=type, lwd=lwd,...)
  difference1 <- which(object$adjusted_pval < alpha1)
  if (length(difference1) > 0) {
    for (j in 1:length(difference1)) {
      min_rect <- abscissa_pval[difference1[j]] - (abscissa_pval[2] - abscissa_pval[1])/2
      max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
      rect(min_rect, par("usr")[3], max_rect, par("usr")[4], col = "gray90", 
           density = -2, border = NA)
    }
    rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], 
         col = NULL, border = "black")
  }
  difference2 <- which(object$adjusted_pval < alpha2)
  if (length(difference2) > 0) {
    for (j in 1:length(difference2)) {
      min_rect <- abscissa_pval[difference2[j]] - (abscissa_pval[2] - abscissa_pval[1])/2
      max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
      rect(min_rect, par("usr")[3], max_rect, par("usr")[4], col = "gray80", 
           density = -2, border = NA)
    }
    rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col = NULL, border = "black")
  }
  for (j in 0:10) {
    abline(h = j / 10, col = 'lightgray', lty = "dotted")
  }
  points(abscissa_pval, object$adjusted_pval, type=type,lwd=2)
  
  devAskNewPage(ask = FALSE)
}

