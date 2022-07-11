################################################################################
# Set up a cluster for parallel computations
clust <- makeCluster(12)

################################################################################
# 1) Test NAP vs AWAKE
# Considering the whole functional data, we run a Interval-Wise Test to check 
# whether there are differences between the patients who took a nap and the 
# subjects who stayed awake.
# N.B. Test one variable at the time

# Prepare data: AWAKE data vs NAP data
variable.dfs.AWK <- lapply(variable.dfs, function(x)
{
  x %>% dplyr::select(contains("A"))
})
# glimpse(variable.dfs.AWK)
variable.dfs.NAP <- lapply(variable.dfs, function(x)
{
  x %>% dplyr::select(contains("N"))
})
# glimpse(variable.dfs.NAP)

# Join in a single data structure to use parallel computation
# n.b. Structure: Variables -> NAP & AWK -> dataframes
variable.dfs.AN <- vector("list", n.vars)
names(variable.dfs.AN) <- names(variable.dfs)
# Initialize the nested lists
variable.dfs.AN <- lapply(variable.dfs.AN, function(x)
{
  x <- vector("list", 2)
  names(x) <- c("AWK", "NAP")
  return(x)
})
# variable.dfs.AN <- mapply(x = variable.dfs.AN, y = variable.dfs.AWK, z = variable.dfs.NAP,
#                            function(x,y,z) {x$AWK <- y; x$NAP <- z; return(x)})
for(i in 1:n.vars)
{
  variable.dfs.AN[[i]]$AWK <- variable.dfs.AWK[[i]] %>% t() %>% as.matrix()
  variable.dfs.AN[[i]]$NAP <- variable.dfs.NAP[[i]] %>% t() %>% as.matrix()
}
# glimpse(variable.dfs.AN)

# Export cluster for parallel computations
clusterExport(cl = clust, varlist = list("ITP2bspline", "IWT2", "is.fd"), 
              envir = globalenv())

# Initialize the container
test.AN.list <- vector("list", n.vars)
test.AN.list <- pbsapply(X = variable.dfs.AN, cl = clust, simplify = F, 
  function(x) 
  {
    ITP2bspline(x$AWK,
                x$NAP,
                mu = 0,
                order = 4,
                nknots = 100,
                B = 10000,
                paired = FALSE)
    # IWT2(x$AWK,
    #      x$NAP,
    #      mu = 0,
    #      B = 1000,
    #      paired = FALSE)
  })

glimpse(test.AN.list)
plot(test.AN.list$BetaMod)

# Plot 
df.plot <- data.frame(p.val = test.AN.list$BetaMod$corrected.pval)
df.plot <- df.plot %>% mutate(is.sign = ifelse(p.val<=0.05,T,F))
df.plot <- df.plot %>% mutate(abscissa = seq(1, data.len, length.out = 102))

gg <- ggplot(data = df.plot, 
             aes(x = abscissa, y = p.val, color = is.sign)) + 
  geom_point() + 
  geom_hline(yintercept=0.05, linetype="dashed", color = "darkgrey") +
  scale_colour_manual(name='Significance',values = c('black','steelblue')) +
  scale_x_continuous(labels=NULL) +
  scale_y_continuous(breaks = c(0,0.05,1)) +
  labs(title = "Test of Beta Modulation", 
       subtitle = "AWAKE subjects vs. NAP subjects",
       x = "Abscissa", y="Adjusted p-values") +
  theme(legend.position = "None")+
  theme(plot.title=element_text(size=24, face="bold")) +
  theme(plot.subtitle=element_text(size=18, face ="italic")) +
  theme(axis.text.y = element_text(size=18)) +
  theme(axis.title.y= element_text(size=18,face ="bold")) +
  theme(plot.title=element_text(hjust=0)) +
  theme(legend.title=element_text(size=15)) +
  theme(legend.text=element_text(size=12)) +
  theme(strip.text.x=element_text(size=13,face ="bold")) +
  theme(strip.text.y=element_text(size=13,face ="bold")) +
  theme(legend.key.size = unit(1,"cm")) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  removeGrid()#ggExtra
gg

stopCluster(clust)

# NO STATISTICAL DIFFERENCE!!!

################################################################################
### PERMUTATIONAL REPEATED MEASURES ANOVA with FUNCTIONAL DATA
################################################################################
### Perform a test to check whether there is a difference between the sets 
### (due to the so-called PRACTICE EFFECT)

# CONTEXT: We have n=23 subjects. For each one we have k=14 repeated measures (sets).
# Each measure consists in a functional datum based on l=56 values.
# TEST: We want to prove that the average functional value differs in at least 
# two sets. Namely:
#           H0: mu1 = mu2 = ... = muk vs H1: (H0)^c

# LIK-INV. TRANSFORMATION: permutation of the individuals between the sets.
# Indeed, under H0 the set-wise means are equivalent!
# This transformation must guarantee that in each group we have one and only one
# record for each subject.
# TEST STAT: MANOVA T-STAT -> (- Wilks' Lambda)
#sum of the squared differences of the group mean and the overall mean
# (multiplied by the # of elements for each group)

# N.B. How to actually build the ANOVA:
# We have K buckets each one containing n_j functions. (in this case K=14, n_j=23)
# Notice that each function is given by a vector of l=56 elements. 
# Clearly, when doing the permutations we cannot separate the elements of a functions.
# Namely, we must permute the functions between the groups without changing their
# inner structure (i.e. the vector of points).
# Therefore, to efficiently handle this situation we would need a dataframe with
# k*sum(n_j) rows and l columns. Then, the permutation scheme is applied by 
# permuting the rows and splitting the df again into the buckets.

# glimpse(AbsDirErr.data)
# 
# res <- perm.func.anova(data = BetaMod.data, niter = 1000,
#                        n.by.group = rep(23, 7), l = 56*2)
# res$p.val
# hist(res$T.stats)
# abline(v = res$T0, col = 'red')

# COMMENT on the results
# Quantities which are significantly different in different sets:
# - PVel (at 10%), AbsDirErr, OnsetT

################################################################################
### PAIRWISE TESTS BETWEEN THE SETS for BetaMod

# N.B. If there is time do also a plot with comparison of the means

# Actually we are only interested in the tests for the BetaMod
# n.b. (7*13 tests per variable) x (5 variables)
# Thus if we only consider BetaMod we have 91 tests
comb.idx <- combn(n.sets,2)
B <- 2000
test.SET.list <- vector("list", ncol(comb.idx))
# Give names to the elements to better identify
names.vec <- rep(NA, ncol(comb.idx))
for(j in 1:ncol(comb.idx))
  names.vec[j] <- paste("Set", comb.idx[1,j], "vs.Set", comb.idx[2,j], sep = ".")
names(test.SET.list) <- names.vec
# Iterate over the combinations of sets
for(j in 1:ncol(comb.idx)) 
{
  invisible(capture.output(
    test <- IWT2(data1 = variable.dfs.byset$BetaMod[[comb.idx[1,j]]] %>% t(),
                 data2 = variable.dfs.byset$BetaMod[[comb.idx[2,j]]] %>% t(),
                 mu = 0, B = B, paired = T)
  ))
  test.SET.list[[j]] <- test
  print(paste("Test nr", j, " DONE", sep = ""))
}

# Check the effect of practice on BETAMOD
# We have 91 tests to analyze -> Too messy to visualize them all
# We focus on initial sets vs final sets
# 8 comparisons: 1 vs (11,12,13,14) & 2 vs (11,12,13,14)
# PLOT OF THE P-VALUE FUNCTIONS for the selected comparisons
# - First we need to create a data.frame object to exploit ggplot2
p.values <- c(test.SET.list$Set.1.vs.Set.11$adjusted_pval,
              test.SET.list$Set.1.vs.Set.12$adjusted_pval,
              test.SET.list$Set.1.vs.Set.13$adjusted_pval,
              test.SET.list$Set.1.vs.Set.14$adjusted_pval,
              test.SET.list$Set.2.vs.Set.11$adjusted_pval,
              test.SET.list$Set.2.vs.Set.12$adjusted_pval,
              test.SET.list$Set.2.vs.Set.13$adjusted_pval,
              test.SET.list$Set.2.vs.Set.14$adjusted_pval)
set.1 <- factor(rep(c("Set.1","Set.2"), c(4*set.len, 4*set.len)))
set.2 <- paste("Set", c(11,12,13,14), sep = ".") %>% rep(rep(set.len,4)) %>% 
  rep(2) %>% as.factor()
abscissa.set <- rep(1:set.len, 8)
p.value.df <- tibble(p.values = p.values, 
                     set.1 = set.1, # labels for set 1 and 2
                     set.2 = set.2, # labels for sets 11 to 14
                     abscissa = abscissa.set,
                     is.sign = p.values<0.05)
gg <- ggplot(data = p.value.df, 
             aes(x = abscissa, y = p.values, color = is.sign)) + 
  geom_point() + 
  geom_hline(yintercept=0.05, linetype="dashed", color = "darkgrey") +
  scale_colour_manual(name='Significance',values = c('black','chartreuse2')) +
  scale_x_continuous(labels=NULL) +
  scale_y_continuous(breaks = c(0.05,1)) +
  facet_grid(rows = vars(set.1), cols = vars(set.2)) +
  labs(title = "BetaMod Set-wise Comparisons", 
       subtitle = "Initial sets vs. Final sets",
       x = NULL, y="Adjusted p-values") +
  theme(legend.position = "bottom")+
  theme(plot.title=element_text(size=24, face="bold")) +
  theme(plot.subtitle=element_text(size=18, face ="italic")) +
  theme(axis.text.y=element_text(size=15)) +
  theme(axis.title.y=element_text(size=18,face ="bold")) +
  theme(plot.title=element_text(hjust=0)) +
  theme(legend.title=element_text(size=15)) +
  theme(legend.text=element_text(size=12)) +
  theme(strip.text.x=element_text(size=13,face ="bold")) +
  theme(strip.text.y=element_text(size=13,face ="bold")) +
  theme(legend.key.size = unit(1,"cm")) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  removeGrid()#ggExtra
gg


# Plot the corresponding curves
x11()
par(mfrow = c(2,4))
s.names <- rbind(rep(c("Set.1", "Set.2"), c(4,4)),
                 rep(c("Set.11", "Set.12", "Set.13", "Set.14"), 2))
for(i in 1:8)
{
  name1 <- s.names[1,i]
  name2 <- s.names[2,i]
  par(mai = c(0.35,0.25,0.65,0.15))
  matplot(variable.dfs.byset$BetaMod[[name1]], 
          ylim = c(0,3.5), lwd=2, col="black", lty=1, type = "l",
          ylab = "Beta Modulation", xlab = "Trial Nr.")
  title(main = paste(name1, "vs", name2), cex = 0.8, line = 0.8)
  matlines(variable.dfs.byset$BetaMod[[name2]], 
           lwd=2, col="chartreuse3", type="l", lty = 1)
}
# Add a general title
mtext("Curves Set-wise Comparisons", side = 3, 
      line = -2.2, outer = TRUE, cex = 2, at = 0.15)

# RESULTS: as expected, BetaMod in the last sets is greater than the one in the 
# initial sets (red curves vs black curves)

### TESTS FOR ALL THE VARIABLES!
# comb.idx <- combn(n.sets,2)
# B <- 200
# test.SET.list <- vector("list", n.vars)
# names(test.SET.list) <- names(variable.dfs)
# test.SET.list <- lapply(test.SET.list,
#                         function(x) {vector("list", ncol(comb.idx))})
# # p.values <- vector('list', n.vars)
# # p.values <- lapply(p.values, function(x) {matrix(NA, nrow = ncol(comb.idx), ncol = 102)})
# for(i in 1:n.vars)
# {
#   curr.data <- variable.dfs.byset[[i]]
#   for(j in 1:ncol(comb.idx)) # iterate over the combinations of sets
#   {
#     # invisible(capture.output(
#     #   test <- ITP2bspline(curr.data[[comb.idx[1,j]]],
#     #                       curr.data[[comb.idx[2,j]]],
#     #                       mu = 0,
#     #                       order = 4,
#     #                       nknots = set.len,
#     #                       B = B,
#     #                       paired = T)))
#     invisible(capture.output(
#       test <- IWT2(data1 = curr.data[[comb.idx[1,j]]] %>% t(),
#                    data2 = curr.data[[comb.idx[2,j]]] %>% t(),
#                    mu = 0, B = B, paired = T)
#     ))
#     test.SET.list[[i]][[j]] <- test
#     # p.values[[i]][j,] <- test$corrected.pval
#     print(paste("Test nr", j, " DONE", sep = ""))
#   }
#   print(paste("Variable nr.", i, " is completed!!!", sep = ""))
# }
# # Give names to the elements to better identify
# for(i in 1:n.vars)
# {
#   names.vec <- rep(NA, ncol(comb.idx))
#   for(j in 1:ncol(comb.idx))
#     names.vec[j] <- paste("Set", comb.idx[1,j],
#                           "vs.Set", comb.idx[2,j], sep = ".")
#   names(test.SET.list[[i]]) <-  names.vec
# }
# names(test.SET.list$AbsDirErr$Set.1.vs.Set.2)
# plot.IWT2(test.SET.list[[1]][[1]])
# glimpse(test.SET.list)

# # Check which p-values are under a certain threshold
# # I am interested in the row (since each row corresponds to a certain comparison)
# signif.comp <- which(p.values[[6]] < 0.05, arr.ind = T) %>% as.data.frame()
# signif.comp[order(signif.comp$row),]
# p.values[[1]][19,3]

################################################################################
# TEST OF THE GLOBAL EFFECT OF PRACTICE (on BetaMod)
# Let c0_i, i = 1, ..., 6 be the average value (across subjects and tasks) in 
# set#1 of each one of the chosen variables. We want to do an Interval-Wise 
# test considering the entire functional data to see whether the average value  
# of the curves is equal (or not) to c0_i.
#           H0: mu_i = c0_i vs H1: mu_i != c0_i, for i = 1,...,6
c0 <- rbind(variable.dfs.byset$BetaMod$Set.1,
        variable.dfs.byset$BetaMod$Set.2) %>% rowMeans() %>% mean()

# Perform the test using a Interval-Wise approach
test.c0 <- IWT1(data = t(variable.dfs$BetaMod), 
                mu = c0, B = 1000)

# Plot the result with ggplot
c0.df <- tibble(values = as.matrix(variable.dfs$BetaMod) %>% as.numeric(),
                abscissas = rep(1:data.len, N),
                pval = rep(test.c0$adjusted_pval, N),
                subj = factor(rep(1:N, rep(data.len, N))))
c0.df <- c0.df %>% mutate(is.signif = pval<0.05)
c0.df

# Create a Data.frame to add text in the plot
text.df <- data.frame(set = paste("Set", 1:n.sets, sep = " "))
x.pos <- lapply(set.delim[-length(set.delim)], function(x)
{
  knots <- seq(1, set.len, length.out = 3)
  # Drop external knots
  knots <- knots[-c(1,3)]
  # Add the set position
  knots <- knots + x
  return(knots)
})
x.pos <- unlist(x.pos)
text.df <- text.df %>% mutate(x.pos = x.pos, y.pos = rep(10,n.sets))

ggg <- ggplot(c0.df, aes(x = abscissas, y = values,
                         group = subj, col = is.signif)) + 
  geom_vline(xintercept = set.delim, col = "darkgrey", size = 0.8) +
  geom_hline(yintercept = c0, col = "red", size = 1.8) +
  geom_line() +
  scale_colour_manual(name='Significance',values = c('black','chartreuse3')) +
  labs(title = "Test on Beta Modulation",
       subtitle = "Beta Modulation vs Mean in 1st and 2nd sets",
       x = "Trial Nr.", y = "Beta Modulation") +
  theme(legend.position = "bottom")+
  theme(plot.title=element_text(size=24, face="bold")) +
  theme(plot.subtitle=element_text(size=18, face ="italic")) +
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.text.x=element_text(size=12)) + 
  theme(axis.title.y=element_text(size=16,face ="bold")) +
  theme(axis.title.x=element_text(size=16,face ="bold")) +
  theme(plot.title=element_text(hjust=0)) +
  theme(legend.title=element_text(size=15)) +
  theme(legend.text=element_text(size=12)) + 
  theme(legend.key.size = unit(1,"cm")) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  removeGrid()
ggg
  
# c0.vect <- lapply(group.means.pw.byset, function(x)
#   {as.numeric(x[1,1])})
# 
# # # Export the cluster for parallel computations
# # clusterExport(cl = clust, varlist = list("c0.vect", "variable.dfs"))
# 
# c0.tests <- vector("list", n.vars)
# c0.tests <- pbmapply(x = variable.dfs, y = c0.vect, function(x,y)
#   {
#     ITP1bspline(data = t(x),
#                 mu = y,
#                 order = 4,
#                 nknots = 100,
#                 B = 1000)
#   }, SIMPLIFY = F, USE.NAMES = T) 
# # N.B. Only "pbsapply" works in parallel! (because it can work column by column)
# 
# c0.tests[[1]]$corrected.pval
# plot(c0.tests[[2]])
# 
# stopCluster(clust)
# 
# # Repeat a similar test considering, instead of the constant function c0_i, the
# # actual sample mean value (across the subjects) of the function in set#1, 
# # repeated over all the sets
# group.mean.functions <- lapply(group.stats, function(x) 
#   {
#     rep(x$Set.1[1,], n.sets) %>% as.numeric()
#   }) 
# 
# plot(group.mean.functions[[1]])
# 
# c1.tests <- vector("list", n.vars)
# c1.tests <- pbmapply(x = variable.dfs, y = group.mean.functions, function(x,y)
# {
#   ITP1bspline(data = t(x),
#               mu = y,
#               order = 4,
#               nknots = 100,
#               B = 1000)
# }, SIMPLIFY = F, USE.NAMES = T) 
# 
# plot(c1.tests$BetaMod)
# c1.tests$BetaMod$corrected.pval
# 
# # To do the test on the values and not on the coefficients of the b-spline
# IWT1(data = t(variable.dfs[[1]]), mu = group.mean.functions[[1]])






