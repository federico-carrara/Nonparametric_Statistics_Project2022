library(roahd)
library(ggplot2)
library(ggrepel)

# In order to exploit the potential of "roahd" package we need to transform the 
# data into "fData" objects.
# To do this we just need data to be stored in a matrix format, with records on
# on the rows (e.g. 1st row 1st subject, ...). Then we call the function "fData",
# specifying the abscissa values.

# We need data organize in the matrix form specified above:
# - Transpose each dataframe in order to have subjects on the rows
# - Cast the data type into matrix type
variable.dfs.fData <- lapply(variable.dfs, function(x) {as.matrix(t(x))})
# Transform into fData object
abscissa <- 1:data.len
variable.dfs.fData <- lapply(variable.dfs.fData, function(x)
  {fData(grid = abscissa, values = x)})

# Consider only the BetaMod
BetaMod.fData <- variable.dfs.fData$BetaMod

# Compute the BD for BetaMod
BetaMod.MBD <- MBD(BetaMod.fData)

# Plot the results
subject.IDs <- names(BetaMod.MBD)
group.IDs <- sapply(subject.IDs, function(x) { grepl(pattern = "A", x = x)})
MBD.df <- data.frame(ID = subject.IDs, 
                    group = as.factor(group.IDs),
                    value = BetaMod.MBD)
glimpse(MBD.df)


g <- ggplot(MBD.df, aes(x = ID, y = value)) +
  theme(legend.position="none") +
  geom_label(data = MBD.df,
             aes(fill = group, label = ID),
             label.padding = unit(0.85, "lines"),
             label.size = 0.75,
             color = "black") +
  labs(title = "Modified Band Depth for Beta Modulation", 
       subtitle = "Awake vs. Nap",
       x = NULL, y="MBD") + 
  scale_fill_discrete(name = "Group", labels = c("Awake", "Nap")) +
  theme(legend.position = "bottom")+
  theme(plot.title=element_text(size=24, face="bold")) +
  theme(plot.subtitle=element_text(size=18, face ="italic")) +
  theme(axis.text.y=element_text(size=8)) +
  theme(axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=18,face ="bold")) +
  theme(plot.title=element_text(hjust=0)) +
  theme(legend.title=element_text(size=18)) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.key.size = unit(1,"cm")) +
  removeGrid()#ggExtra
plot(g)

##################################################################
# Now treat the data as multivariate fData and repeat the analysis
variable.dfs.multifData <- as.mfData(variable.dfs.fData) 
names(variable.dfs.multifData$fDList) <- names(variable.dfs.fData)
plot(variable.dfs.multifData)

# Compute the multivariate MBD and the median subject
multi.MBD <- multiMBD(variable.dfs.multifData)
subject.IDs <- rownames(variable.dfs.fData$AbsArea$values)
group.IDs <- sapply(subject.IDs, function(x) { grepl(pattern = "A", x = x)})
multi.MBD <- data.frame(ID = as.factor(subject.IDs), 
                        group = as.factor(group.IDs),
                        value = multi.MBD)
glimpse(multi.MBD)

subject.IDs[which.max(multi.MBD$value)]
subject.IDs[which.min(multi.MBD$value)]

# Plot the results
g <- ggplot(multi.MBD, aes(x = ID, y = value)) +
            geom_point() +
  theme(legend.position="none") +
  geom_label(data = multi.MBD,
             aes(fill = group),
             label = multi.MBD$ID,
             label.padding = unit(0.55, "lines"),
             label.size = 0.35,
             color = "black")
plot(g)
# N.B. This plot can be used before the test NAP vs AWAKE

###############################
# Compute the Spearman's Matrix
spearman.matrix <- cor_spearman(variable.dfs.multifData)
spearman.df <- expand.grid(X = row.names(spearman.matrix), 
                           Y = row.names(spearman.matrix)) 
spearman.df$value <- as.numeric(spearman.matrix) 
g2 <- ggplot(spearman.df, aes(x = X, y = Y, fill = value)) + 
  geom_tile() + 
  # scale_fill_gradient(low = "white", high = "blue") +
  scale_fill_distiller(name = "Legend") +
  # scale_fill_viridis(option = "C", discrete = F) + 
  labs(title = "Spearman Matrix", 
       subtitle = "Spearman coefficients for Functional Data") +
  theme(legend.position = "right")+
  theme(plot.title=element_text(size=24, face="bold")) +
  theme(plot.subtitle=element_text(size=18, face ="italic")) +
  theme(axis.text.y=element_text(size=15)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x=element_blank()) +
  theme(plot.title=element_text(hjust=0)) +
  theme(legend.title=element_text(size=20)) +
  theme(legend.text=element_text(size=15)) +
  removeGrid()#ggExtra
  
plot(g2)

######################################
# Functional boxplots and outliergrams
# N.B. Here we need to treat the data as univariate fData
for(i in 1:length(variable.dfs.fData))
{
  set.seed(1234)
  x11()
  # fbplot(variable.dfs.fData[[i]], 
  #        main = paste("Functional BoxPlot -", names(variable.dfs.fData)[i]), 
  #        adjust = list(VERBOSE=F)) %>% invisible()
  outliergram(variable.dfs.fData[[i]],
              main = paste("Outliergram -", names(variable.dfs.fData)[i]),
              adjust = list(VERBOSE=F))
  print(paste("Here", i))
}

outliergram(variable.dfs.fData$BetaMod,
            main = paste("Outliergram -", names(variable.dfs.fData)[5]),
            adjust = list(VERBOSE=F))

roahd::fbplot(variable.dfs.fData$AbsArea, 
       main = paste("Functional BoxPlot -", names(variable.dfs.fData)[1]), 
       adjust = list(VERBOSE=F))
