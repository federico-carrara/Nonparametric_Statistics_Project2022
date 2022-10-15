################################################################################
source("Auxiliary_Functions_New.R")
################################################################################
### LOAD DATA
# Set WD
setwd("C:/Users/fede1/OneDrive - Politecnico di Milano/PROJECT - NPS/DATA")

# N.B. Data are stored in an Excel file. Each sheet corresponds to a different
# subject.
# Read data for all the subjects and store them in a list
n.vars <- 6
n.sets <- 14
set.len <- 56 # elements for each set
set.delim <- c(0, set.len*(1:n.sets))
N <- 23 # subjects
data.len <- n.sets * set.len

# Initialize a list. Its elements will be the datasets for each variable
variable.dfs <- vector("list", n.vars) 
for(i in 1:n.vars)
{
  variable.dfs[[i]] <- read_xlsx(path = "Denoise_Data_Block1.xlsx", 
                                 sheet = i)
  print(paste("Data of Variable [", i, "] loaded"))
}

names(variable.dfs) <- excel_sheets(path = "Denoise_Data_Block1.xlsx")

# Drop ERS column (since it is the same as BetaMod)
variable.dfs <- variable.dfs[-5]
n.vars <- 5

################################################################################
### SPLIT THE DATA SET BY SET
# This is useful since we want to do a set-wise smoothing in the following
# Data divided by set.
# N.B. To do this we already implemented ad-hoc the split.sets function.
# N.B. The ouutput of this function is a nested-list structured object.
# The elements of the outer list are the variable. Then, for each variable
# we have an inner list, whose elements are datasets, where on the columns
# there are the N=23 subjects and on the rows there are the n.sets=56 tasks.
variable.dfs.byset <- lapply(variable.dfs, split.sets, set.len)
# glimpse(variable.dfs.byset$AbsArea)


################################################################################
### PLOT OF THE COMPLETE DATASET
# Add for each variable columns for the plot
variable.dfs.plot <- lapply(variable.dfs, function(x)
{
  tibble(values = as.matrix(x) %>% as.numeric(),
         abscissas = rep(1:data.len, N),
         subj = factor(rep(1:N, rep(data.len, N))))
})
names(variable.dfs.plot) <- names(variable.dfs)
glimpse(variable.dfs.plot)
# variable.dfs.plot <- lapply(variable.dfs, function(x)
# {
#   x <- x %>% mutate(abscissas = 1:data.len)
# })
# names(variable.dfs.plot) <- names(variable.dfs)

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

# ggg <- ggplot() + 
#   geom_vline(xintercept = set.delim, col = "darkgrey", size = 0.8)
# for(i in 1:N)
# {
#   curr.subj <-  variable.dfs.plot$AbsArea[[i]]
#   ggg <- ggg + 
#     geom_line(data = variable.dfs.plot$AbsArea, aes(x = abscissas, y = curr.subj))
# }
#   geom_line() +
#   # geom_label(data = text.df,
#   #            aes(x = x.pos, y = y.pos),
#   #            label.padding = unit(0.55, "lines"),
#   #            label.size = 0.35,
#   #            color = "black") + 
#   labs(title = "Raw Data Visualization - Absolute Norm. Area",
#        x = "Trial Nr.", y = "Absolute Area") +
#   scale_color_viridis(option ="C") +
#   theme(legend.position = "None")+
#   theme(plot.title=element_text(size=20, face="bold")) +
#   theme(axis.text.y=element_text(size=12)) +
#   theme(axis.text.x=element_text(size=12)) + 
#   theme(axis.title.y=element_text(size=16,face ="bold")) +
#   theme(axis.title.x=element_text(size=16,face ="bold")) +
#   theme(plot.title=element_text(hjust=0)) +
#   theme(legend.title=element_text(size=15)) +
#   theme(legend.text=element_text(size=12)) + 
#   theme(legend.key.size = unit(1,"cm")) +
#   guides(colour = guide_legend(override.aes = list(size=2))) +
#   removeGrid()
# 
# ggg

ggg <- ggplot(variable.dfs.plot$PVel, aes(x = abscissas, y = values,
                                             group = subj, color = subj)) + 
  geom_vline(xintercept = set.delim, col = "darkgrey", size = 0.8) +
  geom_line() +
  # geom_label(data = text.df,
  #            aes(x = x.pos, y = y.pos),
  #            label.padding = unit(0.55, "lines"),
  #            label.size = 0.35,
  #            color = "black") + 
  labs(title = "Peak Velocity",
       subtitle = "Unit of Measure [dm/s]",
       x = "Trial Nr.", y = "Peak Velocity") +
  scale_color_viridis(option = "H", discrete = T) + 
  theme(legend.position = "None") +
  theme(plot.subtitle=element_text(size=18, face ="italic")) +
  theme(plot.title=element_text(size=20, face="bold")) +
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
