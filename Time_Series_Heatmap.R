library(viridis) # colour blind friendly palette, works in B&W also
library(ggplot2)
library(ggExtra)
# library(Interpol.T) #  will generate a large dataset on initial load
# library(lubridate) # for easy date manipulation
# library(ggExtra) # because remembering ggplot theme options is beyond me
# library(ggplot2)
# library(tidyr)

# Prepare the dataframe
# I need the following variables:
# - y-axis: SUBJECT 
# - x-axis: TASKNO
# - squares grouping: SET
# - numeric variables: 1 value per square ( 1 plot per each variable)

# Put all the numeric data on a column
ts.map.df <- lapply(variable.dfs, function(x) {x %>% as.matrix() %>% as.numeric()})
ts.map.df <- bind_cols(ts.map.df)
# Add the required features for plotting
subjectID <- names(variable.dfs$AbsArea)
ts.map.df <- ts.map.df %>% mutate(set = rep(rep(1:n.sets, rep(56, n.sets)), N),
                                  subject = rep(subjectID, rep(784, N)),
                                  taskNo = rep(1:784, N))

# Let's do the plot for AbsDirErr
p <- ggplot(ts.map.df, aes(taskNo, subject, fill = BetaMod)) +
  geom_tile() +
  # geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis(option = "D") +
  geom_vline(xintercept = set.delim) + 
  # facet_grid(cols = vars(set)) +
  # scale_x_continuous(breaks = c(1,15,30,45,56)) +
  theme_minimal(base_size = 8) +
  labs(title = "Beta Modulation", 
       subtitle = "Unit of Measure [Hz]",
       x = "Task Nr.", y = NULL) + 
  theme(legend.position = "right") +
  theme(plot.title=element_text(size=24, face="bold")) +
  theme(plot.subtitle=element_text(size=18, face ="italic")) +
  theme(plot.title=element_text(hjust=0)) +
  theme(axis.text.y=element_text(size=12, margin=margin(r=0))) +
  theme(axis.text.x=element_text(size=12)) + 
  theme(axis.title.y=element_text(size=16,face ="bold")) +
  theme(axis.title.x=element_text(size=16,face ="bold")) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_blank()) +
  theme(strip.text=element_text(size=12,face ="bold")) +
  removeGrid()#ggExtra
p



################################################################################
### AN EXAMPLE
# data <- data(Trentino_hourly_T,package = "Interpol.T")
# 
# names(h_d_t)[1:5]<- c("stationid","date","hour","temp","flag")
# df <- tbl_df(h_d_t) %>%
#   filter(stationid =="T0001")
# 
# df <- df %>% mutate(year = year(date),
#                     month = month(date, label=TRUE),
#                     day = day(date))
# 
# df$date<-ymd(df$date) # not necessary for plot but 
# #useful if you want to do further work with the data
# 
# #cleanup
# rm(list=c("h_d_t","mo_bias","Tn","Tx",
#           "Th_int_list","calibration_l",
#           "calibration_shape","Tm_list"))
# 
# 
# #create plotting df
# df <- df %>% dplyr::select(stationid,day,hour,month,year,temp)
# df$month %>% length()
# statno <-unique(df$stationid)
# 
# p <- ggplot(df, aes(day,hour,fill=temp))+
#   geom_tile(color= "white",size=0.1) + 
#   scale_fill_viridis(name="Hrly Temps C",option ="C")
# p <-p + facet_grid(year~month)
# p <-p + scale_y_continuous(trans = "reverse", breaks = unique(df$hour))
# p <-p + scale_x_continuous(breaks =c(1,10,20,31))
# p <-p + theme_minimal(base_size = 8)
# p <-p + labs(title= paste("Hourly Temps - Station",statno), x="Day", y="Hour Commencing")
# p <-p + theme(legend.position = "bottom")+
#   theme(plot.title=element_text(size = 14))+
#   theme(axis.text.y=element_text(size=6)) +
#   theme(strip.background = element_rect(colour="white"))+
#   theme(plot.title=element_text(hjust=0))+
#   theme(axis.ticks=element_blank())+
#   theme(axis.text=element_text(size=7))+
#   theme(legend.title=element_text(size=8))+
#   theme(legend.text=element_text(size=6))+
#   removeGrid()#ggExtra

