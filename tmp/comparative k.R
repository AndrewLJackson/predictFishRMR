
infile <- "kplotting.csv"   #save a .csv file under this name, with headers as described in the email (obviously save another one for "Wash_range_test_data.cav"

library(ggplot2)
library(tidyverse)
library(broom)
library(nlme)


Data_table = read.csv(infile , header=TRUE) 
Data_table

#run a mixed model
lme <- lme(K ~ M, random =~ 1 | species, data = Data_table, method="ML")
summary(lme)


p <- ggplot(Data_table, (aes(x=M, y=K, color=species))) + 
  theme_classic(base_size = 14) + 
  ylab("k for cooling (°C/min/°C)") + 
  xlab("Mass (kg)") + 
  geom_point(size=4) +
  geom_abline(slope = -0.6276487, intercept = -0.9730899, colour = "red", size = 0.5)+
  

#+ 
#  geom_line(data = augment(ancova), 
#                           mapping = aes(y = .fitted), 
#            linewidth = 1) + 
  #geom_smooth(method=lm, se=FALSE, fullrange=FALSE) + 
  scale_fill_manual(values=c("blue", "cyan4")) + 
  # scale_colour_manual(values=c("white", "black"))
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9)
  scale_shape(solid = FALSE)
  scale_colour_brewer(palette = "RdYlBu")
  

print(p)

#to remove legend
noleg <-p + theme(legend.position = "none")

#to add image to plot
library("jpeg") 
library("patchwork") 

# defining the x coordinates 
xpos <- 1:5 

# defining the y coordinates 
ypos <- xpos**2 

data_frame = data.frame(xpos = xpos, 
                        ypos = ypos) 

print ("Data points") 
print (data_frame) 

# plotting the data 
graph <- ggplot(data_frame, aes(xpos, ypos)) +     
  geom_point() 





            