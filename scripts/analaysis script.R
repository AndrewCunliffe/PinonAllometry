## Script for analysing  pnion allometry data
## Andrew Cunliffe <andrewmcunliffe@gmail.com> & Cameron McIntire
## Started 2020-10-26

# Establish operating environment ----
rm(list=ls())
#setwd("C:\Users\cmcintire\Desktop")
#path <- getwd()

# Load required packages ----
#   # Download
if(!require(tidyverse)) {install.packages("tidyverse"); require(tidyverse)}
if(!require(viridis)) {install.packages("viridis"); require(viridis)}
if(!require(patchwork)) {install.packages("patchwork"); require(patchwork)}
if(!require(propagate)) {install.packages("propagate"); require(propagate)}
if(!require(nlstools)) {install.packages("nlstools"); require(nlstools)}
if(!require(ggpmisc)) {install.packages("ggpmisc"); require(ggpmisc)}
if(!require(polynom)) {install.packages("polynom"); require(polynom)}
if(!require(gvlma)) {install.packages("gvlma"); require(gvlma)}
if(!require(ggpubr)) {install.packages("ggpubr"); require(ggpubr)}
if(!require(propagate)) {install.packages("propagate"); require(propagate)}
if(!require(ggpmisc)) {install.packages("ggpmisc"); require(ggpmisc)}



# Install
library(tidyverse) # For dplyr and ggplot2
library(viridis) # load friendly colour palette for plotting.
library(patchwork) # this is Andy's favorite for multi-panel plots
library(propagate) # required for predicting confidence intervals on ESD : Biomass plot
library(ggpmisc) # for adding model parameters to plots
library(gvlma) # Global Validation of Linear Models Assumptions
library(polynom)
library(nlstools)
library(ggpubr) 


# Load data ----
pinon_data <- read.csv("data/pinon_data.csv", header = T)[-1,]  # Read in summary  data


# Create plotting theme and colour scheme----
theme_coding <- function(){
  theme_classic()+
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
          axis.title = element_text(size = 10),
          panel.grid = element_blank(),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size = 10, vjust = 1, hjust = 0.5),
          legend.title = element_blank(),
          legend.text = element_text(size = 6),
          legend.key.size = unit(1,"line"),
          legend.position = c(0.9, 0.9))
}

# Colour scheme
Colour.wet <- "#440154FF"
Colour.dry <- "#39568CFF"
Colour.image <- "#20A387FF"
Colours <- c(Colour.wet, Colour.dry, Colour.image)

## Prepare data ----
# Specify numeric types
pinon_data$longitude = as.numeric(pinon_data$longitude)
pinon_data$latitude = as.numeric(pinon_data$latitude)
pinon_data$canopy_area = as.numeric(pinon_data$canopy_area)
pinon_data$canopy_area = as.numeric(pinon_data$canopy_area)
pinon_data$CanDia1 = as.numeric(pinon_data$CanDia1)
pinon_data$CanDia2 = as.numeric(pinon_data$CanDia2)
pinon_data$dieback_pc = as.numeric(pinon_data$dieback_pc)
pinon_data$max_height = as.numeric(pinon_data$max_height)
pinon_data$base_live_canopy = as.numeric(pinon_data$base_live_canopy)
pinon_data$crown_depth = as.numeric(pinon_data$crown_depth)
pinon_data$live_crown_ratio = as.numeric(pinon_data$live_crown_ratio)
pinon_data$diameter_at_base_wet = as.numeric(pinon_data$diameter_at_base)
pinon_data$diameter_at_30_cm = as.numeric(pinon_data$diameter_at_30_cm)
pinon_data$disk_diameter_wet = as.numeric(pinon_data$disk_diameter_wet)
pinon_data$disk_diameter_dry = as.numeric(pinon_data$disk_diameter_dry)
pinon_data$base_ba_wet = as.numeric(pinon_data$base_ba_wet)
pinon_data$disk_ba_wet = as.numeric(pinon_data$disk_ba_wet)
pinon_data$disk_diameter_dry = as.numeric(pinon_data$disk_diameter_dry)
pinon_data$inside_bark_area = as.numeric(pinon_data$inside_bark_area)
pinon_data$inside_bark_diam = as.numeric(pinon_data$inside_bark_diam)
pinon_data$heartwood_area = as.numeric(pinon_data$heartwood_area)
pinon_data$sapwood_area = as.numeric(pinon_data$sapwood_area)
pinon_data$bark_thickness = as.numeric(pinon_data$bark_thickness)
pinon_data$wet_mass_small_partititon = as.numeric(pinon_data$wet_mass_small_partititon)
pinon_data$wet_mass_large_partititon = as.numeric(pinon_data$wet_mass_large_partititon)
pinon_data$wet_mass_total = as.numeric(pinon_data$wet_mass_total)
pinon_data$moisture_content_of_small_partition = as.numeric(pinon_data$moisture_content_of_small_partition)
pinon_data$moisture_content_of_large_partition = as.numeric(pinon_data$moisture_content_of_large_partition)
pinon_data$moisture_content_of_total = as.numeric(pinon_data$moisture_content_of_total)
pinon_data$dry_mass_small_partititon = as.numeric(pinon_data$dry_mass_small_partititon)
pinon_data$dry_mass_large_partititon = as.numeric(pinon_data$dry_mass_large_partititon)
pinon_data$dry_mass_total_g = as.numeric(pinon_data$dry_mass_total_g)
pinon_data$dry_mass_total_kg = as.numeric(pinon_data$dry_mass_total_kg)
pinon_data$proportion_subsampled_and_dried_small_partition = as.numeric(pinon_data$proportion_subsampled_and_dried_small_partition)
pinon_data$proportion_subsampled_and_dried_large_partition = as.numeric(pinon_data$proportion_subsampled_and_dried_large_partition)
pinon_data$proportion_subsampled_and_dried_total = as.numeric(pinon_data$proportion_subsampled_and_dried_total)
pinon_data$drone_canopy_height_min = as.numeric(pinon_data$drone_canopy_height_min)
pinon_data$drone_canopy_height_max = as.numeric(pinon_data$drone_canopy_height_max)
pinon_data$HAG_plotmean_of_cellmax_m = as.numeric(pinon_data$HAG_plotmean_of_cellmax_m)
pinon_data$HAG_plotmedian_of_cellmax_m = as.numeric(pinon_data$HAG_plotmedian_of_cellmax_m)
pinon_data$HAG_plot90percentile_of_cellmax_m = as.numeric(pinon_data$HAG_plot90percentile_of_cellmax_m)


# Calculate canopy area from a and b diameters
pinon_data$canopy_area_from_daim <- pi * (pinon_data$CanDia1)/2 * (pinon_data$CanDia2)/2



# Linear models

# Model for comparing canopy area from polygons (CA1) and field obs (CA2)
Model_CA1_CA2 <- lm(canopy_area_from_daim ~ canopy_area, data = pinon_data)
summary(Model_CA1_CA2)
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.01343    0.16860   -0.08    0.937    
# canopy_area  1.07822    0.02318   46.51   <2e-16 ***
   
# Plot drone derived canopy area as a function of field estimated canopy area
 ggplot(data = pinon_data, aes( x = canopy_area, y = canopy_area_from_daim )) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", formula = "y~x", color = "black") +
    geom_abline(intercept = 0, slope = 1, color="black", 
                linetype="dashed", size= 0.5) +
    theme_classic() +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
          axis.title = element_text(size = 10),
          panel.grid = element_blank(),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          panel.border = element_rect(colour = "black", fill=NA, size = 0.8)) +
    scale_y_continuous(limits=c(-0.5,21), breaks=seq(0,22,5), expand = c(0,0)) +
    scale_x_continuous(limits=c(-0.5,21), breaks=seq(0,22,5), expand = c(0,0)) +
    labs(x = expression(CA["1"]~(m^2)),
         y = expression(CA["2"]~(m^2))) +
    stat_poly_eq(aes(label = paste("atop(", stat(adj.rr.label), ",", stat(eq.label), ")", sep = "")), 
                 formula = "y~x", 
                 parse = TRUE) #insert R2 and equation into the plot
 
 #save the figure to WD
 ggsave("CA1_CA2_regression.tiff", width = 10, height = 10, units = "cm", dpi = 500)

#------
# Model for total tree biomass as a function of the field measured BA
Model_BAwet <- lm(dry_mass_total_kg ~ base_ba_wet, data = pinon_data)
summary(Model_BAwet)
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -7.40437    6.59978  -1.122    0.277    
#base_ba_wet  0.33130    0.03106  10.666 3.28e-09 ***

#------
# Model for total tree biomass as a function of the fresh disk measured BA
Model_disk_BAwet <- lm(dry_mass_total_kg ~ disk_diameter_wet, data = pinon_data)
summary(Model_disk_BAwet)
#Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       -30.6859    11.1847  -2.744   0.0134 *  
#disk_diameter_wet   6.4158     0.8209   7.815 3.41e-07 ***

#------
# Model for total tree biomass as a function of the dry disk measured BA
Model_disk_BAdry <- lm(dry_mass_total_kg ~ disk_diameter_dry, data = pinon_data)
summary(Model_disk_BAdry)
#Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       -29.6225    11.2213  -2.640   0.0166 *  
#disk_diameter_dry   6.7203     0.8725   7.702 4.19e-07 ***

#------
# Model for comparing field measured height to drone derived max height
Model_height <- lm(drone_canopy_height_max ~ max_height, data = pinon_data) 
summary(Model_height)
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -0.11238    0.17395  -0.646    0.526    
#max_height   0.98708    0.04484  22.013 1.83e-14 ***

# Testing difference between field measured height and drone measured height
height_difference <- pinon_data$drone_canopy_height_max - pinon_data$max_height
# Assessing data distribution
length(height_difference)
qqnorm(height_difference)
qqline(height_difference)       
shapiro.test(height_difference)
# W = 0.83306, p-value = 0.00281
# these data are signifcantly non-normal, so a non-parametric test is needed.
# Wilcoxon signed-rank test
wilcox.test(pinon_data$drone_canopy_height_max, pinon_data$max_height,
            paired = TRUE, alternative = "less",
            conf.int = TRUE)
# Result: height derived from drone and field measured are not sig different

# Plot drone derive height as a function of measured height
ggplot(data = pinon_data, aes( x = max_height, y = drone_canopy_height_max )) +
   geom_point(size = 2) +
   geom_smooth(method = "lm", formula = "y~x", color = "black") +
   geom_abline(intercept = 0, slope = 1, color="black", 
               linetype="dashed", size= 0.5) +
   theme_classic() +
   theme(axis.text = element_text(size = 10, color = "black"),
         axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
         axis.title = element_text(size = 10),
         panel.grid = element_blank(),
         plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
         panel.border = element_rect(colour = "black", fill=NA, size = 0.8)) +
   scale_y_continuous(limits=c(-0.5,7), breaks=seq(0,6,2), expand = c(0,0)) +
   scale_x_continuous(limits=c(-0.5,7), breaks=seq(0,6,2), expand = c(0,0)) +
   xlab("Measued height (m)") +
   ylab("UAV estimated height (m)") +
   stat_poly_eq(aes(label = paste("atop(", stat(adj.rr.label), ",", stat(eq.label), ")", sep = "")), 
                formula = "y~x", 
                parse = TRUE) #insert R2 and equation into the plot

#save the figure to WD
ggsave("maxheight_droneheight.tiff", width = 10, height = 10, units = "cm", dpi = 500)


#------
# Model for comparing bark-thickness to RCD
Model_barkthickness <- lm(bark_thickness ~ disk_diameter_wet, data = pinon_data)
summary(Model_barkthickness)
#Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       0.129403   0.063568   2.036   0.0568 .  
#disk_diameter_wet 0.055105   0.004666  11.811 6.52e-10 ***

# Plot drone derive height as a function of measured height
ggplot(data = pinon_data, aes( x = disk_diameter_wet, y = bark_thickness )) +
   geom_point(size = 2) +
   geom_smooth(method = "lm", formula = "y~x", color = "black") +
   theme_classic() +
   theme(axis.text = element_text(size = 10, color = "black"),
         axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
         axis.title = element_text(size = 10),
         panel.grid = element_blank(),
         plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
         panel.border = element_rect(colour = "black", fill=NA, size = 0.8)) +
   scale_y_continuous(limits=c(0,2.2), breaks=seq(0,3,0.5), expand = c(0.01,0)) +
   scale_x_continuous(limits=c(0,24), breaks=seq(0,24,4), expand = c(0.01,0)) +
   xlab("Disk diameter (cm)") +
   ylab("Bark thickness (cm)") +
   stat_poly_eq(aes(label = paste("atop(", stat(adj.rr.label), ",", stat(eq.label), ")", sep = "")), 
                formula = "y~x", 
                parse = TRUE) #insert R2 and equation into the plot

#save the figure to WD
ggsave("barkthickness.tiff", width = 10, height = 10, units = "cm", dpi = 500)

#------
# Model for comparing wet and dry disk diameters
Model_disk_wetdry <- lm(disk_diameter_dry ~ disk_diameter_wet, data = pinon_data)
summary(Model_disk_wetdry)
#Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       -0.112872   0.116967  -0.965    0.347    
#disk_diameter_wet  0.950785   0.008585 110.750   <2e-16 ***

# Reformatting RCD data for plotting
# separate the diameter variables
RCD_wet <- data.frame(pinon_data$diameter_at_base_wet)
Disk_wet <- data.frame(pinon_data$disk_diameter_wet)
Disk_dry <- data.frame(pinon_data$disk_diameter_dry)

#combine RCD variable with disk variables
Disk_wet <- cbind(RCD_wet, Disk_wet)
Disk_dry <- cbind(RCD_wet, Disk_dry)

#assign the wet and dry categorigal variables
Disk_wet$Group <- "Wet"
Disk_dry$Group <- "Dry"

#rename columns
colnames(Disk_wet) <- c("RCD", "Disk_diam", "Group")
colnames(Disk_dry) <- c("RCD", "Disk_diam", "Group")

#coimbine into a single dataframe
RCD_compare <- rbind(Disk_wet, Disk_dry)

ggplot(data = RCD_compare, aes( x = RCD, y = Disk_diam, fill = Group )) +
   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
   geom_point(shape = 21, size = 2) +
   geom_smooth(aes(group = Group), method = "lm", formula = "y~x", color = "black", size = 0.6, fill = "grey20") +
   scale_fill_manual(values = c("black", "white")) +
   theme_classic() +
   theme(axis.text = element_text(size = 10, color = "black"),
         axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
         axis.title = element_text(size = 10),
         panel.grid = element_blank(),
         plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
         panel.border = element_rect(colour = "black", fill=NA, size = 0.8),
         legend.title = element_blank(),
         legend.text = element_text(size = 10),
         legend.key.size = unit(1,"line"),
         legend.position = c(0.4, 0.89),
         legend.background=element_blank()) +
   guides(fill=guide_legend(
      keywidth=0.1,
      keyheight=.6,
      default.unit="cm")
   ) + #this bit allows for adjusting the distance between legend symbols so I can allign with the R2 values
   scale_y_continuous(limits=c(0,24), breaks=seq(0,24,4), expand = c(0.01,0)) +
   scale_x_continuous(limits=c(0,24), breaks=seq(0,24,4), expand = c(0.01,0)) +
   xlab("RCD (cm)") +
   ylab("Disk diameter (cm)") +
   stat_poly_eq(formula = "y~x",
                label.x = "centre",
                eq.with.lhs = "italic(hat(y))~`=`~",
                aes(label = paste(stat(adj.rr.label),sep = "")), rr.digits = 3,
                label.x.npc = "left",
                parse = TRUE)


#save the figure to WD
ggsave("RCD_disk_compare.tiff", width = 10, height = 10, units = "cm", dpi = 500)

#------
# Model for comparing biomass as a function of UAV derived canopy area (CA1)
Model_biomass_CA1 <- lm(dry_mass_total_kg ~ canopy_area, data = pinon_data)
summary(Model_biomass_CA1)
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   22.578     16.110   1.402   0.1781  
#canopy_area    4.083      2.215   1.843   0.0818 ns

# Model for comparing biomass as a function of field measured canopy area (CA2)
Model_biomass_CA2 <- lm(dry_mass_total_kg ~ canopy_area_from_daim, data = pinon_data)
summary(Model_biomass_CA2)
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)  
#(Intercept)             22.791     16.067   1.419    0.173  
#canopy_area_from_daim    3.758      2.047   1.836    0.083 ns

#exploratory plot of this relationship
ggplot(data = pinon_data, aes(x = canopy_area_from_daim, y = dry_mass_total_kg)) +
   geom_point() +
   geom_smooth(method = "lm", formula = "y~x")
#There must be an error with the data for tree P13. It shows the highest canopy area (19.9 m2),
#but is quite small by all other metrics: RCD=6.4cm, height=1.2m, Check w/ Andy. 
#Photos of P13 are unremarkable. 

# Reformatting RCD data for plotting
# separate the diameter variables
RCD_wet <- data.frame(pinon_data$diameter_at_base_wet)
Disk_wet <- data.frame(pinon_data$disk_diameter_wet)
Disk_dry <- data.frame(pinon_data$disk_diameter_dry)

#combine RCD variable with disk variables
Disk_wet <- cbind(RCD_wet, Disk_wet)
Disk_dry <- cbind(RCD_wet, Disk_dry)

#assign the wet and dry categorigal variables
Disk_wet$Group <- "Wet"
Disk_dry$Group <- "Dry"

#rename columns
colnames(Disk_wet) <- c("RCD", "Disk_diam", "Group")
colnames(Disk_dry) <- c("RCD", "Disk_diam", "Group")

#coimbine into a single dataframe
RCD_compare <- rbind(Disk_wet, Disk_dry)

#-------------------------------------------------------------------------------
# wet base diameter as predictor of total tree biomass
# The relationship between diameter and biomass is  nonlinear (power)
# nls = nonlinear least squares, a method in R to directly fit an nonlinear model, and the default algorithm used in nls is a Gauss-Newton algorithm.
# Statistical model
model.rcd.biomass <- nls(dry_mass_total_kg ~ a*diameter_at_base_wet^b,
                               data = pinon_data,
                               start = list(a =1, b =1),
                               na.action=na.exclude)
summary(model.rcd.biomass) # Return model parameters
#Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# a  0.01841    0.03096   0.595    0.559    
# b  2.85091    0.55210   5.164 6.52e-05 ***

# Calculate prediction intervals from function predictNLS in the propagate package
preds.rcd <- data.frame(x = seq(1.0, 22.5, 0.5))  # create dataframe with example values for plotting
colnames(preds.rcd) <- c("diameter_at_base_wet")  # Label column.
prop.rcd <- predictNLS(model.rcd.biomass, newdata = preds.rcd)  # Predict upper and lower values

preds.rcd$dry_mass_total_kg <- prop.rcd$summary[,2]  # Mean prediction.
preds.rcd$Lower.CI <- prop.rcd$summary[,5]
preds.rcd$Upper.CI <- prop.rcd$summary[,6]
#THESE PREDICTION DONT LOOK RIGHT!
# https://stats.stackexchange.com/questions/162691/how-do-i-define-a-confidence-band-for-a-custom-nonlinear-function

#Figure: Total biomass as a function of RCD
ggplot(data = pinon_data, aes(x = diameter_at_base_wet, y = dry_mass_total_kg)) +
   coord_cartesian(ylim = c(0, 200)) + 
   #Chojnacky et al. 2014 (Table 5. Woodland, Pinaceae)
   stat_function(fun = function(x) exp(-2.5356+2.4349*log(x)), aes(colour = "A"), size=1, lty = "dashed") +
   #Jenkins et al., 2003 (Table 1. Softwood, Pine)
   stat_function(fun = function(x) exp(-3.2007+2.5339*log(x)), aes(colour = "B"), size=1, lty = "dashed") +
   #Grier et al. 1992 (Table 2. Pinon young+mature combined total) 
   stat_function(fun = function(x) 10^(-1.468+2.582*log10(x)), aes(colour = "c"), size=1, lty = "dashed") +
   geom_point(na.rm = TRUE, alpha=0.7) +
   stat_function(fun = function(x) (coef(summary(model.rcd.biomass))[, "Estimate"])[1]*(x)^(coef(summary(model.rcd.biomass))[, "Estimate"])[2],
                 aes(colour="D"), size = 1, lty = "solid") +
   scale_color_manual(labels = c(  "Chojnacky et al. 2014",
                                   "Jenkins et al. 2003",
                                   "Grier et al. 1992",
                                   "This study"), 
                        values = c("red",
                                   "blue", 
                                   "orange",
                                   "black")) +
  geom_ribbon(data = preds.rcd, aes( x = diameter_at_base_wet, y = dry_mass_total_kg,
                                     ymin = Lower.CI, ymax = Upper.CI), fill = "grey50", alpha = 0.3) + #plot the 95% CI
  labs(x = expression(RCD~(cm)),
       y = expression(Dry~biomass~(Kg))) +
  theme(legend.position=c(0.2, 0.8))

#-------------------------------------------------------------------------------
# Sapwood area as a function of wet disk diameter
model.disk_sapwoodarea <- nls(sapwood_area ~ a*disk_diameter_wet^b,
                         data = pinon_data,
                         start = list(a =1, b =1),
                         na.action=na.exclude)
summary(model.disk_sapwoodarea) # Return model parameters

#Parameters:
#  Estimate Std. Error t value Pr(>|t|)    
#a   0.8777     0.5605   1.566    0.135    
#b   1.6948     0.2156   7.861 3.15e-07 ***
   
# Calculate prediction intervals from function predictNLS in the propagate package
preds.SWA <- data.frame(x = seq(1.0, 22.5, 0.5))  # create dataframe with example values for plotting
colnames(preds.SWA) <- c("disk_diameter_wet")  # Label column.
prop.SWA <- predictNLS(model.disk_sapwoodarea, newdata = preds.SWA)  # Predict upper and lower values

#compile model predictions into a dataframe
preds.SWA$sapwood_area <- prop.SWA$summary[,1]  # Mean prediction
preds.SWA$Lower.CI <- prop.SWA$summary[,5] # Lower bounds of the CI (2.5%)
preds.SWA$Upper.CI <- prop.SWA$summary[,6] # Upper bounds of the CI (97.5%) 

#Figure: Total biomass as a function of RCD
ggplot(data = pinon_data, aes(x = disk_diameter_wet, y = sapwood_area)) +
   geom_point(na.rm = TRUE, alpha=0.7) +
   geom_ribbon(data = preds.SWA, aes(x = disk_diameter_wet, y = sapwood_area,
                                     ymin = Lower.CI, ymax = Upper.CI), fill = "grey50", alpha = 0.3) + #plot the 95% CI. This looks huge!
   stat_function(fun = function(x) (coef(summary(model.disk_sapwoodarea))[, "Estimate"])[1]*(x)^(coef(summary(model.disk_sapwoodarea))[, "Estimate"])[2],
                 size = 1, lty = "solid") +
   theme(axis.text = element_text(size = 10, color = "black"),
         axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
         axis.title = element_text(size = 10),
         panel.grid = element_blank(),
         plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
         panel.border = element_rect(colour = "black", fill=NA, size = 0.8)) +
    labs(x = expression(Sapwood~area~(cm^2)), 
         y = expression(Dry~biomass~(Kg)))


# examining residuals (ABS)      
model.disk.SWA.res <-  resid(model.disk_sapwoodarea)
# 
# Add residuals to dataframe, and compute normalised residuals (for plotting).
pinon_data$swa_residuals  <- model.disk.SWA.res
pinon_data$normalised_swa_residuals <- (pinon_data$swa_residuals/pinon_data$sapwood_area)

# plot the SWA residuals as a function of wet stem diameter
ggplot(data = pinon_data, aes(x=disk_diameter_wet, y=normalised_swa_residuals)) +
   geom_point() +
   geom_hline(yintercept = 0) +
   theme_classic() +
   theme(axis.text = element_text(size = 10, color = "black"),
         axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
         axis.title = element_text(size = 10),
         panel.grid = element_blank(),
         plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
         panel.border = element_rect(colour = "black", fill=NA, size = 0.8)) +
   scale_y_continuous(limits=c(-1.2,1), breaks=seq(-1.5,1,0.5), expand = c(0.01,0)) +
   scale_x_continuous(limits=c(0,24), breaks=seq(0,24,4), expand = c(0.01,0)) +
   xlab("Wet disk diameter (cm)") +
   ylab("SWA normalized residuals")

#save the figure to WD
ggsave("SWA residuals.tiff", width = 10, height = 10, units = "cm", dpi = 500)


 # ----
 # max_height as predictor of total tree biomass
 # The relationship between maximum height and biomass is  nonlinear (power)
 # nls = nonlinear least squares, a method in R to directly fit an nonlinear model, and the default algorithm used in nls is a Gauss-Newton algorithm.
 # Statistical model
 model.maxheight.biomass <- nls(dry_mass_total_kg ~ a*drone_canopy_height_max^b,
                                data = pinon_data,
                                start = list(a =1, b =1),
                                na.action=na.exclude)
 summary(model.maxheight.biomass) # Return model parameters
 #Parameters:
 #  Estimate Std. Error t value Pr(>|t|)    
 #a   5.8322     4.1856   1.393 0.180468    
 #b   1.6255     0.4129   3.937 0.000967 ***
 
 
 
 
 
###################################################
#### Old script from Juniper to recycle/remove ####
###################################################

# # Reformatting data for plotting
# # Reformat canopy area versus biomass object for plotting.
# temp1 <- rep(juniper_data$dry_mass_total, 2)
# temp2 <- c(juniper_data$CA1,
#            juniper_data$CA2)
# temp3 <- rep(c("CA1", "CA2"), each = 20)
# temp4 <- rep(1:2, each = 20)
# CA.biomass <- data.frame(temp1, temp2, temp3, temp4)
# colnames(CA.biomass) <- c("Biomass", "CA", "Source", "Group")
# remove(temp1, temp2, temp3, temp4)
# CA.biomass$Source <- as.factor(CA.biomass$Source)
# CA.biomass$Group <- as.factor(CA.biomass$Group)


# # Reformat basal area versus biomass object for plotting.
# temp1 <- rep(juniper_data$dry_mass_total, 3)
# temp2 <- c(juniper_data$BA_from_wet_diameter,
#            juniper_data$BA_from_dry_diameter,
#            juniper_data$Scale_corrected_BA_image_cm2)
# temp3 <- rep(c("BAwet", "BAdry", "BAimage"), each = 20)
# temp4 <- rep(1:3, each = 20) # grouping variable to overcome alphabetic sequencing in plots.
# BA.biomass <- data.frame(temp1, temp2, temp3, temp4)
# colnames(BA.biomass) <- c("Biomass", "BA", "Source", "Group")
# remove(temp1, temp2, temp3, temp4)
# BA.biomass$Source <- as.factor(BA.biomass$Source)
# BA.biomass$Group <- as.factor(BA.biomass$Group)



# # Stem diameter and sapwood area (SWA) datasets
# temp1 <- rep(all_stems$Scale_corrected_SWA_cm2, 2)
# temp2 <- c(all_stems$RCD_wet_cm)
# all_stems.SWA.long <- data.frame(temp1, temp2)
# remove(temp1, temp2)
# colnames(all_stems.SWA.long) <- c("Scale_corrected_SWA_cm2", "RCD_cm")
# 

# 
# 
# # Statistical modelling ----
# 
# # Linear models
# #   
# #   Model_CA1_CA2 <- lm(CA2 ~ CA1, data = juniper_data)
# #   summary(Model_CA1_CA2)
# #   
# #   Model_BAwet <- lm(dry_mass_total ~ BA_from_wet_diameter, data = juniper_data)
# #   summary(Model_BAwet)
# #   
# #   #Model_BAdry <- lm(dry_mass_total ~ BA_from_dry_diameter, data = juniper_data)
# #   #summary(Model_BAdry)
# #   
# #   Model_Zmax <- lm(dry_mass_total ~ max_height, data = juniper_data) 
# #   summary(Model_Zmax)
# #   
# #   #Model_SA <- lm(dry_mass_total ~ sapwood_area, data = juniper_data)
# #   #summary(Model_SA)
# #   
# # 
# # Mixed models
# # Model_CAfieldPlusMaxZ <- lm(dry_mass_total ~ CA1 + max_height, data = juniper_data) 
# # summary(Model_CAfieldPlusMaxZ)
# # 
# # Model_CA1 <- lm(dry_mass_total ~ 0+CA1, data = juniper_data) 
# # summary(Model_CA1) 
# 
# # Testing difference between CA1 and CA2
# CA_difference <- juniper_data$CA1 - juniper_data$CA2
# # Assessing data distribution
# length(CA_difference)
# qqnorm(CA_difference)
# qqline(CA_difference)       
# shapiro.test(CA_difference)
# # these data are signifcantly non-normal, so a non-parametric test is needed.
# # these differences are highly skewed.
# # Wilcoxon signed-rank test
# wilcox.test(juniper_data$CA1, juniper_data$CA2, paired = TRUE,
#             alternative = "less",
#             # conf.int = TRUE
# )
# 
# # Testing difference between wet versus dry diameters
# difference <- diameter_data_wo_NA$RCD_wet_cm - diameter_data_wo_NA$RCD_dry_cm
# # Assessing data distribution
# length(difference)
# qqnorm(difference)
# qqline(difference)       
# shapiro.test(difference)
# # This is a large (n=199) sample, with a near-normal distribution of differences (qqnorm plot) and signifcant result from the shapiro test.
# # one-tailed paired t-test, to test whether dry diameters are significantly smaller than wet diameters.
# wet_dry_test <- t.test(diameter_data_wo_NA$RCD_wet_cm, diameter_data_wo_NA$RCD_dry_cm, paired = TRUE, alternative = "greater")
# wet_dry_test
# 
# mean_relative_reduction_in_diameter <- mean(diameter_data_wo_NA$REC_wet_dry_diff_propotional)
# mean_relative_reduction_in_diameter
# 
# 
# # Predict SWA from stem diameter
# powermodel.SWA.diameter <- nls(Scale_corrected_SWA_cm2 ~ a*RCD_wet_cm^b,
#                                data = all_stems,
#                                start = list(a= 1, b= 1),
#                                na.action=na.exclude)
# summary(powermodel.SWA.diameter)
# 

# 
# # SWA residuals
# # examining residuals (ABS)      
# model.SWA.diameter.res <-  resid(powermodel.SWA.diameter)
# 
# # Add residuals to dataframe, and compute normalised residuals (for plotting).
# all_stems$swa_residuals  <- model.SWA.diameter.res
# all_stems$normalised_swa_residuals <- (all_stems$swa_residuals/all_stems$Scale_corrected_SWA_cm2)
# 
# # summary statistics
# mean(model.SWA.diameter.res / all_stems$Scale_corrected_SWA_cm2, na.rm=TRUE)  # Mean error in stems.
# min(model.SWA.diameter.res / all_stems$Scale_corrected_SWA_cm2, na.rm=TRUE)  # min error in stems
# max(model.SWA.diameter.res / all_stems$Scale_corrected_SWA_cm2, na.rm=TRUE)  # max error in stems.
# sd(model.SWA.diameter.res / all_stems$Scale_corrected_SWA_cm2, na.rm=TRUE)  # sd error in stems.
# 
# 
# 
# # Figure 1. Biomass predictors ----
# # Figure 1a - ESDwet as a predictor of biomass
# # Fitting two-parameter logarithmic regression models, of the form: log(biomass) = b0 + b1 log(diameter).
# # However, both the 'predict' function (in propogate) and the # the 'stat_function' (in ggplo2) seems to require a different form:
# # "y = exp(a+b*log(x))", which gives a very slightly different answer. This seems to be just due to a number handelling issue in R,
# # and only affects the visual plotting in the graphs so is not a significant issue.
# 
# # Wet
# # Form 1 (correct)
# model.ESDwet.biomass <- nls(log(dry_mass_total) ~ (a+b*log(ESD_wet)), data = juniper_data,
#                             start = list(a =0, b =1), na.action=na.exclude)
# # Form 2
# model.ESDwet.biomass2 <- nls(dry_mass_total ~ exp(a+b*log(ESD_wet)), data = juniper_data,
#                              start = list(a =0, b =1), na.action=na.exclude)
# 
# 
# # Dry
# # Form 1
# model.ESDdry.biomass <- nls(log(dry_mass_total) ~ (a+b*log(ESD_dry)), data = juniper_data,
#                             start = list(a =0, b =1), na.action=na.exclude)
# # Form 2
# model.ESDdry.biomass2 <- nls(dry_mass_total ~ exp(a+b*log(ESD_dry)), data = juniper_data,
#                              start = list(a =0, b =1), na.action=na.exclude)
# 
# # These different forms have different coefficients and error values, and compute (slightly) different values. For example:
# summary(model.ESDwet.biomass)
# overview(model.ESD.wet.biomass)
# model.ESDwet.biomass
# 
# summary(model.ESDdry.biomass)
# overview(model.ESDdry.biomass)
# model.ESDdry.biomass
# 
# summary(model.ESDwet.biomass2)
# summary(model.ESDdry.biomass2)
# 
# 
# # ESD : Biomass prediction 
# # Fit power models:
# model.ESDwet.biomass.pow <- nls(dry_mass_total ~ a*ESD_wet^b, data = juniper_data,
#                                 start = list(a =1, b =1), na.action=na.exclude)
# model.ESDdry.biomass.pow <- nls(dry_mass_total ~ a*ESD_dry^b, data = juniper_data,
#                                 start = list(a =1, b =1), na.action=na.exclude)
# 
# # Summaries power models
# summary(model.ESDwet.biomass.pow)
# # model.ESDwet.biomass.pow
# summary(model.ESDdry.biomass.pow)
# # model.ESDdry.biomass.pow
# 
# 
# # Compute predicition intervals from function predictNLS in the propagate package.
# preds.wet <- data.frame(x = seq(0.1, 50.1, 1))                                # create dataframe with example values for plotting
# colnames(preds.wet) <- c("ESD_wet")                                           # Rename vector
# preds.dry <- data.frame(x = seq(0.1, 50.1, 1))                                # create dataframe with example values for plotting
# colnames(preds.dry) <- c("ESD_dry")                                           # Rename vector
# 
# prop.wet <- predictNLS(model.ESDwet.biomass.pow, newdata = preds.wet)         # predict upper and lower values for ESD.wet model
# prop.dry <- predictNLS(model.ESDdry.biomass.pow, newdata = preds.dry)         # predict upper and lower values for ESD.dry model
# 
# preds.wet$mean <- prop.wet$summary[,2]                                        # Vectorise prediction
# preds.wet$Upper.CI <- prop.wet$summary[,6]                                    # Vectorise prediction
# preds.wet$Lower.CI <- prop.wet$summary[,5]                                    # Vectorise prediction
# preds.wet$Group <- rep(1, each = 51)                                          # Create numerical grouping variable
# preds.wet$status <- rep("wet", each = 51)                                     # Create numerical grouping variable
# 
# preds.dry$mean <- prop.dry$summary[,2]                                        # Vectorise prediction
# preds.dry$Upper.CI <- prop.dry$summary[,6]                                    # Vectorise prediction
# preds.dry$Lower.CI <- prop.dry$summary[,5]                                    # Vectorise prediction
# preds.dry$Group <- rep(2, each = 51)                                          # Create numerical grouping variable
# preds.dry$Group <- rep(2, each = 51)                                          # Create numerical grouping variable
# preds.dry$status <- rep("dry", each = 51)                                     # Create numerical grouping variable
# 
# colnames(preds.wet) <- c("ESD", "Modelled", "Upper.CI", "Lower.CI", "Group", "Status")  # Rename columns
# colnames(preds.dry) <- c("ESD", "Modelled", "Upper.CI", "Lower.CI", "Group", "Status")  # Rename columns
# 
# ESD.Biomass.preds <- rbind(preds.wet, preds.dry)                              # Combine into long dataframe for plotting prediction intervals
# ESD.Biomass.preds$Group <- as.factor(ESD.Biomass.preds$Group)                 # Specify grouping variable
# 
# # Plot ESD as predictor of biomass
# # NB. the prediciton intervals are not currently plotted, because the prediction envelopes are very large, because 95% uncertatny in the a parameter exceeds the parameter value...
# # P <- ggplot(data = ESD.Biomass.preds, aes(x = ESD, y = Upper.CI)) +
# #   geom_ribbon(aes(ymax = Upper.CI, ymin = Lower.CI, x = ESD, group=Group, colour=Group),
# #               inherit.aes = TRUE, alpha = 0.3, na.rm = TRUE, show.legend = FALSE, linetype = 0) +
# # Plot without prediciton intervals
# P <- ggplot(data = ESD.biomass, aes(x = ESD, y = Biomass)) +
#   #   coord_cartesian(ylim = c(0, 800)) + 
#   stat_function(fun = function(x) (coef(summary(model.ESDwet.biomass.pow))[, "Estimate"])[1]*(x)^(coef(summary(model.ESDwet.biomass.pow))[, "Estimate"])[2],
#                 aes(), size = 1, lty = "solid", colour=Colour.wet) +
#   stat_function(fun = function(x) (coef(summary(model.ESDdry.biomass.pow))[, "Estimate"])[1]*(x)^(coef(summary(model.ESDdry.biomass.pow))[, "Estimate"])[2],
#                 aes(), size = 1, lty = "solid", colour=Colour.dry) +
#   geom_point(data =  ESD.biomass, aes(x = ESD, y = Biomass, colour=Group), na.rm = TRUE, alpha=1, size=1) +
#   scale_color_manual(labels = c("ESDwet", "ESDdry"), values = c(Colour.wet, Colour.dry)) +
#   theme_coding() +
#   labs(y = expression(Biomass~(kg)),
#        x = expression(ESD~(cm)),
#        title = "Biomass predicted from ESD") +
#   theme(legend.position=c(0.16, 0.64)) +
#   annotate("text", x = 18, y = 600, label = "italic(Biomass)==0.002~italic(ESDwet)^3.207", parse = TRUE, color = Colour.wet, size = 3, family = "serif") +
#   annotate("text", x = 18, y = 540, label = "italic(Biomass)==0.005~italic(ESDdry)^3.038", parse = TRUE, color = Colour.dry, size = 3, family = "serif")
# biomass_from_esd <- P 
# 
# 
# # Linear models
# Model.CA1 <- lm(dry_mass_total ~ CA1, data = juniper_data)
# summary(Model.CA1)
# 
# Model.CA1.constrained <- lm(dry_mass_total ~ 0 + CA1, data = juniper_data)
# summary(Model.CA1.constrained)
# 
# Model.CA2 <- lm(dry_mass_total ~ CA2, data = juniper_data)
# summary(Model.CA2)
# 
# Model.CA2.constrained <- lm(dry_mass_total ~ 0 + CA2, data = juniper_data)
# summary(Model.CA2.constrained)
# 
# 
# par(mfrow=c(2,2))
# plot(Model.CA2.constrained, which=1:4)
# 
# 
# # power models
# Model.CA1.biomass.pow <- nls(dry_mass_total ~ a*CA1^b, data = juniper_data,
#                              start = list(a =1, b =1), na.action=na.exclude)
# Model.CA2.biomass.pow <- nls(dry_mass_total ~ a*CA2^b, data = juniper_data,
#                              start = list(a =1, b =1), na.action=na.exclude)
# 
# # Summaries power models
# summary(Model.CA1.biomass.pow)
# # Model.CA1.biomass.pow
# summary(Model.CA2.biomass.pow)
# # Model.CA2.biomass.pow
# 
# 
# # Figure 1B - Canopy Area (CA) predicting dry biomass plot
# P <- ggplot(data = CA.biomass, aes(x = CA, y = Biomass, colour = Source)) +
#   geom_point(na.rm = TRUE) +
#   scale_color_manual(values=c("black", "darkgrey")) +
#   labs(x = expression(Canopy~area~(m^2)), 
#        y = expression(Dry~biomass~(Kg)), 
#        title = "Biomass predicted from Canopy Area") +
#   theme_coding() +
#   theme(legend.position=c(0.12, 0.68)) +
#   geom_smooth(method="lm", formula = y ~ 0+ x, se=TRUE, alpha=0.25, na.rm = TRUE) +
#   coord_cartesian(ylim = c(0, 750)) +
#   annotate("text", x = 19, y = 720, label = "Y==11.06~CA[1]~~~~R^2==0.962",
#            parse = TRUE, color = "black", size = 3, family = "serif") +
#   annotate("text", x = 19, y = 650, label = "Y==9.99~CA[2]~~~~R^2==0.954",
#            parse = TRUE, color = "grey18", size = 3, family = "serif")
# ca_vs_mass <- P
# 
# 
# Model.CA1.biomass.pow <- nls(dry_mass_total ~ a*CA1^b, data = juniper_data,
#                              start = list(a =1, b =1), na.action=na.exclude)
# 
# # Figure 1C - max_height as predictor of biomass
# # The relationship between maximum height and biomass is  nonlinear, and previous studies have reported
# # power relationships (e.g. Ansley et al., 2012). 
# # nls = nonlinear least squares, a method in R to directly fit an nonlinear model, and the default algorithm used in nls is a Gauss-Newton algorithm.
# 
# # Statistical model
# model.maxheight.biomass <- nls(dry_mass_total ~ a*max_height^b,
#                                data = juniper_data,
#                                start = list(a =1, b =1),
#                                na.action=na.exclude)
# summary(model.maxheight.biomass)
# model.maxheight.biomass # Return model summary
# 
# # Calculate prediction intervals from function predictNLS in the propagate package
# preds.height <- data.frame(x = seq(0.05, 4.6, 0.1))  # create dataframe with example values for plotting
# colnames(preds.height) <- c("max_height")  # Label column.
# 
# prop.height <- predictNLS(model.maxheight.biomass, newdata = preds.height)  # Predict upper and lower values.
# 
# preds.height$Modelled_AGB <- prop.height$summary[,2]  # Mean prediction.
# preds.height$Lower.CI <- prop.height$summary[,5]
# preds.height$Upper.CI <- prop.height$summary[,6]
# 
# colnames(preds.height) <- c("Height", "Modelled_AGB", "Lower.CI", "Upper.CI")
# 
# # Create plot. NB. Equations entered manually from above model fits.
# P <- ggplot(data = preds.height, aes(x = Height, y = Upper.CI)) +
#   coord_cartesian(xlim = c(0, 6)) +
#   geom_ribbon(aes(ymax = Upper.CI, ymin = Modelled_AGB, x = Height, fill = "grey80"),
#               inherit.aes = TRUE, alpha = 0.3, na.rm = TRUE, show.legend = FALSE, linetype = 0) +
#   scale_fill_manual("",values="grey70") + 
#   stat_function(fun = function(x) (coef(summary(model.maxheight.biomass))[, "Estimate"])[1]*(x)^(coef(summary(model.maxheight.biomass))[, "Estimate"])[2],
#                 aes(), size = 1, lty = "solid", xlim=c(0,4.6)) +
#   geom_point(data =  juniper_data, aes(x = max_height, y = dry_mass_total), na.rm = TRUE) +  # Add observations.
#   labs(x = expression(Maximum~tree~height~(m)),
#        y = expression(Dry~biomass~(Kg)),
#        title = "Biomass predicted from Maximum Height") +
#   theme_coding() +
#   coord_cartesian(ylim = c(0, 750)) +
#   annotate("text", x = 1.5, y = 700, label = "italic(Y)==0.773~italic(X)^4.347",
#            parse = TRUE, color = "black", size = 3, family = "serif")
# max_height_vs_mass <- P
# 
# 
# # produce combined plot using ggarrange()[in ggpubr]
# Biomass.predictors <- ggarrange(biomass_from_esd, ca_vs_mass, max_height_vs_mass,
#                                 heights = c(10,10),
#                                 labels = c("(a)", "(b)", "(c)"),
#                                 ncol = 3, nrow = 1,
#                                 align = "v")
# 
# # Exporting figure
# png(filename="Figure_1_Biomass_predict.png", width=27, height=10, units="cm", res=500)
# plot(Biomass.predictors)
# dev.off()
# 
# jpeg(filename="Figure_1_Biomass_predict.jpg", width=27, height=10, units="cm", res=500)
# plot(Biomass.predictors)
# dev.off()
# 
# 
# 
# 
# # Figure 2. Comparison of biomass predicitons ----
# # Figure 2a - from ESD
# P <- ggplot(data = ESD.biomass, aes(x = ESD, y = Biomass, colour = Group)) +
#   coord_cartesian(ylim = c(0, 650)) + 
#   # Chojnacky et al. 2014.
#   stat_function(fun = function(x) exp(-2.709+2.1942*log(x)), aes(colour = "A"), size=1, lty = "dashed") +
#   # Jenkins et al., 2003 (single stem)
#   stat_function(fun = function(x) exp(-0.7152+1.7029*log(-6.8180+1.0222*x+1.8879*1)), aes(colour = "B"), size=1, lty = "dashed") +
#   # Jenkins et al., 2003 (multi stem)
#   stat_function(fun = function(x) exp(-0.7152+1.7029*log(-6.8180+1.0222*x+1.8879*0)), aes(colour = "C"), size=1, lty = "dashed") +
#   # Grier et al 1992.
#   stat_function(fun = function(x) 10^(-1.157+2.086*log10(x)), aes(colour = "D"), size=1, lty = "dashed") +
#   geom_point(na.rm = TRUE, alpha=0.7) +
#   stat_function(fun = function(x) (coef(summary(model.ESDwet.biomass.pow))[, "Estimate"])[1]*(x)^(coef(summary(model.ESDwet.biomass.pow))[, "Estimate"])[2],
#                 aes(), size = 1, lty = "solid", colour=Colour.wet) +
#   stat_function(fun = function(x) (coef(summary(model.ESDdry.biomass.pow))[, "Estimate"])[1]*(x)^(coef(summary(model.ESDdry.biomass.pow))[, "Estimate"])[2],
#                 aes(), size = 1, lty = "solid", colour=Colour.dry) +
#   scale_color_manual(labels = c("ESDwet",
#                                 "ESDdry",
#                                 "Chojnacky et al., 2014",
#                                 "Jenkins et al., 2003 (single stem)",
#                                 "Jenkins et al., 2003 (multi stem)",
#                                 "Grier et al., 1992"), 
#                      values = c(Colour.wet,
#                                 Colour.dry,
#                                 "black",
#                                 "darkseagreen4", 
#                                 "chocolate4",
#                                 "darkgrey")) +
#   scale_shape_manual(labels = c("ESDwet", "ESDdry"), values = c(Colour.wet, Colour.dry)) +
#   labs(x = expression(ESD~(cm)),
#        y = expression(Dry~biomass~(Kg)),
#        title = "Equivalent Stem Diameter") +
#   theme_coding() +
#   theme(legend.position=c(0.30, 0.73))
# Mass_from_ESD <- P
# Mass_from_ESD
# 
# 
# 
# 
# # Figure 2a - from ESD
# P <- ggplot(data = ESD.biomass, aes(x = ESD, y = Biomass, colour = Group)) +
#   coord_cartesian(ylim = c(0, 650)) + 
#   # Chojnacky et al. 2014.
#   stat_function(fun = function(x) exp(-2.709+2.1942*log(x)), aes(colour = "A"), size=1, lty = "dashed") +
#   # Jenkins et al., 2003 (single stem)
#   stat_function(fun = function(x) exp(-0.7152+1.7029*log(-6.8180+1.0222*x+1.8879*1)), aes(colour = "B"), size=1, lty = "dashed") +
#   # Jenkins et al., 2003 (multi stem)
#   stat_function(fun = function(x) exp(-0.7152+1.7029*log(-6.8180+1.0222*x+1.8879*0)), aes(colour = "C"), size=1, lty = "dashed") +
#   # Grier et al 1992.
#   stat_function(fun = function(x) 10^(-1.157+2.086*log10(x)), aes(colour = "D"), size=1, lty = "dashed") +
#   stat_function(fun = function(x) (coef(summary(model.ESDwet.biomass.pow))[, "Estimate"])[1]*(x)^(coef(summary(model.ESDwet.biomass.pow))[, "Estimate"])[2],
#                 aes(), size = 1, lty = "solid", colour=Colour.wet) +
#   stat_function(fun = function(x) (coef(summary(model.ESDdry.biomass.pow))[, "Estimate"])[1]*(x)^(coef(summary(model.ESDdry.biomass.pow))[, "Estimate"])[2],
#                 aes(), size = 1, lty = "solid", colour=Colour.dry) +
#   geom_point(na.rm = TRUE, alpha=0.7) +
#   # scale_shape_manual(labels = c("ESDdry", "ESDwet"), values = c(Colour.dry, Colour.wet)) +
#   scale_color_manual(labels = c("ESDdry",
#                                 "ESDwet",
#                                 "Chojnacky et al., 2014",
#                                 "Jenkins et al., 2003 (single stem)",
#                                 "Jenkins et al., 2003 (multi stem)",
#                                 "Grier et al., 1992"), 
#                      values = c(Colour.wet,
#                                 Colour.dry,
#                                 "black",
#                                 "darkseagreen4", 
#                                 "chocolate4",
#                                 "darkgrey")) +
#   labs(x = expression(ESD~(cm)),
#        y = expression(Dry~biomass~(Kg)),
#        title = "Equivalent Stem Diameter") +
#   theme_coding() +
#   theme(legend.position=c(0.30, 0.73))
# Mass_from_ESD <- P
# 
# 
# 
# 
# # Figure 2b - From canopy area
# P <- ggplot(data = CA.biomass, aes(x = CA, y = Biomass, colour = Source)) +
#   coord_cartesian(ylim = c(0, 650)) +
#   labs(x = expression(Canopy~area~(m^2)),
#        y = expression(Dry~biomass~(Kg)),
#        title = "Canopy Area") +
#   theme_coding() +
#   theme(legend.position=c(0.278, 0.780)) +
#   geom_smooth(method="lm", formula = y ~ 0+x, se=FALSE, alpha=0.15, size=1, na.rm = TRUE) +
#   stat_function(fun = function(x) 11.494*(x), aes(colour = "X"), size=1, lty = "longdash") +
#   stat_function(fun = function(x) 9.7164*(x)+37.506, aes(colour = "Y"), size=1, lty = "longdash") +
#   stat_function(fun = function(x) 8.58*(x)-6.81, aes(colour = "Z"), size=1, lty = "longdash") +
#   scale_color_manual(labels = c("J. monosperma (this study CA1)",
#                                 "J. monosperma (this study CA2)", 
#                                 "J. osteosperma (Miller et al., 1981)",
#                                 "J. occidentalis (Sabin, 2008)",
#                                 "J. pinchotii (Ansley et al. 2012)"), 
#                      values = c("black", 
#                                 "darkgrey", 
#                                 "darkorange3",
#                                 "blue",
#                                 "darkgreen")) +
#   geom_point(na.rm = TRUE, alpha=0.7)
# Mass_from_CA <- P  
# 
# # produce combined plot using ggarrange()[in ggpubr]
# Biomass.comparison <- ggarrange(Mass_from_ESD, Mass_from_CA,
#                                 heights = c(10,10),
#                                 labels = c("(a)", "(b)"),
#                                 ncol = 2, nrow = 1,
#                                 align = "v")
# 
# # Exporting figure
# # png
# png(filename="Figure_2_biomass_comparison.png", width=22, height=10.5, units="cm", res=500)
# plot(Biomass.comparison)
# dev.off()
# # jpeg
# jpeg(filename="Figure_2_biomass_comparison.jpg", width=22, height=10.5, units="cm", res=500)
# plot(Biomass.comparison)
# dev.off()
# 
# # Figure 3. Sapwood area prediction ----
# # Confidence intervals from function predictNLS in the propagate package
# # create data frame for prediction
# preds.swa <- data.frame(x = seq(0.1, 32.1, 1))  # create dataframe with example values for plotting
# colnames(preds.swa) <- c("RCD_wet_cm")
# # predict upper and lower values for model
# prop.swa <- predictNLS(powermodel.SWA.diameter, newdata = preds.swa)
# # Vectorise summaries
# preds.swa$mean <- prop.swa$summary[,2]
# preds.swa$Lower.CI <- prop.swa$summary[,5]
# preds.swa$Upper.CI <- prop.swa$summary[,6]
# # Rename columns
# colnames(preds.swa) <- c("RCD_wet_cm", "Modelled",  "Lower.CI", "Upper.CI")
# 
# # Stem diameter as predictor of sapwood area
# # Comparison of models predicting sapwood area from diameter
# # Stem level comparison
# P <- ggplot(data = preds.swa, aes(x = RCD_wet_cm, y = Modelled)) +
#   geom_ribbon(aes(ymax = Upper.CI, ymin = Lower.CI, x = RCD_wet_cm),
#               inherit.aes = TRUE, alpha = 0.3, na.rm = TRUE, show.legend = FALSE, linetype = 0) +
#   theme_coding() + 
#   coord_cartesian(ylim = c(0, 300)) + 
#   # Draw upper bound, the area of a circle
#   stat_function(fun = function(x) ifelse(x<19, pi*(x/2)^2, NA), aes(colour = "C"), size=0.5, lty = "solid") +
#   geom_point(data =  all_stems.SWA.long, aes(x = RCD_cm, y = Scale_corrected_SWA_cm2), na.rm = TRUE, alpha=1, size=1, shape=1) +
#   labs(y = expression(Sapwood~Area~(cm^2)),
#        x = expression(Stem~Diameter~""[(Wet)]~(cm)),
#        title = "SWA predicted from Stem Diameter") +
#   theme(legend.position=c(0.23, 0.76)) +
#   # Draw power model.
#   stat_smooth(aes(colour = "A"), method = 'nls', formula = 'y ~ a*x^b', method.args = list(start = c(a= 0.1, b= 2)),
#               alpha=1, size=0.5, na.rm = TRUE, se=FALSE, show.legend=FALSE) + 
#   #  Pangle et al. 2015.
#   stat_function(fun = function(x) (0.8227*x)^1.3903, aes(colour = "D"), size=1, lty = "dashed", show.legend = TRUE) +   
#   # McDowell et al., 2008.
#   stat_function(fun = function(x) (4.3*x)-9.8, aes(colour = "E"), size=1, lty = "dashed") +
#   scale_color_manual(labels = c("This study",
#                                 "Upper limit",
#                                 "Pangle et al., 2015",
#                                 "McDowell et al., 2008"),
#                      values = c("black",
#                                 "grey",
#                                 "darkgreen",
#                                 "blue")) +
#   guides(colour = guide_legend(override.aes = list(size = 1)))
# swa_from_stem_diameter <- P
# 
# # Exporting figure
# # png
# png(filename="Figure_X_SWA.png", width=11, height=11, units="cm", res=500)
# plot(swa_from_stem_diameter)
# dev.off()
# # jpg
# jpeg(filename="Figure_X_SWA.jpg", width=11, height=11, units="cm", res=500)
# plot(swa_from_stem_diameter)
# dev.off()
# 
# # create residuals plot
# P <- ggplot(data = all_stems, aes(x = RCD_wet_cm, y = normalised_swa_residuals))
# P <- P + theme_coding()
# P <- P + coord_cartesian(xlim = c(0, 34), ylim = c(-10, 2), expand=FALSE)
# P <- P + scale_y_continuous(breaks = seq(-10,2,2))
# P <- P + geom_point(na.rm = TRUE, alpha=1, size=1, shape=1)
# P <- P + geom_abline(intercept = 0, slope = 0, linetype = "solid")
# P <- P + labs(x = expression(Stem~Diameter~""[(Wet)]~(cm)), 
#               y = expression(Normalised~Residual~SWA), 
#               title = "Normalised SWA residuals")
# swa_residual_plot <- P
# swa_residual_plot
# 
# # exporting plot
# png(filename="Figure_X_SWA_residuals.png", width=12, height=8, units="cm", res=500)
# plot(swa_residual_plot)
# dev.off()
# 
# # produce combined plot using ggarrange()[in ggpubr]
# swa_combined <- ggarrange(swa_from_stem_diameter, swa_residual_plot,
#                           heights = c(16,9),
#                           labels = c("(a)", "(b)"),
#                           ncol = 1, nrow = 2,
#                           align = "v")
# 
# # export plot
# # png
# png(filename="Figure_3_SWA.png", width=10, height=14.5, units="cm", res=500)
# plot(swa_combined)
# dev.off()
# # jpg
# jpeg(filename="Figure_3_SWA.jpg", width=10, height=14.5, units="cm", res=500)
# plot(swa_combined)
# dev.off()
# 
# 
# 
# # Figure 4. Method comparisons ----
# # Figure 4a - Canopy area 1 versus canopy area 2
# (CA1_vs_CA2 <- ggplot(data = juniper_data, aes(x = CA1, y = CA2)) +
#    geom_point(na.rm = TRUE) +
#    labs(x = expression(CA["1"]~(m^2)),
#         y = expression(CA["2"]~(m^2)),
#         title = "Canopy Area") +
#    theme_coding() +
#    coord_cartesian(ylim = c(0, 62),xlim = c(0, 62), expand=FALSE) +
#    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#    geom_line(method=lm, stat="smooth", se=TRUE, color=1, size=0.5, na.rm = TRUE) +
#    stat_poly_eq(aes(label = paste("atop(", ..eq.label.., ",", ..rr.label.., ")", sep="")),
#                 formula = y ~ x, na.rm = TRUE, rr.digits = 3, size = 2.5, parse = TRUE,
#                 label.x.npc = 0.1, label.y.npc = 0.9))  # Specify relative position of equation
# 
# 
# 
# 
# # Figure 4b - wet versus dry stem diameters
# (Dwet_vs_Ddry <- ggplot(data = diameter_data_wo_NA, aes(x = RCD_wet_cm, y = RCD_dry_cm)) +
#     geom_point(na.rm = TRUE, shape = 1, size = 0.8) +
#     labs(x = expression(Wet~Stem~Diameter~(cm)),
#          y=expression(Dry~Stem~Diameter~(cm)),
#          title = "Stem Diameter") + 
#     theme_coding() +
#     coord_cartesian(ylim = c(0, 33),xlim = c(0, 33), expand=FALSE) +
#     geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#     geom_line(method=lm, stat="smooth", se=TRUE, color=1, size=0.5,  na.rm = TRUE) +
#     stat_poly_eq(aes(label = paste("atop(", ..eq.label.., ",", ..rr.label.., ")", sep="")),
#                  formula = y ~ x, na.rm = TRUE, rr.digits = 3, size = 2.5, parse = TRUE,
#                  label.x.npc = 0.1, label.y.npc = 0.9))  # Specify relative position of equation
# 
# 
# # Figure 4c - BA_dry versus BA_image
# # Creating basal area dataframe for subsequent plotting
# Tree_ID <- c(rep(juniper_data$Tree_ID,3))
# Tree_ID <- as.factor(Tree_ID)
# BA_source <- c(rep("BAwet" , 20), rep("BAdry" , 20), rep("BAimage" , 20))
# BA <- c(juniper_data$BA_from_wet_diameter, juniper_data$BA_from_dry_diameter, juniper_data$Scale_corrected_BA_image_cm2)
# BAlong <- data.frame(Tree_ID, BA_source, BA)
# #  Vector for ordering plots.
# BA_level_order <- c('BAwet', 'BAdry', 'BAimage')
# 
# # Actual plot
# (BAcomparison <- ggplot(drop_na(BAlong), aes(x = reorder(Tree_ID, BA), y = BA,
#                                              fill = factor(BA_source, level = BA_level_order))) +
#     geom_col(position = "dodge", na.rm = TRUE) +
#     theme_coding() +
#     theme(legend.text.align=0, 
#           legend.position=c(0.29, 0.67)) +
#     scale_color_manual(values = c(Colour.wet, Colour.dry, Colour.image)) +
#     scale_fill_manual(values = c(Colour.wet, Colour.dry, Colour.image), 
#                       labels = c(expression(BA[Wet]),
#                                  expression(BA[Dry]),
#                                  expression(~BA[Image]))) +
#     labs(x = expression(Tree~ID), 
#          y = expression(Basal~Area~(cm^2)),
#          title = "Basal Area",
#          fill = "Method") +
#     scale_y_continuous(lim = c(0, 2000)))
# 
# 
# # produce combined plot using ggarrange()[in ggpubr]
# method.comparisons <- ggarrange(CA1_vs_CA2, Dwet_vs_Ddry, BAcomparison,
#                                 heights = c(10,10),
#                                 labels = c("(a)", "(b)", "(c)"),
#                                 ncol = 3, nrow = 1,
#                                 align = "v")
# 
# # Exporting figure
# # png
# png(filename="Figure_4_method_comparisons.png", width=18, height=6, units="cm", res=500)
# plot(method.comparisons)
# dev.off()
# # jpeg
# jpeg(filename="Figure_4_method_comparisons.jpg", width=18, height=6, units="cm", res=500)
# plot(method.comparisons)
# dev.off()
# 
# 
# 
# # Figure S1. Tree morphology ----
# # create plot A
# # Fit model
# x <- juniper_data$CA1
# y <- juniper_data$max_height
# 
# model.CA1.Zmax <- nls(y ~ a*x^b, data = juniper_data, start = list(a =1, b =1), na.action=na.exclude)
# 
# summary(model.CA1.Zmax)
# 
# # Create plot  
# (Tree_morphology_1 <- ggplot(data = juniper_data, aes(x = CA1, y = max_height)) +
#     stat_function(fun = function(x) (coef(summary(model.CA1.Zmax))[, "Estimate"])[1]*(x)^(coef(summary(model.CA1.Zmax))[, "Estimate"])[2],
#                   aes(), size = 1, lty = "solid") +
#     geom_point(na.rm = TRUE) +  # Add observations.
#     labs(x = expression(Canopy~area~(m^2)),
#          y = expression(Maximum~height~(m)),
#          title = "Canopy Area") +
#     theme_coding() +
#     annotate("text", x=12, y=5.5, label="italic(Y)==1.038~italic(X)^0.401",
#              parse=TRUE, color="black", size=3, family = "serif")
# )
# 
# 
# # create plot B.  mean canopy diameter as predictor of max height.
# # create mean canopy diameter variable
# juniper_data$CanDiaMean <- (juniper_data$CanDia1 + juniper_data$CanDia2) / 2
# 
# # Fit model
# x <- juniper_data$CanDiaMean
# y <- juniper_data$max_height
# 
# Model.CanDiaMean.Zmax <- nls(y ~ a*x^b, data = juniper_data, start = list(a =1, b =1), na.action=na.exclude)
# 
# summary(Model.CanDiaMean.Zmax)
# 
# # Create plot  
# (Tree_morphology_2 <- ggplot(data = juniper_data, aes(x = CA1, y = max_height)) +
#     geom_point(na.rm = TRUE) +  # Add observations.
#     labs(x = expression(Mean~canopy~diameter~(m)),
#          y = expression(Maximum~height~(m)),
#          title = "Canopy Diameter") +
#     theme_coding() 
# )    
# 
# 
# 
# # produce combined plot using ggarrange()[in ggpubr]
# FigureS1 <- ggarrange(Tree_morphology_1, Tree_morphology_2,
#                       heights = c(10,10),
#                       labels = c("(a)", "(b)"),
#                       ncol = 2, nrow = 1,
#                       align = "v")
# 
# # Exporting figure
# png(filename="Figure_S1.png", width=16, height=9, units="cm", res=500)
# plot(FigureS1)
# dev.off()
# 
# 
# # Figure S2. Proportion of mass in <3 cm diameter ----
# # create plot
# (Prop_small_mass <- ggplot(data = juniper_data, aes(x = dry_mass_total, y = proportion_small_dry)) +
#    geom_point (na.rm = TRUE) +
#    labs(x=expression(Total~dry~biomass~(kg)),
#         y=expression(Proportion~"<"~3~cm~diameter),
#         title = "Proportion of biomass in < 3 cm component") +
#    theme_coding() +
#    coord_cartesian(ylim = c(0, 1)))
# 
# # Exporting figure
# png(filename="Figure_S2.png", width=12, height=12, units="cm", res=500)
# plot(Prop_small_mass)
# dev.off() 
# 
# 
# 
# # Statistical analsysis of the < 3 cm component (not the most elegant code, but it works!).
# MassThreshold <- 40 # Kg. NB. Breakpoint sleected from visual analysis of plot.
# 
# juniper_data %>% 
#   filter(dry_mass_total > MassThreshold) %>%
#   select(proportion_small_dry) %>%
#   summarise(Prop_small_above_threhold_mean = mean(proportion_small_dry, na.rm = TRUE))
# 
# juniper_data %>% 
#   filter(dry_mass_total > MassThreshold) %>%
#   select(proportion_small_dry) %>%
#   summarise(Prop_small_above_threhold_SD = sd(proportion_small_dry, na.rm = TRUE))
# 
# juniper_data %>% 
#   filter(dry_mass_total < MassThreshold) %>%
#   select(proportion_small_dry) %>%
#   summarise(Prop_small_below_threhold_mean = mean(proportion_small_dry, na.rm = TRUE))
# 
# juniper_data %>% 
#   filter(dry_mass_total < MassThreshold) %>%
#   select(proportion_small_dry) %>%
#   summarise(Prop_small_below_threhold_SD = sd(proportion_small_dry, na.rm = TRUE))
# 
# 
# 
# 
# # Figure XXX. ESD : SWA   ----  
# # ESD as predictor of total sapwood area
# model.ESD.SWA <- nls(Scale_corrected_SWA_cm2 ~ a*ESD_wet^b,
#                      data = juniper_data,
#                      start = list(a= 0.1, b= 2),
#                      na.action=na.exclude)
# summary(model.ESD.SWA)
# 
# # Tree-level (ESD) plot comparison
# P <- ggplot(data = juniper_data, aes(x = ESD_wet, y = Scale_corrected_SWA_cm2))
# P <- P + theme_coding()
# P <- P + coord_cartesian(ylim = c(0, 300))
# # Draw upper bound, the area of a circle
# P <- P + stat_function(fun = function(x) ifelse(x<19, pi*(x/2)^2, NA), aes(colour = "C"), size=0.5, lty = "solid")
# P <- P + geom_point(na.rm = TRUE, alpha=1, size=1, shape=1)
# P <- P + labs(x = expression(Stem~Diameter~""[(Wet)]~(cm)), y = expression(Sapwood~Area~(cm^2)), title = "SWA from Equivalent Stem Diameter")
# P <- P + theme(legend.position=c(0.22, 0.80))
# # Draw power model.
# P <- P + stat_smooth(aes(colour = "A"), method = 'nls', formula = 'y ~ a*x^b', method.args = list(start = c(a= 0.1, b= 2)),
#                      alpha=1, size=0.5, na.rm = TRUE, se=FALSE, show.legend=FALSE)
# #  Pangle et al. 2015.
# P <- P + stat_function(fun = function(x) (0.8227*x)^1.3903, aes(colour = "D"), size=1, lty = "dashed", show.legend = TRUE)   
# # McDowell et al., 2008.
# P <- P + stat_function(fun = function(x) (4.3*x)-9.8, aes(colour = "E"), size=1, lty = "dashed")  
# P <- P + scale_color_manual(labels = c("This Study (power)",
#                                        "Upper limit",
#                                        "Pangle et al., 2015",
#                                        "McDowell et al., 2008"),
#                             values = c("black",
#                                        "grey",
#                                        "darkgreen",
#                                        "blue"))
# P <- P + guides(colour = guide_legend(override.aes = list(size = 1)))
# SWA_from_stem_diameter <- P
# 
# # Exporting figure
# # PLOT NO LONGER EXPORTED, BECAUSE WE HAVE DECIDED THAT ESD-LEVEL RESULTS ARE UN-INFORMATIVE
# # png(filename="Figure_S5_SWA-ESD.png", width=11, height=11, units="cm", res=500)
# # plot(SWA_from_stem_diameter)
# # dev.off()
# 
# 
# # Figure X. Other SWA relationships. ----
# # sapwood_area versus biomass
# (SWA_vs_mass <- ggplot(data = juniper_data,
#                        aes(x = Scale_corrected_SWA_cm2, y = dry_mass_total)) +
#    geom_point( na.rm = TRUE) +
#    labs(x = expression(SWA~(cm^2)),
#         y = expression(Dry~biomass~(Kg)),
#         title = "(a) Biomass") +
#    theme_coding()
# )
# 
# # CA as predictor of SA 
# (SWA_vs_CA1 <- ggplot(data = juniper_data, aes(x = Scale_corrected_SWA_cm2, y = CA1)) +
#     geom_point(na.rm = TRUE) +
#     labs(x = expression(SWA~(cm^2)),
#          y = expression(Canopy~area~(m^2)),
#          title = "(b) Canopy Area") +
#     theme_coding()
# )
# 
# # SA with <3 cm component 
# (SWA_vs_SmallMass <- ggplot(data = juniper_data, aes(x = Scale_corrected_SWA_cm2, y = dry_mass_small)) +
#     geom_point(na.rm = TRUE) +
#     labs(x = expression(SWA~(cm^2)),
#          y = expression(Mass~of~"<"~3~cm~component~(Kg)),
#          title = "(c) mass of < 3 cm component") +
#     theme_coding()
# )
# 
# # Mass with number of stems 
# (Stems_vs_Mass <- ggplot(data = juniper_data, aes(x = stem_number, y = dry_mass_total)) +
#     geom_point(na.rm = TRUE) +
#     labs(x = expression(Number~of~stems),
#          y = expression(Biomass~(Kg)),
#          title = "(d) Stems versus mass") +
#     theme_coding()
# )
# 
# 
# 
# # Combining into multiplot figure 
# FigureS6 <- grid.arrange(SWA_vs_mass, SWA_vs_CA1, SWA_vs_SmallMass, Stems_vs_Mass, nrow = 2)
# 
# # Exporting figure
# png(filename="NOT_FOR_PUBLICATION_FIG1.png", width=16, height=16, units="cm", res=500)
# plot(FigureS6)
# dev.off()
# 
# 
# 
# 
# CA.biomass2 <- rbind(CA.biomass, Miller_CA_mass)
# 
# # plot
# (Mass_from_CA_all <- ggplot(data = CA.biomass2, aes(x = CA, y = Biomass, colour = Source)) +
#     coord_cartesian(ylim = c(0, 1000)) +
#     labs(x = expression(Canopy~area~(m^2)),
#          y = expression(Dry~biomass~(Kg)),
#          title = "Canopy Area") +
#     theme_coding() +
#     theme(legend.position=c(0.18, 0.84)) +
#     geom_smooth(method="lm", formula = y ~ 0+ x, se=FALSE, alpha=0.15,
#                 size=1, na.rm = TRUE) +
#     #stat_function(fun = function(x) 11.494*(x), aes(colour = "Miller et al., 1981 (constrained)"), size=1, lty = "longdash") +
#     scale_color_manual(labels = c("CA1", 
#                                   "CA2", 
#                                   "Miller et al., 1981"), 
#                        values = c("black", 
#                                   "darkgrey", 
#                                   "darkorange3")) +
#     geom_point(na.rm = TRUE, alpha=0.7)
# )
# 
# # Combining into multiplot figure 
# Figure_X_CanopyArea <- grid.arrange(Mass_from_CA_all, nrow = 1)
# 
# # Exporting figure
# png(filename="Figure_X_CanopyArea.png", width=16, height=16, units="cm", res=500)
# plot(Figure_X_CanopyArea)
# dev.off()