# Leaf litter weevil richness increases with altitude in a tropical-temperate
# transitional forest in El Cielo Biosphere Reserve, northeastern Mexico
# Damián Villaseñor-Amador
# Milan Janda
# Madai Rosas Mejia
# Fatima Magdalena Sandoval Becerra
# Juan Jose Morrone Lupi*
# corresponding author: morrone@ciencias.unam.mx

#________________________________________________________
#### Loading libraries and data ####
#________________________________________________________

# Please ensure this R script is in the same directory that the database: El_Cielo_weevils_data.xlsx
# Github repository of this project: https://github.com/Damian-VA/CURBC_github

# Set working directory
setwd(".")
dir()

# Set scientific notation off
options(scipen=999)

# Libraries needed:
library(readxl) #to read excel files (the mother database is in .xlsx)
library(vegan) #for species accumulation curves, diversity and CCA
library(fossil) #to calculate Chao1 and Chao2
library(iNEXT) #to compute diversity estimates using Hill numbers
library(ggplot2) #to visualize iNEXT plots
library(ggrepel) #plot pretty CCA plot
library(betapart) #to compute beta diversity
library(raster) #to extract climatic values for each coordinate
library(usdm) #to stack tif or asc WorldClim layers
library(PerformanceAnalytics) #to plot Pearson or Spearman colinearity among all variables
library(MASS) #to run stepwise backwards simplification of models
library(car) #for the Anova() function to partition explained deviance of model
library(picante) #for randomizeMatrix function (community data matrix randomization previous to beta analyses)
library(mgcv) #for generalized additive models (gam) to plot fitted line to beta scatterplot
library(MuMIn) #for model averaging (deltaAICc) and model dredging
library(dplyr) #for easy data manipulation
library(tidyr) #for easy data manipulation
library(plyr) #for ddply function to easily restructre dataframes
library(ggplotify) #combine plots + ggplots
library(patchwork) #easily combine multiple ggplots

# Database with all the processed samples of El Cielo Biosphere Reserve
# (it contains samples with leaf litter weevils + non-leaf litter weevils
# +samples without weevils)
all_data <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                              sheet = "El_Cielo_weevils_data", col_names = T)

# Database that contains leaf litter weevil samples
weevil_data <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                              sheet = "Leaf_litter_weevils_data", col_names = T)

# Transform database object into the data.frame class
weevil_data = as.data.frame(weevil_data)

#____________________________________________________
#### Alpha diversity (figure 1) ####
#____________________________________________________

# Load community abundance matrix
weevil_abundance_matrix <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                               sheet = "Community_abundance_matrix", col_names = T)

# Set the first column (elevational floors) as row names and the
# elevational floors (200:2000 masl) as column names
weevil_abundance_matrix = as.data.frame(weevil_abundance_matrix)
weevil_abundance <- weevil_abundance_matrix[,-1]
rownames(weevil_abundance) <- weevil_abundance_matrix[,1]
names(weevil_abundance) <- seq(200,2000,200)

# iNEXT inputs 
# - according to Hsieh et al., 2016 -
# When there are N assemblages, input data consist of an S by N 
# abundance matrix (S = species as rows, N = assemblages/sites
# as columns) or N lists of species abundances. In iNEXT, 
# this type of abundance data (from 1 to N assemblages) is 
# specified by an argument datatype = "abundance".
# An abundance matrix was used as input, and not a sample
# occurrence matrix (datatype="incidence_freq") because
# Jones et al. (2008, 2012) have used abundance matrices as
# input in their Chiapas' leaf-litter weevils and El Cielo
# phytophagous Apionidae studies. Moctezuma et al. (2016) 
# used an abundance matrix as input to calculate the sample 
# coverage of Carabidae in the MTZ. Ohwake et al. (2021) used
# an abundance matrix to estimate diversity in a beetle-spider
# assamblage in Mt. Fuji, Japan.
head(weevil_abundance)

weevdiv = iNEXT(weevil_abundance, #input matrix
                           q=c(0, #calculate richness (q=0)
                               1, #Shannon diversity (q=1)
                               2  #Simpson diversity (q=2)
                               ),
                           datatype="abundance", #abundance matrix is input
                           se = T #bootstrap -> 95% CI for each estimate
                           )
weevdiv$DataInfo

#For explained iNEXT output check: 
#https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html

# Non-standardized diversity estimates along with 95% CI (qD, qD.LCL, qD.UCL) 
# for each site (each elevational floor). First step is to create the non-standardized 
# richness dataframe. For that, we extract data of the observed abundance of each site.
s200 = weevdiv$iNextEst$size_based[weevdiv$iNextEst$size_based$Assemblage=="200" &
                                   weevdiv$iNextEst$size_based$Method=="Observed" &
                                   weevdiv$iNextEst$size_based$Order=="0",]
s400 = weevdiv$iNextEst$size_based[weevdiv$iNextEst$size_based$Assemblage=="400" &
                                   weevdiv$iNextEst$size_based$Method=="Observed" &
                                   weevdiv$iNextEst$size_based$Order=="0",]
s600 = weevdiv$iNextEst$size_based[weevdiv$iNextEst$size_based$Assemblage=="600" &
                                   weevdiv$iNextEst$size_based$Method=="Observed" &
                                   weevdiv$iNextEst$size_based$Order=="0",]
s800 = weevdiv$iNextEst$size_based[weevdiv$iNextEst$size_based$Assemblage=="800" &
                                   weevdiv$iNextEst$size_based$Method=="Observed" &
                                   weevdiv$iNextEst$size_based$Order=="0",]
s1000 = weevdiv$iNextEst$size_based[weevdiv$iNextEst$size_based$Assemblage=="1000" &
                                   weevdiv$iNextEst$size_based$Method=="Observed" &
                                   weevdiv$iNextEst$size_based$Order=="0",]
s1200 = weevdiv$iNextEst$size_based[weevdiv$iNextEst$size_based$Assemblage=="1200" &
                                   weevdiv$iNextEst$size_based$Method=="Observed" &
                                   weevdiv$iNextEst$size_based$Order=="0",]
s1400 = weevdiv$iNextEst$size_based[weevdiv$iNextEst$size_based$Assemblage=="1400" &
                                   weevdiv$iNextEst$size_based$Method=="Observed" &
                                   weevdiv$iNextEst$size_based$Order=="0",]
s1600 = weevdiv$iNextEst$size_based[weevdiv$iNextEst$size_based$Assemblage=="1600" &
                                   weevdiv$iNextEst$size_based$Method=="Observed" &
                                   weevdiv$iNextEst$size_based$Order=="0",]
s1800 = weevdiv$iNextEst$size_based[weevdiv$iNextEst$size_based$Assemblage=="1800" &
                                   weevdiv$iNextEst$size_based$Method=="Observed" &
                                   weevdiv$iNextEst$size_based$Order=="0",]
s2000 = weevdiv$iNextEst$size_based[weevdiv$iNextEst$size_based$Assemblage=="2000" &
                                   weevdiv$iNextEst$size_based$Method=="Observed" &
                                   weevdiv$iNextEst$size_based$Order=="0",]

#Bind all site rows into a single dataframe
nonstanrich = rbind(s200,s400,s600,s800,s1000,s1200,s1400,s1600,s1800,s2000)
#Add a column with the sites' elevation
nonstanrich$site = seq(200,2000,200)
rownames(nonstanrich) = nonstanrich$site
#Drop method (observed, interpolated, extrapolated) and order (Hill numbers) columns
nonstanrich = nonstanrich[,!(names(nonstanrich) %in% c("Method","Order.q","Assemblage"))]

#Non-standardized  richness dataframe
#m = size (site weevil abundance in this example)
#qD = observed richness
#qD.LCL = lower 95% CI
#qD.UCL = upper 95% CI
#SC = sample coverage
nonstanrich
#        m qD    qD.LCL    qD.UCL        SC    SC.LCL    SC.UCL site
# 200   38  7  4.703359  9.296641 0.8965999 0.8285826 0.9646172  200
# 400   54 12  9.144602 14.855398 0.9266182 0.8735366 0.9796998  400
# 600  126 15 10.527789 19.472211 0.9445711 0.9046113 0.9845310  600
# 800   23  6  4.401820  7.598180 0.9632107 0.8829618 1.0000000  800
# 1000 298 24 20.971651 27.028349 0.9900225 0.9800316 1.0000000 1000
# 1200 338 16 12.516117 19.483883 0.9852246 0.9745464 0.9959029 1200
# 1400 557 25 20.869530 29.130470 0.9892730 0.9820380 0.9965080 1400
# 1600 966 33 29.480519 36.519481 0.9927643 0.9895875 0.9959412 1600
# 1800 790 26 22.832220 29.167780 0.9936773 0.9891006 0.9982540 1800
# 2000 814 22 19.461757 24.538243 0.9950950 0.9916990 0.9984911 2000

# Figure 1A: alpha diversity, non-standardized with 95% CI
plot(nonstanrich$qD ~ nonstanrich$site,
     ylim=c(0,max(nonstanrich$qD.UCL)),
     ylab="Leaf litter weevil richness",
     xlab="Elevation (m a.s.l)",
     xaxt="n",
     bty="l",
     cex.axis=1.5,
     cex.lab=1.5,
     bg="white",
     las=1)
axis(1, at = seq(200,2000,200),labels = seq(200,2000,200),
     cex.axis=1.5, cex.lab=1.5)
arrows(x0=seq(200,2000,200), x1=seq(200,2000,200), 
       y0=nonstanrich$qD.LCL, y1=nonstanrich$qD.UCL, 
       col="black", lwd=1, code=3, angle = 90)
alpha.non.stan <- recordPlot()

# Load leaf litter weevil richness per elevational floor database
weevil_richness_elevation <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                             sheet = "richness_per_elevation", col_names = T)

# Figure 1B: Richness curve + abundance given by dot size
rich.curve = ggplot(weevil_richness_elevation, #data
       aes(x=plot, y=richness)) + #x variable, y variable
       geom_smooth(colour = "black",      # black line
              se = T,                     # turn off confidence band
              method = glm,               # use lm, loess, gam, glm
              n = 1000,                   # the larger the n, the smoother the curve
              size = 1.5,                 # line width
              formula = y~I(x^1)+I(x^2)+I(x^3)) + # specify formula -> best deltaAICc = quadratic
       geom_point(aes(size=Abundance), alpha = 1) +
       scale_size(range = c(2,8)) + #variation between point sizes
       labs(x = "Elevation (m a.s.l)", 
            y = "Litter weevil richness") +
       scale_y_continuous(name ="",
                  limits = c(0,max(nonstanrich$qD.UCL)),
                  breaks = seq(0,30,10)) +
       scale_x_continuous(name ="Elevation (m a.s.l)",
                  limits = c(200,2000),
                  breaks = seq(200,2000,200)) +
       theme_classic() + 
       theme(axis.text=element_text(size=18,colour="black"),
             axis.title.x=element_text(size=18,colour="black",
                                       margin = margin(t = 5, r = 0, b = 0, l = 0,
                                                    unit = "mm")),
             legend.title=element_text(face="bold",size=18),
             legend.text=element_text(size=15),
             legend.position = c(0.15,0.9),
             plot.margin = unit(c(2,0,0.55,0),"cm"),
             axis.ticks.length.x = unit(0.25, "cm")
             )

#Combine both richness plots
figure_richness = cowplot::plot_grid(alpha.non.stan, rich.curve, labels = c("A","B"), 
                   hjust = -3, vjust = 2, label_size = 30)

#Save richness plots as figure
ggsave("./Figures/Figure_1.tiff", figure_richness, 
       width=12000, height=6000, units="px", dpi=600, compression="lzw")

#___________________________________________________________________
#### Models to test relationship between elevation and richness ####
#___________________________________________________________________

# Following methodology of Joaqui et al. 2021:

# (a) linear, (b) quadratic, (c) cubic, (d) exponential, (e) null
# They only tested exponential for neotropical affinity

# Models with ΔAICc < 2 and a normal distribution of residuals
# (Shapiro– Wilk test, p > 0.05) were selected and considered 
# to be models with good fit (Hurvich & Tsai, 1989).

# All weevil richness
linear_all = glm(all$richness ~ all$plot, family = poisson(link=log)) #linear
quad_all = glm(all$richness ~ all$plot+I(all$plot^2), family = poisson(link=log)) #quadratic
cubic_all = glm(all$richness ~ all$plot+I(all$plot^2)+I(all$plot^3), family = poisson(link=log)) #cubic
expo_all = glm(all$richness ~ I(all$plot^2), family = poisson(link=log)) #exponential
null_all = glm(all$richness ~ 1, family = poisson(link=log)) #null

# Model average (ΔAICc, weighted AICc)
round(summary(MuMIn::model.avg(linear_all,quad_all,cubic_all,expo_all,null_all))$msTable,2)
#        df logLik   AICc delta weight
# 12      3 -57.32 122.13  0.00   0.61
# 123     4 -56.35 123.37  1.24   0.33
# 1       2 -61.09 126.89  4.76   0.06
# 2       2 -65.85 136.40 14.27   0.00
# (Null)  1 -86.94 176.09 53.96   0.00

# Distribution of residuals (normal if S-W p > 0.05)
shapiro.test(linear_all$residuals) #p-value = 0.3929
shapiro.test(quad_all$residuals) #p-value = 0.5046
shapiro.test(cubic_all$residuals) #p-value = 0.2735
shapiro.test(expo_all$residuals) #p-value = 0.5464
shapiro.test(null_all$residuals) #p-value = 0.851

#______________________________________________
#### Elevation ranges per genus (figure 2) ####
#______________________________________________

#Load leaf litter weevil database ordered by elevational range
elev_ranges <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                               sheet = "weevils_ordered_by_elevation", col_names = T)
elev_ranges = as.data.frame(elev_ranges)

#Add weevil genera names
#Subset a single unique row for each morphospecies
genera.unique = litter[!duplicated(litter[,c("SpeciesCode")]),]

#Replace "TribeNewGenus" with shorter label
genera.unique$Genus =  gsub('CryptorhynchiniNewGenus', 'Crypt. gen. nov. ', genera.unique$Genus)
genera.unique$Genus =  gsub('LymantiniNewGenus', 'Lyman. gen. nov. ', genera.unique$Genus)
genera.unique$Genus =  gsub('ConotracheliniNewGenus', 'Conot. gen. nov. ', genera.unique$Genus)

#Add column with morphospecies taxonomical names
gen_elev_abun = merge(elev_ranges, genera.unique[,c("SpeciesCode", "Genus")],
                      by="SpeciesCode", all=T, sort=T)

#Create column of elevational floor median for each morph
agg_elev <- plyr::ddply(gen_elev_abun, ~ Genus + Plot, function(x){c(Abundance=nrow(x))})
elev_med = tapply(agg_elev$Plot, agg_elev$Genus, median)
elev_med = data.frame(elev_med)
elev_med$Genus = rownames(elev_med)
elev_dist = merge(gen_elev_abun, elev_med[,c("Genus","elev_med")],
                  by = "Genus", all = F, sort = T)
elev_dist = elev_dist[with(elev_dist, order(elev_med)),]

#Set morphospecies level order as the ascending median of their
#elevational range
elev_dist$Genus = factor(elev_dist$Genus,
                               levels = unique(elev_dist$Genus[
                                 order(elev_dist$elev_med)
                               ]), ordered = T)

#Italicize scientific names
elev_dist$GenusItalics = paste0("italic('",elev_dist$Genus,"')")
elev_dist$GenusItalics = gsub("italic('Crypt. gen. nov. 1')",
                                "'Crypt. gen. nov. 1'",
                                 elev_dist$GenusItalics,
                                 fixed = T) 
elev_dist$GenusItalics = gsub("italic('Crypt. gen. nov. 2')",
                                "'Crypt. gen. nov. 2'",
                                 elev_dist$GenusItalics,
                                 fixed = T) 
elev_dist$GenusItalics = gsub("italic('Crypt. gen. nov. 3')",
                                "'Crypt. gen. nov. 3'",
                                 elev_dist$GenusItalics,
                                 fixed = T) 
elev_dist$GenusItalics = gsub("italic('Crypt. gen. nov. 4')",
                                "'Crypt. gen. nov. 4'",
                                 elev_dist$GenusItalics,
                                 fixed = T) 
elev_dist$GenusItalics = gsub("italic('Lyman. gen. nov. 1')",
                                "'Lyman. gen. nov. 1'",
                                 elev_dist$GenusItalics,
                                 fixed = T) 
elev_dist$GenusItalics = gsub("italic('Conot. gen. nov. 1')",
                                "'Conot. gen. nov. 1'",
                                 elev_dist$GenusItalics,
                                 fixed = T) 

#Set morphospecies level order as the ascending median of their
#elevational range -> for italicized scientific names
elev_dist$GenusItalics = factor(elev_dist$GenusItalics,
                               levels = unique(elev_dist$GenusItalics[
                                 order(elev_dist$elev_med)
                               ]), ordered = T)

# ggplot elevation distribution of weevil genus (violin plot)
# - modified from Perez-Rosales et al. (2022) - 
Elevation_distribution <-
  ggplot(elev_dist, aes(x=Genus, y=Plot)) +
  # geom_path(aes(colour = "black"), size = 1, lineend = "round") +
  geom_violin(aes(fill = Affinity), scale = "width", position = "dodge", trim = T, 
              adjust=0.8, show.legend = T) +
  # geom_dotplot(binaxis = "y", binwidth = 2, stackdir = "center", show.legend = F, dotsize = 1) +
  stat_summary(data=agg_elev, fun = "median", geom = "point", colour = "black", size = 2) +
  # scale_y_reverse(name ="Elevation (m a.s.l)",lim=c(2000,200), breaks = seq(200,2000,200)) +
  scale_fill_manual(values = c("#56B4E9", "#009E73")) +
  scale_y_continuous(name ="Elevation (m a.s.l)",
                  limits = c(200,2000),
                  breaks = seq(200,2000,200)) +
  scale_x_discrete(position = "top", name ="Weevil genus", 
                   labels = parse(text = levels(elev_dist$GenusItalics))) +
  theme_classic() + 
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                          axis.text.x = element_text(angle=70, size=10, hjust=0),
                          axis.text.y = element_text(size=10, colour="black"),
                          axis.title = element_text(size=11,face="bold",colour="black"),
        legend.position = c(0.1,0.9), 
        legend.text = element_text(size=11),
        legend.title = element_text(size=11, face = "bold"))

ggsave("./Figures/Figure_2.tiff", Elevation_distribution, 
       width=6000, height=4500, units="px", dpi=600, compression="lzw")

#_________________________________________________________
#### Multiple beta-diversity distributions (figure 3) ####
#_________________________________________________________

# - modified from Wayman et al. 2022 -

# Load community abundance matrix
weevil_abundance_matrix <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                               sheet = "Community_abundance_matrix", col_names = T)

# Set the first column (elevational floors) as row names and the
# elevational floors (200:2000 masl) as column names
weevil_abundance_matrix = as.data.frame(weevil_abundance_matrix)
am <- weevil_abundance_matrix[,-1]
rownames(am) <- weevil_abundance_matrix[,1]
names(am) <- seq(200,2000,200)

# Transform abundance-based matrix to incidence-based matrix
im = am
im[im > 0] <- 1

#Load transect 1 abundance matrix (sampled in 2019)
am <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                        sheet = "abundance_matrix_T1", col_names = T)
am2019 = am[,-1]
rownames(am2019) = am[,1]
im2019 = am2019
im2019[im2019 > 0] <- 1

#Load transect 2 abundance matrix (sampled in 2020)
am <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                        sheet = "abundance_matrix_T2", col_names = T)
am2020 = am[,-1]
rownames(am2020) = am[,1]
im2020 = am2020
im2020[im2020 > 0] <- 1

#Input databases requiered: presence-absence (or abundance-based) matrices
#were columns are species (or any studied taxa level) and rows are sites.
#In our work columns are leaf litter weevil morphospecies (55 columns)
#and rows are the elevational floors sampled (10 rows)
#Abundance-based matrices
View(t(am)) #both transects
View(t(am2019)) #transect 1
View(t(am2020)) #transect 2
#Incidence-based matrices
View(t(im)) #both transects
View(t(im2019)) #transect 1
View(t(im2020)) #transect 2

# BRAY-CURTIS index distributions #
# (using abundance-based matrix)

#Sample the distribution (1000) using a subset of sites (1/2 = 5 elevation floors)
samp_19 <- betapart::beta.sample.abund(t(am2019),index.family = "bray",
                                       samples = 1000,sites = round(nrow(t(am2019))/2))
samp_20 <- betapart::beta.sample.abund(t(am2020),index.family = "bray",
                                       samples = 1000,sites = round(nrow(t(am2020))/2))

#Bray's balanced variation = Jaccard turnover = Sorensen sim
bal_19 = samp_19[[1]][1]
bal_20 = samp_20[[1]][1]

#Bray's abundance gradient = Jaccard nestedness = Sorensen sne
gra_19 = samp_19[[1]][2]
gra_20 = samp_20[[1]][2]

#Bray's total dissimilarity = Jaccard total dissimilarity = Sorensen total dissimilarity
bray_19 = samp_19[[1]][3]
bray_20 = samp_20[[1]][3]

#The below code "mded" measures the difference between two independent or 
#non-independent empirical distributions and returns a significance level of 
#difference on the basis of the methods proposed by Poe et al. (1997, 2005). 
out_bal <- mded::mded(bal_19$beta.BRAY.BAL, bal_20$beta.BRAY.BAL, independent = T)
out_gra <- mded::mded(gra_19$beta.BRAY.GRA, gra_20$beta.BRAY.GRA, independent = T)
out_bray<- mded::mded(bray_19$beta.BRAY, bray_20$beta.BRAY, independent = T)

#Get p-value
p_bal <- round(out_bal$stat, 3);p_bal
p_gra <- round(out_gra$stat, 3);p_gra
p_bray <- round(out_bray$stat, 3);p_bray

#Add transect number as a column
bal_19$Transect = 1
bal_20$Transect = 2

gra_19$Transect = 1
gra_20$Transect = 2

bray_19$Transect = 1
bray_20$Transect = 2

#Get the balanced variation data ready for plotting
mu_bal <- data.frame(mean = numeric(0), Transect = character(0))
mu_bal[1,1] <- mean(bal_19$beta.BRAY.BAL)
mu_bal[1,2] <- 1
mu_bal[2,1] <- mean(bal_20$beta.BRAY.BAL)
mu_bal[2,2] <- 2
bal_df <- rbind(bal_19, bal_20)
bal_df$Transect <- as.factor(bal_df$Transect)

#Get the abundance gradient data ready for plotting
mu_gra <- data.frame(mean = numeric(0), Transect = character(0))
mu_gra[1,1] <- mean(gra_19$beta.BRAY.GRA)
mu_gra[1,2] <- 1
mu_gra[2,1] <- mean(gra_20$beta.BRAY.GRA)
mu_gra[2,2] <- 2
gra_df <- rbind(gra_19, gra_20)
gra_df$Transect <- as.factor(gra_df$Transect)

#Get the bray dissimilarity data ready for plotting
mu_bray <- data.frame(mean = numeric(0), Transect = character(0))
mu_bray[1,1] <- mean(bray_19$beta.BRAY)
mu_bray[1,2] <- 1
mu_bray[2,1] <- mean(bray_20$beta.BRAY)
mu_bray[2,2] <- 2
bray_df <- rbind(bray_19, bray_20)
bray_df$Transect <- as.factor(bray_df$Transect)

#Plot multiple beta-diversity distributions

#Bray total dissimilarity
plot_bray <- ggplot(bray_df, aes(x = beta.BRAY, group = Transect, fill = Transect)) + 
  geom_density(alpha = 0.4) +
  geom_vline(data=mu_bray, aes(xintercept = mean),
             linetype="dashed") + 
  scale_fill_manual(#values = c("#56B4E9","#D55E00")
                    values = c("#56B4E9","#D55E00")) + 
  labs(x = expression("MβD"[BRAY]), y = "Density") + 
  theme_classic() +
  theme(axis.title.x = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.margin = margin(10,0,10,0)) +   # top, right, bottom, left 
  annotate("text", y = max(density(bray_df$beta.BRAY)$y), 
           x = max(bray_df$beta.BRAY), 
           label = bquote(italic(p)==.(p_bray)), 
           hjust = 0, size = 5)

#Balanced variation
plot_bal <- ggplot(bal_df, aes(x = beta.BRAY.BAL, group = Transect, fill = Transect)) + 
  geom_density(alpha = 0.4) +
  geom_vline(data=mu_bal, aes(xintercept = mean),
             linetype="dashed") + 
  scale_fill_manual(#values = c("#56B4E9","#D55E00")
                    values = c("#56B4E9","#D55E00")) + 
  labs(x = expression("MβD"[BRAY.BAL]), y = "") + 
  theme_classic() +
  theme(axis.title.x = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.margin = margin(10,0,10,0)) +   # top, right, bottom, left 
  annotate("text", y = 4, 
           x = 0.3, 
           label = bquote(italic(p)==.(p_bal)), 
           hjust = 0, size = 5)

#Abundance gradient
plot_gra <- ggplot(gra_df, aes(x = beta.BRAY.GRA, group = Transect, fill = Transect)) + 
  geom_density(alpha = 0.4) +
  geom_vline(data=mu_gra, aes(xintercept = mean),
             linetype="dashed") + 
  scale_fill_manual(#values = c("#56B4E9","#D55E00")
                    values = c("#56B4E9","#D55E00")) + 
  labs(x = expression("MβD"[BRAY.GRA]), y = "") + 
  theme_classic() +
  theme(axis.title.x = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.margin = margin(10,0,10,0)) +   # top, right, bottom, left 
  annotate("text", y = max(density(gra_df$beta.BRAY.GRA)$y), 
           x = max(gra_df$beta.BRAY.GRA), 
           label = bquote(italic(p)==.(p_gra)), 
           hjust = 1, size = 5)


# JACCARD index distributions #

#Sample the distribution (1000) using a subset of sites (1/2 = 5 elevation floors)
samp_19 <- betapart::beta.sample(t(im2019),index.family = "jaccard",
                                       samples = 1000,sites = round(nrow(t(am2019))/2))
samp_20 <- betapart::beta.sample(t(im2020),index.family = "jaccard",
                                       samples = 1000,sites = round(nrow(t(am2020))/2))

#Bray's balanced variation = Jaccard turnover = Sorensen sim
turn_19 = samp_19[[1]][1]
turn_20 = samp_20[[1]][1]

#Bray's abundance gradient = Jaccard nestedness = Sorensen sne
nes_19 = samp_19[[1]][2]
nes_20 = samp_20[[1]][2]

#Bray's total dissimilarity = Jaccard total dissimilarity = Sorensen total dissimilarity
jac_19 = samp_19[[1]][3]
jac_20 = samp_20[[1]][3]

#The below code "mded" measures the difference between two independent or 
#non-independent empirical distributions and returns a significance level of 
#difference on the basis of the methods proposed by Poe et al. (1997, 2005). 
out_turn <- mded::mded(turn_19$beta.JTU, turn_20$beta.JTU, independent = T)
out_nes <- mded::mded(nes_19$beta.JNE, nes_20$beta.JNE, independent = T)
out_jac<- mded::mded(jac_19$beta.JAC, jac_20$beta.JAC, independent = T)

#Get p-value
p_turn <- round(out_turn$stat, 3);p_turn
p_nes <- round(out_nes$stat, 3);p_nes
p_jac <- round(out_jac$stat, 3);p_jac

#Add transect number as a column
turn_19$Transect = 1
turn_20$Transect = 2

nes_19$Transect = 1
nes_20$Transect = 2

jac_19$Transect = 1
jac_20$Transect = 2

#Get the turnover data ready for plotting
mu_turn <- data.frame(mean = numeric(0), Transect = character(0))
mu_turn[1,1] <- mean(turn_19$beta.JTU)
mu_turn[1,2] <- 1
mu_turn[2,1] <- mean(turn_20$beta.JTU)
mu_turn[2,2] <- 2
turn_df <- rbind(turn_19, turn_20)
turn_df$Transect <- as.factor(turn_df$Transect)

#Get the nestedness data ready for plotting
mu_nes <- data.frame(mean = numeric(0), Transect = character(0))
mu_nes[1,1] <- mean(nes_19$beta.JNE)
mu_nes[1,2] <- 1
mu_nes[2,1] <- mean(nes_20$beta.JNE)
mu_nes[2,2] <- 2
nes_df <- rbind(nes_19, nes_20)
nes_df$Transect <- as.factor(nes_df$Transect)

#Get the jaccard dissimilarity data ready for plotting
mu_jac <- data.frame(mean = numeric(0), Transect = character(0))
mu_jac[1,1] <- mean(jac_19$beta.JAC)
mu_jac[1,2] <- 1
mu_jac[2,1] <- mean(jac_20$beta.JAC)
mu_jac[2,2] <- 2
jac_df <- rbind(jac_19, jac_20)
jac_df$Transect <- as.factor(jac_df$Transect)

#Plot multiple beta-diversity distributions

#Jaccard total dissimilarity
plot_jac <- ggplot(jac_df, aes(x = beta.JAC, group = Transect, fill = Transect)) + 
  geom_density(alpha = 0.4) +
  geom_vline(data=mu_jac, aes(xintercept = mean),
             linetype="dashed") + 
  scale_fill_manual(#values = c("#009E73","#CC79A7")
                    values = c("#56B4E9","#D55E00")) + 
  labs(x = expression("MβD"[JAC]), y = "Density") + 
  theme_classic() +
  theme(axis.title.x = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.margin = margin(10,0,10,0)) +   # top, right, bottom, left 
  annotate("text", y = max(density(jac_df$beta.JAC)$y), 
           x = 0.7, 
           label = bquote(italic(p)==.(p_jac)), 
           hjust = 0, size = 5)

#Jaccard turnover
plot_turn <- ggplot(turn_df, aes(x = beta.JTU, group = Transect, fill = Transect)) + 
  geom_density(alpha = 0.4) +
  geom_vline(data=mu_turn, aes(xintercept = mean),
             linetype="dashed") + 
  scale_fill_manual(#values = c("#009E73","#CC79A7")
                    values = c("#56B4E9","#D55E00")) + 
  labs(x = expression("MβD"[JAC.TURN]), y = "") + 
  theme_classic() +
  theme(axis.title.x = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.margin = margin(10,0,10,0)) +   # top, right, bottom, left 
  annotate("text", y = 9, 
           x = 0.4, 
           label = bquote(italic(p)==.(p_turn)), 
           hjust = 0, size = 5)

#Jaccard nestedness
plot_nes <- ggplot(nes_df, aes(x = beta.JNE, group = Transect, fill = Transect)) + 
  geom_density(alpha = 0.4) +
  geom_vline(data=mu_nes, aes(xintercept = mean),
             linetype="dashed") + 
  scale_fill_manual(#values = c("#009E73","#CC79A7")
                    values = c("#56B4E9","#D55E00")) + 
  labs(x = expression("MβD"[JAC.NES]), y = "") + 
  theme_classic() +
  theme(axis.title.x = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.margin = margin(10,0,10,0)) +   # top, right, bottom, left 
  annotate("text", y = max(density(nes_df$beta.JNE)$y), 
           x = max(nes_df$beta.JNE), 
           label = bquote(italic(p)==.(p_nes)), 
           hjust = 1, size = 5)

# SORENSEN index distributions #
# (using incidence-based matrix)

#Sample the distribution (1000) using a subset of sites (1/2 = 5 elevation floors)
samp_19 <- betapart::beta.sample(t(im2019),index.family = "sorensen",
                                       samples = 1000,sites = round(nrow(t(am2019))/2))
samp_20 <- betapart::beta.sample(t(im2020),index.family = "sorensen",
                                       samples = 1000,sites = round(nrow(t(am2020))/2))

#Bray's balanced variation = Jaccard turnover = Sorensen sim
sim_19 = samp_19[[1]][1]
sim_20 = samp_20[[1]][1]

#Bray's abundance gradient = Jaccard nestedness = Sorensen sne
sne_19 = samp_19[[1]][2]
sne_20 = samp_20[[1]][2]

#Bray's total dissimilarity = Jaccard total dissimilarity = Sorensen total dissimilarity
sor_19 = samp_19[[1]][3]
sor_20 = samp_20[[1]][3]

#The below code "mded" measures the difference between two independent or 
#non-independent empirical distributions and returns a significance level of 
#difference on the basis of the methods proposed by Poe et al. (1997, 2005). 
out_sim <- mded::mded(sim_19$beta.SIM, sim_20$beta.SIM, independent = T)
out_sne <- mded::mded(sne_19$beta.SNE, sne_20$beta.SNE, independent = T)
out_sor<- mded::mded(sor_19$beta.SOR, sor_20$beta.SOR, independent = T)

#Get p-value
p_sim <- round(out_sim$stat, 3);p_sim
p_sne <- round(out_sne$stat, 3);p_sne
p_sor <- round(out_sor$stat, 3);p_sor

#Add transect number as a column
sim_19$Transect = 1
sim_20$Transect = 2

sne_19$Transect = 1
sne_20$Transect = 2

sor_19$Transect = 1
sor_20$Transect = 2

#Get the sorensen turnover data ready for plotting
mu_sim <- data.frame(mean = numeric(0), Transect = character(0))
mu_sim[1,1] <- mean(sim_19$beta.SIM)
mu_sim[1,2] <- 1
mu_sim[2,1] <- mean(sim_20$beta.SIM)
mu_sim[2,2] <- 2
sim_df <- rbind(sim_19, sim_20)
sim_df$Transect <- as.factor(sim_df$Transect)

#Get the sorensen nestedness data ready for plotting
mu_sne <- data.frame(mean = numeric(0), Transect = character(0))
mu_sne[1,1] <- mean(sne_19$beta.SNE)
mu_sne[1,2] <- 1
mu_sne[2,1] <- mean(sne_20$beta.SNE)
mu_sne[2,2] <- 2
sne_df <- rbind(sne_19, sne_20)
sne_df$Transect <- as.factor(sne_df$Transect)

#Get the sorensen dissimilarity data ready for plotting
mu_sor <- data.frame(mean = numeric(0), Transect = character(0))
mu_sor[1,1] <- mean(sor_19$beta.SOR)
mu_sor[1,2] <- 1
mu_sor[2,1] <- mean(sor_20$beta.SOR)
mu_sor[2,2] <- 2
sor_df <- rbind(sor_19, sor_20)
sor_df$Transect <- as.factor(sor_df$Transect)

#Plot multiple beta-diversity distributions

#Sorensen total dissimilarity
plot_sor <- ggplot(sor_df, aes(x = beta.SOR, group = Transect, fill = Transect)) + 
  geom_density(alpha = 0.4) +
  geom_vline(data=mu_sor, aes(xintercept = mean),
             linetype="dashed") + 
  scale_fill_manual(#values = c("#661100","#44AA99")
                    values = c("#56B4E9","#D55E00")) + 
  labs(x = expression("MβD"[SOR]), y = "Density") + 
  theme_classic() +
  theme(axis.title.x = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.margin = margin(10,0,10,0)) +   # top, right, bottom, left 
  annotate("text", y = 8.5, 
           x = 0.55, 
           label = bquote(italic(p)==.(p_sor)), 
           hjust = 0, size = 5)

#Sorensen turnover
plot_sim <- ggplot(sim_df, aes(x = beta.SIM, group = Transect, fill = Transect)) + 
  geom_density(alpha = 0.4) +
  geom_vline(data=mu_sim, aes(xintercept = mean),
             linetype="dashed") + 
  scale_fill_manual(#values = c("#661100","#44AA99")
                    values = c("#56B4E9","#D55E00")) + 
  labs(x = expression("MβD"[SOR.SIM]), y = "") + 
  theme_classic() +
  theme(axis.title.x = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.margin = margin(10,0,10,0)) +   # top, right, bottom, left 
  annotate("text", y = 7, 
           x = 0.2, 
           label = bquote(italic(p)==.(p_sim)), 
           hjust = 0, size = 5)

#Sorensen nestedness
plot_sne <- ggplot(sne_df, aes(x = beta.SNE, group = Transect, fill = Transect)) + 
  geom_density(alpha = 0.4) +
  geom_vline(data=mu_sne, aes(xintercept = mean),
             linetype="dashed") + 
  scale_fill_manual(#values = c("#661100","#44AA99")
                    values = c("#56B4E9","#D55E00")) + 
  labs(x = expression("MβD"[SOR.SNE]), y = "") + 
  theme_classic() +
  theme(axis.title.x = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.margin = margin(10,0,10,0)) +   # top, right, bottom, left 
  annotate("text", y = max(density(sne_df$beta.SNE)$y), 
           x = max(sne_df$beta.SNE), 
           label = bquote(italic(p)==.(p_sne)), 
           hjust = 1, size = 5)

#Save all 9 multiple-beta-diversity plots as one figure
MBD_dist_CURBC = 
  plot_jac + plot_turn + plot_nes +
  plot_sor + plot_sim + plot_sne +
  plot_bray + plot_bal + plot_gra +
  patchwork::plot_layout(ncol = 3, nrow = 3, guides = 'collect')
ggsave(MBD_dist_CURBC, filename = "./Figures_3.tiff",
       height = 6300, width = 6300, units = "px", dpi = 600, compression = "lzw")

#___________________________________________
#### Beta diversity heatmaps (figure 4) ####
#___________________________________________

# Load community abundance matrix
weevil_abundance_matrix <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                               sheet = "Community_abundance_matrix", col_names = T)

# Set the first column (elevational floors) as row names and the
# elevational floors (200:2000 masl) as column names
weevil_abundance_matrix = as.data.frame(weevil_abundance_matrix)
am <- weevil_abundance_matrix[,-1]
rownames(am) <- weevil_abundance_matrix[,1]
names(am) <- seq(200,2000,200)

# Transform abundance-based matrix to incidence-based matrix
im = am
im[im > 0] <- 1

#Multiple site measures (Baselga & Orme, 2012) incidence-based
#Species (columns) and sites (rows)
core = betapart::betapart.core(t(im)) #to get data description
multi = betapart::beta.multi(core, index.family = "sorensen")
multi
# $beta.SIM (turnover)
# 0.6640827
# $beta.SNE (nestedness)
# 0.127584
# $beta.SOR (total dissimilarity)
# 0.7916667

#Abundance-based multiple sites dissimilarity (Baselga, 2017)
#Species (columns) and sites (rows)
core.ab = betapart::betapart.core.abund(t(am))
multi.ab = betapart::beta.multi.abund(core.ab, index.family = "bray")
multi.ab
# $beta.BRAY.BAL
# [1] 0.5974458
# 
# $beta.BRAY.GRA
# [1] 0.2509932
# 
# $beta.BRAY
# [1] 0.8484389

# Abundance-based assemblage dissimilarity can be partitioned into 
# components accounting for (i) balanced variation in abundance, 
# whereby the individuals of some species in one site are substituted
# by the same number of individuals of different species in another site,
# and (ii) abundance gradients, whereby some individuals are lost 
# from one site to the other.

#β-diversity patterns within each elevational site:
#total dissimilarity among all points with weevil 
#morphospecies recorded.
#Sampling across equal sites (Baselga & Orme, 2012)
samp.sor = betapart::beta.sample(core, sites = 9, samples = 999, index.family = "sorensen")
samp.sor
# $mean.values
#  beta.SIM  beta.SNE  beta.SOR 
# 0.5874736 0.1940410 0.7815146 
# 
# $sd.values
#   beta.SIM   beta.SNE   beta.SOR 
# 0.11901078 0.08776046 0.04172784 

# JACCARD
samp.jac = betapart::beta.sample(core, sites = 9, samples = 999, index.family = "jaccard")
samp.jac
# $mean.values
#  beta.JTU  beta.JNE  beta.JAC 
# 0.7281220 0.1473452 0.8754671 
# 
# $sd.values
#   beta.JTU   beta.JNE   beta.JAC 
# 0.10528481 0.08471461 0.02714653 

# BRAY-CURTIS
samp.bray = betapart::beta.sample.abund(am, sites = 9, samples = 999, index.family = "bray")
samp.bray
# $mean.values
# beta.BRAY.BAL beta.BRAY.GRA     beta.BRAY 
#     0.5427082     0.3859278     0.9286360 
# 
# $sd.values
# beta.BRAY.BAL beta.BRAY.GRA     beta.BRAY 
#    0.17549712    0.17486360    0.03057849  

# Concatenate into a single object all nine beta diversity values
beta.values = round(c(samp.jac$mean.values,samp.sor$mean.values,samp.bray$mean.values),2)

# - modified from Perez-Rosales et al. 2022 -

#Measure pairwise dissimilarity with betapart package
#Parwise dissimilarity analysed with jaccard index
pair_jac = betapart::beta.pair(t(im), index.family = "jaccard")
#Parwise dissimilarity analysed with sorensen index
pair_sor = betapart::beta.pair(t(im), index.family = "sorensen")
#Parwise dissimilarity analysed with bray-curtis index
pair_bray = betapart::beta.pair.abund(t(am), index.family = "bray")

#Create total dissimilarity, turnover and nestedness matrices
#Total dissimilarity: beta.jac = beta.sor = beta.bray
#Turnover: beta.jtu = beta.sim = beta.bray.bal
#Nestedness: beta.jne = beta.sne = beta.bray.gra
total_jac = pair_jac$beta.jac; total_jac_mat = as.matrix(total_jac)
total_sor = pair_sor$beta.sor; total_sor_mat = as.matrix(total_sor)
total_bray = pair_bray$beta.bray; total_bray_mat = as.matrix(total_bray)

turn_jac = pair_jac$beta.jtu; turn_jac_mat = as.matrix(turn_jac)
turn_sor = pair_sor$beta.sim; turn_sor_mat = as.matrix(turn_sor)
turn_bray = pair_bray$beta.bray.bal; turn_bray_mat = as.matrix(turn_bray)

nes_jac = pair_jac$beta.jne; nes_jac_mat = as.matrix(nes_jac)
nes_sor = pair_sor$beta.sne; nes_sor_mat = as.matrix(nes_sor)
nes_bray = pair_bray$beta.bray.gra; nes_bray_mat = as.matrix(nes_bray)

#Eliminate duplicate lower triangle matrix values
total_jac_mat[upper.tri(total_jac_mat)] <- NA
total_sor_mat[upper.tri(total_sor_mat)] <- NA
total_bray_mat[upper.tri(total_bray_mat)] <- NA

turn_jac_mat[upper.tri(turn_jac_mat)] <- NA
turn_sor_mat[upper.tri(turn_sor_mat)] <- NA
turn_bray_mat[upper.tri(turn_bray_mat)] <- NA

nes_jac_mat[upper.tri(nes_jac_mat)] <- NA
nes_sor_mat[upper.tri(nes_sor_mat)] <- NA
nes_bray_mat[upper.tri(nes_bray_mat)] <- NA

#Reshape beta matrices into three column dataframes, with colnames:
#elev1, elev2, beta-diversity values
#Total dissimilarity dataframes
total_jac_mat = reshape2::melt(total_jac_mat, 
                               measure.vars="nobserv", 
                               na.rm=F); colnames(total_jac_mat) = c(
                                 "elev1","elev2","total_jac")
total_sor_mat = reshape2::melt(total_sor_mat, 
                               measure.vars="nobserv", 
                               na.rm=F); colnames(total_sor_mat) = c(
                                 "elev1","elev2","total_sor")
total_bray_mat = reshape2::melt(total_bray_mat, 
                               measure.vars="nobserv", 
                               na.rm=F); colnames(total_bray_mat) = c(
                                 "elev1","elev2","total_bray")
#Turnover dataframes
turn_jac_mat = reshape2::melt(turn_jac_mat, 
                               measure.vars="nobserv", 
                               na.rm=F); colnames(turn_jac_mat) = c(
                                 "elev1","elev2","turn_jac")
turn_sor_mat = reshape2::melt(turn_sor_mat, 
                               measure.vars="nobserv", 
                               na.rm=F); colnames(turn_sor_mat) = c(
                                 "elev1","elev2","turn_sor")
turn_bray_mat = reshape2::melt(turn_bray_mat, 
                               measure.vars="nobserv", 
                               na.rm=F); colnames(turn_bray_mat) = c(
                                 "elev1","elev2","turn_bray")
#Nestedness dataframes
nes_jac_mat = reshape2::melt(nes_jac_mat, 
                               measure.vars="nobserv", 
                               na.rm=F); colnames(nes_jac_mat) = c(
                                 "elev1","elev2","nes_jac")
nes_sor_mat = reshape2::melt(nes_sor_mat, 
                               measure.vars="nobserv", 
                               na.rm=F); colnames(nes_sor_mat) = c(
                                 "elev1","elev2","nes_sor")
nes_bray_mat = reshape2::melt(nes_bray_mat, 
                               measure.vars="nobserv", 
                               na.rm=F); colnames(nes_bray_mat) = c(
                                 "elev1","elev2","nes_bray")
#Eliminate NAs
total_jac_mat = na.omit(total_jac_mat)
total_sor_mat = na.omit(total_sor_mat)
total_bray_mat = na.omit(total_bray_mat)

turn_jac_mat = na.omit(turn_jac_mat)
turn_sor_mat = na.omit(turn_sor_mat)
turn_bray_mat = na.omit(turn_bray_mat)

nes_jac_mat = na.omit(nes_jac_mat)
nes_sor_mat = na.omit(nes_sor_mat)
nes_bray_mat = na.omit(nes_bray_mat)

#Eliminate all beta-diversity zero-values
total_jac_mat = total_jac_mat[!(total_jac_mat$total_jac==0),]
total_sor_mat = total_sor_mat[!(total_sor_mat$total_sor==0),]
total_bray_mat = total_bray_mat[!(total_bray_mat$total_bray==0),]
 
turn_jac_mat = turn_jac_mat[!(turn_jac_mat$turn_jac==0),]
turn_sor_mat = turn_sor_mat[!(turn_sor_mat$turn_sor==0),]
turn_bray_mat = turn_bray_mat[!(turn_bray_mat$turn_bray==0),]

nes_jac_mat = nes_jac_mat[!(nes_jac_mat$nes_jac==0),]
nes_sor_mat = nes_sor_mat[!(nes_sor_mat$nes_sor==0),]
nes_bray_mat = nes_bray_mat[!(nes_bray_mat$nes_bray==0),]

#Merge total dissimilarity, turnover and nestedness dataframes
#according to their respective index
jac_mat = merge(total_jac_mat, turn_jac_mat, by=c("elev1","elev2"), all.x=T, all.y=F)
jac_mat = merge(jac_mat, nes_jac_mat, by=c("elev1","elev2"), all.x=T, all.y=F)

sor_mat = merge(total_sor_mat, turn_sor_mat, by=c("elev1","elev2"), all.x=T, all.y=F)
sor_mat = merge(sor_mat, nes_sor_mat, by=c("elev1","elev2"), all.x=T, all.y=F)

bray_mat = merge(total_bray_mat, turn_bray_mat, by=c("elev1","elev2"), all.x=T, all.y=F)
bray_mat = merge(bray_mat, nes_bray_mat, by=c("elev1","elev2"), all.x=T, all.y=F)

weevil_beta_pair_mat = merge(jac_mat, sor_mat, by=c("elev1","elev2"), all=T)
weevil_beta_pair_mat = merge(weevil_beta_pair_mat, bray_mat, by=c("elev1","elev2"), all=T)

#Set elevation columns as factors
weevil_beta_pair_mat$elev1 = factor(weevil_beta_pair_mat$elev1, 
                                    levels=c("200", "400", "600", "800", "1000", 
                                             "1200", "1400", "1600", "1800", "2000"))
weevil_beta_pair_mat$elev2 = factor(weevil_beta_pair_mat$elev2, 
                                    levels=c("200", "400", "600", "800", "1000", 
                                             "1200", "1400", "1600", "1800", "2000"))

#JACCARD DISSIMILARITY
total_jac_plot <- 
  ggplot(weevil_beta_pair_mat, aes(elev1, elev2, fill=total_jac))+
  geom_tile(color="white")+
  scale_fill_gradient2(#low="lemonchiffon", mid="orange", high="red",
                       low="#ffffb2", mid="#fecc5c", high="#bd0026",  
                       midpoint = 0.5, limit=c(0,1), 
                       breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), 
                       space="Lab", name="Dissimilarity")+
  scale_x_discrete(position="top")+
  scale_y_discrete(limits=rev)+
  theme_bw()+
  coord_fixed()+
  labs(x="Total dissimilarity",y="Elevation (m a.s.l)", title = "") +
  annotate(geom="text", x=2, y=2, size=6, label=expression("β"[JAC]))+
  annotate(geom="text", x=4, y=2, size=6, label=paste0("= ",beta.values["beta.JAC"][[1]]))+
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=11, face="bold", colour="black"), 
        strip.text = element_text(size=10, colour="black"),
        axis.text.x = element_text(angle=20, size=8, hjust = 0),
        axis.text.y = element_text(size=8),
        legend.position = "none")

#JACCARD TURNOVER
turn_jac_plot <- 
  ggplot(weevil_beta_pair_mat, aes(elev1, elev2, fill=turn_jac))+
  geom_tile(color="white")+
  scale_fill_gradient2(#low="lemonchiffon", mid="orange", high="red",
                       low="#ffffb2", mid="#fecc5c", high="#bd0026",  
                       midpoint = 0.5, limit=c(0,1), 
                       breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), 
                       space="Lab", name="Turnover")+
  scale_x_discrete(position="top")+
  scale_y_discrete(limits=rev)+
  theme_bw()+
  coord_fixed()+
  labs(x="Turnover",y="") +
  annotate(geom="text", x=2, y=2, size=6, label=expression("β"[JTU]))+
  annotate(geom="text", x=4, y=2, size=6, label=paste0("= ",beta.values["beta.JTU"][[1]]))+
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=11, face="bold", colour="black"), 
        strip.text = element_text(size=10, colour="black"),
        axis.text.x = element_text(angle=20, size=8, hjust = 0),
        axis.text.y = element_text(size=8),
        legend.position = "none")

#JACCARD NESTEDNESS
nes_jac_plot <- 
  ggplot(weevil_beta_pair_mat, aes(elev1, elev2, fill=nes_jac))+
  geom_tile(color="white")+
  scale_fill_gradient2(#low="lemonchiffon", mid="orange", high="red",
                       low="#ffffb2", mid="#fecc5c", high="#bd0026",  
                       midpoint = 0.5, limit=c(0,1), 
                       breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), 
                       space="Lab", name="JACCARD")+
  scale_x_discrete(position="top")+
  scale_y_discrete(limits=rev)+
  theme_bw()+
  coord_fixed()+
  labs(x="Nestedness",y="") +
  annotate(geom="text", x=2, y=2, size=6, label=expression("β"[JNE]))+
  annotate(geom="text", x=4, y=2, size=6, label=paste0("= ",beta.values["beta.JNE"][[1]]))+
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=11, face="bold", colour="black"), 
        strip.text = element_text(size=10, colour="black"),
        axis.text.x = element_text(angle=20, size=8, hjust = 0),
        axis.text.y = element_text(size=8))

#SORENSEN DISSIMILARITY
total_sor_plot <- 
  ggplot(weevil_beta_pair_mat, aes(elev1, elev2, fill=total_sor))+
  geom_tile(color="white")+
  scale_fill_gradient2(#low="lemonchiffon", mid="orange", high="red",
                       low="#ffffb2", mid="#fecc5c", high="#bd0026",  
                       midpoint = 0.5, limit=c(0,1), 
                       breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), 
                       space="Lab", name="Dissimilarity")+
  scale_x_discrete(position="top")+
  scale_y_discrete(limits=rev)+
  theme_bw()+
  coord_fixed()+
  labs(x="",y="Elevation (m a.s.l)") +
  annotate(geom="text", x=2, y=2, size=6, label=expression("β"[SØR]))+
  annotate(geom="text", x=4, y=2, size=6, label=paste0("= ",beta.values["beta.SOR"][[1]]))+
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=11, face="bold", colour="black"), 
        strip.text = element_text(size=10, colour="black"),
        axis.text.x = element_text(angle=20, size=8, hjust = 0),
        axis.text.y = element_text(size=8),
        legend.position = "none")

#SORENSEN TURNOVER
turn_sor_plot <- 
  ggplot(weevil_beta_pair_mat, aes(elev1, elev2, fill=turn_sor))+
  geom_tile(color="white")+
  scale_fill_gradient2(#low="lemonchiffon", mid="orange", high="red",
                       low="#ffffb2", mid="#fecc5c", high="#bd0026",  
                       midpoint = 0.5, limit=c(0,1), 
                       breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), 
                       space="Lab", name="Turnover")+
  scale_x_discrete(position="top")+
  scale_y_discrete(limits=rev)+
  theme_bw()+
  coord_fixed()+
  labs(#x=expression("Turnover"[(SOR.SIM)]),
       x="",
       y="") +
  annotate(geom="text", x=2, y=2, size=6, label=expression("β"[SIM]))+
  annotate(geom="text", x=4, y=2, size=6, label=paste0("= ",beta.values["beta.SIM"][[1]]))+
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=11, colour="black"), 
        strip.text = element_text(size=10, colour="black"),
        axis.text.x = element_text(angle=20, size=8, hjust = 0),
        axis.text.y = element_text(size=8),
        legend.position = "none")

#SORENSEN NESTEDNESS
nes_sor_plot <- 
  ggplot(weevil_beta_pair_mat, aes(elev1, elev2, fill=nes_sor))+
  geom_tile(color="white")+
  scale_fill_gradient2(#low="lemonchiffon", mid="orange", high="red",
                       low="#ffffb2", mid="#fecc5c", high="#bd0026",  
                       midpoint = 0.5, limit=c(0,1), 
                       breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), 
                       space="Lab", name="SORENSEN")+
  scale_x_discrete(position="top")+
  scale_y_discrete(limits=rev)+
  theme_bw()+
  coord_fixed()+
  labs(#x=expression("Nestedness"[(SOR.SNE)]),
       x="",
       y="") +
  annotate(geom="text", x=2, y=2, size=6, label=expression("β"[SNE]))+
  annotate(geom="text", x=4, y=2, size=6, label=paste0("= ",beta.values["beta.SNE"][[1]]))+
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=11, colour="black"), 
        strip.text = element_text(size=10, colour="black"),
        axis.text.x = element_text(angle=20, size=8, hjust = 0),
        axis.text.y = element_text(size=8))

#BRAY-CURTIS DISSIMILARITY
total_bray_plot <- 
  ggplot(weevil_beta_pair_mat, aes(elev1, elev2, fill=total_bray))+
  geom_tile(color="white")+
  scale_fill_gradient2(#low="lemonchiffon", mid="orange", high="red",
                       low="#ffffb2", mid="#fecc5c", high="#bd0026",  
                       midpoint = 0.5, limit=c(0,1), 
                       breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), 
                       space="Lab", name="Dissimilarity")+
  scale_x_discrete(position="top")+
  scale_y_discrete(limits=rev)+
  theme_bw()+
  coord_fixed()+
  labs(x="",y="Elevation (m a.s.l)") +
  annotate(geom="text", x=2, y=2, size=6, label=expression("β"[BRAY]))+
  annotate(geom="text", x=4.2, y=2, size=6, label=paste0("= ",beta.values["beta.BRAY"][[1]]))+
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=11, face="bold", colour="black"), 
        strip.text = element_text(size=10, colour="black"),
        axis.text.x = element_text(angle=20, size=8, hjust = 0),
        axis.text.y = element_text(size=8),
        legend.position = "none")

#BRAY-CURTIS TURNOVER
turn_bray_plot <- 
  ggplot(weevil_beta_pair_mat, aes(elev1, elev2, fill=turn_bray))+
  geom_tile(color="white")+
  scale_fill_gradient2(#low="lemonchiffon", mid="orange", high="red",
                       low="#ffffb2", mid="#fecc5c", high="#bd0026",  
                       midpoint = 0.5, limit=c(0,1), 
                       breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), 
                       space="Lab", name="Turnover")+
  scale_x_discrete(position="top")+
  scale_y_discrete(limits=rev)+
  theme_bw()+
  coord_fixed()+
  labs(#x=expression("Turnover"[(BRAY.BAL)]),
       x="",
       y="") +
  annotate(geom="text", x=2, y=2, size=6, label=expression("β"[BAL]))+
  annotate(geom="text", x=4, y=2, size=6, label=paste0("= ",beta.values["beta.BRAY.BAL"][[1]]))+
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=11, colour="black"), 
        strip.text = element_text(size=10, colour="black"),
        axis.text.x = element_text(angle=20, size=8, hjust = 0),
        axis.text.y = element_text(size=8),
        legend.position = "none")

#BRAY-CURTIS NESTEDNESS
nes_bray_plot <- 
  ggplot(weevil_beta_pair_mat, aes(elev1, elev2, fill=nes_bray))+
  geom_tile(color="white")+
  scale_fill_gradient2(#low="lemonchiffon", mid="orange", high="red",
                       low="#ffffb2", mid="#fecc5c", high="#bd0026",
                       midpoint = 0.5, limit=c(0,1), 
                       breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), 
                       space="Lab", name="BRAY-CURTIS")+
  scale_x_discrete(position="top")+
  scale_y_discrete(limits=rev)+
  theme_bw()+
  coord_fixed()+
  labs(#x=expression("Nestedness"[(BRAY.GRA)]),
       x="",
       y="") +
  annotate(geom="text", x=2, y=2, size=6, label=expression("β"[GRA]))+
  annotate(geom="text", x=4, y=2, size=6, label=paste0("= ",beta.values["beta.BRAY.GRA"][[1]]))+
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=11, colour="black"), 
        strip.text = element_text(size=10, colour="black"),
        axis.text.x = element_text(angle=20, size=8, hjust = 0),
        axis.text.y = element_text(size=8))

#Save all 9 heatmaps as a single figure
beta_pairwise_heatmaps_CURBC = 
  total_jac_plot + turn_jac_plot + nes_jac_plot +
  total_sor_plot + turn_sor_plot + nes_sor_plot +
  total_bray_plot + turn_bray_plot + nes_bray_plot +
  patchwork::plot_layout(ncol = 3, nrow = 3, guides = 'keep')
ggsave(beta_pairwise_heatmaps_CURBC, filename = "./Figure_4.tiff",
       height = 6000, width = 6800, units = "px", dpi = 600, compression = "lzw")

#_____________________________________________________
#### Canonical Correspondence Analysis (figure 5) ####
#_____________________________________________________

# Differences between RDA and CCA: RDA is based in PCA, and PCA assumes normal
# distributions for its variables, so to run an RDA you need to transform your
# species abundance data. CCA is based in CA that can handle untransformed species
# abundance data, nevertheless the problem with CCA is rare species: 
# "CCA use should be limited to situations where rare species are well sampled and
# are seen as potential indicators of particular characteristics of an ecosystem; 
# the alternative is to eliminate rare species from the data table before CCA."
# Brocard et al., 2018. Numerical Ecology with R.

CCA_response_matrix = readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                               sheet = "CCA_response_matrix", col_names = T)
CCA_response_matrix = as.data.frame(CCA_response_matrix)
responsedf <- CCA_response_matrix[,-1]
rownames(responsedf) <- CCA_response_matrix[,1]

explanatorydf = enviro_predictors <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                                       sheet = "Environmental_predictors", col_names = T)
explanatorydf = data.frame(explanatorydf)

# Extract response variables (abundances of each weevil morph per plot-transect) and
# remove rare species (singletons and doubletons)
spe = responsedf
spe.sin = responsedf[,colSums(responsedf)>1]
spe.sindou = responsedf[,colSums(responsedf)>2]

# Explanatory variables (with VIF<1.5 altitude, litter depth,
# bio15(precip. seasonality), dbh_median, tree_rich, pinus_prop)
env = explanatorydf[c("litter_depth",
                    "bio15",
                    "dbh_median",
                    "tree_rich",
                    "pinus_prop",
                    "quercus_prop",
                    "ternstroemia_prop",
                    #"tropical_veg_prop",
                    #"boreal_veg_prop",
                    "altitude")]
names(env) = c("LITTER", 
               "RAIN",
               "DBH",
               "TREE",
               "PINUS",
               "QUERCUS",
               "TERNSTROEMIA",
               #"TROPICAL_TREES",
               #"BOREAL_TREES",
               "ELEVATION")

# Run DCA (detrended constrained Analysis) to decide whether a linear (RDA)
# or unimodal (CCA) ordination method should be used. According to Leps & Smilauer, 2003
# the axis lengths < 3 show a linear trend (RDA), while axis lengths > 3 unimodal (CCA)
# Use vegan::decorana() function to run a DCA with the response matrix:
(vegan::decorana(spe)) # 4.5590 axis length
(vegan::decorana(vegan::decostand(spe, "hellinger"))) # 4.3857 axis length
(vegan::decorana(log1p(spe))) # 3.5446 axis length
# CCA is chosen

#Use raw species abundance data (Brocard et al., 2018)
spe.cca <- vegan::cca(spe.sindou ~ ., env)
summary(spe.cca) # Scaling 2 (default)

# As mentioned by Brocard et al. (2018), the R2 of a RDA or CCA is biased like
# the ordinary R2 of multiple regression, because, on the one hand, any
# variable included in an explanatory matrix X increases the R2, irrespective 
# of it being related, or not, to the response data. On the other hand, the 
# accumulation of explanatory variables inflates the apparent amount of explained
# variance because of random correlations. This problem can be cured by adjusting
# the R2 with vegan's function RsquareAdj()
# Unadjusted and adjusted R^2 - like statistics
# Unadjusted and adjusted R^2 - like statistics
vegan::RsquareAdj(spe.cca) #0.3188059

# Forwards selection and parsimonious CCA
# CCA-based forward selection using vegan's ordistep()
# This function allows the use of factors (categorical 
# explanatory variables)
cca.step.forward <-
  vegan::ordistep(vegan::cca(spe.sindou ~ 1, data = env),
           scope = formula(spe.cca),
           direction = "forward",
           permutations = how(nperm = 999))
cca.step.forward$anova
#             Df    AIC      F Pr(>F)    
# + ELEVATION  1 123.22 4.3255  0.001 ***
# + QUERCUS    1 121.39 3.5872  0.001 ***
# + DBH        1 121.01 2.0222  0.023 *  
RsquareAdj(cca.step.forward) #0.3062024

# Parsimonious CCA using ELEVATION + DBH + QUERCUS 
spe.cca.pars <- cca(spe.sindou ~ ELEVATION + DBH + QUERCUS, data = env)
anova(spe.cca.pars, permutations = how(nperm = 999))
#          Df ChiSquare      F Pr(>F)    
# Model     3    0.9872 3.6899  0.001 ***
# Residual 16    1.4269       
RsquareAdj(spe.cca.pars) #0.3055038

# Compare variance inflation factors
vif.cca(spe.cca)
# ELEVATION       LITTER         RAIN          DBH         TREE        PINUS      QUERCUS TERNSTROEMIA 
#  3.843103     3.557135     2.154848     1.371694     3.650465     3.266333     5.642100     3.198570 
vif.cca(spe.cca.pars)
# ELEVATION       DBH   QUERCUS 
#  1.709390  1.189045  1.901465 

# Partitioning R2-adj for each of the CCA factors:
full <- cca(spe.sindou ~ ELEVATION + DBH + QUERCUS, data = env)
no.elev <- cca(spe.sindou ~ DBH + QUERCUS, data = env)
no.dbh <- cca(spe.sindou ~ ELEVATION + QUERCUS, data = env)
no.quercus <- cca(spe.sindou ~ ELEVATION + DBH, data = env)
round(RsquareAdj(full)[[2]] - RsquareAdj(no.elev)[[2]],3) #Elevation R2: 0.117
round(RsquareAdj(full)[[2]] - RsquareAdj(no.dbh)[[2]],3) #DBH R2: 0.042
round(RsquareAdj(full)[[2]] - RsquareAdj(no.quercus)[[2]],3) #Quercus R2: 0.098

#Full CCA is "spe.cca", most parsimonious CCA is "spe.cca.pars"
#(the one with only elevation, dbh and quercus)

weevil.cca <- plot(spe.cca,
                   #plot(spe.cca.pars,
                   display = c("sp","cn","lc"),
                   scaling = 2,
                   main = "CCA without singletons and doubletons",
                   type="points"
                   )

# Explained variance (adjusted R2)
RsquareAdj(spe.cca) #0.3137338
RsquareAdj(spe.cca.pars) #0.2992039

#Efficiency of the variability possible to explain by any predictors
#in each canonical axis for full CCA:
round(anova(spe.cca, permutations = how(nperm = 999), by = "axis")[,2][1:2],3)
# 0.584 0.384

#For most parsimonious CCA:
round(anova(spe.cca.pars, permutations = how(nperm = 999), by = "axis")[,2][1:2],3)
#0.516 0.336

#Extract CCA coordinates (dimensions) for each score type: predictor scores (biplot
#scores), species scores and site scores (linear combinations of constraining variables)
sp.scor = weevil.cca$species
pred.scor = weevil.cca$biplot
site.scor = weevil.cca$constraints

#Change site scores "2019" to "T1" and "2020" to "T2"
rownames(site.scor) = gsub('_2019', '_T1', rownames(site.scor))
rownames(site.scor) = gsub('_2020', '_T2', rownames(site.scor))

#Subset the coordinates (dimensions) of the first two axis (CCA1, CCA2)
sp.scor = sp.scor[,c(1,2)]
pred.scor = pred.scor[,c(1,2)]
site.scor = site.scor[,c(1,2)]

#Create dataframe with coordinates (dimensions) as the first two columns, then the
#score type (species, predictor, site) and the label (weevil morphospecies, environmental
#predictor name and elevation floor per year)
weevil.cca.df = rbind(pred.scor, sp.scor, site.scor)

#Add score type
scores = c(rep("predictor",dim(pred.scor)[1]),
                        rep("species",dim(sp.scor)[1]),
                        rep("site",dim(site.scor)[1]))
weevil.cca.df = cbind(weevil.cca.df, scores)

#Add labels
labels = rownames(weevil.cca.df)
weevil.cca.df = cbind(weevil.cca.df, labels)

#Ensure the data is stored as a dataframe
weevil.cca.df = as.data.frame(weevil.cca.df)
weevil.cca.df$CCA1 = as.numeric(weevil.cca.df$CCA1)
weevil.cca.df$CCA2 = as.numeric(weevil.cca.df$CCA2)
weevil.cca.df$scores = as.factor(weevil.cca.df$scores)

#Add genera names
#Subset a single unique row for each morphospecies
genera.unique = litter[!duplicated(litter[,c("SpeciesCode")]),]
# Add column with morphospecies taxonomical names
weevil.cca.genera = merge(weevil.cca.df, genera.unique[,c("SpeciesCode", "Genus")],
                          by.x="labels",
                          by.y="SpeciesCode", 
                          all=F, 
                          sort=T)
names(weevil.cca.genera) = c("ID","CCA1","CCA2","scores","labels")

#Replace "NewGenus" for ".gen.nov."
weevil.cca.genera$labels = gsub("NewGenus", ".gen.nov.", weevil.cca.genera$labels)

#Italicize scientific names
weevil.cca.genera$labels = paste0("italic('",weevil.cca.genera$labels,"')")
weevil.cca.genera$labels = gsub("italic('Cryptorhynchini.gen.nov.1')",
                                "'Cryptorhynchini gen. nov. 1'",
                                 weevil.cca.genera$labels,
                                 fixed = T) 
weevil.cca.genera$labels = gsub("italic('Cryptorhynchini.gen.nov.2')",
                                "'Cryptorhynchini gen. nov. 2'",
                                 weevil.cca.genera$labels,
                                 fixed = T) 
weevil.cca.genera$labels = gsub("italic('Cryptorhynchini.gen.nov.3')",
                                "'Cryptorhynchini gen. nov. 3'",
                                 weevil.cca.genera$labels,
                                 fixed = T) 
weevil.cca.genera$labels = gsub("italic('Cryptorhynchini.gen.nov.4')",
                                "'Cryptorhynchini gen. nov. 4'",
                                 weevil.cca.genera$labels,
                                 fixed = T) 
weevil.cca.genera$labels = gsub("italic('Lymantini.gen.nov.1')",
                                "'Lymantini gen. nov. 1'",
                                 weevil.cca.genera$labels,
                                 fixed = T) 
weevil.cca.genera$labels = gsub("italic('Conotrachelini.gen.nov.1')",
                                "'Conotrachelini gen. nov. 1'",
                                 weevil.cca.genera$labels,
                                 fixed = T) 

#Wrap quotes around the labels of scores, sites and species
#so that the PARSE=T argument runs (if not it will stop at
#any label that begins with a number or has a space)
weevil.cca.df$labels = paste0("'",weevil.cca.df$labels,"'")

#Plot CCA with ggplot
CCA_CURBC <- ggplot(weevil.cca.df, 
       aes(x = CCA1, y = CCA2, col = scores, shape = scores,label = labels)) +
  geom_vline(xintercept = 0, lty = "dashed", alpha = .5) +
  geom_hline(yintercept = 0, lty = "dashed", alpha = .5) +
  geom_point(size=2, aes(colour=scores)) +
  labs(x = "CCA 1 (58.4%)", y = "CCA 2 (38.4%)") + #full CCA
  #labs(x = "CCA 1 (51.6%)", y = "CCA 2 (33.6%)") + #most parsimonious CCA
  geom_segment(data=weevil.cca.df[weevil.cca.df$scores=="predictor",],
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               size=2,
               arrow = arrow(length = unit(0.5, "cm")),show.legend=FALSE) +
  #for bind_rows function you need dplyr library 
  #for geom_label_repel function you need ggrepel library
  geom_label_repel(data=bind_rows(weevil.cca.genera,
                                  weevil.cca.df[weevil.cca.df$scores=="site",],
                                  #weevil.cca.genera[grepl("Eurhoptus", weevil.cca.genera$labels),],
                                  weevil.cca.df[weevil.cca.df$scores=="predictor",]),
                   show.legend = F,
                   segment.alpha = 0.5,
                   point.padding = unit(0.5, "points"),
                   max.overlaps = 100,
                   size = 4,
                   parse = T,
                   box.padding = unit(0.5, "lines")) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  theme_classic() +
  theme(axis.text=element_text(size=18), #axis tickmarks size
        axis.title=element_text(size=20), #axis title size
        legend.title=element_blank(), #no legend title
        legend.text=element_text(size=30),
        #legend.position=c(0.1,0.9),
        legend.position="none" #remove legend
        )
ggsave("./Figures/Figure_5.tiff", CCA_CURBC, 
       width=8000, height=6300, units="px", dpi=600, compression="lzw")

#_________________________________
#### χ2 goodness-of-fit tests ####
#_________________________________

# Database that contains leaf litter weevil richness per biogeographic affinity
chi_table_richness <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                              sheet = "chi-squared-table_richness", col_names = T)
# Transform the database object into a table object
chi_richness = tapply(chi_table_richness$Freq, 
                      list(chi_table_richness$Affinity, chi_table_richness$Plot), 
                      sum)

# Database that contains leaf litter weevil abundance per biogeographic affinity
chi_table_abundance <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                              sheet = "chi-squared-table_abundance", col_names = T)
# Transform the database object into a table object
chi_abundance = tapply(chi_table_abundance$Freq, 
                      list(chi_table_abundance$Affinity, chi_table_abundance$Plot), 
                      sum)

# Database that contains leaf litter weevil sample occurrence per biogeographic affinity
chi_table_sampleocc <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                              sheet = "chi-squared-table_sampleocc", col_names = T)
# Transform the database object into a table object
chi_sampleocc = tapply(chi_table_sampleocc$Freq, 
                      list(chi_table_sampleocc$Affinity, chi_table_sampleocc$Plot), 
                      sum)

#Run the Chi2 likelihood test with richness
chisq.test(chi_richness[1:2,]) #X-squared = 3.6297, df = 9, p-value = 0.9341

#With abundance
chisq.test(chi_abundance[1:2,]) #X-squared = 557.6, df = 9, p-value < 0.00000000000000022

#With sample occurrence
chisq.test(chi_sampleocc[1:2,]) # X-squared = 30.687, df = 9, p-value = 0.0003349

#_________________________________________________
#### Colinearity of environmental predictors ####
#_________________________________________________

# Test for microclimatic and climatic environmental predictors colinearity
enviro_predictors <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                                       sheet = "Environmental_predictors", col_names = T)
enviro_predictors = data.frame(enviro_predictors)
names(enviro_predictors)[names(enviro_predictors) == 'altitude'] <- "elevation"

# Variance Inflation Factor (VIF) to test for colinearity
usdm::vifstep(enviro_predictors[
  , names(enviro_predictors) %in% 
    c(                    # Uncomment variables which you want to test for colinearity
       #"T",               # microclimatic -> VIF = 1.5 no-colinear with other microclimatic and WorldClim
       #"RH",              # microclimatic -> VIF = 10 no-colinear with other microclimatic and WorldClim
       "litter_depth",    # microclimatic used by Jones et al. 2012 -> VIF = 10 or 1.5 no-colinear with other microclimatic and WorldClim
       "litter_humidity", # microclimatic used by Jones et al. 2012 -> VIF = 10 no-colinear with other microclimatic and WorldClim
       "bio1",   # WorldClim 2.5arcmin
       "bio2",   # WorldClim 2.5arcmin -> VIF = 10 no-colinear with other WorldClim or microclimatic
       "bio3",   # WorldClim 2.5arcmin
       "bio4",   # WorldClim 2.5arcmin
       "bio5",   # WorldClim 2.5arcmin
       "bio6",   # WorldClim 2.5arcmin
       "bio7",   # WorldClim 2.5arcmin Temperature Annual Range (BIO5-BIO6) -> VIF = 10 or 1.5 no-colinear with other WorldClim
       "bio10",  # WorldClim 2.5arcmin
       "bio11",  # WorldClim 2.5arcmin
       "bio12",  # WorldClim 2.5arcmin
       "bio13",  # WorldClim 2.5arcmin
       "bio14",  # WorldClim 2.5arcmin
       "bio15",  # WorldClim 2.5arcmin Precipitation Seasonality (Coefficient of Variation) -> VIF = 10 or 1.5 no-colinear with other WorldClim and microclimatic
       "bio16",  # WorldClim 2.5arcmin
      #"bio17",  # WorldClim 2.5arcmin Precipitation of Driest Quarter -> VIF = 10 no-colinear with other microclimatic
      #"bio8",   # WorldClim 2.5arcmin -> Escobar et al. 2014 excluded it since they show odd discontinuities between neighbouring pixels
      #"bio9",   # WorldClim 2.5arcmin -> Escobar et al. 2014 excluded it since they show odd discontinuities between neighbouring pixels
      #"bio18",  # WorldClim 2.5arcmin -> Escobar et al. 2014 excluded it since they show odd discontinuities between neighbouring pixels
      #"bio19",  # WorldClim 2.5arcmin -> Escobar et al. 2014 excluded it since they show odd discontinuities between neighbouring pixels
      #"litter_weevil_richness",
      #"neotropical",
      #"nearctic",
      #"unknwon_affinity",
      #"ID",
      #"plot",
      #"year",
      #"lon",
      #"lat",
      "dbh_mean",
      "dbh_median",
      "tree_rich",
      "tree_abun",
      "pinus_prop",
      "quercus_prop",
      "ternstroemia_prop",
      #"boreal_veg_prop",
      #"tropical_veg_prop",
      "elevation"
      )
  ], 
  th=1.5 #Montgomery and Peck (1992) threshold VIF < 10. Zuur et al. (2009) threshold VIF < 1.5
  ) 
#With 25/04/22 configuration (VIF < 1.5)
#      Variables      VIF
# 1 litter_depth 1.192442
# 2     dbh_mean 1.455153
# 3    tree_rich 1.366024
# 4   pinus_prop 1.150687
# 5        bio15 1.385891

# Plot Spearman colinearity between chosen variables
x11()
PerformanceAnalytics::chart.Correlation(enviro_predictors
                  [,names(enviro_predictors) %in%
                      c(                     # Uncomment variables which you want to plot colinearity scatterplots
                        #"T",               # microclimatic -> VIF = 1.5 no-colinear with other microclimatic and WorldClim
                        #"RH",              # microclimatic -> VIF = 10 no-colinear with other microclimatic and WorldClim
                        "litter_depth",    # microclimatic used by Jones et al. 2012 -> VIF = 10 or 1.5 no-colinear with other microclimatic and WorldClim
                        #"litter_humidity", # microclimatic used by Jones et al. 2012 -> VIF = 10 no-colinear with other microclimatic and WorldClim
                        #"bio1",   # WorldClim 2.5arcmin
                        #"bio2",   # WorldClim 2.5arcmin -> VIF = 10 no-colinear with other WorldClim or microclimatic
                        #"bio3",   # WorldClim 2.5arcmin
                        #"bio4",   # WorldClim 2.5arcmin
                        #"bio5",   # WorldClim 2.5arcmin
                        #"bio6",   # WorldClim 2.5arcmin
                        #"bio7",   # WorldClim 2.5arcmin Temperature Annual Range (BIO5-BIO6) -> VIF = 10 or 1.5 no-colinear with other WorldClim
                        #"bio10",  # WorldClim 2.5arcmin
                        #"bio11",  # WorldClim 2.5arcmin
                        #"bio12",  # WorldClim 2.5arcmin
                        #"bio13",  # WorldClim 2.5arcmin
                        #"bio14",  # WorldClim 2.5arcmin
                        "bio15",  # WorldClim 2.5arcmin Precipitation Seasonality (Coefficient of Variation) -> VIF = 10 or 1.5 no-colinear with other WorldClim and microclimatic
                        #"bio16",  # WorldClim 2.5arcmin
                        #"bio17",  # WorldClim 2.5arcmin Precipitation of Driest Quarter -> VIF = 10 no-colinear with other microclimatic
                        #"bio8",   # WorldClim 2.5arcmin -> Escobar et al. 2014 excluded it since they show odd discontinuities between neighbouring pixels
                        #"bio9",   # WorldClim 2.5arcmin -> Escobar et al. 2014 excluded it since they show odd discontinuities between neighbouring pixels
                        #"bio18",  # WorldClim 2.5arcmin -> Escobar et al. 2014 excluded it since they show odd discontinuities between neighbouring pixels
                        #"bio19",  # WorldClim 2.5arcmin -> Escobar et al. 2014 excluded it since they show odd discontinuities between neighbouring pixels
                        #"abundance",
                        #"richness",
                        #"neotropical",
                        #"nearctic",
                        #"unknwon_affinity",
                        #"ID",
                        #"plot",
                        #"year",
                        #"lon",
                        #"lat",
                        #"dbh_mean",
                        "dbh_median",
                        "tree_rich",
                        #"tree_abun",
                        #"pinus_prop",
                        #"quercus_prop",
                        #"ternstroemia_prop",
                        #"boreal_veg_prop",
                        #"tropical_veg_prop",
                        "elevation"
                      )],
                  histogram=F,
                  pch=21,
                  method = "spearman",
                  cex=5)

# 25/04/22: substitute dbh_mean for dbh_median because no colinearity.
#           pinus proportion was removed due to colinearity with altitude

#________________________________________________
#### Richness ~ environmental predictors GLM ####
#________________________________________________

# Load leaf litter weevil richness standardized predictors table.
# Standardization was carried out by subtracting the average of 
# each value and dividing the result by its standard deviation.
# Test for microclimatic and climatic environmental predictors colinearity
enviro_predictors <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                                       sheet = "Environmental_predictors", col_names = T)
enviro_predictors = data.frame(enviro_predictors)

# Testing normality of litter weevil richness
shapiro.test(enviro_predictors$richness) #p-value = 0.851 normal!
qqnorm(enviro_predictors$richness); qqline(enviro_predictors$richness)
hist(enviro_predictors$richness)

# Testing normality of litter weevil abundance
shapiro.test(enviro_predictors$abundance) #p-value = 0.02205 not-normal!
qqnorm(enviro_predictors$abundance); qqline(enviro_predictors$abundance)
hist(enviro_predictors$abundance)

# #Optional step (not currently used)
# #Normalizing variables = scale independent numeric variables because their units
# #are in different orders of magnitude, so in order to compare variation between
# #the variables, they should all go from 0 to 1
# normalize<-function(x){
#   (x-min(x))/(max(x)-min(x))
# }

# Run the model with non-colinear predictors (see section Colinearity of environmental predictors)
m <- glm(richness ~ (altitude + litter_depth + dbh_median + tree_rich + bio15)^2,
         family = poisson(link=log),
         data = enviro_predictors)

# Standardized predictors
m_stan <- glm(richness ~ (altitude_stan + litter_depth_stan + dbh_median_stan + tree_rich_stan + bio15_stan)^2,
         family = poisson(link=log),
         data = enviro_predictors)

# Compare which model has best fit (lowest AIC)
AIC(m, m_stan)
#        df      AIC
# m      16 120.7175
# m_stan 16 120.7175
# Same AIC value, best use non-standardized predictors for easier biological interpretation

# Stepwise backwards simplification with Crawley (2013) method using Chi-squared tests:
m2 <- update(m, ~.- altitude:tree_rich)
anova(m2, m, test = "Chisq") #comparison between m2 and m3 non-significant? drop the variable
                             #comparison between m2 and m3 significant? keep the variable
                             #P = 0.5986, non-significant, drop the variable
m3 <- update(m2, ~.- litter_depth:tree_rich)
anova(m3, m2, test = "Chisq") #P = 0.7317
m4 <- update(m3, ~.- tree_rich:bio15)
anova(m4, m3, test = "Chisq") #P = 0.3239
m5 <- update(m4, ~.- altitude:litter_depth)
anova(m5, m4, test = "Chisq") #P = 0.2495
m6 <- update(m5, ~.- altitude:bio15)
anova(m6, m5, test = "Chisq") #P = 0.2144
m7 <- update(m6, ~.- litter_depth:bio15)
anova(m7, m6, test = "Chisq") #P = 0.001356

# Minimum model chosen:
mm <- m6
summary(mm)

# Check for overdispersion by dividing residual deviance / residual degrees of freedom
mm$deviance/mm$df.residual # = 1.042931, dispersion parameter for poisson family taken
                           # to be 1

# Correct for overdispersion by using quasipoisson(link=log) error family
mm <- glm(formula = richness ~ altitude + litter_depth + dbh_median + 
    tree_rich + bio15 + altitude:dbh_median + litter_depth:dbh_median + 
    litter_depth:bio15 + dbh_median:tree_rich + dbh_median:bio15, 
    family = quasipoisson(link = log), data = enviro_predictors)
summary(mm) # dispersion parameter of quasipoisson taken to be 1.045859

#                            Estimate  Std. Error t value Pr(>|t|)   
# (Intercept)              33.1109466  18.2768838   1.812  0.10347   
# altitude                  0.0039357   0.0009746   4.038  0.00294 **
# litter_depth            -28.4389331   9.1223803  -3.117  0.01237 * 
# dbh_median              335.5521584  84.2075506   3.985  0.00318 **
# tree_rich                 0.4658854   0.1239381   3.759  0.00449 **
# bio15                    -0.4510877   0.2387698  -1.889  0.09145 . 
# altitude:dbh_median      -0.0340815   0.0100214  -3.401  0.00786 **
# litter_depth:dbh_median  28.1040087   7.4874720   3.753  0.00453 **
# litter_depth:bio15        0.3369767   0.1108087   3.041  0.01399 * 
# dbh_median:tree_rich     -4.8582935   1.2319600  -3.944  0.00339 **
# dbh_median:bio15         -3.9608610   1.0078510  -3.930  0.00346 **

# Total variation explained by the model (using McFadden pseudo-R2)
# https://web.archive.org/web/20130701052120/http://www.ats.ucla.edu:80/stat/mult_pkg/faq/general/Psuedo_RSquareds.htm
round((mm$null.deviance-mm$deviance)/mm$null.deviance,3) # = 0.896
round(1-(mm$deviance/mm$null.deviance),3) # = 0.896

# Visualize partitioned deviance of minimum model
pvalue = summary(mm)$coef[,4] #extract p-values of minimum model
pseudoR2 = anova(mm)[2]/mm$null.deviance #partitioned McFadden pseudo R2
names(pseudoR2) = "pseudoR2"
mm.R.Pvalue = round(cbind(pseudoR2,pvalue),3)[-1,]; mm.R.Pvalue

# GLM results table
#                         pseudoR2 pvalue
# altitude                   0.602  0.003
# litter_depth               0.029  0.012
# dbh_median                 0.003  0.003
# tree_rich                  0.009  0.004
# bio15                      0.009  0.091
# altitude:dbh_median        0.002  0.008
# litter_depth:dbh_median    0.001  0.005
# litter_depth:bio15         0.003  0.014
# dbh_median:tree_rich       0.030  0.003
# dbh_median:bio15           0.209  0.003

#____________________________________________
#### Beta ~ environmental predictors GLM ####
#____________________________________________

# Load leaf litter weevil richness standardized predictors table.
# Standardization was carried out by subtracting the average of 
# each value and dividing the result by its standard deviation.
enviro_predictors <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                                       sheet = "Environmental_predictors", col_names = T)
enviro_predictors = data.frame(enviro_predictors)

# Testing normality of litter weevil mean βsor
shapiro.test(enviro_predictors$mean_Bsor) #p-value = 0.05656 normal!
qqnorm(enviro_predictors$mean_Bsor); qqline(enviro_predictors$mean_Bsor)
hist(enviro_predictors$mean_Bsor)

# GLM of beta taxonomical diversity (Sorensen index) given by
# the chosen environmental predictors:
m <- glm(mean_Bsor ~ (altitude * litter_depth * dbh_median * tree_rich * bio15 * pinus_prop),
         family = binomial(link = logit), #binomial family for proportion data
         data = enviro_predictors)
summary(m) #ALL P VALUES > 0.9 (NOT SIGNIFICANT IN THE SLIGHTEST!)

#### Abundance, richness, density, Chao 2 and Jack 1 ####
#________________________________________________________

# Load community abundance matrix
weevil_abundance_matrix <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                               sheet = "Community_abundance_matrix", col_names = T)

# Set the first column (elevational floors) as row names and the
# elevational floors (200:2000 masl) as column names
weevil_abundance_matrix = as.data.frame(weevil_abundance_matrix)
weevil_abundance <- weevil_abundance_matrix[,-1]
rownames(weevil_abundance) <- weevil_abundance_matrix[,1]
names(weevil_abundance) <- seq(200,2000,200)

# Abundance measurements for first results paragraph

sum(weevil_abundance) #4,004 leaf litter weevil individuals
dim(weevil_abundance) #55 morphospecies in 10 elevational floors

# Most abundant morphospecies
sort(rowSums(weevil_abundance), decreasing = T)
#CURC002 1113
#CURC010 760
#CURC001 378

# Doubletons
length(rowSums(weevil_abundance)[rowSums(weevil_abundance)==2]) #5

# Singletons
length(rowSums(weevil_abundance)[rowSums(weevil_abundance)==1]) #8

# Richness and density measurements per 1 m2.
# Read weevils per sample database
weevils_per_sample <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                               sheet = "weevils_per_sample", col_names = T)

# Set the first column (morphospecies names) as row names and the
# sample names as column names
weevils_per_sample = as.data.frame(weevils_per_sample)
weevils_sample <- weevils_per_sample[,-1]
rownames(weevils_sample) <- weevils_per_sample[,1]
names(weevils_sample) <- awm[1,]

#Mean and SD density per sample
round(mean(colSums(weevils_sample)),2) #13.76
round(sd(colSums(weevils_sample)),2) #16.38

#Mean and SD richness per sample
rich.per.sample = weevils_sample
rich.per.sample[rich.per.sample > 0] <- 1
round(mean(colSums(rich.per.sample)),2) #3.63
round(sd(colSums(rich.per.sample)),2) #3.02

#Richness estimators
ace <- vegan::estimateR(rowSums(weevils_sample)); ace
#     S.obs   S.chao1  se.chao1     S.ACE    se.ACE 
# 55.000000 59.666667  4.489587 59.292163  3.733364  

chao_jack <- vegan::specpool(t(weevils_sample)); chao_jack
#     Species     chao chao.se    jack1 jack1.se    jack2     boot  boot.se   n
# All      55 72.93814 14.3402 66.95876 3.729718 74.91748 60.25303 2.046085 291

c(round(chao_jack[1,"chao"]),round(1.96 * chao_jack[1,"chao.se"])) #73+-28
c(round(chao_jack[1,"jack1"]),round(1.96 * chao_jack[1,"jack1.se"])) #67+-7

#Chao 2 % completeness
round((chao_jack[1,"Species"]*100)/chao_jack[1,"chao"]) #75%

#Jackknife 1 % completeness
round((chao_jack[1,"Species"]*100)/chao_jack[1,"jack1"]) #82%

#________________________________________________________
#### Inventory completeness (supplementary figure 2) ####
#________________________________________________________

# Create input file for EstimateS v.9

# Subset morphospecies and Plot_Transect_Line_Sample to obtain samples dataframe
sdf = as.data.frame(weevil_data[,c("SpeciesCode","PTLS","Abundance","Plot","Transect")])
sdf = na.omit(sdf) # remove NAs

# Rename columns
colnames(sdf) = c("Morphospecies","Sample","Abundance","Plot","Year")

# Eliminate the "_1" or "_2" from Transect column, so as
# to only have the year (2019 or 2020)
sdf$Year =  gsub('_1', '', sdf$Year)
sdf$Year =  gsub('_2', '', sdf$Year)
sdf$Year = as.numeric(sdf$Year) #set Year as numeric
sdf$Sample = as.factor(sdf$Sample) #set sample as factor

# Sort dataframe in ascending year and plot order
sdf = sdf[with(sdf, order(Year,Plot)), ]

# Drop Plot and Year columns
sdf = sdf[,c("Morphospecies","Sample","Abundance")]

# Reshape three column dataframe to an abundance weevil matrix
sdf$Sample = factor(sdf$Sample, levels=unique(as.character(sdf$Sample)))
awm = xtabs(Abundance~Morphospecies+Sample, data=sdf)

# Save it in format required by EstimateS v.9
write("RBC litter weevil abundance data (2019-2020)", #file title
      file="./leaf_litter_weevil_abundance_EstimateS_input.txt")
write(paste0(dim(awm)[1], #morphospecies number (rows)
             "\t",
             dim(awm)[2]), #sample number (columns)
      file = "./leaf_litter_weevil_abundance_EstimateS_input.txt",
      append=T)
write.table(awm, file = "./leaf_litter_weevil_abundance_EstimateS_input.txt",
            sep="\t",
            row.names = T,
            col.names = T,
            append=T)

# As of 11 of March 2022 Colwell's ESTIMATES v9 software was ran with the following indications:
# input file was "leaf_litter_weevil_abundance_EstimateS_input.txt", 1000 runs, non-biased 
# corrected Chao, sampling with replacement, each knot at every sample (sample size + extrapolated),
# extrapolation up to double the sample size (see Colwell & Chao, 2014). 
# Read the output file:
EstimateS <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                               sheet = "EstimateS_output", col_names = T)
EstimateS = as.data.frame(EstimateS)

# Supplementary figure 2A: ESTIMATES plot
par(mar=c(2.5, 2.5, 0.1, 0.1))
plot(EstimateS$S.est.[1:231], 
     type = "l", 
     ylim=c(0,60),
     xlim=c(0,230),
     xlab = "",
     ylab = "",
     yaxt = "n",
     xaxt = "n",
     col = "#661100",
     lwd = 5,
     cex.axis = 1,
     bg = "white",
     las = 1)
mtext("Species diversity", side = 2, line = 1.5, cex = 1.2)
mtext("Number of samples", side = 1, line = 1.5, cex = 1.2)
axis(2, at=c(0,10,20,30,40,50,60), labels = c(0,10,20,30,40,50,60), las=1, hadj = 0.3, tck = -0.01)
axis(1, at=c(0,50,100,150,200), labels = c(0,50,100,150,200), las=1, padj = -1.2, tck = -0.01)
points(EstimateS$ICE_Mean[1:231], type = "l", col="#88CCEE", lwd = 3)#ICE
points(EstimateS$ACE_Mean[1:231], type = "l", col="#CC6677", lwd = 3)#ACE
points(EstimateS$Chao_2_Mean[1:231], type = "l", col="#DDCC77", lwd = 3)#Chao2
points(EstimateS$Jack_1_Mean[1:231], type = "l", col="#117733", lwd = 3)#Jack1
points(EstimateS$Singletons_Mean[1:231], type = "l", col="#44AA99", lwd = 3)#singletons
points(EstimateS$Doubletons_Mean[1:231], type = "l", col="#AA4499", lwd = 3)#doubletons
legend(x=120, y=40,
       legend=c("Jackknife 1","Doubletons","Singletons","Richness","Chao 2","ACE","ICE"),
       col = c("#117733","#AA4499","#44AA99","#661100","#DDCC77","#CC6677","#88CCEE"),
       cex = 1,
       lty = 1, 
       lwd = 2,
       y.intersp = 1.2,
       bty = "n")
ESTIMATESplot <- recordPlot()

# Create leaf litter weevil community abundance matrix where:
# columns are the 10 elevational sites (200:2000 masl) and 
# rows are the 55 leaf litter weevil morphospecies

# Subset morphospecies, Plot and Abundance
elevplot = as.data.frame(weevil_data[,c("SpeciesCode","Plot","Abundance","Transect")])
elevplot = na.omit(elevplot) # remove NAs

# Turn transect into year (delete "_1", "_2")
elevplot$Transect =  gsub('_1', '', elevplot$Transect)
elevplot$Transect =  gsub('_2', '', elevplot$Transect)

# Rename columns
colnames(elevplot) = c("Morphospecies","Plot","Abundance","Year")

# Sort dataframe in ascending year and plot order
elevplot = elevplot[with(elevplot, order(Year,Plot)), ] #sort by descending order
elevplot$Year = as.factor(elevplot$Year)
elevplot$Plot = as.factor(elevplot$Plot)

# Abundance matrix per elevation plot
abunelevplot = xtabs(Abundance~Morphospecies+Plot, data=elevplot)

# Save the matrix
write.csv(abunelevplot, file = "./Community_abundance_matrix.csv")

# Load community abundance matrix
weevil_abundance_matrix <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                               sheet = "Community_abundance_matrix", col_names = T)

# Set the first column (elevational floors) as row names and the
# elevational floors (200:2000 masl) as column names
weevil_abundance_matrix = as.data.frame(weevil_abundance_matrix)
weevil_abundance <- weevil_abundance_matrix[,-1]
rownames(weevil_abundance) <- weevil_abundance_matrix[,1]
names(weevil_abundance) <- seq(200,2000,200)

# Create list of selected sites
weevils_Sites_Abun = list("0200 m" = apply(weevil_abundance[1,],2,as.integer), 
                         "0400 m" = apply(weevil_abundance[2,],2,as.integer),
                         "0600 m" = apply(weevil_abundance[3,],2,as.integer),
                         "0800 m" = apply(weevil_abundance[4,],2,as.integer),
                         "1000 m" = apply(weevil_abundance[5,],2,as.integer),
                         "1200 m" = apply(weevil_abundance[6,],2,as.integer), 
                         "1400 m" = apply(weevil_abundance[7,],2,as.integer), 
                         "1600 m" = apply(weevil_abundance[8,],2,as.integer),
                         "1800 m" = apply(weevil_abundance[9,],2,as.integer),
                         "2000 m" = apply(weevil_abundance[10,],2,as.integer))

# Apply iNEXt() on data
# - modified code from Moses et al. 2021 -
out_weevilSites_Abun_q0 <- iNEXT(weevils_Sites_Abun, q=0, datatype="abundance",
                                    se=TRUE, 
                                    conf=0.95, 
                                    nboot=999, 
                                    endpoint=30, #max of 30 units 
                                    knots = 40)

out_weevilSites_Abun_q0_type2 <- iNEXT(weevils_Sites_Abun, q=0, datatype="abundance",
                                    se=TRUE,
                                    conf=0.95,
                                    nboot=999,
                                    endpoint=100, #max of 100 units
                                    knots = 40)

# Create input data for iNEXT ggplots
# - modified code from Moses et al. 2021 -

# For accumulation curve (iNEXT curve type 1): 
# number of samples by species diversity
df_outSITES_q0type1_data <- ggiNEXT(out_weevilSites_Abun_q0, type = 1)$data %>% 
  dplyr::mutate(id = rownames(.)) %>% 
  dplyr::select(id, everything()) %>% 
  tidyr::as_tibble()

# For coverage-based rarefaction curve (iNEXT curve type 3): 
# sample coverage by species diversity
df_outSITES_q0type3_data <- ggiNEXT(out_weevilSites_Abun_q0, type = 3)$data %>%
  dplyr::mutate(id = rownames(.)) %>%
  dplyr::select(id, everything()) %>%
  tidyr::as_tibble()

# For sample completeness rarefaction curve (iNEXT curve type 2): 
# number of individuals by sample coverage
df_outSITES_q0type2_data <- ggiNEXT(out_weevilSites_Abun_q0_type2, type = 2)$data %>%
  dplyr::mutate(id = rownames(.)) %>%
  dplyr::select(id, everything()) %>%
  tidyr::as_tibble()

# Supplementary figure 2B: Accumulation curve of number of samples by species diversity
numbersamples_speciesdiversity <- ggplot2::ggplot(df_outSITES_q0type1_data, aes(x = x, y = y)) +
  ggplot2::geom_line(aes(x = x, y = y, col = Assemblage), lty = 1, lwd = 2) +
  ggplot2::geom_ribbon(aes(ymin = y.lwr, ymax = y.upr,
                           fill = Assemblage, colour = NULL), alpha = 0.25, show.legend = F) +
  #Numbers after id %in% represent the last row of values of each elevational floor
  #(ex. for elevational floor 0200 m there're 30 rows of values, therefore its number
  #will be 30, and so the first number between parentheses after id %in% c() is "30")
  ggplot2::geom_point(data = subset(df_outSITES_q0type1_data, x==15), #for type 1 all numbers are
                                                                      #set for 15, because 15 samples
                                                                      #were processed for each floor
             aes(x = x, y = y, col = Assemblage, shape = Assemblage), size = 5) +
  ggplot2::scale_colour_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288",
                                          "#AA4499", "#44AA99", "#999933", "#882255", "#661100")) + 
  ggplot2::scale_fill_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288",
                                          "#AA4499", "#44AA99", "#999933", "#882255", "#661100")) + 
  ggplot2::scale_shape_manual(values = c(22, 20, 21, 19, 18, 17, 16, 15, 14, 13)) + 
  ggplot2::xlim(0,15) + ggplot2::ylim(0,10) +
  ggplot2::labs(x = "Number of samples", y = "Species diversity") +
  ggplot2::theme_bw(base_size = 15, base_family = "sans") +
  ggplot2::theme(
    axis.text = element_text(color = "black"),
    legend.justification = c(0,1), 
    legend.position = c(0.001, 0.9),
    legend.direction = "vertical",
    legend.title = element_blank(),
    strip.background = element_blank(), 
    strip.text.x = element_blank(),
    legend.background = element_rect(fill = "transparent", 
                                     size = 0.5, 
                                     linetype = 1,
                                     colour = "transparent"),
    legend.text = element_text(size = 10), 
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"))

# Supplementary figure 2C: Coverage based rarefaction curve of sample coverage by species diversity
samplecoverage_speciesdiversity <- ggplot2::ggplot(df_outSITES_q0type3_data, aes(x = x, y = y)) +
  ggplot2::geom_line(aes(x = x, y = y, col = Assemblage), lty = 1, lwd = 2) +
  ggplot2::geom_ribbon(aes(ymin = y.lwr, ymax = y.upr,
                  fill = Assemblage, colour = NULL), alpha = 0.25, show.legend = F) +
  #Numbers after id %in% represent the last row of values of each elevational floor
  #(ex. for elevational floor 0200 m there're 30 rows of values, therefore its number
  #will be 30, and so the first number between parentheses after id %in% c() is "30")
  ggplot2::geom_point(data = subset(df_outSITES_q0type3_data, id %in% c("30","60","90","110","147",
                                                                        "177","207","237","267","297")), 
             aes(x = x, y = y, col = Assemblage, shape = Assemblage), size = 5) +
  ggplot2::scale_colour_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288",
                                          "#AA4499", "#44AA99", "#999933", "#882255", "#661100")) + 
  ggplot2::scale_fill_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288",
                                          "#AA4499", "#44AA99", "#999933", "#882255", "#661100")) + 
  ggplot2::scale_shape_manual(values = c(22, 20, 21, 19, 18, 17, 16, 15, 14, 13)) + 
  ggplot2::xlim(0,1) + ggplot2::ylim(0,15) +
  ggplot2::labs(x = "Sample coverage", y = "Species diversity") +
  ggplot2::theme_bw(base_size = 15, base_family = "sans") +
  ggplot2::theme(
    axis.text = element_text(color = "black"),
    legend.justification = c(0,1), 
    legend.position = c(0.001, 0.8),
    legend.direction = "vertical",
    legend.title = element_blank(),
    strip.background = element_blank(), 
    strip.text.x = element_blank(),
    legend.background = element_rect(fill = "transparent", 
                                     size = 0.5, 
                                     linetype = 1,
                                     colour = "transparent"),
    legend.text = element_text(size = 10), 
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"))

# Supplementary figure 2D: Sample-completeness rarefaction curve of individuals by sample coverage
individuals_samplecoverage <- ggplot2::ggplot(df_outSITES_q0type2_data, aes(x = x, y = y)) +
  ggplot2::geom_line(aes(x = x, y = y, col = Assemblage), lty = 1, lwd = 2) +
  ggplot2::geom_ribbon(aes(ymin = y.lwr, ymax = y.upr,
                  fill = Assemblage, colour = NULL), alpha = 0.25, show.legend = F) +
  #Numbers after id %in% represent the last row of values of each elevational floor
  #(ex. for elevational floor 0200 m there're 40 rows of values, therefore its number
  #will be 40, and so the first number between parentheses after id %in% c() is "40")
  ggplot2::geom_point(data = subset(df_outSITES_q0type2_data, id %in% c("40","80","120","140","200",
                                                                        "240","280","320","360","400")), 
             aes(x = x, y = y, col = Assemblage, shape = Assemblage), size = 5) +
  ggplot2::scale_colour_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288",
                                          "#AA4499", "#44AA99", "#999933", "#882255", "#661100")) + 
  ggplot2::scale_fill_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288",
                                          "#AA4499", "#44AA99", "#999933", "#882255", "#661100")) + 
  ggplot2::scale_shape_manual(values = c(22, 20, 21, 19, 18, 17, 16, 15, 14, 13)) + 
  ggplot2::xlim(0,100) + ggplot2::ylim(0,1) +
  ggplot2::labs(x = "Number of individuals", y = "Sample coverage") +
  ggplot2::theme_bw(base_size = 15, base_family = "sans") +
  ggplot2::theme(
    axis.text = element_text(color = "black"),
    legend.justification = c(0,1), 
    legend.position = c(0.7, 0.7),
    legend.direction = "vertical",
    legend.title = element_blank(),
    strip.background = element_blank(), 
    strip.text.x = element_blank(),
    legend.background = element_rect(fill = "transparent", 
                                     size = 0.5, 
                                     linetype = 1,
                                     colour = "transparent"),
    legend.text = element_text(size = 10), 
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"))

#Combine the four inventory completeness subplots
figure_inventory = cowplot::plot_grid(ESTIMATESplot, numbersamples_speciesdiversity,
                                      samplecoverage_speciesdiversity , individuals_samplecoverage,
                                     labels = c("A","B","C","D"), hjust = -2.5,
                                     vjust = 1.5, label_size = 30)

#Save richness plots as figure
ggsave("./Figures/Suppl_Figure_2.tiff", figure_inventory, 
       width=6000, height=6000, units="px", dpi=600, compression="lzw")

#_____________________________________________________________________
#### Richness per biogeographic affinity (supplementary figure 3) ####
#_____________________________________________________________________

# Load leaf litter weevil richness per elevational floor database
weevil_richness_elevation <- readxl::read_xlsx("./El_Cielo_weevils_data.xlsx",
                             sheet = "richness_per_elevation", col_names = T)

# Plot weevil richnes by biogeographical affinity
rich.per.affinity <- ggplot(data=affinity_rich, aes(x=plot)) + 
  
  #Nearctic affinity points and regression line
  geom_point(aes(y=nearctic, color="#56B4E9"), size=3) + 
  geom_smooth(aes(plot,nearctic), method=lm, se=T, size=1.5,
              alpha=0.2, colour="#56B4E9", fill="#56B4E9",
              formula = y ~ x) + #best deltaAICc = linear (although S-W P=0.027)
  
  #Neotropical affinity points and regression line
  geom_point(aes(y=neotropical, color="#009E73"), size=3) + 
  geom_smooth(aes(plot,neotropical), method=lm, se=T, size=1.5,
              alpha=0.2, colour="#009E73", fill="#009E73",
              formula = y ~ I(x^1)+I(x^2)) + #best deltaAICc = quadratic
  
  #Legend
  scale_color_manual(name = "Affinity",
                     values = c("#009E73","#56B4E9"),
                     labels = c("Neotropical", "Nearctic")) +
  
  #X-axis and Y-axis title labels  
  labs(x = "Elevation (m)", y = "Leaf litter weevil richness", color="Affinity") + 
  
  #Axis limits
  ylim(1, 30) +
  
  #Plot theme (classic eliminates background gridlines and the full plot box)
  theme_classic() + 
  
  #Set axis and legend characteristics
  theme(
    axis.text = element_text(color = "black"),
    legend.position = c(0.15,0.9), 
    legend.text = element_text(size=11),
    legend.title = element_text(size=11, face = "bold")
  )

ggsave("./Figures/Suppl_Figure_3.tiff", rich.per.affinity, 
       width=4000, height=2500, units="px", dpi=600, compression="lzw")

# Models to test relationship between elevation and richness 
# following methodology of Joaqui et al. 2021:
# (a) linear, (b) quadratic, (c) cubic, (d) exponential, (e) null
# They only tested exponential for neotropical affinity
# Models with ΔAICc < 2 and a normal distribution of residuals
# (Shapiro– Wilk test, p > 0.05) were selected and considered 
# to be models with good fit (Hurvich & Tsai, 1989).

# Neotropical weevil richness
linear_neo = glm(neo$richness ~ neo$plot, family = poisson(link=log)) #linear
quad_neo = glm(neo$richness ~ neo$plot+I(neo$plot^2), family = poisson(link=log)) #quadratic
cubic_neo = glm(neo$richness ~ neo$plot+I(neo$plot^2)+I(neo$plot^3), family = poisson(link=log)) #cubic
expo_neo = glm(neo$richness ~ I(neo$plot^2), family = poisson(link=log)) #exponential
null_neo = glm(neo$richness ~ 1, family = poisson(link=log)) #null

# Model average (ΔAICc, weighted AICc)
round(summary(model.avg(linear_neo,quad_neo,cubic_neo,expo_neo,null_neo))$msTable,2)
#        df logLik   AICc delta weight
# 13      3 -56.46 120.42  0.00   0.57 -> quadratic
# 123     4 -55.32 121.30  0.88   0.37
# 3       2 -60.07 124.86  4.43   0.06
# 1       2 -64.32 133.34 12.91   0.00
# (Null)  1 -81.98 166.18 45.76   0.00

# Distribution of residuals (normal if S-W p > 0.05)
shapiro.test(linear_neo$residuals) #p-value = 0.4035
shapiro.test(quad_neo$residuals) #p-value = 0.6724
shapiro.test(cubic_neo$residuals) #p-value = 0.4546
shapiro.test(expo_neo$residuals) #p-value = 0.5581
shapiro.test(null_neo$residuals) #p-value = 0.7268

# Nearctic weevil richness
linear_arc = glm(arc$richness ~ arc$plot, family = poisson(link=log)) #linear
quad_arc = glm(arc$richness ~ arc$plot+I(arc$plot^2), family = poisson(link=log)) #quadratic
cubic_arc = glm(arc$richness ~ arc$plot+I(arc$plot^2)+I(arc$plot^3), family = poisson(link=log)) #cubic
expo_arc = glm(arc$richness ~ I(arc$plot^2), family = poisson(link=log)) #exponential
null_arc = glm(arc$richness ~ 1, family = poisson(link=log)) #null

# Model average (ΔAICc, weighted AICc)
round(summary(model.avg(linear_arc,quad_arc,cubic_arc,expo_arc,null_arc))$msTable,2)
#        df logLik  AICc delta weight
# 1       2 -20.76 46.23  0.00   0.50 -> linear
# 2       2 -21.31 47.32  1.09   0.29
# 12      3 -20.50 48.49  2.27   0.16
# 123     4 -20.46 51.59  5.36   0.03
# (Null)  1 -25.23 52.69  6.46   0.02

# Distribution of residuals (normal if S-W p > 0.05)
shapiro.test(linear_arc$residuals) #p-value = 0.02705
shapiro.test(quad_arc$residuals) #p-value = 0.004291
shapiro.test(cubic_arc$residuals) #p-value = 0.002465
shapiro.test(expo_arc$residuals) #p-value = 0.07844
shapiro.test(null_arc$residuals) #p-value = 0.009588