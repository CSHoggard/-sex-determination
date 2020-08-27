### ARTICLE TITLE: 
### Investigating sexual dimorphism in the human mandible 
### through a geometric morphometric approach: an Indian population case study

### AUTHORS: 
### Abraham Johnson, Christian Steven Hoggard, Sraddha Singh and Anju Thomas

### ABSTRACT:
### The determination of sex through morphological studies of the mandible 
### is of particular significance in forensic anthropology and bioarchaeological 
### cases, where fragments often remain. The present research utilizes recent 
### developments in geometric morphometrics to examine the degree to which the mandible 
### can be a marker of sex. Using a modern Indian population, and orthopantomographic 
### images, two-dimensional landmark and semi-landmark analyses were performed. 
### Through a multivariate analytical and exploratory framework, incorporating 
### ordination-based methods and supervised classificatory techniques, morphological 
### changes were identified to be consistent with sex. The findings of this research 
### were highly optimistic, supporting previous studies, and suggests that this 
### approach may be used to assess sexual dimorphism. Thus, a new reference for sex 
### determination by orthopantomography is created that will be of great benefit to 
### the entire forensic science community. 

### SCRIPT AUTHOR: 
### Dr Christian Steven Hoggard

### SCRIPT CONTACT: 
### C.S.Hoggard@soton.ac.uk 

### LAST EDITED: 
### 20/08/2020 

### PACKAGE REQUIREMENTS ###

if(!require("geomorph")) install.packages('geomorph', repos='http://cran.us.r-project.org')
if(!require("tidyverse")) install.packages('tidyverse', repos='http://cran.us.r-project.org')
if(!require("Momocs")) install.packages('Momocs', repos='http://cran.us.r-project.org')
if(!require("extrafont")) install.packages('extrafont', repos='http://cran.us.r-project.org')

library(geomorph)
library(tidyverse)
library(Momocs)
library(extrafont)

### DATA IMPORTING ###

shape.table <- read.csv("Johnson_et_al_2020.csv", header = T, row.names = 1)
shape.table$sex <- as.factor(shape.table$sex)
shape.table <- as_tibble(shape.table)

shape.data <- import_tps("Johnson_et_al_2020.tps")

### CLASS CREATION ###

shape.ldk  <- Ldk(shape.data$coo, fac = shape.table, links = NULL)

### PANEL VISUALISATION ###

panel(shape.ldk)

### PROCRUSTES (FG) ###

shape.gpa <- fgProcrustes(shape.ldk)
stack(shape.gpa, 
      title = "", 
      centroid = FALSE, 
      xy.axis = FALSE,
      ldk_cex = 0.75, 
      meanshape = TRUE)

### PCA ###

shape.pca <- PCA(shape.gpa)
scree(shape.pca)

shape.pca.data <- as_tibble(shape.pca$x)
shape.data.new <- cbind(shape.table, shape.pca.data)

ggplot(shape.data.new, aes(PC1, PC2, colour = sex)) + 
  geom_point(alpha = 0.9) + 
  lims(x = c(-0.1,0.1), y = c(-0.1,0.1)) + 
  labs(x = "Principal Component 1 (43.0%)",
       y = "Principal Component 2 (18.2%)",
       colour = "Sex") +
  scale_colour_manual(values=c("orange2","royalblue4"), labels = c("Female", "Male")) +
  stat_ellipse(level = 0.66) +
  theme_minimal() +
  coord_fixed() +
  theme(
    text = element_text(family = "Open Sans"),
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
    plot.margin = unit(c(5, 5, 5, 5), "mm")
  )

ggsave("figure_2.tiff", plot = last_plot(), dpi = 300)

### DISCRIMINATION ANALYSIS (PC SCORES AND PROCRUSTES COORDINATES) ###

shape.da.pca <- LDA(shape.pca, ~sex)
shape.da.gpa <- LDA(shape.gpa, ~sex)

classification_metrics(shape.da.pca)
classification_metrics(shape.da.gpa)

### PROCRUSTES ANOVA ###

shape.data.2 <- readland.tps("Johnson_et_al_2020.tps")

shape.gpa2 <- gpagen(shape.data.2)
plot(shape.gpa2)

shape.gdf <- geomorph.data.frame(shape.gpa2, sex = shape.table$sex)
anova(procD.lm(coords ~ sex, data = shape.gdf))

### SHAPE CHANGES ###

shape.pca.2 <- gm.prcomp(shape.gpa2$coords)


plotRefToTarget(shape.pca.2$shapes$shapes.comp1$min, 
                shape.pca.2$shapes$shapes.comp1$max, 
                method = "vector",
                label = FALSE,
                gridPars = gridPar(pt.size = 1.25))

title("First Principal Component (43.0%)", cex.main=1.5, line = 0)

plotRefToTarget(shape.pca.2$shapes$shapes.comp2$min, 
                shape.pca.2$shapes$shapes.comp2$max, 
                method = "vector",
                label = FALSE,
                gridPars = gridPar(pt.size = 1.25))

title("Second Principal Component (18.2%)", cex.main=1.5, line = 0)

plotRefToTarget(shape.pca.2$shapes$shapes.comp3$min, 
                shape.pca.2$shapes$shapes.comp3$max, 
                method = "vector",
                label = FALSE,
                gridPars = gridPar(pt.size = 1.25))

title("Third Principal Component (11.3%)", cex.main=1.5, line = 0)

### MEAN SHAPES ###

mean.shape.sex <- MSHAPES(shape.gpa, ~sex)
coo_plot(mean.shape.sex$shp$male, poly = FALSE, pch = 16, cex = 1.25, border = "royalblue4", centroid = FALSE, xy.axis = FALSE, first.point = FALSE)
coo_draw(mean.shape.sex$shp$female, poly = FALSE, pch = 16, cex = 1.25, border = "orange2", centroid = FALSE, xy.axis = FALSE, first.point = FALSE)
legend("topright", lwd=1,
       col=c("royalblue4", "orange2"), legend=c("Male", "Female"))
