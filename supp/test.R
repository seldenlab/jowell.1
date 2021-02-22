# load packages

devtools::install_github("MomX/Momocs")
library(here)
library(wesanderson)
library(ggplot2)
library(Momocs)

# read images + attribute data
jpg.list <- list.files(here("./jpegs"), full.names = TRUE)
att.data <- read.csv("att.data.csv", header = TRUE, as.is = TRUE)

# attributes to factors
att.data$site <- as.factor(att.data$site)

# generate outlines
outlines <- jpg.list %>%
  import_jpg()

# add attributes
data.out <- Out(outlines, 
                fac = att.data)

# center, scale, align, and rotate specimens
norm.outlines <- data.out %>% 
  coo_center %>%
  coo_scale() %>%
  coo_align() %>% 
  coo_rotate()

## 41an13
an13 <- filter(norm.outlines, 
               site %in% "41AN13")

an13 <- an13 %>% 
  coo_scale() %>%
  coo_align() %>%
  coo_rotate() %>% 
  coo_center()

## 41hs261
hs261 <- filter(norm.outlines, 
                site %in% "41HS261")

hs261 <- hs261 %>% 
  coo_scale() %>%
  coo_align() %>%
  coo_rotate() %>% 
  coo_center()

# render figure
#par(mfrow=c(2, 1))
stack(an13, title = "41AN13", xy.axis = TRUE, centroid = FALSE, points =  FALSE)
stack(hs261, title = "41HS261", xy.axis = TRUE, centroid = FALSE, points = FALSE)

# calibrate how many harmonics needed
calibrate_harmonicpower_efourier(norm.outlines, nb.h = 30)

# 9 harmonics needed to capture 99 percent of variation
calibrate_reconstructions_efourier(norm.outlines, range = 1:9)

# generate efa outlines with 9 harmonics
efa.outlines <- efourier(norm.outlines, nb.h = 9, norm = TRUE)

# use efa.outlines for pca
pca.outlines <- PCA(efa.outlines)

# pca 
scree_plot(pca.outlines)
##############################################################################
# plot pca by site
plot_PCA(pca.outlines,
         morphospace_position = "range_axes",
         ~site, 
         zoom = 1.5)

# plot pca by context
plot(pca.outlines,
     pos.shp = "range_axes",
     ~site,
     chull = TRUE,
     morphospace = TRUE,
     labelsgroups = TRUE,
     cex.labelsgroups = 0.5,
     rect.labelsgroups = TRUE,
     rug = TRUE,
     grid = TRUE,
     zoom = 0.95)

# contribution of each pc
boxplot(pca.outlines, ~site, nax = 1:5)

# mean shape + 2sd for the first 10 pcs
PCcontrib(pca.outlines, nax = 1:5)

# manova
# does shape differ between sites?
MANOVA(pca.outlines, 'site')

# mean shapes
ms.1 <- MSHAPES(efa.outlines, ~site)
plot_MSHAPES(ms.1, size = 0.75)
