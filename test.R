# Elliptic Fourier Analysis

## Load packages + data

# load packages
library(here)
library(wesanderson)
library(Momocs)

# read images + attribute data
jpg.list <- list.files(here("/jpegs"), full.names = TRUE)
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

# oriented based upon side with highest shoulder height
norm.outlines$coo[which(att.data$sh.up=="dn")] <- 
  lapply(att.data$coo[which(att.data$sh.up=="dn")], coo_flipy)

norm.outlines.out$coo[which(att.data$sh.up=="dn")] <- 
  lapply(att.data$coo[which(att.data$sh.up=="dn")], coo_rev)

# outline pile
pile(norm.outlines)

# mosaic of individual specimens from the different sites
mosaic(norm.outlines, ~site)

# calibrate how many harmonics needed
calibrate_harmonicpower_efourier(norm.outlines, nb.h = 30)

# 11 harmonics needed to capture 99 percent of variation
calibrate_reconstructions_efourier(norm.outlines, range = 1:11)

# generate efa outlines with 11 harmonics
efa.outlines <- efourier(norm.outlines, nb.h = 11, norm = TRUE)

# use efa.outlines for pca
pca.outlines <- PCA(efa.outlines)

# pca 
scree_plot(pca.outlines)

# plot pca by site
plot_PCA(pca.outlines, morphospace_position = "range",
         ~site, zoom = 1.25, palette = pal_qual_Dark2)

# contribution of each pc
boxplot(pca.outlines, ~site, nax = 1:5, palette = pal_qual_Dark2)

# mean shape + 2sd for the first 10 pcs
PCcontrib(pca.outlines, nax = 1:10)

# manova + manova_pw
MANOVA(pca.outlines, 'site')
MANOVA_PW(pca.outlines, 'site')

# mean shapes
ms.1 <- MSHAPES(efa.outlines, ~site)
plot_MSHAPES(ms.1, size = 0.75)
