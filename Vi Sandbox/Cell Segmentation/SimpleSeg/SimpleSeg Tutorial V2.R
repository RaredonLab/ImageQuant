# SET WD
setwd("Z:/")

# PACKAGES
library(simpleSeg)
library(ggplot2)
library(EBImage)
library(cytomapper)

# Get path to image directory
pathToImages <- "Z:/Raredon_Lab_Personal_Folders/Satoshi/EVOS/2D_BASC_Cytospin_P0+P1_#1-3_052524"

# Get directories of images
imageDirs <- dir(pathToImages, "P0_#3", full.names = TRUE)
names(imageDirs) <- dir(pathToImages, "P0_#3", full.names = FALSE)

# Get files in each directory
files <- files <- lapply(
  imageDirs,
  list.files,
  pattern = "TR",
  full.names = TRUE
)

# Read files with readImage from EBImage
images <- lapply(files, EBImage::readImage, as.is = TRUE)

# Convert to cytoImageList
images <- cytomapper::CytoImageList(images)
mcols(images)$imageID <- names(images)