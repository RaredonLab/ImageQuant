# SET WD
setwd("C:/Users/vl325/OneDrive - Yale University/Lab/ImageQuant/Vi Sandbox/Cell Segmentation/SimpleSeg")

# PACKAGES
library(simpleSeg)
library(ggplot2)
library(EBImage)
library(cytomapper)

### LOAD IMAGES ###

# Get directories of images
pathToImages <- "C:/Users/vl325/OneDrive - Yale University/Lab/ImageQuant/Vi Sandbox/Cell Segmentation/SimpleSeg/Test Images"
imageDirs <- dir(pathToImages, "P", full.names = TRUE)
names(imageDirs) <- dir(pathToImages, "P", full.names = FALSE)

# Get files in each directory
files <- files <- list.files(imageDirs, full.names = TRUE)

# Read files with readImage from EBImage
images <- lapply(files, EBImage::readImage, as.is = TRUE)

# Convert to cytoImageList
images <- cytomapper::CytoImageList(images)
mcols(images)$imageID <- names(images)

### SEGMENTATION ###
masks <- simpleSeg::simpleSeg(images,
                              nucleus = "DAPI",
                              transform = "sqrt"
)
# Visualise segmentation performance one way.
EBImage::display(colorLabels(masks[[1]]))