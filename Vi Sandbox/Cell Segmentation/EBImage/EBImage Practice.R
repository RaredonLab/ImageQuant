# SET WD
setwd("C:/Users/vl325/OneDrive - Yale University/Lab/ImageQuant/Vi Sandbox/Cell Segmentation")

# PACKAGES
library(EBImage)

# CELL SEGMENTATION EXAMPLE # 
DAPI = readImage("C:/Users/vl325/OneDrive - Yale University/Lab/ImageQuant/Vi Sandbox/Cell Segmentation/SimpleSeg/Test Images/P0/DAPI_Raw.TIF")
RED = readImage("C:/Users/vl325/OneDrive - Yale University/Lab/ImageQuant/Vi Sandbox/Cell Segmentation/SimpleSeg/Test Images/P0/Vimentin_Raw.TIF")
GREEN = readImage("C:/Users/vl325/OneDrive - Yale University/Lab/ImageQuant/Vi Sandbox/Cell Segmentation/SimpleSeg/Test Images/P0/EPCAM_Raw.TIF")
CYAN = readImage("C:/Users/vl325/OneDrive - Yale University/Lab/ImageQuant/Vi Sandbox/Cell Segmentation/SimpleSeg/Test Images/P0/Sox9_Raw.TIF")

#nuc = readImage('/gpfs/gibbs/project/raredon/vl325/CytospinQuant/Cytospin Images/BSL.BDL2/BASC#6_cytospin/HOPX_BSL.BDL_Cytospin_032224_BASC#6 pre_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TR_p00_0_A01f00d0.TIF')
#cel = readImage('/gpfs/gibbs/project/raredon/vl325/CytospinQuant/Cytospin Images/BSL.BDL2/BASC#6_cytospin/DAPI_BSL.BDL_Cytospin_032224_BASC#6 pre_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TR_p00_0_A01f00d1.TIF')

#creating composite image
CELLS = rgbImage(blue = 2*DAPI, red = 1.5*RED, green = 2*GREEN) #assign channels to colors of an overlay (1.5 to make green brighter)
#CELLS = rgbImage(green=1.5*EPCAM, blue=4*DAPI, red=2*CD31) #assign channels to colors of an overlay (1.5 to make green brighter)
display(CELLS, all = TRUE)

#Segment nuclei
Dmask = thresh(DAPI, w=10, h=10, offset=0.0005)
Dmask = opening(Dmask, makeBrush(2, shape='disc'))
Dmask = fillHull(Dmask)
Dmask = bwlabel(Dmask)
display(Dmask, all=TRUE)

#Voronai Tesselation using segmented nuclei as seeds for cytoplasm segmentation
#Red Mask
ctmask = opening(RED>0.0035, makeBrush(3, shape='disc'))
Rmask = propagate(RED, seeds=Dmask, mask=ctmask)
display(ctmask, all=TRUE)
#Green Mask
ctmask2 = opening(GREEN>0.003, makeBrush(7, shape='disc'))
Gmask = propagate(GREEN, seeds=Dmask, mask=ctmask2)
#display(ctmask2, all=TRUE)
#Cyan Mask
ctmask3 = opening(CYAN>0.0025, makeBrush(3, shape='disc'))
Cmask = propagate(CYAN, seeds=Dmask, mask=ctmask3)
#display(ctmask3, all=TRUE)

#display segmentation outlines
segmented = paintObjects(Rmask, CELLS, col='red')
segmented = paintObjects(Dmask, segmented, col='blue')
segmented = paintObjects(Gmask, segmented, col='green')
segmented = paintObjects(Cmask, segmented, col='cyan')
display(segmented, all=TRUE)


writeImage(segmented, "TEST Segmentation.jpeg", quality = 100)
