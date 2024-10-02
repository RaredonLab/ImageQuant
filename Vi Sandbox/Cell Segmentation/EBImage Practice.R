# SET WD
setwd("C:/Users/vl325/OneDrive - Yale University/Lab/ImageQuant/Vi Sandbox/Cell Segmentation")

# PACKAGES
library(EBImage)

# CELL SEGMENTATION EXAMPLE # 
DAPI = readImage('BLUE.tif')
RED = readImage('RED.tif')
GREEN = readImage('GREEN.tif')
#CYAN = readImage('CYAN.tif')

#nuc = readImage('/gpfs/gibbs/project/raredon/vl325/CytospinQuant/Cytospin Images/BSL.BDL2/BASC#6_cytospin/HOPX_BSL.BDL_Cytospin_032224_BASC#6 pre_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TR_p00_0_A01f00d0.TIF')
#cel = readImage('/gpfs/gibbs/project/raredon/vl325/CytospinQuant/Cytospin Images/BSL.BDL2/BASC#6_cytospin/DAPI_BSL.BDL_Cytospin_032224_BASC#6 pre_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TR_p00_0_A01f00d1.TIF')

#creating composite image
CELLS = rgbImage(blue=2*DAPI, red=1.5*RED, green = 2*GREEN) #assign channels to colors of an overlay (1.5 to make green brighter)
#CELLS = rgbImage(green=1.5*EPCAM, blue=4*DAPI, red=2*CD31) #assign channels to colors of an overlay (1.5 to make green brighter)
display(CELLS, all = TRUE)
#segment nuclei
Dmask = thresh(DAPI, w=10, h=10, offset=0.00001)
Dmask = opening(Dmask, makeBrush(5, shape='disc'))
Dmask = fillHull(Dmask)
Dmask = bwlabel(Dmask)
display(Dmask, all=TRUE)
#voronai tesselation using segmented nuclei as seeds for cytoplasm segmentation
ctmask = opening(RED>0.005, makeBrush(7, shape='disc'))
Rmask = propagate(RED, seeds=Dmask, mask=ctmask)
display(ctmask, all=TRUE)

ctmask2 = opening(GREEN>0.00005, makeBrush(7, shape='disc'))
Gmask = propagate(GREEN, seeds=Dmask, mask=ctmask2)
display(ctmask2, all=TRUE)
#display segmentation outlines
segmented = paintObjects(Rmask, CELLS, col='red')
segmented = paintObjects(Dmask, segmented, col='blue')
segmented = paintObjects(Gmask, segmented, col='yellow')
display(segmented, all=TRUE)


writeImage(segmented, "TEST.jpeg", quality = 100)
