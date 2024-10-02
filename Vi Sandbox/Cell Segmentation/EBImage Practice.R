# SET WD
setwd("~/project/CytospinQuant/ImageQuant.git/Vi Sandbox/EBImage Method")

# PACKAGES
library(EBImage)

# CELL SEGMENTATION EXAMPLE # 
nuc = readImage('/gpfs/gibbs/project/raredon/vl325/CytospinQuant/ImageQuant.git/Vi Sandbox/EBImage Method/111623_Col IV_Coated_SOX9m_ABCA3.510_Vimentin.690/111623_Col IV_Coated_SOX9m_ABCA3.510_Vimentin.690_TTF1.490_Top Slide_TD_p00_0_A01f00d1.TIF')
cel = readImage('/gpfs/gibbs/project/raredon/vl325/CytospinQuant/ImageQuant.git/Vi Sandbox/EBImage Method/111623_Col IV_Coated_SOX9m_ABCA3.510_Vimentin.690/111623_Col IV_Coated_SOX9m_ABCA3.510_Vimentin.690_TTF1.490_Top Slide_TD_p00_0_A01f00d0.TIF')


nuc = readImage('/gpfs/gibbs/project/raredon/vl325/CytospinQuant/Cytospin Images/BSL.BDL2/BASC#6_cytospin/HOPX_BSL.BDL_Cytospin_032224_BASC#6 pre_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TR_p00_0_A01f00d0.TIF')
cel = readImage('/gpfs/gibbs/project/raredon/vl325/CytospinQuant/Cytospin Images/BSL.BDL2/BASC#6_cytospin/DAPI_BSL.BDL_Cytospin_032224_BASC#6 pre_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TR_p00_0_A01f00d1.TIF')

#creating composite image
cells = rgbImage(green=1.5*cel, blue=nuc) #assign channels to colors of an overlay (1.5 to make green brighter)
display(cells, all = TRUE)
#segment nuclei
nmask = thresh(nuc, w=10, h=10, offset=0.05)
nmask = opening(nmask, makeBrush(5, shape='disc'))
nmask = fillHull(nmask)
nmask = bwlabel(nmask)
display(nmask, all=TRUE)
#voronai tesselation using segmented nuclei as seeds for cytoplasm segmentation
ctmask = opening(cel>0.1, makeBrush(5, shape='disc'))
cmask = propagate(cel, seeds=nmask, mask=ctmask)
display(ctmask, all=TRUE)
#display segmentation outlines
segmented = paintObjects(cmask, cells, col='#ff00ff')
segmented = paintObjects(nmask, segmented, col='#ffff00')
display(segmented, all=TRUE)
