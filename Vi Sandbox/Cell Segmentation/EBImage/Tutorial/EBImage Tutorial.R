# SET WD
  setwd("~/project/CytospinQuant/ImageQuant.git/Vi Sandbox/EBImage Method")

# PACKAGES
  install.packages("BiocManager")
  BiocManager::install("EBImage")
  library(EBImage)

#LOAD IMAGE(S)
  f = system.file("images", "sample.png", package="EBImage")
  img = readImage(f)
  display(img, method="browser") #JavaScript Viewer

#EDIT/ANNOTATE
  display(img, method="raster") #R plot Viewer to edit/add text
  text(x = 20, y = 20, label = "Parrots", adj = c(0,1), col = "orange", cex = 2) #Adding text

#SAVE
  filename = "parrots.jpg"
  dev.print(jpeg, filename = filename , width = dim(img)[1], height = dim(img)[2])
  file.info(filename)$size #Verify file size in (MB)

#EXAMPLE COLOR PHOTO
  imgcol = readImage(system.file("images", "sample-color.png", package="EBImage"))
  display(imgcol)

#LOAD EXAMPLE MULTI-FRAME CELL FILES
  nuc = readImage(system.file("images", "nuclei.tif", package="EBImage"))
  display(nuc, method = "raster", all = TRUE) #all = TRUE displays all frames in a grid

#CONVERT FILE TYPE
  writeImage(imgcol, "sample.jpeg", quality = 85)
  display(imgcol, method = "raster", all=TRUE)

#ANALYZING THE IMAGE DATA FRAME
  str(img) #recall data about the data composing the image by pixel intensities
  dim(img) #recall the size of the image in pixel dimensions
  imageData(img)[1:3, 1:6] #access individual pixel data
  as.array(img) #convert image to array
  hist(img) #produce histogram of pixel intensities
  img #summary
  print(img, short=TRUE) #short summary
  numberOfFrames(imgcol, type = "render")
  numberOfFrames(imgcol, type = "total")

#COLOR MANAGEMENT
  ##color to grayscale
  colorMode(imgcol) = Grayscale #changes a color image to grayscale frames of each color channel, does not change file properties
  display(imgcol, all=TRUE)
  #creating a color matrix image
  colorMat = matrix(rep(c("red","green", "#0000ff"), 25), 5, 5)
  colorImg = Image(colorMat)
  colorImg
  display(colorImg, interpolate=FALSE)
  colorMode(colorImg) = Grayscale

#IMAGE MANIPULATION
  ##inverting
  img_neg = max(img)-img #inverting the pixel values
  display(img_neg)
  ##brightness adjustment
  img_comb = combine(
    img,
    img + 0.3, #+0.3 pixel intensity
    img * 2, #double pixel intensity
    img ^ 0.5 #sqrt pixel intensity
  )
  display(img_comb, all = TRUE)
  ##crop and threshold
  img_crop = img[366:749, 58:441] #crop to specific pixels
  img_thresh = img_crop > .5 #binarizes pixel values according to the set boundary condition of 0.5 intensity
  display(img_thresh)
  ##transpose image
  img_t = transpose(img)
  display(img)
  display(img_t)
  
#LINEAR FILTERING
  #gaussian smoothing/low pass filtering, manual
  w = makeBrush(size = 31, shape = 'gaussian', sigma = 5) #filter brush parameters
  plot(w[(nrow(w)+1)/2, ], ylab = "w", xlab = "", cex = 0.7) #visualization of brush parameters
  img_flo = filter2(img, w) #filter according to the brush w
  display(img_flo)
  #guassian smoothing, automatic
  nuc_gblur = gblur(nuc, sigma = 5)
  display(nuc_gblur, all=TRUE )
  #laplacian sharpening/high pass filtering
  fhi = matrix(1, nrow = 3, ncol = 3)
  fhi[2, 2] = -8
  img_fhi = filter2(img, fhi)
  display(img_fhi)
  #noise reduction
    #making a noisy image for demonstration
    l = length(img)
    n = l/10
    pixels = sample(l, n)
    img_noisy = img
    img_noisy[pixels] = runif(n, min=0, max=1)
    display(img_noisy)
  ##applying a median filter which replaces pixels with a median of their neighbors
  img_median = medianFilter(img_noisy, 1)
  display(img_median)

#THRESHOLDING
  #adaptive manual, thresholds relative to the background area, good for images with differences in background lighting
  disc = makeBrush(31, "disc")
  disc = disc / sum(disc)
  offset = 0.05
  nuc_bg = filter2( nuc, disc )
  nuc_th = nuc > nuc_bg + offset
  display(nuc_th, all=TRUE)
  #adaptive using thresh
  display( thresh(nuc, w=15, h=15, offset=0.05), all=TRUE )
  
#IMAGE SEGMENTATION
  #non-touching objects using bwlabel
  #*making logo*#
  shapes = readImage(system.file('images', 'shapes.png', package='EBImage'))
  logo = shapes[110:512,1:130]
  display(logo)
  kern = makeBrush(5, shape='diamond')
  display(kern, interpolate=FALSE)
  logo_erode= erode(logo, kern)
  logo_dilate = dilate(logo, kern)
  display(combine(logo_erode, logo_dilate), all=TRUE)
  logo_label = bwlabel(logo) #segment groups of touching pixels into integer-labelled groups
  table(logo_label)
  display( normalize(logo_label) ) #normalize to 0-1 pixel values for display
  display( colorLabels(logo_label) ) #display segmented groups in distinct colors
  #touching objects using watershed
  nmask = watershed( distmap(nuc_th), 2 )
  display(colorLabels(nmask), all=TRUE)
  display(nuc)
  #voronoi tesselation to use seed points to define surrounding regions
  voronoiExamp = propagate(seeds = nmask, x = nmask, lambda = 100)
  voronoiPaint = colorLabels (voronoiExamp)
  display(voronoiPaint)

#FLOODFILL
  rgblogo = toRGB(logo)
  points = rbind(c(50, 50), c(100, 50), c(150, 50))
  colors = c("red", "green", "blue")
  rgblogo = floodFill(rgblogo, points, colors)
  display( rgblogo )
  
  
  
# CELL SEGMENTATION EXAMPLE # 
  nuc = readImage(system.file('images', 'nuclei.tif', package='EBImage'))
  cel = readImage(system.file('images', 'cells.tif', package='EBImage'))
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
  