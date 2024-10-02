# SET WD
setwd("~/project/CytospinQuant/ImageQuant.git/Vi Sandbox")

# PACKAGES
require(EBImage)
require(reshape2)
require(ggplot2)
require(Seurat)
require(ggpubr)

# LOAD FUNCTIONS
source("~/project/CytospinQuant/ImageQuant.git/R Functions/ChunkIt.R")
source("~/project/CytospinQuant/ImageQuant.git/R Functions/LoadImageChannels.R")
source("~/project/CytospinQuant/ImageQuant.git/R Functions/MeltIntoDataframe.R")
source("~/project/CytospinQuant/ImageQuant.git/R Functions/NormalizeImages.R")
source("~/project/CytospinQuant/ImageQuant.git/R Functions/PreProcessImage.R")
source("~/project/CytospinQuant/ImageQuant.git/R Functions/QuantFunction.R")
source("~/project/CytospinQuant/ImageQuant.git/R Functions/QuantifyOneLabel.R")


--------------------------------------------------------------------------------
  
  #***PRE-SEEDING***
  # LOAD SINGLE CHANNEL IMAGES
  images.raw <- LoadImageChannels(channel.filepaths = 
                                    c("~/project/CytospinQuant/Cytospin Images/BSL.BDL1/BASC#3_cytospin/DAPI_012424_BDL1 pre_Cytospin_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TR_p00_0_A01f00d1.TIF",
                                      "~/project/CytospinQuant/Cytospin Images/BSL.BDL1/BASC#3_cytospin/HOPX_012424_BDL1 pre_Cytospin_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TR_p00_0_A01f00d0.TIF",
                                      "~/project/CytospinQuant/Cytospin Images/BSL.BDL1/BASC#3_cytospin/VIM_012424_BDL1 pre_Cytospin_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TR_p00_0_A01f00d2.TIF",
                                      "~/project/CytospinQuant/Cytospin Images/BSL.BDL1/BASC#3_cytospin/SOX9_012424_BDL1 pre_Cytospin_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TR_p00_0_A01f00d3.TIF"
                                    ),
                                  
                                  channel.names = 
                                    c('DAPI','HOPX','Vim','Sox9')
  )
  
  # display(images.raw$DAPI)
  # display(images.raw$HOPX)
  # display(images.raw$Vim)
  # display(images.raw$Sox9)

  # Normalize each channel
  images.norm <- NormalizeImages(images.raw)
  # display(images.norm$DAPI)
  # display(images.norm$Sox9)
  
  # PreProcess image channels
  images.processed <- PreProcessImage(images.norm = images.norm,
                                      box.dim = 1000,
                                      offset = 0.01) # These parameters will make a huge difference in your outputs and should be tuned + checked
  
  # display(images.processed$DAPI)
  # display(images.processed$HOPX)
  # display(images.processed$Vim)
  # display(images.processed$Sox9)
  
  # Melt into dataframe
  image.object <- MeltIntoDataframe(images.processed)
  
  # Define number of pixels within each chunk
  num.pixels.per.chunk = 400000 # experiment with chunk size to yield a good number of samples
  num.pix <- nrow(image.object) # total number of pixels
  n.chunks <- floor(num.pix/num.pixels.per.chunk)
  message(paste(num.pixels.per.chunk,'pixels per chunk will yield',n.chunks,'total chunks'))
  
  chunked.data <- ChunkIt(image.object = image.object,
                          n.chunks = n.chunks)
  
  # Build the entire dataset needed for plotting
  output.data <- data.frame('All Sox9-Positive' = QuantifyOneLabel(chunked.data,
                                                                    label.name = 'All Sox9-Positive',
                                                                    base.criteria.positive = c('DAPI'),
                                                                    label.criteria.positive = c('Sox9'),
                                                                    label.criteria.negative = NULL),
                            'Sox9-Positive Only' = QuantifyOneLabel(chunked.data,
                                                                    label.name = 'Sox9-Positive Only',
                                                                    base.criteria.positive = c('DAPI'),
                                                                    label.criteria.positive = c('Sox9'),
                                                                    label.criteria.negative = c('Vim','HOPX')),
                            'ATI-like Transitional' = QuantifyOneLabel(chunked.data,
                                                                       label.name = 'ATI-like Transitional',
                                                                       base.criteria.positive = c('DAPI'),
                                                                       label.criteria.positive = c('Sox9','HOPX'),
                                                                       label.criteria.negative = c('Vim')),
                            'ATI-like End Stage' = QuantifyOneLabel(chunked.data,
                                                                    label.name = 'ATI-like End Stage',
                                                                    base.criteria.positive = c('DAPI'),
                                                                    label.criteria.positive = c('HOPX'),
                                                                    label.criteria.negative = c('Vim','Sox9')),
                            'EMT-like' = QuantifyOneLabel(chunked.data,
                                                          label.name = 'EMT-like',
                                                          base.criteria.positive = c('DAPI'),
                                                          label.criteria.positive = c('Sox9','Vim'),
                                                          label.criteria.negative = NULL),
                            # label.criteria.negative = c('HOPX'))
                            'Mesenchyme' = QuantifyOneLabel(chunked.data,
                                                            label.name = 'Mesenchyme',
                                                            base.criteria.positive = c('DAPI'),
                                                            label.criteria.positive = c('Vim'),
                                                            label.criteria.negative = c('Sox9','HOPX'))
                            
  )
  
  
  # Label with Sample
  output.data$Sample <- 'Pre-Seeding 1' #*************UPDATE NUMBER*************
  
  # Stash so we don't overwrite
  PreSeed1Data <- melt(output.data) #*************UPDATE NUMBER*************
  
  
  
  
  --------------------------------------------------------------------------------
    
    #***BSL***
    # LOAD SINGLE CHANNEL IMAGES
    images.raw <- LoadImageChannels(channel.filepaths = 
                                      c("~/project/CytospinQuant/Cytospin Images/BSL.BDL1/BSL1_cytospin/DAPI_01.25.24_BSL1_Post-Culture_Cytsopin_HOP_690_Vimentin_510_SOX9_490_20x_Top Slide_TR_p00_0_A01f00d1.TIF",
                                        "~/project/CytospinQuant/Cytospin Images/BSL.BDL1/BSL1_cytospin/HOPX_01.25.24_BSL1_Post-Culture_Cytsopin_HOP_690_Vimentin_510_SOX9_490_20x_Top Slide_TR_p00_0_A01f00d0.TIF",
                                        "~/project/CytospinQuant/Cytospin Images/BSL.BDL1/BSL1_cytospin/VIM_01.25.24_BSL1_Post-Culture_Cytsopin_HOP_690_Vimentin_510_SOX9_490_20x_Top Slide_TR_p00_0_A01f00d2.TIF",
                                        "~/project/CytospinQuant/Cytospin Images/BSL.BDL1/BSL1_cytospin/SOX9_01.25.24_BSL1_Post-Culture_Cytsopin_HOP_690_Vimentin_510_SOX9_490_20x_Top Slide_TR_p00_0_A01f00d3.TIF"
                                      ),
                                    
                                    channel.names = 
                                      c('DAPI','HOPX','Vim','Sox9')
    )
  
  # display(images.raw$DAPI)
  # display(images.raw$HOPX)
  # display(images.raw$Vim)
  # display(images.raw$Sox9)
  
  # Normalize each channel
  images.norm <- NormalizeImages(images.raw)
  # display(images.norm$DAPI)
  # display(images.norm$Sox9)
  
  # PreProcess image channels
  images.processed <- PreProcessImage(images.norm = images.norm,
                                      box.dim = 1000,
                                      offset = 0.01) # These parameters will make a huge difference in your outputs and should be tuned + checked
  
  # display(images.processed$DAPI)
  # display(images.processed$HOPX)
  # display(images.processed$Vim)
  # display(images.processed$Sox9)
  
  # Melt into dataframe
  image.object <- MeltIntoDataframe(images.processed)
  
  # Define number of pixels within each chunk
  num.pixels.per.chunk = 600000 # experiment with chunk size to yield a good number of samples
  num.pix <- nrow(image.object) # total number of pixels
  n.chunks <- floor(num.pix/num.pixels.per.chunk)
  message(paste(num.pixels.per.chunk,'pixels per chunk will yield',n.chunks,'total chunks'))
  
  chunked.data <- ChunkIt(image.object = image.object,
                          n.chunks = n.chunks)
  
  # Build the entire dataset needed for plotting
  output.data <- data.frame('All Sox9-Positive' = QuantifyOneLabel(chunked.data,
                                                                   label.name = 'All Sox9-Positive',
                                                                   base.criteria.positive = c('DAPI'),
                                                                   label.criteria.positive = c('Sox9'),
                                                                   label.criteria.negative = NULL),
                            'Sox9-Positive Only' = QuantifyOneLabel(chunked.data,
                                                                    label.name = 'Sox9-Positive Only',
                                                                    base.criteria.positive = c('DAPI'),
                                                                    label.criteria.positive = c('Sox9'),
                                                                    label.criteria.negative = c('Vim','HOPX')),
                            'ATI-like Transitional' = QuantifyOneLabel(chunked.data,
                                                                       label.name = 'ATI-like Transitional',
                                                                       base.criteria.positive = c('DAPI'),
                                                                       label.criteria.positive = c('Sox9','HOPX'),
                                                                       label.criteria.negative = c('Vim')),
                            'ATI-like End Stage' = QuantifyOneLabel(chunked.data,
                                                                    label.name = 'ATI-like End Stage',
                                                                    base.criteria.positive = c('DAPI'),
                                                                    label.criteria.positive = c('HOPX'),
                                                                    label.criteria.negative = c('Vim','Sox9')),
                            'EMT-like' = QuantifyOneLabel(chunked.data,
                                                          label.name = 'EMT-like',
                                                          base.criteria.positive = c('DAPI'),
                                                          label.criteria.positive = c('Sox9','Vim'),
                                                          label.criteria.negative = NULL),
                            # label.criteria.negative = c('HOPX'))
                            'Mesenchyme' = QuantifyOneLabel(chunked.data,
                                                            label.name = 'Mesenchyme',
                                                            base.criteria.positive = c('DAPI'),
                                                            label.criteria.positive = c('Vim'),
                                                            label.criteria.negative = c('Sox9','HOPX'))
                            
  )
  
  
  # Label with Sample
  output.data$Sample <- 'BSL1' #*************UPDATE NUMBER*************
  
  # Stash so we don't overwrite
  BSL1Data <- melt(output.data) #*************UPDATE NUMBER*************
  
  
  
  
  --------------------------------------------------------------------------------
    
    #***BDL***
    # LOAD SINGLE CHANNEL IMAGES
    images.raw <- LoadImageChannels(channel.filepaths = 
                                      c("~/project/CytospinQuant/Cytospin Images/BSL.BDL1/BDL1_cytospin/DAPI_012424_BDL1 post_Cytospin_HOP.488_Vimentin.555_SOX9.647_ver2_Bottom Slide_TR_p00_0_A01f00d1.TIF",
                                        "~/project/CytospinQuant/Cytospin Images/BSL.BDL1/BDL1_cytospin/HOP_012424_BDL1 post_Cytospin_HOP.488_Vimentin.555_SOX9.647_ver2_Bottom Slide_TR_p00_0_A01f00d0.TIF",
                                        "~/project/CytospinQuant/Cytospin Images/BSL.BDL1/BDL1_cytospin/VIM_012424_BDL1 post_Cytospin_HOP.488_Vimentin.555_SOX9.647_ver2_Bottom Slide_TR_p00_0_A01f00d2.TIF",
                                        "~/project/CytospinQuant/Cytospin Images/BSL.BDL1/BDL1_cytospin/SOX9_012424_BDL1 post_Cytospin_HOP.488_Vimentin.555_SOX9.647_ver2_Bottom Slide_TR_p00_0_A01f00d3.TIF"
                                      ),
                                    
                                    channel.names = 
                                      c('DAPI','HOPX','Vim','Sox9')
    )
  
  # display(images.raw$DAPI)
  
  # Normalize each channel
  images.norm <- NormalizeImages(images.raw)
  # display(images.norm$DAPI)
  # display(images.norm$Sox9)
  
  # PreProcess image channels
  images.processed <- PreProcessImage(images.norm = images.norm,
                                      box.dim = 1000,
                                      offset = 0.01) # These parameters will make a huge difference in your outputs and should be tuned + checked
  
  # display(images.processed$DAPI)
  # display(images.processed$HOPX)
  # display(images.processed$Vim)
  # display(images.processed$Sox9)
  
  # Melt into dataframe
  image.object <- MeltIntoDataframe(images.processed)
  
  # Define number of pixels within each chunk
  num.pixels.per.chunk = 440000 # experiment with chunk size to yield a good number of samples
  num.pix <- nrow(image.object) # total number of pixels
  n.chunks <- floor(num.pix/num.pixels.per.chunk)
  message(paste(num.pixels.per.chunk,'pixels per chunk will yield',n.chunks,'total chunks'))
  
  chunked.data <- ChunkIt(image.object = image.object,
                          n.chunks = n.chunks)
  
  # Build the entire dataset needed for plotting
  output.data <- data.frame('All Sox9-Positive' = QuantifyOneLabel(chunked.data,
                                                                   label.name = 'All Sox9-Positive',
                                                                   base.criteria.positive = c('DAPI'),
                                                                   label.criteria.positive = c('Sox9'),
                                                                   label.criteria.negative = NULL),
                            'Sox9-Positive Only' = QuantifyOneLabel(chunked.data,
                                                                    label.name = 'Sox9-Positive Only',
                                                                    base.criteria.positive = c('DAPI'),
                                                                    label.criteria.positive = c('Sox9'),
                                                                    label.criteria.negative = c('Vim','HOPX')),
                            'ATI-like Transitional' = QuantifyOneLabel(chunked.data,
                                                                       label.name = 'ATI-like Transitional',
                                                                       base.criteria.positive = c('DAPI'),
                                                                       label.criteria.positive = c('Sox9','HOPX'),
                                                                       label.criteria.negative = c('Vim')),
                            'ATI-like End Stage' = QuantifyOneLabel(chunked.data,
                                                                    label.name = 'ATI-like End Stage',
                                                                    base.criteria.positive = c('DAPI'),
                                                                    label.criteria.positive = c('HOPX'),
                                                                    label.criteria.negative = c('Vim','Sox9')),
                            'EMT-like' = QuantifyOneLabel(chunked.data,
                                                          label.name = 'EMT-like',
                                                          base.criteria.positive = c('DAPI'),
                                                          label.criteria.positive = c('Sox9','Vim'),
                                                          label.criteria.negative = NULL),
                            # label.criteria.negative = c('HOPX'))
                            'Mesenchyme' = QuantifyOneLabel(chunked.data,
                                                            label.name = 'Mesenchyme',
                                                            base.criteria.positive = c('DAPI'),
                                                            label.criteria.positive = c('Vim'),
                                                            label.criteria.negative = c('Sox9','HOPX'))
                            
  )
  
  
  # Label with Sample
  output.data$Sample <- 'BDL1' #*************UPDATE NUMBER*************
  
  # Stash so we don't overwrite
  BDL1Data <- melt(output.data) #*************UPDATE NUMBER*************
  
  
  
  
  --------------------------------------------------------------------------------
    
    #********* Plotting Exploration ********#
    
    # Bind two datasets together for plotting
  data <- rbind(PreSeed1Data,BSL1Data,BDL1Data) #*************UPDATE NUMBER*************   
  data_sorted <- data
  data_sorted$Sample <- factor(data_sorted$Sample,
                               c("Pre-Seeding 1","BSL1","BDL1")) #*************UPDATE NUMBER*************
  
  png('BSL1.BDL1.BASC3 Violin Plot.png',width = 15,height = 8,units='in',res=1000)
  ggplot(data_sorted, aes(x = variable,y=value,fill=Sample))+
    geom_violin(
      mapping = NULL,
      data = data,
      stat = "ydensity",
      position = "dodge",
      draw_quantiles = NULL,
      trim = FALSE,
      scale = "count",
      na.rm = FALSE,
      orientation = NA,
      show.legend = NA,
      inherit.aes = TRUE,
    )+
  stat_ydensity(
    mapping = NULL,
    data = NULL,
    geom = "violin",
    position = "dodge",
    bw = "nrd0",
    adjust = 1,
    kernel = "gaussian",
    trim = FALSE,
    scale = "width",
    na.rm = FALSE,
    orientation = NA,
    show.legend = NA,
    inherit.aes = TRUE,
  )+
    geom_point(position = position_jitterdodge(dodge.width = 0.9,jitter.width=0.25),size=0.08,color='black')+
    theme_classic(base_size=20)+
    theme(axis.text.x = element_text(size=12))+
    ggtitle('Effect of Dynamic Ventilation on BASC Character in 3D Bioengineered Rat Lungs')+
    ylab('Fraction of DAPI+ Pixels with Character')+
    xlab('Character')+
   #compare_means(Sample ~ value,  data = data_sorted, ref.group = data_sorted$Sample("Pre-Seeding 3"),
  #                method = "t.test")+
    scale_fill_discrete(breaks=c("Pre-Seeding 1","BSL1","BDL1"))    #*************UPDATE NUMBER************* 
      #scale_fill_manual(values = c('lightpink','purple','cyan'))+

    #scale_color_manual(values = c('#93B7BE','#F19A3E','#3D3B8E','#E072A4','lightgreen'))
  
  dev.off()
