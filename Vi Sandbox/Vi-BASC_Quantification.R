# SET WD
setwd("~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/Vi Sandbox")

# PACKAGES
require(EBImage)
require(reshape2)
require(ggplot2)
require(Seurat)
require(ggpubr)

# LOAD FUNCTIONS
source("~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/R Functions/ChunkIt.R")
source("~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/R Functions/LoadImageChannels.R")
source("~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/R Functions/MeltIntoDataframe.R")
source("~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/R Functions/NormalizeImages.R")
source("~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/R Functions/PreProcessImage.R")
source("~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/R Functions/QuantFunction.R")
source("~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/R Functions/QuantifyOneLabel.R")


--------------------------------------------------------------------------------

#***PRE-SEEDING***
    # LOAD SINGLE CHANNEL IMAGES
    images.raw <- LoadImageChannels(channel.filepaths = 
                                      c("~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/Vi Sandbox/Cytospin Images/BSL.BDL3/BSL.BDL3.Pre.Cytospin/DAPI_BSL.BDL_Cytospin_032224_BASC#5 pre_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TD_p00_0_A01f00d1.TIF",
                                        "~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/Vi Sandbox/Cytospin Images/BSL.BDL3/BSL.BDL3.Pre.Cytospin/HOP_BSL.BDL_Cytospin_032224_BASC#5 pre_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TD_p00_0_A01f00d3.TIF",
                                        "~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/Vi Sandbox/Cytospin Images/BSL.BDL3/BSL.BDL3.Pre.Cytospin/VIM_BSL.BDL_Cytospin_032224_BASC#5 pre_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TD_p00_0_A01f00d2.TIF",
                                        "~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/Vi Sandbox/Cytospin Images/BSL.BDL3/BSL.BDL3.Pre.Cytospin/SOX9_BSL.BDL_Cytospin_032224_BASC#5 pre_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TR_p00_0_A01f00d3.TIF"
                                        ),
                                                      
                                      channel.names = 
                                        c('DAPI','HOP','Vim','Sox9')
                                    )
    
    display(images.raw$DAPI)
    
    # Normalize each channel
    images.norm <- NormalizeImages(images.raw)
    # display(images.norm$DAPI)
    # display(images.norm$Sox9)
    
    # PreProcess image channels
    images.processed <- PreProcessImage(images.norm = images.norm,
                                        box.dim = 1000,
                                        offset = 0.01) # These parameters will make a huge difference in your outputs and should be tuned + checked
    
    # display(images.processed$DAPI)
    # display(images.processed$HOP)
    # display(images.processed$Vim)
    # display(images.processed$Sox9)
    
    # Melt into dataframe
    image.object <- MeltIntoDataframe(images.processed)
    
    # Define number of pixels within each chunk
    num.pixels.per.chunk = 500000 # experiment with chunk size to yield a good number of samples
    num.pix <- nrow(image.object) # total number of pixels
    n.chunks <- floor(num.pix/num.pixels.per.chunk)
    message(paste(num.pixels.per.chunk,'pixels per chunk will yield',n.chunks,'total chunks'))
    
    chunked.data <- ChunkIt(image.object = image.object,
                            n.chunks = n.chunks)
    
    # Build the entire dataset needed for plotting
    output.data <- data.frame('ATI-like Transitional' = QuantifyOneLabel(chunked.data,
                                                                         label.name = 'ATI-like Transitional',
                                                                         base.criteria.positive = c('DAPI'),
                                                                         label.criteria.positive = c('Sox9','HOP'),
                                                                         label.criteria.negative = c('Vim')),
                              'ATI-like End Stage' = QuantifyOneLabel(chunked.data,
                                                                      label.name = 'ATI-like End Stage',
                                                                      base.criteria.positive = c('DAPI'),
                                                                      label.criteria.positive = c('HOP'),
                                                                      label.criteria.negative = c('Vim','Sox9')),
                              'Sox9-Positive Only' = QuantifyOneLabel(chunked.data,
                                                                      label.name = 'Sox9-Positive Only',
                                                                      base.criteria.positive = c('DAPI'),
                                                                      label.criteria.positive = c('Sox9'),
                                                                      label.criteria.negative = c('Vim','HOP')),
                              'Mesenchyme' = QuantifyOneLabel(chunked.data,
                                                              label.name = 'Mesenchyme',
                                                              base.criteria.positive = c('DAPI'),
                                                              label.criteria.positive = c('Vim'),
                                                              label.criteria.negative = c('Sox9','HOP')),
                              'EMT-like' = QuantifyOneLabel(chunked.data,
                                                            label.name = 'EMT-like',
                                                            base.criteria.positive = c('DAPI'),
                                                            label.criteria.positive = c('Sox9','Vim'),
                                                            # label.criteria.negative = c('HOP'))
                              )
    )
    
    
    
    # Label with condition
    output.data$Condition <- 'Pre-Seeding 3' #*************UPDATE NUMBER*************
    
    # Stash so we don't overwrite
    PreSeed3Data <- melt(output.data) #*************UPDATE NUMBER*************




--------------------------------------------------------------------------------

#***BSL***
    # LOAD SINGLE CHANNEL IMAGES
    images.raw <- LoadImageChannels(channel.filepaths = 
                                      c("~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/Vi Sandbox/Cytospin Images/BSL.BDL3/BSL3.Post.Cytospin/DAPI_BSL.BDL_Cytospin_032224_BSL3_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TD_p00_0_A01f00d1.TIF",
                                        "~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/Vi Sandbox/Cytospin Images/BSL.BDL3/BSL3.Post.Cytospin/HOP_BSL.BDL_Cytospin_032224_BSL3_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TD_p00_0_A01f00d0.TIF",
                                        "~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/Vi Sandbox/Cytospin Images/BSL.BDL3/BSL3.Post.Cytospin/VIM_BSL.BDL_Cytospin_032224_BSL3_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TD_p00_0_A01f00d2.TIF",
                                        "~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/Vi Sandbox/Cytospin Images/BSL.BDL3/BSL3.Post.Cytospin/SOX9_BSL.BDL_Cytospin_032224_BSL3_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TD_p00_0_A01f00d3.TIF"
                                        ),
                                        
                                      channel.names = 
                                        c('DAPI','HOP','Vim','Sox9')
                                    )
    
      # display(images.raw$DAPI)
      # display(images.raw$HOP)
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
      # display(images.processed$HOP)
      # display(images.processed$Vim)
      # display(images.processed$Sox9)
        
    # Melt into dataframe
    image.object <- MeltIntoDataframe(images.processed)
    
    # Define number of pixels within each chunk
    num.pixels.per.chunk = 500000 # experiment with chunk size to yield a good number of samples
    num.pix <- nrow(image.object) # total number of pixels
    n.chunks <- floor(num.pix/num.pixels.per.chunk)
    message(paste(num.pixels.per.chunk,'pixels per chunk will yield',n.chunks,'total chunks'))
    
    chunked.data <- ChunkIt(image.object = image.object,
                            n.chunks = n.chunks)
    
    # Build the entire dataset needed for plotting
    output.data <- data.frame('ATI-like Transitional' = QuantifyOneLabel(chunked.data,
                                                          label.name = 'ATI-like Transitional',
                                                          base.criteria.positive = c('DAPI'),
                                                          label.criteria.positive = c('Sox9','HOP'),
                                                          label.criteria.negative = c('Vim')),
                              'ATI-like End Stage' = QuantifyOneLabel(chunked.data,
                                                          label.name = 'ATI-like End Stage',
                                                          base.criteria.positive = c('DAPI'),
                                                          label.criteria.positive = c('HOP'),
                                                          label.criteria.negative = c('Vim','Sox9')),
                              'Sox9-Positive Only' = QuantifyOneLabel(chunked.data,
                                                          label.name = 'Sox9-Positive Only',
                                                          base.criteria.positive = c('DAPI'),
                                                          label.criteria.positive = c('Sox9'),
                                                          label.criteria.negative = c('Vim','HOP')),
                              'Mesenchyme' = QuantifyOneLabel(chunked.data,
                                                          label.name = 'Mesenchyme',
                                                          base.criteria.positive = c('DAPI'),
                                                          label.criteria.positive = c('Vim'),
                                                          label.criteria.negative = c('Sox9','HOP')),
                              'EMT-like' = QuantifyOneLabel(chunked.data,
                                                          label.name = 'EMT-like',
                                                          base.criteria.positive = c('DAPI'),
                                                          label.criteria.positive = c('Sox9','Vim'),
                                                        # label.criteria.negative = c('HOP'))
                              )
    )
    
    
    # Label with condition
    output.data$Condition <- 'BSL3' #*************UPDATE NUMBER*************
    
    # Stash so we don't overwrite
    BSL3Data <- melt(output.data) #*************UPDATE NUMBER*************




--------------------------------------------------------------------------------

#***BDL***
  # LOAD SINGLE CHANNEL IMAGES
  images.raw <- LoadImageChannels(channel.filepaths = 
                                    c("~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/Vi Sandbox/Cytospin Images/BSL.BDL3/BDL3.Post.Cytospin/DAPI_BSL.BDL_Cytospin_032224_BDL3_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TD_p00_0_A01f00d1.TIF",
                                      "~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/Vi Sandbox/Cytospin Images/BSL.BDL3/BDL3.Post.Cytospin/HOP_BSL.BDL_Cytospin_032224_BDL3_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TD_p00_0_A01f00d0.TIF",
                                      "~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/Vi Sandbox/Cytospin Images/BSL.BDL3/BDL3.Post.Cytospin/VIM_BSL.BDL_Cytospin_032224_BDL3_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TD_p00_0_A01f00d2.TIF",
                                      "~/project/BSL.BDL3.CytospinQuant/ImageQuant.git/Vi Sandbox/Cytospin Images/BSL.BDL3/BDL3.Post.Cytospin/SOX9_BSL.BDL_Cytospin_032224_BDL3_HOP.488_Vimentin.555_SOX9.647_Bottom Slide_TD_p00_0_A01f00d3.TIF"
                                    ),
                                  
                                  channel.names = 
                                    c('DAPI','HOP','Vim','Sox9')
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
    # display(images.processed$HOP)
    # display(images.processed$Vim)
    # display(images.processed$Sox9)
  
  # Melt into dataframe
  image.object <- MeltIntoDataframe(images.processed)
  
  # Define number of pixels within each chunk
  num.pixels.per.chunk = 500000 # experiment with chunk size to yield a good number of samples
  num.pix <- nrow(image.object) # total number of pixels
  n.chunks <- floor(num.pix/num.pixels.per.chunk)
  message(paste(num.pixels.per.chunk,'pixels per chunk will yield',n.chunks,'total chunks'))
  
  chunked.data <- ChunkIt(image.object = image.object,
                          n.chunks = n.chunks)
  
  # Build the entire dataset needed for plotting
  output.data <- data.frame('ATI-like Transitional' = QuantifyOneLabel(chunked.data,
                                                        label.name = 'ATI-like Transitional',
                                                        base.criteria.positive = c('DAPI'),
                                                        label.criteria.positive = c('Sox9','HOP'),
                                                        label.criteria.negative = c('Vim')),
                               'ATI-like End Stage' = QuantifyOneLabel(chunked.data,
                                                        label.name = 'ATI-like End Stage',
                                                        base.criteria.positive = c('DAPI'),
                                                        label.criteria.positive = c('HOP'),
                                                        label.criteria.negative = c('Vim','Sox9')),
                               'Sox9-Positive Only' = QuantifyOneLabel(chunked.data,
                                                        label.name = 'Sox9-Positive Only',
                                                        base.criteria.positive = c('DAPI'),
                                                        label.criteria.positive = c('Sox9'),
                                                        label.criteria.negative = c('Vim','HOP')),
                                       'Mesenchyme' = QuantifyOneLabel(chunked.data,
                                                        label.name = 'Mesenchyme',
                                                        base.criteria.positive = c('DAPI'),
                                                        label.criteria.positive = c('Vim'),
                                                        label.criteria.negative = c('Sox9','HOP')),
                                         'EMT-like' = QuantifyOneLabel(chunked.data,
                                                        label.name = 'EMT-like',
                                                        base.criteria.positive = c('DAPI'),
                                                        label.criteria.positive = c('Sox9','Vim'),
                                                      # label.criteria.negative = c('HOP'))  
                                                      ))
  
  
  # Label with condition
  output.data$Condition <- 'BDL3' #*************UPDATE NUMBER*************
  
  # Stash so we don't overwrite
  BDL3Data <- melt(output.data) #*************UPDATE NUMBER*************
  
  


--------------------------------------------------------------------------------

#********* Plotting Exploration ********#

# Bind two datasets together for plotting
data <- rbind(BSL3Data,BDL3Data)

png('Test.Output.png',width = 7,height = 5,units='in',res=300)
ggplot(data=data,
       aes(x=Condition, y=value,fill=variable,color=variable))+
  geom_violin()+
  geom_point(color='black',
             size=0.25,
             position = position_jitterdodge(dodge.width = 0.9,jitter.width=0.25))+
  theme_classic()+
  ylab('Fraction of DAPI Pixels')
dev.off()

png('Test.Output.Demo.png',width = 6,height = 5,units='in',res=300)
ggplot(data = data,
       aes(x = Condition,y=value,fill=variable,color=variable))+
  geom_violin()+
  geom_point(position = position_jitterdodge(dodge.width = 0.9,jitter.width=0.25),size=0.1,color='black')+
  theme_classic()+
  ggtitle('Effect of Dynamic Ventilation on BASC Character in 3D')+
  ylab('Fraction of DAPI+ Pixels')+
  scale_fill_manual(values = c('#93B7BE','#F19A3E','#3D3B8E','#E072A4','lightgreen'))+
  scale_color_manual(values = c('#93B7BE','#F19A3E','#3D3B8E','#E072A4','lightgreen'))
dev.off()
