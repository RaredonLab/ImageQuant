# Set WD
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Lab Members/Satoshi")


# Packages
require(EBImage)
require(reshape2)
require(ggplot2)
require(Seurat)
require(ggpubr)

# Load local quantification functions
setwd("/Users/msbr/GitHub/ImageQuant")
source('R Functions/LoadImageChannels.R')
source('R Functions/PreProcessImage.R')
source('R Functions/QuantFunction.R')
source('R Functions/QuantifyOneLabel.R')
source('R Functions/ChunkIt.R')
source('R Functions/MeltIntoDataframe.R')
source('R Functions/NormalizeImages.R')

#### Matrigel - ####
# Load image channels
images.raw <- LoadImageChannels(channel.filepaths = 
                    c("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Lab Members/Satoshi/EVOS images_101323/BASC isolation#1_P0(P1)_SOX9m/101323_BASC isolation #1_P0(P1)_SOX9m_ABCA3.488_Vimentin.555_SOX9.647_Bottom Slide_TR_p00_0_A01f00d0.TIF",
                      "~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Lab Members/Satoshi/EVOS images_101323/BASC isolation#1_P0(P1)_SOX9m/101323_BASC isolation #1_P0(P1)_SOX9m_ABCA3.488_Vimentin.555_SOX9.647_Bottom Slide_TR_p00_0_A01f00d1.TIF",
                      "~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Lab Members/Satoshi/EVOS images_101323/BASC isolation#1_P0(P1)_SOX9m/101323_BASC isolation #1_P0(P1)_SOX9m_ABCA3.488_Vimentin.555_SOX9.647_Bottom Slide_TR_p00_0_A01f00d2.TIF",
                      "~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Lab Members/Satoshi/EVOS images_101323/BASC isolation#1_P0(P1)_SOX9m/101323_BASC isolation #1_P0(P1)_SOX9m_ABCA3.488_Vimentin.555_SOX9.647_Bottom Slide_TR_p00_0_A01f00d3.TIF"),
                    channel.names = 
                      c('Abca3','DAPI','Vim','Sox9'))

display(images.raw$DAPI)

# Normalize each channel
images.norm <- NormalizeImages(images.raw)
display(images.norm$DAPI)
display(images.norm$Sox9)

# PreProcess image channels
images.processed <- PreProcessImage(images.norm = images.norm,
                                    box.dim = 1000,
                                    offset = 0.01) # These parameters will make a huge difference in your outputs and should be tuned + checked
display(images.processed$DAPI)
display(images.processed$Sox9)
display(images.processed$Abca3)
display(images.processed$Vim)

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
output.data <- data.frame('BASC-like Epithelium' = QuantifyOneLabel(chunked.data,
                                                                    label.name = 'BASC-like Epithelium',
                                                                    base.criteria.positive = c('DAPI'),
                                                                    label.criteria.positive = c('Sox9','Abca3'),
                                                                    label.criteria.negative = c('Vim')),
                          'Tuft-like Epithelium' = QuantifyOneLabel(chunked.data,
                                                                    label.name = 'Tuft-like Epithelium',
                                                                    base.criteria.positive = c('DAPI'),
                                                                    label.criteria.positive = c('Sox9'),
                                                                    label.criteria.negative = c('Vim','Abca3')),
                          'Sox9-Positive Epithelium' = QuantifyOneLabel(chunked.data,
                                                                        label.name = 'Sox9-Positive Epithelium',
                                                                        base.criteria.positive = c('DAPI'),
                                                                        label.criteria.positive = c('Sox9'),
                                                                        label.criteria.negative = c('Vim')),
                          'Mesenchyme' = QuantifyOneLabel(chunked.data,
                                                                 label.name = 'Mesenchyme',
                                                                 base.criteria.positive = c('DAPI'),
                                                                 label.criteria.positive = c('Vim'),
                                                                 label.criteria.negative = c('Sox9','Abca3')),
                          'Epi.Tuft-like' = QuantifyOneLabel(chunked.data,
                                                          label.name = 'Mesenchyme',
                                                          base.criteria.positive = c('DAPI'),
                                                          base.criteria.negative = c('Vim'),
                                                          label.criteria.positive = c('Sox9'),
                                                          label.criteria.negative = c('Abca3')))
# Label with condition
output.data$Condition <- 'Matrigel-'

# Stash so we don't overwrite
data1 <- melt(output.data)

#### Matrigel + ####
# Load image channels
images.raw <- LoadImageChannels(channel.filepaths = 
                                  c("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Lab Members/Satoshi/EVOS images_101323/BASC isolation#1_P0(P1)_SOX9m+Matrigel/101323_BASC isolation #1_P0(P1)_SOX9m+Matrigel_ABCA3.488_Vimentin.555_SOX9.647_Bottom Slide_TR_p00_0_A01f00d0.TIF",
                                    "~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Lab Members/Satoshi/EVOS images_101323/BASC isolation#1_P0(P1)_SOX9m+Matrigel/101323_BASC isolation #1_P0(P1)_SOX9m+Matrigel_ABCA3.488_Vimentin.555_SOX9.647_Bottom Slide_TR_p00_0_A01f00d1.TIF",
                                    "~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Lab Members/Satoshi/EVOS images_101323/BASC isolation#1_P0(P1)_SOX9m+Matrigel/101323_BASC isolation #1_P0(P1)_SOX9m+Matrigel_ABCA3.488_Vimentin.555_SOX9.647_Bottom Slide_TR_p00_0_A01f00d2.TIF",
                                    "~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Raredon_Lab_Administration/Lab Members/Satoshi/EVOS images_101323/BASC isolation#1_P0(P1)_SOX9m+Matrigel/101323_BASC isolation #1_P0(P1)_SOX9m+Matrigel_ABCA3.488_Vimentin.555_SOX9.647_Bottom Slide_TR_p00_0_A01f00d3.TIF"),
                                    channel.names = c('Abca3','DAPI','Vim','Sox9'))

display(images.raw$DAPI)

# Normalize each channel
images.norm <- NormalizeImages(images.raw)
display(images.norm$DAPI)

# PreProcess image channels
images.processed <- PreProcessImage(images.norm = images.norm,
                                    box.dim = 1000,
                                    offset = 0.01)
display(images.processed$DAPI)
display(images.processed$Sox9)
display(images.processed$Abca3)
display(images.processed$Vim)

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
output.data <- data.frame('BASC-like Epithelium' = QuantifyOneLabel(chunked.data,
                                                                    label.name = 'BASC-like Epithelium',
                                                                    base.criteria.positive = c('DAPI'),
                                                                    label.criteria.positive = c('Sox9','Abca3'),
                                                                    label.criteria.negative = c('Vim')),
                          'Tuft-like Epithelium' = QuantifyOneLabel(chunked.data,
                                                                    label.name = 'Tuft-like Epithelium',
                                                                    base.criteria.positive = c('DAPI'),
                                                                    label.criteria.positive = c('Sox9'),
                                                                    label.criteria.negative = c('Vim','Abca3')),
                          'Sox9-Positive Epithelium' = QuantifyOneLabel(chunked.data,
                                                                        label.name = 'Sox9-Positive Epithelium',
                                                                        base.criteria.positive = c('DAPI'),
                                                                        label.criteria.positive = c('Sox9'),
                                                                        label.criteria.negative = c('Vim')),
                          'Mesenchyme' = QuantifyOneLabel(chunked.data,
                                                          label.name = 'Mesenchyme',
                                                          base.criteria.positive = c('DAPI'),
                                                          label.criteria.positive = c('Vim'),
                                                          label.criteria.negative = c('Sox9','Abca3')))# Label with condition
output.data$Condition <- 'Matrigel+'

# Stash so we don't overwrite
data2 <- melt(output.data)

#### Plotting Exploration ####

# Bind two datasets together for plotting
data <- rbind(data1,data2)

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
  ggtitle('Effect of Matrigel on Epithelial Stemness in 2D Culture')+
  ylab('Fraction of DAPI+ Pixels')+
  scale_fill_manual(values = c('#93B7BE','#F19A3E','#3D3B8E','#E072A4'))+
  scale_color_manual(values = c('#93B7BE','#F19A3E','#3D3B8E','#E072A4'))
dev.off()
