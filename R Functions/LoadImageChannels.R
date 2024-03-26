LoadImageChannels <- function(channel.filepaths = 
                                c("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Experiments/3. Pneumonectomy/EdU - Pneumonectomy/2021.11.30 PNEdU ABCA3-g EdU-r Sox9-c/PN0/2021.11.30 Sox9-c ABCA3-g EdU-r tiled/Control  tiled_DAPI.TIF",
                                  "/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Experiments/3. Pneumonectomy/EdU - Pneumonectomy/2021.11.30 PNEdU ABCA3-g EdU-r Sox9-c/PN0/2021.11.30 Sox9-c ABCA3-g EdU-r tiled/Control  tiled Sox9 .TIF",
                                  "/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Experiments/3. Pneumonectomy/EdU - Pneumonectomy/2021.11.30 PNEdU ABCA3-g EdU-r Sox9-c/PN0/2021.11.30 Sox9-c ABCA3-g EdU-r tiled/Control  tiled EdU .TIF",
                                  "/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Experiments/3. Pneumonectomy/EdU - Pneumonectomy/2021.11.30 PNEdU ABCA3-g EdU-r Sox9-c/PN0/2021.11.30 Sox9-c ABCA3-g EdU-r tiled/Control  tiled  ABCA3 .TIF"),
                              channel.names = 
                                c('DAPI','Sox9','EdU','Abca3')
){
  # Check that length of filepaths equals length of names
  if(length(channel.filepaths) != length(channel.names)){stop('length of filepaths does not equal length of names')}else{}
  
  # Load each channel of an image and stash each as a separate grayscale matrix
  message('Loading raw image data...')
  images.raw <- list()
  for(i in 1:length(channel.names)){
    images.raw[[channel.names[i]]] <- EBImage::channel(EBImage::readImage(channel.filepaths[i]),'gray')
  }
  return(images.raw)
}
