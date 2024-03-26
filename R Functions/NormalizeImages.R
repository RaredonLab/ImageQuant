NormalizeImages <- function(images.raw){
  
        # Identify channel names
      channel.names <- names(images.raw)
      message(paste('Using',paste(channel.names,collapse=" & "),'as channels...'))
      
      # Normalize each channel
      message('Normalizing data...')
      images.norm <- list()
      for(i in 1:length(channel.names)){
        images.norm[[channel.names[i]]] <- EBImage::normalize(images.raw[[i]])

      }
      return(images.norm)
}
