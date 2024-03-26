MeltIntoDataframe <- function(images.processed){
  
  # Identify channel names
  channel.names <- names(images.processed)
  message(paste('Using',paste(channel.names,collapse=" & "),'as channels...'))
  
  # Merge all information from image into a single dataframe
  message('Melting into dataframe...')
  images.melt <- list()
  for(i in 1:length(channel.names)){
    images.melt[[channel.names[i]]] <- reshape2::melt(images.processed[[channel.names[i]]])
    colnames(images.melt[[channel.names[i]]]) <- c('x','y',channel.names[i])
  }
  
  # Combine into an 'image.object' dataframe
  image.object <- data.frame('x' = images.melt[[1]]$x,
                             'y' = images.melt[[1]]$y)
  for(i in 1:length(channel.names)){
    image.object[[channel.names[i]]] <- images.melt[[channel.names[i]]][[channel.names[i]]]
  }
  return(image.object)
}
