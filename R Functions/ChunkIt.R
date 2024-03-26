ChunkIt <- function(image.object,n.chunks){
  # Split image.object into chunks
  message(paste('Splitting into',n.chunks,'chunks...'))
  chunked.data <- split(image.object, (seq(nrow(image.object))-1) %/% (nrow(image.object)/n.chunks)) # break into chunks
  return(chunked.data)
  message('Done.')
  }
