QuantifyOneLabel <- function(chunked.data,
                             label.name = 'BASC-like Quiescent',
                             base.criteria.positive = c('DAPI'),
                             base.criteria.negative = NULL,
                             label.criteria.positive = c('Sox9','Abca3'),
                             label.criteria.negative = c('EdU')){
  # Load QuantFunction
  #source('QuantFunction.R')
  
  n.chunks <- length(chunked.data)
  
  # Apply QuantFunction
  message(paste('Applying QuantFunction to',n.chunks,'chunks...'))
  test.output <- lapply(X = chunked.data,
                        FUN = QuantFunction,
                        label.name = label.name,
                        base.criteria.positive = base.criteria.positive,
                        base.criteria.negative = base.criteria.negative,
                        label.criteria.positive = label.criteria.positive,
                        label.criteria.negative = label.criteria.negative)
  
  # Bind into datafram
  message('Binding into dataframe structure...')
  test.output <- dplyr::bind_rows(test.output)
  test.output$chunk <- c(1:nrow(test.output))
  colnames(test.output)[which(names(test.output) == 'frac.positive')] <- label.name
  
  # Convert to percentages
  message('Converting to percentages...')
  test.output[[label.name]] <- test.output[[label.name]]*100 # convert to percentages
  
  # Output measurements
  message('Outputting binned measurements...')
  return(test.output[[label.name]])
}
