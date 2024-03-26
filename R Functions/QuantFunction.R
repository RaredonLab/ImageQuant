
QuantFunction <- function(pixels,
                          label.name = 'BASC-like Quiescent',
                          base.criteria.positive = c('DAPI'),
                          base.criteria.negative = NULL,
                          label.criteria.positive = c('Sox9','Abca3'),
                          label.criteria.negative = c('EdU')){
  
  # temp for testing
  #pixels <- test[[1]]
  
  # What fraction of pixels satisfy the base.criteria.positive?
  if(!is.null(base.criteria.positive)){
    sub <- pixels
      for(i in 1:length(base.criteria.positive)){
        sub <- sub[sub[[base.criteria.positive[i]]]>0,] # All base.criteria.positive conditions must be positive
      }
      }else{
        sub <- pixels # If no base criteria are provided, then consider all pixels
      }
  
  # What fraction of the above satisfy the base.criteria.negative?
  if(!is.null(base.criteria.negative)){
    #sub <- sub
    for(i in 1:length(base.criteria.negative)){
      sub <- sub[sub[[base.criteria.negative[i]]]==0,] # All base.criteria.negative conditions must be negative
    }
  }else{
    #sub <- sub # If no base.criteria.negative are provided, then consider all pixels meeting base.criteria.negative
  }
  
  # What fraction of the above base pixels satisfy label.criteria.positive?
  if(!is.null(label.criteria.positive)){
    sub.2 <- sub
      for(i in 1:length(label.criteria.positive)){
        sub.2 <- sub.2[sub.2[[label.criteria.positive[i]]]>0,] # How many of the remaining pixels meet the positive criteria?
      }
    }else{
      sub.2 <- sub # If no label.criteria.positive are provided, then consider all pixels meeting base.criteria
    }
  # What fraction of the above pixels satisfy label.criteria.negative?
  if(!is.null(label.criteria.negative)){
    sub.3 <- sub.2
    for(i in 1:length(label.criteria.negative)){
      sub.3 <- sub.3[sub.3[[label.criteria.negative[i]]]==0,] # How many of the remaining pixels meet the negative criteria?
    }
    }else{
      sub.3 <- sub.2 # If no label.criteria.negative are provided, then consider all pixels meeting base.criteria
    }
  
  # Only count as true measure if we measured some base pixels, otherwise, output NA
  if(nrow(sub)>0){
    frac.positive <- nrow(sub.3)/nrow(sub)
    }else{
    frac.positive <- NA
  }
  return(data.frame(frac.positive))
}

# 
# 
#     dapi.basc <- sum(sub$abca3>0 & sub$sox9>0)/nrow(sub) # 0.6% BASC-like
#     dapi.tuft <- sum(sub$abca3==0 & sub$sox9>0)/nrow(sub) # 1.9% Tuft-like
#     dapi.atii <- sum(sub$abca3>0 & sub$sox9==0)/nrow(sub) # 16.6% ATII-like
#     dapi.other <- sum(sub$abca3==0 & sub$sox9==0)/nrow(sub) # 80.7% Other
#     
#     dapi.prolif <- sum(sub$edu>0)/nrow(sub) # 0.2% Proliferating
#     dapi.quiesc <- sum(sub$edu==0)/nrow(sub) # 99.7% Quiescent
#     
#     dapi.prolif.basc <- sum(sub$abca3>0 & sub$sox9>0 & sub$edu>0)/nrow(sub) # 0% BASC-prolif
#     dapi.prolif.tuft <- sum(sub$abca3==0 & sub$sox9>0 & sub$edu>0)/nrow(sub) # 0.18% Tuft-prolif
#     dapi.prolif.atii <- sum(sub$abca3>0 & sub$sox9==0 & sub$edu>0)/nrow(sub) # 0% ATII-prolif
#     dapi.prolif.other <- sum(sub$abca3==0 & sub$sox9==0 & sub$edu>0)/nrow(sub) # 0.09% Other-prolif
#     
#     dapi.quiesc.basc <- sum(sub$abca3>0 & sub$sox9>0 & sub$edu==0)/nrow(sub) # 0.6% BASC-quiesc
#     dapi.quiesc.tuft <- sum(sub$abca3==0 & sub$sox9>0 & sub$edu==0)/nrow(sub) # 1.7% Tuft-quiesc
#     dapi.quiesc.atii <- sum(sub$abca3>0 & sub$sox9==0 & sub$edu==0)/nrow(sub) # 16.6% ATII-quiesc
#     dapi.quiesc.other <- sum(sub$abca3==0 & sub$sox9==0 & sub$edu==0)/nrow(sub) # 80.6% Other-quiesc
#   }else{
#     dapi.basc <- NA
#     dapi.tuft <- NA
#     dapi.atii <- NA
#     dapi.other <- NA
#     dapi.prolif <- NA
#     dapi.quiesc <- NA
#     dapi.prolif.basc <- NA
#     dapi.prolif.tuft <- NA
#     dapi.prolif.atii <- NA
#     dapi.prolif.other <- NA
#     dapi.quiesc.basc <- NA
#     dapi.quiesc.tuft <- NA
#     dapi.quiesc.atii <- NA
#     dapi.quiesc.other <- NA
#   }
#   # What fraction of Tuft-like pixels are also ... ?
#   sub <- X[X$sox9>0 & X$abca3==0,]
#   if(sum(sub$sox9)>0){
#     tuft.prolif <- sum(sub$edu>0)/nrow(sub) # 13.6% Proliferating
#     tuft.quiesc <- sum(sub$edu==0)/nrow(sub) # 86.3% Quiescent
#   }else{
#     tuft.prolif <- NA
#     tuft.quiesc <- NA}
#   # What fraction of BASC-like pixels are also ... ?
#   sub <- X[X$abca3>0 & X$sox9>0,]
#   if(sum(sub$sox9)>0 & sum(sub$abca3)>0 ){
#     basc.prolif <- sum(sub$edu>0)/nrow(sub) # 68.5% Proliferating
#     basc.quiesc <- sum(sub$edu==0)/nrow(sub) # 31.4% Quiescent
#   }else{
#     basc.prolif <- NA
#     basc.quiesc <- NA}
#   # What fraction of ATII-like pixels are also ... ?
#   sub <- X[X$abca3>0 & X$sox9==0,]
#   if(sum(sub$abca3)>0){
#     atii.prolif <- sum(sub$edu>0)/nrow(sub) # 0% Proliferating
#     atii.quiesc <- sum(sub$edu==0)/nrow(sub) # 100% Quiescent
#   }else{
#     atii.prolif <- NA
#     atii.quiesc <- NA}
#   # What fraction of Other pixels are also ... ?
#   sub <- X[X$abca3==0 & X$sox9==0,]
#   if(sum(sub$dapi)>0){
#     other.prolif <- sum(sub$edu>0)/nrow(sub) # Proliferating
#     other.quiesc <- sum(sub$edu==0)/nrow(sub) # Quiescent
#   }else{
#     other.prolif <- NA
#     other.quiesc <- NA}
#   # What fraction of proliferating pixels are also ... ?
#   sub <- X[X$edu>0,]
#   if(sum(sub$edu)>0){
#     prolif.basc <- sum(sub$abca3>0 & sub$sox9>0)/nrow(sub) # 17.7% BASC-like
#     prolif.tuft <- sum(sub$abca3==0 & sub$sox9>0)/nrow(sub) # 23.2% Tuft-like
#     prolif.atii <- sum(sub$abca3>0 & sub$sox9==0)/nrow(sub) # 0.3% ATII-like
#     prolif.other <- sum(sub$abca3==0 & sub$sox9==0)/nrow(sub) # 40.8% Other
#   }else{
#     prolif.basc <- NA
#     prolif.tuft <- NA
#     prolif.atii <- NA
#     prolif.other <- NA
#   }
#   return(data.frame(tuft.prolif,tuft.quiesc,
#                     basc.prolif,basc.quiesc,
#                     atii.prolif,atii.quiesc,
#                     other.prolif,other.quiesc,
#                     prolif.basc,prolif.tuft,prolif.atii,prolif.other,
#                     dapi.basc,
#                     dapi.tuft,
#                     dapi.atii,
#                     dapi.other,
#                     dapi.prolif,
#                     dapi.quiesc,
#                     dapi.prolif.basc,
#                     dapi.prolif.tuft,
#                     dapi.prolif.atii,
#                     dapi.prolif.other,
#                     dapi.quiesc.basc,
#                     dapi.quiesc.tuft,
#                     dapi.quiesc.atii,
#                     dapi.quiesc.other))
# }