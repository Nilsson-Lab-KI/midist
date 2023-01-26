# 
# gold standard functions
#

#' @export
# function to convolute: x + z -> y
convolute_gs <- function(z, x, y, gs){
  return((1-gs[z,x])*gs[x,y] + (1-gs[x,z])*gs[z,y])
}

#' @export
get_convoluted_gs <- function(mmm, gs_raw, input, symmetrize_by = max){
  # compare this combined and filtered pairwise matrix to the fractional gold standard
  unique_z <- unique(mmm[which(is.na(mmm) == F)])
  
  # convolute the gold standard for UNWEIGHTED
  gs <- gs_raw
  for (r2 in 1:nrow(gs_raw)){
    for (d2 in r2:nrow(gs_raw)){
      if (mmm[r2,d2] %in% unique_z){
        if (length(get_avg_mid(input$midata, r2, 1)) < length(get_avg_mid(input$midata, d2, 1))) 
          gs[r2,d2] <- gs[d2,r2] <- convolute_gs(mmm[r2,d2], r2, d2, gs_raw) else
            gs[r2,d2] <- gs[d2,r2] <- convolute_gs(mmm[r2,d2], d2, r2, gs_raw)
      }
    }
  }
  # now symmetric
  for (r2 in 1:nrow(gs)){
    for (d2 in r2:ncol(gs)){
      gs[r2,d2] <- gs[d2,r2] <- symmetrize_by(gs[r2,d2], gs[d2,r2])
    }
  }
  
  return(gs)
}
  



