# Benchmark-only alternative; production remains unchanged.
subset_tag_positions_interval <- function(tag, positions){
  if(!inherits(tag, "Rle")) return(tag[positions])
  if(is.numeric(positions) && !anyNA(positions) &&
     all(is.finite(positions) & positions >= 1 & positions <= length(tag) & positions == trunc(positions))){
    ends <- cumsum(as.double(runLength(tag)))
    return(as.vector(runValue(tag))[findInterval(positions - 1, ends) + 1L])
  }
  as.vector(tag)[positions]
}
