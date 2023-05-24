#' Convert grange list to dataframe
#' 
#' @param grange.list
#' @param anno.id
#'
#' @keywords internal
grList_to_df = function(grange.list, anno.id="cell.population"){
  df.list = list()
  for(name in names(grange.list)){
    df.list[[name]] = grange.list[[name]] %>% as.data.frame()
    df.list[[name]][[anno.id]] = name
  }
  return(Reduce(rbind, df.list))
}

#' Normalize ranks
#' 
#' @param x
#'
#' @keywords internal
range01 = function(x){
    return((x-min(x))/(max(x)-min(x)))
}
