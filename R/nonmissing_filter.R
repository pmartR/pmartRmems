#' Title
#'
#' @param rRNA_obj
#'
#' @return
#' @export
#'
#' @examples
nonmissing_filter <- function(rRNA_obj) {

  if (class(rRNA_obj) != "rRNAdata") stop("rRNA_obj must be of class 'rRNAdata'")
  gdf <- attr(rRNA_obj, "group_DF")
  if (is.null(gdf)) stop("group_designation() must be run prior to nonmissing_filter()")

  nonmiss_per_group <- nonmissing_per_group(rRNA_obj)
  output <- nonmiss_per_group$nonmiss_totals

  orig_class <- class(output)
  class(output) <- c("nonmissFilt", orig_class)

  attr(output, "group_sizes") <- nonmiss_per_group$group_sizes

  return(output)
}
