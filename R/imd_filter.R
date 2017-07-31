#' Title
#'
#' @param rRNAdata
#'
#' @return
#' @export
#'
#' @examples
imd_filter <- function(rRNAdata) {

  if (class(rRNAdata) != "rRNAdata") stop("rRNAdata must be of class 'rRNAdata'")
  gdf <- attr(rRNAdata, "group_DF")
  if (is.null(gdf)) stop("group_designation() must be run prior to nonmissing_filter()")

  nonmiss_per_group <- nonmissing_per_group(rRNAdata)
  output <- nonmiss_per_group$nonmiss_totals

  orig_class <- class(output)
  class(output) <- c("imdFilt", orig_class)

  attr(output, "group_sizes") <- nonmiss_per_group$group_sizes

  return(output)
}
