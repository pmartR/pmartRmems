#' Nonmissing per group filter
#'
#' @param rRNAdata
#'
#' @return a filter object that can be used with applyFilt to filter out missing OTUs.
#' @export
#'
#' @examples
imd_filter <- function(rRNAdata) {

  if (class(rRNAdata) != "rRNAdata") stop("rRNAdata must be of class 'rRNAdata'")
  if (is.null(attr(rRNAdata, "group_DF")))
    stop("group_designation() must be run prior to nonmissing_filter()")

  nonmiss_per_group <- nonmissing_per_group(rRNAdata)
  output <- nonmiss_per_group$nonmiss_totals

  orig_class <- class(output)
  class(output) <- c("imdFilt", orig_class)

  attr(output, "group_sizes") <- nonmiss_per_group$group_sizes

  return(output)
}
