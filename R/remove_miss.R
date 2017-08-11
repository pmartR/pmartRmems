#' Remove Missing Rows
#'
#' Typically run within applyFilt to make sure there are no rows with all missing values
#' after a filter is applied.
#'
#' @param rRNAdata
#'
#' @return filtered rRNAdata object
#' @export
#'
#' @examples
remove_miss <- function(rRNAdata) {

  edata_id <- attr(rRNAdata, "cnames")$edata_cname

  # check e_data for missing rows
  not_missing <- rowSums(rRNAdata$e_data[, names(rRNAdata$e_data) != edata_id]) > 0
  rRNAdata$e_data <- rRNAdata$e_data[not_missing, ]

  # if OTUs were missing, remove them from e_meta
  if (sum(!not_missing) > 0 && !is.null(rRNAdata$e_meta)) {
    rRNAdata$e_meta <- rRNAdata$e_meta[not_missing, ]
  }

  return(rRNAdata)
}
