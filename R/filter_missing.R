#' Title
#'
#' @param rRNA_obj
#' @param min_nonmiss
#'
#' @return
#' @export
#'
#' @examples

filter_missing <- function(rRNA_obj, min_nonmiss = 2) {

  # error checks
  if (is.null(rRNA_obj)) stop("rRNA_obj cannot be NULL")
  if (class(rRNA_obj) != "rRNAdata") stop("rRNA_obj must be of class 'rRNAdata'")

  # get name attributes
  cname <- attr(rRNA_obj, "cnames")$edata_cname
  sampid <- attr(rRNA_obj, "cnames")$fdata_cname
  main_effects <- attr(rRNA_obj, "main_effects")
  gdf <- attr(rRNA_obj, "group_DF")
  if (is.null(gdf)) stop("group_designation() must be run prior to filter_nonmiss_per_group()")

  # initalize new rRNA_obj
  ret <- rRNA_obj
  edata <- ret$e_data
# browser()
  # put in format to count missing data
  edata_gather <- tidyr::gather_(edata, sampid, "Value", names(edata)[-which(names(edata) == cname)])
  comb <- merge(edata_gather, gdf, by = sampid)

  # get OTUs to keep
  countdf <- dplyr::summarise(dplyr::group_by_(comb, main_effects, cname), n = sum(Value > 0))
  tmp <- tidyr::spread_(countdf, cname, "n")
  which_keep <- apply(tmp[,-1], 2, function(x) !any(x < min_nonmiss))
  keep_otus <- names(which_keep)[which_keep]
  drop_otus <- names(which_keep)[!which_keep]

  # subset edata
  lvec <- edata[[cname]] %in% keep_otus
  new_edata <- edata[lvec, ]

  ret$e_data <- new_edata

  # record OTUs removed
  attr(ret, "OTUs_removed") <- drop_otus
  attr(ret, "num_OTUs_removed") <- length(drop_otus)
  message(paste("Number of OTUs removed:", length(drop_otus)))
  message(paste("OTUs remaining:", nrow(ret$e_data)))

  return(ret)
}


# filter_samples <- function(rRNA_obj, group = NULL, keep_levels = NULL,
#                            keep_samples = NULL, drop_samples = NULL) {
#
#   # get name attributes
#   cname <- attr(rRNA_obj, "cnames")$edata_cname
#   sampid <- attr(rRNA_obj, "cnames")$fdata_cname
#   main_effects <- attr(rRNA_obj, "main_effects")
#   gdf <- attr(rRNA_obj, "group_DF")
#
#   # initalize new rRNA_obj
#   ret <- rRNA_obj
#   edata <- ret$e_data
#
#   if (!is.null(group)) {
#     if (is.null(gdf)) stop("group_designation() must be applied when 'group' is specified")
#     if (is.null(keep_levels)) {
#       message("No filter applied. Please indicate which levels in the group to keep.")
#       return(rRNA_obj)
#     }
#
#
#   }
#
#
#
#   return(ret)
# }


