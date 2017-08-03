#' Title
#'
#' @param rRNAdata
#'
#' @return
#' @export
#'
#' @examples
fdata_filter <- function(rRNAdata, grp_level_rmv = NULL, grp_level_keep = NULL,
                         cname = NULL, level_rmv = NULL, level_keep = NULL) {

  #### Initial Error Checks ####
  if (class(rRNAdata) != "rRNAdata") stop("rRNAdata must be of class 'rRNAdata'")
  # count how many "level" arguments are specified
  sum_args <- sum(sapply(list(grp_level_rmv, grp_level_keep, level_rmv, level_keep), is.null))
  # only one can be specified, so three must be NULL
  if (sum_args != 3) {
    if (sum_args < 4)
      stop("Only one of grp_level_rmv, grp_level_keep, level_rmv, or level_keep can be specified")
    if (sum_args == 4)
      stop("Please specify one of grp_level_rmv, grp_level_keep, level_rmv, or level_keep")
  }
  if (!is.null(level_rmv) || !is.null(level_keep)) {
    if (is.null(cname))
      stop("cname must be specified when either level_rmv or level_keep is used")
    if (length(cname) > 1) {
      warning("only the first element of cname will be used")
      cname <- cname[1]
    }
    if (!(cname %in% names(rRNAdata$f_data)))
      stop("cname must be a column name in f_data")
  }
  #### End Initial Error Checks ####

  fdata_cname <- attr(rRNAdata, "cnames")$fdata_cname

  ### removal begins in group_DF
  if (!is.null(grp_level_rmv) || !is.null(grp_level_keep)) {

    gdf <- attr(rRNAdata, "group_DF")
    if (is.null(gdf)) stop("group_designation() must be run when filtering group levels")

    # remove
    if (!is.null(grp_level_rmv)) {
      if (!all(grp_level_rmv %in% gdf$Group))
        stop("grp_level_rmv must contain levels of 'Group' within group_DF")
      filter_levels <- gdf$Group %in% grp_level_rmv

    # keep
    } else {
      if (!all(grp_level_keep %in% gdf$Group))
        stop("grp_level_keep must contain levels of 'Group' within group_DF")
      filter_levels <- !(gdf$Group %in% grp_level_keep)
    }

    # samples to remove
    output <- gdf[[fdata_cname]][filter_levels]

  ### removal begins in f_data
  } else {
    f_data <- rRNAdata$f_data

    # remove
    if (!is.null(level_rmv)) {
      if (!all(level_rmv %in% f_data[[cname]]))
        stop("level_rmv must contain levels of cname within f_data")
      filter_levels <- f_data[[cname]] %in% level_rmv

    # keep
    } else {
      if (!all(level_keep %in% f_data[[cname]]))
        stop("level_keep must contain levels of cname within f_data")
      filter_levels <- !(f_data[[cname]] %in% level_keep)
    }

    # samples to remove
    output <- f_data[[fdata_cname]][filter_levels]
  }

  output <- as.character(output)
  orig_class <- class(output)
  class(output) <- c("fdataFilt", orig_class)

  return(output)
}
