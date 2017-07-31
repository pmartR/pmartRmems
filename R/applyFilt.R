#' Title
#'
#' @param filter_object
#' @param rRNAdata
#'
#' @return
#' @export
#'
#' @examples
applyFilt <- function(filter_object, rRNAdata, ...) {
  UseMethod("applyFilt")
}

#' @export
#' @name applyFilt
#' @rdname applyFilt
applyFilt.imdFilt <- function(filter_object, rRNAdata, min_nonmiss = 2) {

  group_sizes <- attr(filter_object, "group_sizes")
  nonmiss_per_group <- list(nonmiss_totals = filter_object, group_sizes = group_sizes)

  groupDF <- attributes(rRNAdata)$group_DF
  e_data <- rRNAdata$e_data
  edata_cname <- attr(rRNAdata, "cnames")$edata_cname
  samp_cname <- attr(rRNAdata, "cnames")$fdata_cname
  emeta_cname <- attributes(rRNAdata)$cnames$emeta_cname

  filter.edata <- anova_filter(nonmiss_per_group = nonmiss_per_group,
                               min_nonmiss_anova = min_nonmiss, cname_id = edata_cname)

  #checking that filter.edata does not specify all of e_data in rRNAdata
  if (all(rRNAdata$e_data[[edata_cname]] %in% filter.edata)) {
    stop("filter.edata specifies all of e_data in rRNAdata")
  }

  filter_object_new <- list(edata_filt = filter.edata, emeta_filt = NULL, samples_filt = NULL)

  # call the function that does the filter application
  results_pieces <- MSomics_filter_worker(filter_object = filter_object_new, omicsData = rRNAdata)

  # return filtered data object #
  rRNAdata$e_data <- results_pieces$temp.pep2
  rRNAdata$f_data <- results_pieces$temp.samp2
  rRNAdata$e_meta <- results_pieces$temp.meta1
  results <- rRNAdata

  # if group attribute is present, re-run group_designation in case filtering any items
  # impacted the group structure
  if (!is.null(attr(results, "group_DF"))){
    results <- group_designation(results,
                                 main_effects = attr(attr(rRNAdata, "group_DF"), "main_effects"))
  } else {
    # Update attributes (7/11/2016 by KS) - this is being done already in group_designation
    attributes(results)$data_info$num_edata <- length(unique(results$e_data[, edata_cname]))
    attributes(results)$data_info$num_miss_obs <-
      sum(is.na(results$e_data[,-which(names(results$e_data) == edata_cname)]))
    attributes(results)$data_info$num_prop_missing <-
      mean(is.na(results$e_data[,-which(names(results$e_data) == edata_cname)]))
    attributes(results)$data_info$num_samps <- ncol(results$e_data) - 1

    if (!is.null(results$e_meta)) {
      # number of unique proteins that map to a peptide in e_data #
      if (!is.null(emeta_cname)) {
        num_emeta <-
          length(unique(results$e_meta[which(as.character(results$e_meta[, edata_cname]) %in%
                                               as.character(results$e_data[, edata_cname])),
                                       emeta_cname]))
      } else {
        num_emeta <- NULL
      }
    } else {
      num_emeta <- NULL
    }
    attr(results, "data_info")$num_emeta <- num_emeta
  }

}




anova_filter <- function(nonmiss_per_group, min_nonmiss_anova=2, cname_id){

  # min_nonmiss_anova must be >=2
  if(!is.null(min_nonmiss_anova)){
    if(min_nonmiss_anova<2){
      stop("min_nonmiss_anova must be >=2")
    }
  }

  # check that nonmiss_per_group is a list of length 2 #
  if(class(nonmiss_per_group) != "list" | length(nonmiss_per_group)!=2) stop("nonmiss_per_group must be a list of length 2.")

  # column names of nonmiss_per_group$nonmiss_totals
  my.names <- names(nonmiss_per_group$nonmiss_totals)
  inds.names <- my.names %in% c(cname_id, "NA")
  inds.names <- which(inds.names==TRUE)
  my.names <- my.names[-inds.names] # remove cname_id and "NA" (if there's a column named "NA")

  # need at least n=2 per group in at least 2 groups in order to keep the biomolecule
  temp <- (nonmiss_per_group$nonmiss_totals[,which(names(nonmiss_per_group$nonmiss_totals) %in% my.names)] >= min_nonmiss_anova)

  # sum the number of groups that meet the nonmissing per group requirement
  temp2 <- rowSums(temp)

  # append the vector of biomolecules
  temp3 <- data.frame(nonmiss_per_group$nonmiss_totals[, 1], temp2)
  names(temp3)[1] <- names(nonmiss_per_group$nonmiss_totals)[1]

  # create indicator for which rows do not meet the nonmissing requirement for at least 2 groups
  inds.rm <- which(temp2 < 2) # these are the rows to remove since they do not have at least 2 groups not meeting nonmissing requirements

  # get names of biomolecules to be filtered
  filter.ids <- as.vector(temp3[inds.rm, 1])

  # output names of peptides/proteins/genes to be filtered
  return(filter.ids)
}
