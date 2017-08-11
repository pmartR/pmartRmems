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
applyFilt.customFilt <- function(filter_object, rRNAdata){


  edata_cname <- attributes(rRNAdata)$cnames$edata_cname
  fdata_cname <- attributes(rRNAdata)$cnames$fdata_cname
  emeta_cname <- attributes(rRNAdata)$cnames$emeta_cname

  #if filter_object contains removes
  if(!is.null(filter_object$e_data_remove)||!is.null(filter_object$f_data_remove)||!is.null(filter_object$e_meta_remove)){

    filter_object_new = list(edata_filt = filter_object$e_data_remove, emeta_filt = filter_object$e_meta_remove, samples_filt = filter_object$f_data_remove)

    # check that edata_filt doesn't specify ALL the items in rRNAdata #
    if(all(rRNAdata$e_data[, edata_cname] %in% filter_object_new$edata_filt)){stop("edata_filt specifies all the items in the data")}

    # check that samples_filt doesn't specify ALL the items in rRNAdata #
    if(all(rRNAdata$f_data[, fdata_cname] %in% filter_object_new$samples_filt)){stop("samples_filt specifies all the items in the data")}

    # check that emeta_filt doesn't specify ALL the items in rRNAdata, emeta_filt is present #
    if(!is.null(rRNAdata$e_meta[, emeta_cname])){
      if(all(rRNAdata$e_meta[, emeta_cname] %in% filter_object_new$emeta_filt)){stop("emeta_filt specifies all the items in the data")}
    }


  }

  else{
    filter_object_new = list(edata_keep = filter_object$e_data_keep, emeta_keep = filter_object$e_meta_keep, samples_keep = filter_object$f_data_keep)

    # check that edata_keep doesn't specify ALL the items in rRNAdata #
    if(all(rRNAdata$e_data[, edata_cname] %in% filter_object_new$edata_keep)){stop("edata_keep specifies all the items in the data")}

    # check that samples_keep doesn't specify ALL the items in rRNAdata #
    if(all(rRNAdata$f_data[, fdata_cname] %in% filter_object_new$samples_keep)){stop("samples_keep specifies all the items in the data")}

    # check that emeta_keep doesn't specify ALL the items in rRNAdata #
    if(!is.null(rRNAdata$e_meta[, emeta_cname])){
      if(all(rRNAdata$e_meta[, emeta_cname] %in% filter_object_new$emeta_keep)){stop("emeta_keep specifies all the items in the data")}

    }
  }

  # call the function that does the filter application
  results_pieces <- MSomics_filter_worker(rRNAdata = rRNAdata, filter_object = filter_object_new)

  # return filtered data object #
  results <- rRNAdata
  results$e_data <- results_pieces$temp.pep2
  results$f_data <- results_pieces$temp.samp2
  if(!is.null(rRNAdata$e_meta)){ # if-statement added by Kelly 3/24/2017 #
    results$e_meta <- data.frame(results_pieces$temp.meta1)
    names(results$e_meta)[which(names(rRNAdata$e_meta) == emeta_cname)] <- emeta_cname
  }else{
    # e_meta is null
    results$e_meta <- NULL
  }

  # remove any OTUs with all missing data
  results <- remove_miss(results)


  # if group attribute is present, re-run group_designation in case filtering any items impacted the group structure #
  if(!is.null(attr(results, "group_DF"))){
    results <- group_designation(rRNAdata = results, main_effects = attr(attr(rRNAdata, "group_DF"), "main_effects"), covariates = attr(attr(rRNAdata, "group_DF"), "covariates"), time_course = attr(attr(rRNAdata, "group_DF"), "time_course"))
  }else{
    # Update attributes (7/11/2016 by KS) - this is being done already in group_designation
    attributes(results)$data_info$num_edata = length(unique(results$e_data[, edata_cname]))
    attributes(results)$data_info$num_miss_obs = sum(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
    attributes(results)$data_info$num_prop_missing = mean(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
    attributes(results)$data_info$num_samps = ncol(results$e_data) - 1

    if(!is.null(results$e_meta)){
      # number of unique proteins that map to a peptide in e_data #
      if(!is.null(emeta_cname)){
        num_emeta = length(unique(results$e_meta[which(as.character(results$e_meta[, edata_cname]) %in% as.character(results$e_data[, edata_cname])), emeta_cname]))
      }else{num_emeta = NULL}
    }else{
      num_emeta = NULL
    }
    attr(results, "data_info")$num_emeta = num_emeta
    ## end of update attributes (7/11/2016 by KS)
  }


  # check to see whether a customFilt has already been run on rRNAdata #
  if("customFilt" %in% names(attributes(rRNAdata)$filters)){
    # yes, so be sure to keep track of this customFilt in a separate $filters attribute #

    # get number of previous customFilts applied #
    n_prev_filts <- length(grep("customFilt", names(attributes(rRNAdata)$filters)))
    # will need to append this number + 1 to the end of the $filters$customFilt attribute #


    # set attributes for which filters were run
    attr(results, "filters")[[(paste("customFilt", n_prev_filts + 1, sep=""))]] <- list(report_text = "", threshold = c(), filtered = c())

    n_edata_filtered <- nrow(rRNAdata$e_data) - nrow(results$e_data)
    n_fdata_filtered <- nrow(rRNAdata$f_data) - nrow(results$f_data)
    if(!is.null(rRNAdata$e_meta)){
      n_emeta_filtered <- nrow(rRNAdata$e_meta) - nrow(results$e_meta)
    }else{
      n_emeta_filtered = NA
    }
    if(!is.null(rRNAdata$e_meta)){
      attr(results, "filters")[[(paste("customFilt", n_prev_filts + 1, sep=""))]]$report_text <- paste("A custom filter was applied to the data, removing ", n_edata_filtered, " ", edata_cname, "s from e_data, ", n_emeta_filtered, " ", emeta_cname, "s from e_meta, and ", n_fdata_filtered, " ", fdata_cname, "s from f_data.", sep="")
    }else{
      attr(results, "filters")[[(paste("customFilt", n_prev_filts + 1, sep=""))]]$report_text <- paste("A custom filter was applied to the data, removing ", n_edata_filtered, " ", edata_cname, "s from e_data and ", n_fdata_filtered, " ", fdata_cname, "s from f_data.", sep="")
    }
    attr(results, "filters")[[(paste("customFilt", n_prev_filts + 1, sep=""))]]$filtered <- filter_object

  }else{ # no previous customFilt, so go ahead and name the attribute like normal #

    # set attributes for which filters were run
    attr(results, "filters")$customFilt <- list(report_text = "", threshold = c(), filtered = c())
    n_edata_filtered <- nrow(rRNAdata$e_data) - nrow(results$e_data)
    n_fdata_filtered <- nrow(rRNAdata$f_data) - nrow(results$f_data)
    if(!is.null(rRNAdata$e_meta)){
      n_emeta_filtered <- nrow(rRNAdata$e_meta) - nrow(results$e_meta)
    }else{
      n_emeta_filtered = NA
    }
    if(!is.null(rRNAdata$e_meta)){
      attr(results, "filters")$customFilt$report_text <- paste("A custom filter was applied to the data, removing ", n_edata_filtered, " ", edata_cname, "s from e_data, ", n_emeta_filtered, " ", emeta_cname, "s from e_meta, and ", n_fdata_filtered, " ", fdata_cname, "s from f_data.", sep="")
    }else{
      attr(results, "filters")$customFilt$report_text <- paste("A custom filter was applied to the data, removing ", n_edata_filtered, " ", edata_cname, "s from e_data and ", n_fdata_filtered, " ", fdata_cname, "s from f_data.", sep="")
    }
    attr(results, "filters")$customFilt$filtered <- filter_object
  }
  return(results)
}





#' @export
#' @name applyFilt
#' @rdname applyFilt
applyFilt.subsetFilt <- function(filter_object, rRNAdata) {

  fdata_cname <- attr(rRNAdata, "cnames")$fdata_cname
  edata_cname <- attr(rRNAdata, "cnames")$edata_cname

  # filter samples from f_data
  samples_keep_fdata <- !(rRNAdata$f_data[[fdata_cname]] %in% filter_object)
  new_fdata <- rRNAdata$f_data[samples_keep_fdata, ]

  # filter samples from e_data
  samples_keep_edata <- !(names(rRNAdata$e_data) %in% filter_object)
  new_edata <- rRNAdata$e_data[, samples_keep_edata]

  results <- rRNAdata
  results$e_data <- new_edata
  results$f_data <- droplevels(new_fdata)

  # remove any OTUs with all missing data
  results <- remove_miss(results)

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
  }

  return(results)
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

  return(results)
}


anova_filter <- function(nonmiss_per_group, min_nonmiss_anova=2, cname_id){

  # min_nonmiss_anova must be >=2
  if(!is.null(min_nonmiss_anova)){
    if(min_nonmiss_anova<2){
      stop("min_nonmiss_anova must be >=2")
    }
  }

  # check that nonmiss_per_group is a list of length 2 #
  if(class(nonmiss_per_group) != "list" | length(nonmiss_per_group)!=2)
    stop("nonmiss_per_group must be a list of length 2.")

  # column names of nonmiss_per_group$nonmiss_totals
  my.names <- names(nonmiss_per_group$nonmiss_totals)
  inds.names <- my.names %in% c(cname_id, "NA")
  inds.names <- which(inds.names==TRUE)
  my.names <- my.names[-inds.names] # remove cname_id and "NA" (if there's a column named "NA")

  # need at least n=2 per group in at least 2 groups in order to keep the biomolecule
  temp <- (nonmiss_per_group$nonmiss_totals[,which(names(nonmiss_per_group$nonmiss_totals) %in%
                                                     my.names)] >= min_nonmiss_anova)

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
