
#' Apply a transformation to rRNA data
#'
#' @param rRNA_obj
#' @param method
#'
#' @return
#' @export
#'
#' @examples
transform_data <- function(rRNA_obj, method = "clr", shift = 0.5, basis_ind = NULL,
                           sampleMC = FALSE, ...) {
  .transform_data(rRNA_obj, method, shift, basis_ind, sampleMC, ...)
}

.transform_data <- function(rRNA_obj, method = "clr", shift = 0.5, basis_ind = NULL,
                            sampleMC = FALSE, MCdistn = "dirichlet", niter = 200,
                            include_zero_rows = FALSE, quiet = FALSE, track_seed = FALSE) {

  if (class(rRNA_obj) != "rRNAdata")
    stop("rRNA_obj must be of class 'rRNAdata.' See ?as.rRNAdata for more info.")
  if(is.null(method))
    stop("method must be one of clr, alr, iqlr, zero, or user")

  # ignore case
  method <- tolower(method)

  if (!(method %in% c("clr", "alr", "iqlr", "zero", "user")))
    stop("method must be one of clr, alr, iqlr, zero, or user")
  if (!(class(shift) %in% c("numeric", "integer")))
    stop("shift must be a number")


  ret <- rRNA_obj
  transdata <- ret$e_data
  edata_cname <- attributes(rRNA_obj)$cnames$edata_cname
  cind <- which(colnames(transdata) == edata_cname)

  # remove any OTUs containing all zeros
  if (!include_zero_rows) {
    Not_zero <- transdata[,-cind] != 0
    ind_allzero <- which(rowSums(Not_zero) == 0)
    if (length(ind_allzero > 0)) {
      if (method %in% c("user", "alr")) {
        if (any(basis_ind %in% ind_allzero))
          stop("At least one of selected OTUs for basis_ind is all zeros")

      }
      transdata <- transdata[-ind_allzero, ]
      if(!quiet) {
        message(paste0("Removing the following ", length(ind_allzero), " row(s) containing all zeros:\n",
                      paste(ind_allzero, collapse = ", ")))
        attributes(ret)$OTUs_removed <-
          as.character(unlist(ret$e_data[ind_allzero, cind], use.names = FALSE))
      }
    }
  }

  # shift data
  transdata[,-cind] <- transdata[,-cind] + shift


  # # If user supplies a distribution a function will be called to do the sampling,
  # # something like this...
  # if (sampleMC) {
  #   transdata <- get_MC_samples_fun(transdata, MCdistn, niter, track_seed)
  # }

  # log transform
  transdata[,-cind] <- log2(transdata[,-cind])


  # CLR transformation
  if (method == "clr") {
    transdata[,-cind] <- apply(transdata[,-cind], 2, function(x){x - mean(x)})


  # ALR transformation
  } else if (method == "alr") {

    # check for errors
    if (is.null(basis_ind))
      stop("Please supply an integer for basis_ind to specify which OTU will be used as the basis")
    if (length(basis_ind) != 1)
      stop("basis_ind must be of length 1 when 'alr' method is used")

    # use an OTU as a basis to scale the others
    transdata[,-cind] <- apply(transdata[,-cind], 2, function(x){x - x[basis_ind]})

    # OTU used for basis is taken out
    transdata <- transdata[-basis_ind,]


  # IQLR transformation
  } else if (method == "iqlr") {
    # first run CLR
    tempdata <- apply(transdata[,-cind], 2, function(x){x - mean(x)})

    # calculate variance of OTUs
    otu_var <- apply(tempdata, 1, var)

    # get IQR
    q <- quantile(otu_var, c(0.25, 0.75), names = FALSE)

    # get indices for the data with which to calculate the geometric mean
    geom_ind <- which(otu_var >= q[1] & otu_var <= q[2])

    # get means
    new_means <- colMeans(transdata[geom_ind, -cind])

    # scale data
    transdata[,-cind] <- transdata[,-cind] - matrix(new_means, nrow = nrow(transdata),
                                              ncol = length(new_means), byrow = TRUE)


  # Zero transformation
  } else if (method == "zero") {
    if (is.null(attributes(ret)$group_DF))
      stop("No group designated for rRNA_obj. See ?group_designation for details.")

    # get group designation and sample ID
    gdf <- attributes(ret)$group_DF
    samp_id <- attributes(ret)$cname$fdata_cname

    # define function for finding nonzero OTUs within each group and calculating sample means
    nonzero_trans <- function(group) {

      # subset by the group provided
      temp <- subset(gdf, Group == group)
      samps <- as.character(temp[, samp_id])

      # find the rows within that group (of original data) that aren't zero
      # ind <- which(rowSums(ret$e_data[, samps]) != 0)

      if (length(ind_allzero) == 0) {
        Not_zero <- ret$e_data[, samps] != 0
      } else {
        Not_zero <- ret$e_data[-ind_allzero, samps] != 0
      }

      # calculate the means of each sample with only the nonzero rows
      cmeans <- colMeans(transdata[rowSums(Not_zero) > 0, samps])

      return(cmeans)
    }


    # get all means and order them properly
    temp_means <- unlist(lapply(unique(gdf$Group), nonzero_trans))
    nm <- data.frame(Names = names(temp_means), Values = temp_means, row.names = NULL)
    new_means <- nm$Values[match(colnames(transdata[,-cind]), nm$Names)]

    # scale data
    transdata[,-cind] <- transdata[,-cind] - matrix(new_means, nrow = nrow(transdata),
                                                    ncol = length(new_means), byrow = TRUE)


  # User transformation
  } else if (method == "user") {
    if (is.null(basis_ind))
      stop("basis_ind cannot be NULL for this method, please supply a vector of indices")

    # get means
    new_means <- colMeans(transdata[basis_ind, -cind])

    # scale data
    transdata[,-cind] <- transdata[,-cind] - matrix(new_means, nrow = nrow(transdata),
                                                    ncol = length(new_means), byrow = TRUE)
  }


  ret$e_data <- transdata

  # add attribute
  attributes(ret)$data_info$transformation_method <- method

  return(ret)
}




