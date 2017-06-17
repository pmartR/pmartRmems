
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
                           sampleMC = FALSE, MCdistn = "dirichlet", niter = 200) {

  if (class(rRNA_obj) != "rRNAdata")
    stop("rRNA_obj must be of class 'rRNAdata.' See ?as.rRNAdata for more info.")
  if(is.null(method))
    stop("method must be one of clr, alr, or iqlr")
  if (!(tolower(method) %in% c("clr", "alr", "iqlr")))
    stop("method must be one of clr, alr, or iqlr")
  if (!(class(shift) %in% c("numeric", "integer")))
    stop("shift must be a number")
  if (!is.null(basis_ind)) {
    if (!(class(basis_ind) %in% c("numeric", "integer")))
      stop("basis_ind must be a positive integer no greater than the number of rows in e_data")
    if(basis_ind < 1 | basis_ind > nrow(rRNA_obj$e_data) | basis_ind %% 1 != 0)
      stop("basis_ind must be a positive integer no greater than the number of rows in e_data")
  }


  ret <- rRNA_obj
  transdata <- ret$e_data

  # # Here a separate function will be called to get the shifted data
  # transdata <- get_shift_fun(ret$e_data, shift)

  # for now we'll do this...
  edata_cname <- attributes(rRNA_obj)$cnames$edata_cname
  ind <- which(colnames(transdata) == edata_cname)
  transdata[,-ind] <- log2(transdata[,-ind] + shift)

  # # If user supplies a distribution, a function will be called to do the sampling
  # if (sampleMC) {
  #   transdata <- get_MC_samples_fun(transdata, MCdistn, niter)
  # }
  # # probably do the log transformation at this point once the other functions are in
  # transdata[,-ind] <- log2(transdata[,-ind])


  # CLR transformation
  if (tolower(method) == "clr") {
    transdata[,-ind] <- apply(transdata[,-ind], 2, function(x){x - mean(x)})
  }


  # ALR transformation
  if (tolower(method) == "alr") {
    # use an OTU as a basis to scale the others
    denom <- ifelse(is.null(basis_ind), nrow(transdata), basis_ind)
    transdata[,-ind] <- apply(transdata[,-ind], 2, function(x){x - x[denom]})

    # OTU used for basis is taken out
    transdata <- transdata[-denom,]
  }


  # IQLR transformation
  if (tolower(method) == "iqlr") {
    # first run CLR
    tempdata <- apply(transdata[,-ind], 2, function(x){x - mean(x)})

    # calculate variance of OTUs
    otu_var <- apply(tempdata, 1, var)

    # get IQR
    q <- quantile(otu_var, c(0.25, 0.75), names = FALSE)

    # get indices for the data with which to calculate the geometric mean
    geom_ind <- which(otu_var >= q[1] & otu_var <= q[2])

    # get means
    new_means <- colMeans(transdata[geom_ind, -ind])

    # scale data
    transdata[,-ind] <- transdata[,-ind] - matrix(new_means, nrow = nrow(transdata[,-1]),
                                              ncol = length(new_means), byrow = TRUE)
  }

  ret$e_data <- transdata

  # add attribute
  attributes(ret)$data_info$transformation_method <- method

  return(ret)
}


