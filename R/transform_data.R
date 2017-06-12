
#' Apply a transformation to rRNA data
#'
#' @param rRNA_obj
#' @param method
#'
#' @return
#' @export
#'
#' @examples
transform_data <- function(rRNA_obj, method = "clr", shift = 0.5) {

  if (class(rRNA_obj) != "rRNAdata")
    stop("rRNA_obj must be of class 'rRNAdata.' Use ?as.rRNAdata for more info.")

  dat <- rRNA_obj$e_data

  # CLR transformation
  if(method == "clr") {
    # get ID variable
    edata_cname <- attributes(rRNA_obj)$cnames$edata_cname
    ind <- which(colnames(dat) == edata_cname)

    # log transform all but ID
    ## shift data to get rid of zeros
    num_data <- dat[,-ind] + shift
    log2dat <- log2(num_data)

    # centered log ratio
    temp <- apply(log2dat, 2, function(x){x - mean(x, na.rm = TRUE)})

    # reattach and name ID column
    transdata <- as.data.frame(cbind(dat[,ind], temp), stringsAsFactors = FALSE)
    colnames(transdata)[1] <- edata_cname

    # make sure vectors are numeric
    transdata[,-1] <- apply(transdata[,-1], 2, as.numeric)

    # return to original order of samples in untransformed data frame
    transdata <- dplyr::select(transdata, dplyr::one_of(colnames(dat)))
  }


  # ALR transformation
  if(method == "alr") {
    transdata <- dat
    for(i in 2:ncol(dat)) {
      vec <- dat[,i]
      # get last nonzero entry for a basis
      denom <- rev(which(vec > 0, arr.ind = TRUE))[1]
      # transform data
      transdata[,i] <- log2(vec/denom)
    }
  }

  rRNA_obj$e_data <- transdata
  ret <- rRNA_obj

  # add attribute
  attributes(ret)$data_info$transformation_method <- method

  return(ret)
}


