#' Remove items that need to be filtered out
#'
#' This function removes
#'
#' @param omicsData an object of the class \code{pepData}, \code{proData}, \code{lipidData}, or \code{metabData}, usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.lipidData}}, or \code{\link{as.metabData}}, respectively.
#' @param filter_object a list created by the functions above
#' @return list
#' @author Kelly Stratton, Lisa Bramer
#'
MSomics_filter_worker <- function(filter_object, omicsData){
  # pull column names from omicR_data attributes #
  col_nms = attr(omicsData, "cnames")
  samp_cname = col_nms$fdata_cname
  edata_cname = col_nms$edata_cname
  emeta_cname = col_nms$emeta_cname

  # pull group_DF attribute #
  group_DF = attr(omicsData, "group_DF")

  #check if filter object contains remove arguments
  if(!is.null(filter_object$edata_filt)||!is.null(filter_object$emeta_filt)||!is.null(filter_object$samples_filt)){

    ## check to see if e_meta is provided ##
    # if not provided we only need to worry about e_data and f_data #
    if(attr(omicsData, "meta_info") == FALSE){

      ## remove entries in edata ##
      if(!is.null(filter_object$edata_filt) & !is.null(edata_cname)){

        temp.pep = omicsData$e_data

        # have to check that at least one of the items is present in the data #
        edat_ids = which(temp.pep[,edata_cname] %in% filter_object$edata_filt)

        if(length(edat_ids) > 0){
          # identify which peptides in e_data match filter list and remove #
          temp.pep1 = temp.pep[-which(temp.pep[,edata_cname] %in% filter_object$edata_filt),]
        }else{temp.pep1 = temp.pep}

      }else{ # no entries in edata need to be removed
        temp.pep1 = omicsData$e_data
      }

      ## remove samples ##
      if(!is.null(filter_object$samples_filt) & !is.null(samp_cname)){
        # identify which samples in f_data match filter list #
        temp.samp = omicsData$f_data

        # check that at least one sample is in f_data and e_data #
        fdat_ids = which(temp.samp[,samp_cname] %in% filter_object$samples_filt)
        edat_ids2 = which(names(temp.pep1) %in% filter_object$samples_filt)

        if(length(fdat_ids) > 0){
          temp.samp2 = temp.samp[-which(temp.samp[,samp_cname] %in% filter_object$samples_filt),]
        }else{temp.samp2 = temp.samp}

        # identify which samples in e_data match filter list and remove #
        if(length(edat_ids2) > 0){
          temp.pep2 = temp.pep1[, -which(names(temp.pep1) %in% filter_object$samples_filt)]
        }else{temp.pep2 = temp.pep1}

      }else{ # no entries in f_data need to be removed
        temp.samp2 = omicsData$f_data
        temp.pep2 = temp.pep1
      }

      temp.meta2 = NULL


    }else{ # e_meta is present, so we need to work with it as well
      ## remove entries in edata ##
      if(!is.null(filter_object$edata_filt) & !is.null(edata_cname)){

        temp.pep = omicsData$e_data

        # have to check that at least one of the items is present in the data #
        edat_ids = which(temp.pep[,edata_cname] %in% filter_object$edata_filt)

        if(length(edat_ids) > 0){
          # identify which peptides in e_data and e_meta match filter list and remove#
          temp.pep1 = temp.pep[-which(temp.pep[,edata_cname] %in% filter_object$edata_filt),]
        }else{temp.pep1 = temp.pep}

        temp.meta = omicsData$e_meta

        # check that at least one of the peptides is present in e_meta #
        emeta_ids = which(temp.meta[,edata_cname] %in% filter_object$edata_filt)

        if(length(emeta_ids) > 0){
          temp.meta1 = temp.meta[-which(temp.meta[,edata_cname] %in% filter_object$edata_filt),]
        }else{temp.meta1 = temp.meta}

      }else{
        temp.pep1 = omicsData$e_data
        temp.meta1 = omicsData$e_meta
      }

      ## remove samples ##
      if(!is.null(filter_object$samples_filt) & !is.null(samp_cname)){
        # identify which samples in f_data match filter list #
        temp.samp = omicsData$f_data

        # check that at least one sample is in f_data and e_data #
        fdat_ids = which(temp.samp[,samp_cname] %in% filter_object$samples_filt)
        edat_ids2 = which(names(temp.pep1) %in% filter_object$samples_filt)

        if(length(fdat_ids) > 0){
          temp.samp2 = temp.samp[-which(temp.samp[,samp_cname] %in% filter_object$samples_filt),]
        }else{temp.samp2 = temp.samp}

        # identify which samples in e_data match filter list and remove #
        if(length(edat_ids2) > 0){
          inds = which(names(temp.pep1) %in% filter_object$samples_filt)
          temp.pep2 = temp.pep1[, -inds]
        }else{temp.pep2 = temp.pep1}

      }else{
        temp.samp2 = omicsData$f_data
        temp.pep2 = temp.pep1
      }

      ## remove entries in emeta ##
      if(!is.null(filter_object$emeta_filt) & !is.null(emeta_cname)){
        # identify which proteins in data match filter list and remove from e_meta #
        temp.meta = temp.meta1

        # check that at least one of the proteins is in e_meta #
        if(!is.null(ncol(temp.meta))){
          emeta_ids2 = which(as.character(temp.meta[,emeta_cname]) %in% filter_object$emeta_filt)
        }else{
          emeta_ids2 = which(temp.meta %in% filter_object$emeta_filt)
        }


        if(length(emeta_ids2) > 0){
          if(!is.null(ncol(temp.meta))){
            temp.meta2 = temp.meta[-which(temp.meta[,emeta_cname] %in% filter_object$emeta_filt),]
          }else{
            temp.meta2 = temp.meta[-which(temp.meta %in% filter_object$emeta_filt)]
          }
        }else{temp.meta2 = temp.meta}
      }else{
        temp.meta2 = temp.meta1
      }


      # check for rogue entries in edata #
      if(!is.null(ncol(temp.meta2))){
        edat_ids2 = which(!(temp.pep2[,edata_cname] %in% temp.meta2[,edata_cname]))
      }else{
        edat_ids2 = which(!(temp.pep2[,edata_cname] %in% temp.meta2))
      }


      # filter out edata entries which no longer have mappings to emeta entries #
      if(length(edat_ids2) > 0){
        #temp.pep2 = temp.pep2[-which(!(temp.pep2[,edata_cname] %in% temp.meta2[,edata_cname])),]
        temp.pep2 = temp.pep2[-edat_ids2,]
      }


      # filter out fdata entries which no longer have mappings to edata entries #
      temp.fdata <- omicsData$f_data
      fdat_ids <- which(!(temp.fdata[,samp_cname] %in% names(temp.pep2)[names(temp.pep2) != edata_cname]))
      temp.samp2 <- temp.samp2[-fdat_ids,]

    }

    output <- list(temp.pep2 = temp.pep2, temp.samp2 = temp.samp2, temp.meta1 = temp.meta2, edata_cname = edata_cname, emeta_cname = emeta_cname, samp_cname = samp_cname)
  }


  #if filter object contains keep arguments
  else{

    ## check to see if e_meta is provided ##
    # if not provided we only need to worry about e_data and f_data #
    if(attr(omicsData, "meta_info") == FALSE){

      ## keep entries in edata ##
      if(!is.null(filter_object$edata_keep) & !is.null(edata_cname)){

        temp.pep = omicsData$e_data

        # have to check that at least one of the items is present in the data #
        edat_ids = which(temp.pep[,edata_cname] %in% filter_object$edata_keep)

        if(length(edat_ids) > 0){
          # identify which peptides in e_data match filter list and keep #
          temp.pep1 = temp.pep[which(temp.pep[,edata_cname] %in% filter_object$edata_keep),]
        }else{temp.pep1 = temp.pep}

      }else{ # no entries in edata need to be removed
        temp.pep1 = omicsData$e_data
      }

      ## keep samples ##
      if(!is.null(filter_object$samples_keep) & !is.null(samp_cname)){
        # identify which samples in f_data match filter list #
        temp.samp = omicsData$f_data

        # check that at least one sample is in f_data and e_data #
        fdat_ids = which(temp.samp[,samp_cname] %in% filter_object$samples_keep)
        edat_ids2 = which(names(temp.pep1) %in% filter_object$samples_keep)

        if(length(fdat_ids) > 0){
          temp.samp2 = temp.samp[which(temp.samp[,samp_cname] %in% filter_object$samples_keep),]
        }else{temp.samp2 = temp.samp}

        # identify which samples in e_data match filter list and keep #
        if(length(edat_ids2) > 0){
          edata_cname_id = which(names(temp.pep1) == edata_cname)
          temp.pep2 = temp.pep1[,c(edata_cname_id,(which(names(temp.pep1) %in% filter_object$samples_keep)))]
        }else{temp.pep2 = temp.pep1}

      }else{ # no entries in f_data need to be removed
        temp.samp2 = omicsData$f_data
        temp.pep2 = temp.pep1
      }

      temp.meta2 = NULL




    }else{ # e_meta is present, so we need to work with it as well
      ## keep entries in edata ##
      if(!is.null(filter_object$edata_keep) & !is.null(edata_cname)){

        temp.pep = omicsData$e_data

        # have to check that at least one of the items is present in the data #
        edat_ids = which(temp.pep[,edata_cname] %in% filter_object$edata_keep)

        if(length(edat_ids) > 0){
          # identify which peptides in e_data and e_meta match filter list and keep#
          temp.pep1 = temp.pep[which(temp.pep[,edata_cname] %in% filter_object$edata_keep),]
        }else{temp.pep1 = temp.pep}

        temp.meta = omicsData$e_meta

        # check that at least one of the peptides is present in e_meta #
        emeta_ids = which(temp.meta[,edata_cname] %in% filter_object$edata_keep)

        if(length(emeta_ids) > 0){
          temp.meta1 = temp.meta[which(temp.meta[,edata_cname] %in% filter_object$edata_keep),]
        }else{temp.meta1 = temp.meta}

      }else{
        temp.pep1 = omicsData$e_data
        temp.meta1 = omicsData$e_meta
      }

      ## keep samples ##
      if(!is.null(filter_object$samples_keep) & !is.null(samp_cname)){
        # identify which samples in f_data match filter list #
        temp.samp = omicsData$f_data

        # check that at least one sample is in f_data and e_data #
        fdat_ids = which(temp.samp[,samp_cname] %in% filter_object$samples_keep)
        edat_ids2 = which(names(temp.pep1) %in% filter_object$samples_keep)

        if(length(fdat_ids) > 0){
          temp.samp2 = temp.samp[which(temp.samp[,samp_cname] %in% filter_object$samples_keep),]
        }else{temp.samp2 = temp.samp}

        # identify which samples in e_data match filter list and keep #
        if(length(edat_ids2) > 0){
          edata_cname_id = which(names(temp.pep1) == edata_cname)

          inds = which(names(temp.pep1) %in% filter_object$samples_keep)
          temp.pep2 = temp.pep1[,c(edata_cname_id,inds)]
        }else{temp.pep2 = temp.pep1}

      }else{
        temp.samp2 = omicsData$f_data
        temp.pep2 = temp.pep1
      }

      ## keep entries in emeta ##
      if(!is.null(filter_object$emeta_keep) & !is.null(emeta_cname)){
        # identify which proteins in data match filter list and keep in e_meta #
        temp.meta = temp.meta1
        temp.meta_not_kept = omicsData$e_meta[-which(omicsData$e_meta$Mass_Tag_ID %in% temp.meta$Mass_Tag_ID),]

        # check that at least one of the proteins is in e_meta (this is e_meta after e_data_keep has been applied) #
        if(!is.null(ncol(temp.meta))){
          emeta_ids2 = which(as.character(temp.meta[,emeta_cname]) %in% filter_object$emeta_keep)
        }else{
          emeta_ids2 = which(temp.meta %in% filter_object$emeta_keep)
        }

        # check that at least one of the proteins is in temp_meta_not_kept #
        if(!is.null(ncol(temp.meta_not_kept))){
          emeta_ids_not_kept = which(as.character(temp.meta_not_kept[,emeta_cname]) %in% filter_object$emeta_keep)
        }else{
          emeta_ids_not_kept = which(temp.meta_not_kept %in% filter_object$emeta_keep)
        }

        #if there are e_meta_keep (proteins) in the part of e_meta that we previously kept, then these proteins have already been kept
        if(length(emeta_ids2) > 0){
          if(!is.null(ncol(temp.meta))){
            temp.meta2 = temp.meta
          }else{
            temp.meta2 = temp.meta
          }
        }else{temp.meta2 = temp.meta}

        #if there are e_meta_keep(proteins) outside of e_meta that we previously kept we will keep these
        if(length(emeta_ids_not_kept) > 0){
          if(!is.null(ncol(temp.meta_not_kept))){
            temp.meta3 = temp.meta_not_kept[which(temp.meta_not_kept[,emeta_cname] %in% filter_object$emeta_keep),]
            temp.meta2 = rbind(temp.meta2,temp.meta3)
          }else{
            temp.meta3 = temp.meta_not_kept[which(temp.meta_not_kept %in% filter_object$emeta_keep)]
            temp.meta2 = rbind(temp.meta2, temp.meta3)
          }
        }else{temp.meta2 = temp.meta}

      }else{
        temp.meta2 = temp.meta1
      }


      # check for entries in e_meta[,emeta_cname] that are not in e_data[,edata_cname]#
      if(!is.null(ncol(temp.meta2))){
        edat_ids2 = which(!temp.meta2[,edata_cname] %in% (temp.pep2[,edata_cname]))
      }else{
        edat_ids2 = which(!(temp.meta2 %in% temp.pep2[,edata_cname]))
      }


      # add edata entries which were present in emeta but not edata #
      if(length(edat_ids2) > 0){
        additional_peps<- temp.meta2$Mass_Tag_ID[edat_ids2]
        edata_cname_id = which(names(temp.pep1) == edata_cname)

        if(is.null(filter_object$samples_keep)){
          inds = which(names(temp.pep1) %in% temp.samp2[,samp_cname])
        }
        else inds = which(names(temp.pep1) %in% filter_object$samples_keep)

        temp.pep2 = rbind(temp.pep2, omicsData$e_data[which(omicsData$e_data$Mass_Tag_ID %in% additional_peps) ,c(edata_cname_id,inds)])
      }



    }

    output <- list(temp.pep2 = temp.pep2, temp.samp2 = temp.samp2, temp.meta1 = temp.meta2, edata_cname = edata_cname, emeta_cname = emeta_cname, samp_cname = samp_cname)

  }


  # return the pieces needed to assemble a proData/pepData/lipidData/metabData object
  return(output)
}
