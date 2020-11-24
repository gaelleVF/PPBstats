# 0. help -----------------------------------------------------------------
#' Analysis of mixtures performances 
#' 
#' @param res_model output from \code{\link{check_model.fit_model_bh_intra_location}} function
#' 
#' @param data_mixture output from \code{\link{format_data_PPBstats.data_mixture.R}} function
#' A list of 3 elements
#' \itemize{
#'     \item Mixture_query
#'     \item Mixture_classic
#'     \item Mixture_S
#' }
#' 
#' 
#' @return results from the model including weighted mean of components for each mixture
#' 
#' @author Gaelle Van Frank
#' 
#' @seealso \code{\link{PPBstats::analyse.outputs}}, \code{\link{shinemas2R::get.data}}
#' 
#' @import qdapRegex
#' 
#' @export
#' 
analyse_mixtures = function(res_model,
                            data_mixture,
                            type = "comparison"   # comparison ou modalities, à retirer pour le package
)
{
  

  #1. Error message
  
  
  #2. define parameters
  
  if("PPBstats" %in% class(res_model)){
    if("check_model_bh_intra_location" %in% class(res_model)){param <- "mu"}
    if("check_model_bh_variance_intra" %in% class(res_model)){param <- "sigma"}
  }else{
    stop("res_model must come from PPBstats::check_model_bh_intra_location or PPBstats::check_model_bh_variance_intra")
  }
  
  # Select analysis to make
  if(type == "comparison"){get_distrib <- TRUE}else{get_distrib <- FALSE}
  
  ## fonction à intégrer dans package ==> retirer tout ce qui est fait avec get_modalities == TRUE
  if(type == "modalities"){get_modalities <- TRUE}else{get_modalities <- FALSE}
  
  # 3. Analysis
  
  #if(get_distrib){
    
    # 1. Récupérer les données pour chaque mélange
    if(get_distrib){d_env <- plyr:::splitter_d(data_mixture$data_all, .(germplasm))}
    if(get_modalities){d_env <- plyr:::splitter_d(data_mixture$data_all, .(expe_melange))} 
    
    distrib <- mclapply(d_env, function(d_mix){
      if(length(which(!is.na(d_mix$expe_melange))) == 0){return(NULL)}
        
        #For each mixture, keep only 1st year of cultivation
        d_env_yr <- plyr:::splitter_d(d_mix, .(year))
        
        if(length(d_env_yr) < 2){return(NULL)}
        
        if(get_distrib){d_yr <- list(d_env_yr[[2]])}   # First year is mixture event year : for get_distrib
        if(get_modalities){d_yr <- d_env_yr}
        
        ## Attention : ne pas oublier de faire aussi pour les environnements où ça n'a pas convergé !!
        res_yr <- lapply(d_yr, function(x){
          x$entry <- unlist(lapply(as.character(x$seed_lot), function(x) strsplit(x,"_")[[1]][1]))
          x$ID <- paste(param,"[", x$entry, ",", x$location, ":", x$year, "]", sep="")
          res_loc <-res_model$MCMC[, grep(unique(paste(x$location, x$year,sep=":")),names(res_model$MCMC))]
          
          # Select results for the mixture and its components
          res_mix <- res_loc[,colnames(res_loc) %in% x$ID]
          if(class(res_mix) == "numeric"){
            res_mix <- as.data.frame(res_mix)
            colnames(res_mix) <- unique(colnames(res_loc)[which(colnames(res_loc) %in% x$ID)])
          }
          if(ncol(res_mix) == 0){return(NULL)}
            
          if(get_distrib){
            comp <- data_mixture$data_mix[data_mixture$data_mix$germplasm %in% x$germplasm,]
            comp <- unique(comp[which(!is.na(comp$father_germplasm)),])
            res_comp <- res_loc[unlist(rm_between(colnames(res_loc), "[", ",", extract=TRUE)) %in% comp$entry]
            
            # A retirer pour le package
            if(ncol(res_comp) < length(unique(comp$entry))){
              res_comp <- res_loc[grep(paste(comp$sl_father_mod2,collapse="|"), names(res_loc))]
            }
            if(ncol(res_comp) < length(unique(comp$entry))){
              res_comp <- res_loc[grep(paste(comp$sl_father_mod1,collapse="|"), names(res_loc))]
            }
            
            # Pb si on a plusieurs sélections différentes --> à retirer pour le package !!
            for (x in comp$father_germplasm){
              b <- grep(x, names(res_comp))
              if(length(b) > 1){
                to_add <- apply(res_comp[,b], 1, mean)
                res_comp <- cbind(res_comp, to_add)
                nom <- names(res_comp)[b[1]]
                res_comp <- res_comp[,-b]
                names(res_comp)[ncol(res_comp)] = nom
              }
            }
            
            if(ncol(res_comp) == length(unique(comp$entry))){missingcomp = FALSE}else{missingcomp = TRUE}
            
            if(!missingcomp){
              
              if(ncol(res_comp) == 0){return(NULL)}
              
              if("proportion" %in% colnames(comp)){
                prop <- na.omit(comp[,c("father_germplasm", "proportion")])
                noms <- names(res_comp)
                noms <- unlist(lapply(noms, function(x) strsplit(strsplit(strsplit(x,"[[]")[[1]][2],",")[[1]][1],"#")[[1]][1]))
                vec <- prop[match(prop$father_germplasm, noms),"proportion"]
                
                MeanComp = apply(res_comp, 1, function(x){return(sum(vec*x))})
              }else{
                MeanComp = apply(res_comp, 1, mean)
              }
              M <- cbind(res_mix,MeanComp,res_comp)
              colnames(M)[colnames(M) %in% "MeanComp"] = paste(param,"[", "MeanComp",",",unique(na.omit(x$location)),":",unique(na.omit(x$year)),"]",sep="")
              
            }else{
              M = cbind(res_mix, res_comp)
            }
            
            if(ncol(M) < 2){return(NULL)}
            
            # Get results from model
            res_check = list("MCMC" = M, "MCMC_conv_not_ok" = NULL, "data_env_with_no_controls" = NULL, "data_env_whose_param_did_not_converge" = NULL, "data_ggplot" = NULL)
            class(res_check) <- c("PPBstats", "check_model_bh_intra_location")
            comp_mu <- mean_comparisons(res_check, parameter = param)
            
            return(comp_mu$data_mean_comparisons[[1]])
            
          } #end if get_distrib
          
          if(get_modalities){
            # Retirer les sélections de l'année en cours
            a <- data_mixture$data_selection[data_mixture$data_selection$son %in% x$seed_lot,]
            if(nrow(a) > 0){
              b <- as.character(a[grep("bouquet",a$sl_statut),"son"])
              if(length(b) > 0){
                d <- x[as.character(x$seed_lot) %in% b,"ID"]
                x <- x[!(as.character(x$seed_lot) %in% b),]
                res_mix <- res_mix[!(names(res_mix) %in% d)]
              }
            }
            if(ncol(res_mix) < 2){return(NULL)}
            
            res_check = list("MCMC" = res_mix, "MCMC_conv_not_ok" = NULL, "data_env_with_no_controls" = NULL, "data_env_whose_param_did_not_converge" = NULL, "data_ggplot" = NULL)
            class(res_check) <- c("PPBstats", "check_model_bh_intra_location")
            comp_mu <- mean_comparisons(res_check, parameter = param)$data_mean_comparisons[[1]]
            
            comp_mu$mean.comparisons$modalite = unlist(lapply(as.character(comp_mu$mean.comparisons$entry),function(y){
              if(length(grep("[.]2",y)) == 1){
                if(length(grep("#JB",y)) == 1){
                  return("Mélange issu 1 année sélection \n  dans composantes puis 1 année sélection\n  dans mélange (Mod2)")
                }else{
                  return("Mélange issu 1 année sélection \n dans composantes (Mod2)")
                }
              }
              if(length(grep("#B",y)) == 1){
                if(length(grep("#BB",y)) == 1){
                  return("Mélange sélectionné 2 années (Mod3)")
                }else{
                  return("Mélange sélectionné 1 année (Mod3)")
                }
              }
              if(length(grep("[.]3",y)) == 1){return("Mélange issu 2 années sélection \n dans composantes (Mod1)")}
              if(length(grep("[.]2",y)) == 0 & length(grep("#B",y)) == 0 &  length(grep("[.]3",y)) == 0){return("Mélange non sélectionné (Mod4)")}
            }))
            return(comp_mu)
            
            
          } # end if get_modalities
          
        }) # end lapply(d_yr)

  }, mc.cores = 3) # end lapply

    distrib <- distrib[!sapply(distrib, is.null)]
    return(distrib)

  
  
} # end function
