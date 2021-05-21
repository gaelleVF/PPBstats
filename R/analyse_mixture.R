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
#' @param type type
#' 
#' @return results from the model including weighted mean of components for each mixture
#' 
#' @author Gaelle Van Frank
#' 
#' @seealso \code{\link{PPBstats::analyse.outputs}}, \code{\link{shinemas2R::get.data}}
#' 
#' @import qdapRegex
#' @import tidyverse
#' 
#' @export
#' 
analyse_mixtures = function(res_model,
                            data_mixture,
                            type = "comparison"   # comparison ou modalities, à retirer pour le package
)
{
  

  ## STILL TO BE DONE
  # Gérer les cas où pas de résultats dans res_model
  # Gérer les résultats varintra
  
  
  
  #1. Functions -----------------
  weightedvar <- function(x, prop, na.rm = FALSE) {
    # Get variance of components
    sum.prop <- sum(prop)
    sum.prop2 <- sum(prop^2)
    mean.prop <- sum(x * prop) / sum(prop)
    (sum.prop / (sum.prop^2 - sum.prop2)) * sum(prop * (x - mean.prop)^2)
  }
  
  get_comp_modalities <- function(p1, p2, data){
    # Get gain and pvalue of mixtures selection modalities
    x1 <- as.character(data$mean.comparisons[grep(p1,data$mean.comparisons$modalite), "parameter"])
    x2 <- as.character(data$mean.comparisons[grep(p2,data$mean.comparisons$modalite), "parameter"])
    if(length(x1) + length(x2) == 2){
      pval <- max(data$Mpvalue[which(rownames(data$Mpvalue) %in% c(x1, x2)),
                               which(colnames(data$Mpvalue) %in% c(x1, x2))])
      gain <- 100*(data$mean.comparisons[data$mean.comparisons$parameter == x1, "median"]/data$mean.comparisons[data$mean.comparisons$parameter == x2, "median"] -1)
      return(data.frame("gain" = gain, "pval" = pval))
    }else{
      return(data.frame("gain" = NA, "pval" = NA))
    }
  }
  
  
  #2. define parameters -----------------
  
  if("PPBstats" %in% class(res_model)){
    if("check_model_bh_intra_location" %in% class(res_model)){param <- "mu"}
    if("check_model_bh_variance_intra" %in% class(res_model)){param <- "sigma"}
  }else{
    stop("res_model must come from PPBstats::check_model_bh_intra_location or PPBstats::check_model_bh_variance_intra")
  }
  
  # Select analysis to make
  if(type == "comparison"){get_distrib <- TRUE}else{get_distrib <- FALSE}
  
  ## fonction à integrer dans package ==> retirer tout ce qui est fait avec get_modalities == TRUE
  if(type == "modalities"){get_modalities <- TRUE}else{get_modalities <- FALSE}
  
  # 3. Analysis --------------------
  ## 3.1 Per mixture ------------------------
  if(get_distrib){d_env <- plyr:::splitter_d(data_mixture$data_all, .(germplasm))}
  if(get_modalities){d_env <- plyr:::splitter_d(data_mixture$data_all, .(expe_melange))} 
  
  distrib <- mclapply(d_env, function(d_mix){
     if(length(which(!is.na(d_mix$expe_melange))) == 0){return(NULL)}
     
     d_env_yr <- plyr:::splitter_d(d_mix, .(year))
     if(length(d_env_yr) < 2){return(NULL)}
     if(get_distrib){d_yr <- list(d_env_yr[[2]])}   # First year is mixture event year : for get_distrib
     if(get_modalities){d_yr <- d_env_yr}
     
     ## Attention : ne pas oublier de faire aussi pour les environnements où ça n'a pas converge !!
     res_yr <- lapply(d_yr, function(x){
       x$entry <- unlist(lapply(as.character(x$seed_lot), function(x) strsplit(x,"_")[[1]][1]))
       x$ID <- paste(param,"[", x$entry, ",", x$location, ":", x$year, "]", sep="")
       res_loc <- res_model$MCMC[, grep(unique(paste(x$location, x$year,sep=":")),names(res_model$MCMC))]
       
       # Delete year's selections
       to_delete <- NULL
       to_delete <- data_mixture$data_selection[data_mixture$data_selection$son %in% x$seed_lot & str_detect(as.character(data_mixture$data_selection$sl_statut),"bouquet"),"son"]
       to_delete <- unlist(lapply(as.character(to_delete), function(y) strsplit(y,"_")[[1]][1]))
       # Delete what is not supposed to be here ... !
       if(length(grep("#JB", x$entry)) == 1 & length(grep("#JB", to_delete)) == 0 & unique(x$location) !="JSG"){
         a <- as.character(x[grep("[.]2", x$entry),"entry"])
         if(length(a) > 1){
           to_delete <- c(to_delete, a[-grep("#JB",a)])
         }
       }
       # Ajouter les sélections de 2019 qui ne sont pas encore dans le fichier mixture_selection... =_=
       if(length(grep("MélangeC.2-JSG#JB_JSG_2019", x$seed_lot)) >0){
         to_delete <- c(to_delete, as.character(x[grep("MélangeC.2-JSG#JB_JSG_2019",x$seed_lot),"entry"]))
       }
       if(length(grep("MélangeC-JSG#BB_JSG_2019", x$seed_lot)) >0){
         to_delete <- c(to_delete, as.character(x[grep("MélangeC-JSG#BB_JSG_2019",x$seed_lot),"entry"]))
       }
#       if(unique(x$location) %in% "JSG" & unique(x$year) %in% 2019){
#         to_delete <- c(to_delete, as.character(x[grep("#BB",x$seed_lot),"seed_lot"]), as.character(x[grep("#JB",x$seed_lot),"seed_lot"]))
#       }
       if(unique(x$location) %in% "RAB" & unique(x$year) %in% 2019){
         to_delete <- c(to_delete, as.character(x[grep("#BC",x$seed_lot),"entry"]))
       }
       # Retirer #VA pour le mélange antiverse qui est une composante
       if(length(grep("#VA",x$entry)) > 0){
         to_delete <- c(to_delete, x[grep("#VA",x$entry), "entry"])
       }
       
       if(length(to_delete) > 0){
         a <- which(x$entry %in% to_delete)
         x <- x[-a,]
        }
       
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
         
         # Pb si on a plusieurs selections differentes --> à retirer pour le package !!
         nom <- NULL
         for (pop in comp$father_germplasm){
           b <- grep(pop, names(res_comp))
           if(length(b) > 1){
             to_add <- apply(res_comp[,b], 1, mean)
             res_comp <- cbind(res_comp, to_add)
             nom <- c(nom, names(res_comp)[b[1]])
             res_comp <- res_comp[,-b]
             names(res_comp)[ncol(res_comp)] = nom[length(nom)]
             nom[length(nom)] <- unlist(rm_between(nom[length(nom)], "[", ",", extract=TRUE))
           }else{
             nom <- c(nom, ifelse(length(names(res_comp)[grep(pop, names(res_comp))]) > 0, unlist(rm_between(names(res_comp)[grep(pop, names(res_comp))], "[", ",", extract=TRUE)), 
                                  comp[grep(pop, comp$father_germplasm),"entry"]))
           }
         }
         comp$entry <- nom
         
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
         comp_mu <- mean_comparisons(res_check, parameter = param)$data_mean_comparisons[[1]]
         comp_mu$mean.comparisons <- merge(comp_mu$mean.comparisons, comp[,c("entry","group","proportion")], by="entry", all=T)
         comp_mu$mean.comparisons[comp_mu$mean.comparisons$entry == "MeanComp","group"] = "MC"
         comp_mu$mean.comparisons[is.na(comp_mu$mean.comparisons$group),"group"] = "M"
         comp_mu$mean.comparisons$germplasm <- unlist(lapply(as.character(comp_mu$mean.comparisons$entry), function(x) strsplit(x,"#")[[1]][1]))
         
         return(comp_mu)
         
       } #end if get_distrib
       
       if(get_modalities){
         if(FALSE){
           # Déjà fait plus haut !
           # Retirer les selections de l'annee en cours
           a <- data_mixture$data_selection[data_mixture$data_selection$son %in% x$seed_lot,]
           if(!is.null(a)){
             b <- as.character(a[grep("bouquet",a$sl_statut),"son"])
             if(length(b) > 0){
               d <- x[as.character(x$seed_lot) %in% b,"ID"]
               x <- x[!(as.character(x$seed_lot) %in% b),]
               res_mix <- res_mix[!(names(res_mix) %in% d)]
             }
           }
         }
        
         if(ncol(res_mix) < 2){return(NULL)}
         
         # Régler pb FRC 2 sélections dans mélange --> moyenne des sélections
         if(length(grep("#BA", names(res_mix)))>1){
           a <- data.frame(apply(res_mix[grep("#BA", names(res_mix))], 1, mean))
           names(a) = names(res_mix)[grep("#BA", names(res_mix))][1]
           res_mix <- res_mix[-grep("#BA", names(res_mix))]
           res_mix <- cbind(res_mix,a)
         }
         
         res_check = list("MCMC" = res_mix, "MCMC_conv_not_ok" = NULL, "data_env_with_no_controls" = NULL, "data_env_whose_param_did_not_converge" = NULL, "data_ggplot" = NULL)
         if(param == "mu"){class(res_check) <- c("PPBstats", "check_model_bh_intra_location")}
         if(param == "sigma"){class(res_check) <- c("PPBstats", "check_model_bh_variance_intra")}
         comp_mu <- mean_comparisons(res_check, parameter = param)$data_mean_comparisons[[1]]
         
         comp_mu$mean.comparisons$modalite = unlist(lapply(as.character(comp_mu$mean.comparisons$entry),function(y){
           if(length(grep("[.]2",y)) == 1){
             if(length(grep("#JB",y)) == 1){
               return("Melange issu 1 annee selection \n  dans composantes puis 1 annee selection\n  dans melange (Mod2)")
             }else{
               return("Melange issu 1 annee selection \n dans composantes (Mod2)")
             }
           }
           if(length(grep("#B",y)) == 1){
             if(length(grep("#BB",y)) == 1){
               return("Melange selectionne 2 annees (Mod3)")
             }else{
               return("Melange selectionne 1 annee (Mod3)")
             }
           }
           if(length(grep("[.]3",y)) == 1){return("Melange issu 2 annees selection \n dans composantes (Mod1)")}
           if(length(grep("[.]2",y)) == 0 & length(grep("#B",y)) == 0 &  length(grep("[.]3",y)) == 0){return("Melange non selectionne (Mod4)")}
         }))
         return(comp_mu)
         
         
       } # end if get_modalities
       
     }) # end lapply(d_yr)
}, mc.cores = 3) # end lapply

  distrib <- distrib[!sapply(distrib, is.null)]
    
  ## 3.2 Over all mixtures ---------------------------
  
  Distrib_all <- lapply(distrib, function(mel){

    if(get_distrib){
      tab <- mel[[1]]$mean.comparisons
      if(is.null(tab) | length(grep("MC", tab$group)) == 0){return(NA)}
      tab <- tab[order(tab$median),]
      
      pval <- unlist(lapply(as.character(mel[[1]]$mean.comparisons$parameter), function(x){
        mix <- as.character(tab[grep("^M$", tab$group), "parameter"])
        b <- mel[[1]]$Mpvalue[rownames(mel[[1]]$Mpvalue) %in% c(x,mix), colnames(mel[[1]]$Mpvalue) %in% c(x,mix)]
        return(max(b))
      }))
      tab$pval <- pval
      gain <-  100*(as.numeric(as.character(tab[grep("^M$", tab$group), "median"])) - as.numeric(as.character(tab[grep("^MC$", tab$group), "median"])))/as.numeric(as.character(tab[grep("^MC$", tab$group), "median"]))
      vect <- data.frame("location" = as.character(unique(tab$location)),
                         "year" = as.character(unique(tab$year)),
                         "mixture" = tab[grep("^M$", tab$group), "entry"],
                         "modality" = as.character(tab[grep("^M$", tab$group), "group"]),
                         "value_mix" = tab[grep("^M$", tab$group), "median"],
                         "value_mean_comp" = tab[grep("^MC$", tab$group), "median"],
                         "overyielding" = gain,
                         "pvalue" = tab[grep("^MC$", tab$group), "pval"],
                         "low_comp" =  min(tab[grep("^C$", tab$group), "median"]),
                         "high_comp" =  max(tab[grep("^C$", tab$group), "median"]),
                         "nb_comp" = length(grep("^C$", tab$group)),
                        "weighted_var_comp" = weightedvar(tab[tab$group == "C","median"], tab[tab$group == "C","proportion"])   
                          )
    } #end if get_distrib
    
    if(get_modalities){
      mel <- mel[!sapply(mel, is.null)]
      if(length(mel) > 0){
        vect <- lapply(1:length(mel), function(yr){
          mel_yr <- mel[[yr]]
          p14 <- get_comp_modalities("Mod1", "Mod4", mel_yr)
          p24 <- get_comp_modalities("Mod2", "Mod4", mel_yr)
          p34 <- get_comp_modalities("Mod3", "Mod4", mel_yr)
          p12 <- get_comp_modalities("Mod1", "Mod2", mel_yr)
          p13 <- get_comp_modalities("Mod1", "Mod3", mel_yr)
          p23 <- get_comp_modalities("Mod2", "Mod3", mel_yr)
          
          comp_mod <- data.frame("melange" = unique(gsub("[.]2|[.]3|#JA|#JAa|#JAb|#JB|#JBa|#JBb|#BA|#BAa|#BAb|#BB|#BBa|#BBb", "", mel_yr$mean.comparisons$entry)),
                                 "location" = unique(mel_yr$mean.comparisons$location),
                                 "year_mix" = unique(mel_yr$mean.comparisons$year),
                                 "mod_year" = yr,
                                 "gain_M1" = p14$gain,
                                 "pval_M1" = p14$pval,
                                 "gain_M2" = p24$gain,
                                 "pval_M2" = p24$pval,
                                 "gain_M3" = p34$gain,
                                 "pval_M3" = p34$pval,
                                 "gain_M12" = p12$gain,
                                 "pval_M12" = p12$pval,
                                 "gain_M23" = p23$gain,
                                 "pval_M23" = p23$pval,
                                 "gain_M13" = p13$gain,
                                 "pval_M13" = p13$pval,
                                 "M1" = ifelse(length(grep("Mod1", mel_yr$mean.comparisons$modalite)) > 0, mel_yr$mean.comparisons[grep("Mod1", mel_yr$mean.comparisons$modalite), "median"], NA),
                                 "M2" = ifelse(length(grep("Mod2", mel_yr$mean.comparisons$modalite)) > 0, mel_yr$mean.comparisons[grep("Mod2", mel_yr$mean.comparisons$modalite), "median"], NA),
                                 "M3" = ifelse(length(grep("Mod3", mel_yr$mean.comparisons$modalite)) > 0, mel_yr$mean.comparisons[grep("Mod3", mel_yr$mean.comparisons$modalite), "median"], NA),
                                 "M4" = ifelse(length(grep("Mod4", mel_yr$mean.comparisons$modalite)) > 0, mel_yr$mean.comparisons[grep("Mod4", mel_yr$mean.comparisons$modalite), "median"], NA)
          )
          return(comp_mod)
        })
        vect <- do.call(rbind, vect)
      }else{vect <- NULL}
    }#end if get_modalities
    return(vect)
  })

  Distrib_all <- do.call(rbind, Distrib_all)
  
  out <- list("indiv" = distrib,
              "global" = Distrib_all)
  class(out) <- c("PPBstats", "analyse_mixture")
  if(get_distrib){class(out) = c(class(out), "get_distrib")}
  if(get_modalities){class(out) = c(class(out), "get_modalities")}
  
  return(out)

} # end function
