# 0. help -----------------------------------------------------------------
#' Summary of of mixtures performances from \code{\link{analyse_mixture}} object
#' 
#' @param data output from \code{\link{analyse_mixture}} function. 
#' If \code{correl_between_traits} is TRUE then data must be a list of 2 elements, each containing the output from  \code{\link{analyse_mixture}} function for the traits the correlation is to be done for.
#' 
#' @param type The type of table to be returned. It can be :
#' \itemize{
#'  \item "distribution" : contains the number of mixtures, the propotion of mixtures higher than the weighted mean of their components, the number of significant positive and negative overyieldings,
#'  the proportion of mixtures lower than their lowest components and higher than their highest components, the mean value of mixtures and the mean value of the weighted mean of their components,
#'  the mean overyielding and pvalue of the comparison to 0, and overyieldings Std. Dev.
#'  \item "correlation" : contains the correlations of overyielding with the number of components, with the weighted variance of components, and with the difference between the highest and lowest component. 
#'  If correl_between_traits is TRUE, the correlation of the overyieldings of the 2 traits, and the correlation between the overyielding of trait 1 with the weighted variance of components of trait 2
#'  \item "selection_modalities" : a supprimer du package, uniquement pour expe mélange !
#'  \item "DS_RS" : à supprimer du package, uniquement pour expe melange !
#' }
#' 
#' @return Summary of results from \code{\link{analyse_mixture}}
#' 
#' @author Gaelle van Frank
#' 
#' @seealso \code{\link{PPBstats::analyse_mixture}}
#' 
#' @import Hmisc
#' 
#' @export
#' 

PPBsummary.analyse_mixture = function(data, type = "distribution"){

  #1. Check arguments
  if("analyse_mixture" %in% class(data) & type == "correlation"){
    correl_between_traits = FALSE
  }else{
    correl_between_traits = TRUE
  }
 
  if(type == "distribution"){
    if(shapiro.test(as.numeric(as.character(data$global$overyielding)))$p.value <0.05){
      Signif = wilcox.test(as.numeric(as.character(data$global$overyielding)),mu=0)$p.value
    }else{
      Signif = t.test(as.numeric(as.character(data$global$overyielding)),mu=0)$p.value
    }
    
    a = data.frame("number_mixtures" = nrow(data$global),
                   "Proportion_blend_sup_components_mean" = round(sum(as.numeric(as.character(data$global$value_mix))>as.numeric(as.character(data$global$value_mean_comp)))*100/nrow(data$global),3), 
                   "Number_significant_positive_overyieldings" = sum(as.numeric(as.character(data$global$pvalue))<=0.05 & as.numeric(as.character(data$global$overyielding))>0),
                   "Number_significant_negative_overyieldings" = sum(as.numeric(as.character(data$global$pvalue))<=0.05 & as.numeric(as.character(data$global$overyielding))<0),
                   "Proportion_blend_inf_lowest_component" = round(sum(as.numeric(as.character(data$global$value_mix))<as.numeric(as.character(data$global$low_comp)))*100/nrow(data$global),3),
                   "Proportion_blend_highest_component" = round(sum(as.numeric(as.character(data$global$value_mix))>as.numeric(as.character(data$global$high_comp)))*100/nrow(data$global),3),
                   "Weighted_mean_components_value" = mean(as.numeric(as.character(data$global$value_mean_comp))),
                   "Mean_Mixtures_value" = mean(as.numeric(as.character(data$global$value_mix))), 
                   "Overyielding" = mean(as.numeric(data$global$overyielding)),
                   "Std_Dev_overyielding" = sd(as.numeric(data$global$overyielding)),
                   "pvalue_overyielding" = Signif,
                   "stars" = get_stars(Signif)
    )
    return(a)
  } #end if table.type == distribution
  
  if(type == "correlation"){
    if(!correl_between_traits){
      # 1. overyielding ~ number of components
      cor_nb_comp <- data.frame("R" = min(rcorr(data$global$overyielding, data$global$nb_comp)$r), 
                                "pvalue" = rcorr(data$global$overyielding, data$global$nb_comp)$P[1,2], 
                                "n" = max(rcorr(data$global$overyielding, data$global$nb_comp)$n)
      )
      # 2. overyielding ~ variabilité composantes
      cor_var_comp <- data.frame("R" = min(rcorr(data$global$overyielding, data$global$weighted_var_comp)$r), 
                                 "pvalue" = rcorr(data$global$overyielding, data$global$weighted_var_comp)$P[1,2], 
                                 "n" = max(rcorr(data$global$overyielding, data$global$weighted_var_comp)$n)
      )
      
      
      # 3. overyielding ~ différence entre composantes extrêmes
      cor <- rcorr(data$global$overyielding, as.numeric(as.character(data$global$high_comp))-as.numeric(as.character(data$global$low_comp)))
      cor_diff_comp <- data.frame("R" = min(cor$r), 
                                  "pvalue" = cor$P[1,2], 
                                  "n" = max(cor$n)
      )
      
      correl <- rbind(cor_nb_comp, cor_var_comp, cor_diff_comp)
      rownames(correl) <- c("cor_nb_comp","cor_weightedvar_comp","cor_diffextr_comp")
    }# end if !correl_between_traits
    
    if(correl_between_traits){
      D <- merge(data[[1]]$global[,c("location","year","mixture","overyielding")], data[[2]]$global[,c("location","year","mixture","overyielding","weighted_var_comp")], by=c("mixture","location","year"))
      
      # 4. overyieldings entre caractères
      cor_overyiel_traits <- data.frame("R" = min(rcorr(D$overyielding.x, D$overyielding.y)$r), 
                                        "pvalue" = rcorr(D$overyielding.x, D$overyielding.y)$P[1,2], 
                                        "n" = max(rcorr(D$overyielding.x, D$overyielding.y)$n)
      )
      # 5. overyielding caractère 1 ~ variabilité composantes caractère 2
      cor_traits_comp <- data.frame("R" = min(rcorr(D$overyielding.x, D$weighted_var_comp)$r), 
                                    "pvalue" = rcorr(D$overyielding.x, D$weighted_var_comp)$P[1,2], 
                                    "n" = max(rcorr(D$overyielding.x, D$weighted_var_comp)$n)
      )
      
      correl <- rbind(cor_overyiel_traits, cor_traits_comp)
      rownames(correl) <- c("cor_overyied_traits","cor_overyield1_weightedvar_comp2")
    } # end if correl_between_traits
    
    return(correl)
  } #end if table.type == corrélation
  
  
  if(type == "selection_modalities"){
    Tab <- lapply(c("all","all_but_1",unique(data$global$mod_year)), function(yr){
      if(yr == "all"){
        D <- data$global
      }else if(yr == "all_but_1"){
        D <- data$global[data$global$mod_year != 1,]
      }else{
        D <- data$global[data$global$mod_year == yr,]
      }
      m1 <- ifelse(length(na.omit(D$gain_M1))>2, ifelse(shapiro.test(D$gain_M1)$p.value<=0.05,t.test(D$gain_M1,mu=0)$p.value,wilcox.test(D$gain_M1,mu=0)$p.value),NA)
      m2 <- ifelse(length(na.omit(D$gain_M2))>2, ifelse(shapiro.test(D$gain_M2)$p.value<=0.05,t.test(D$gain_M2,mu=0)$p.value,wilcox.test(D$gain_M2,mu=0)$p.value),NA)
      m3 <- ifelse(length(na.omit(D$gain_M3))>2, ifelse(shapiro.test(D$gain_M3)$p.value<=0.05,t.test(D$gain_M3,mu=0)$p.value,wilcox.test(D$gain_M3,mu=0)$p.value),NA)
      
      m12 <- ifelse(length(na.omit(D$gain_M12))>2, ifelse(shapiro.test(D$gain_M12)$p.value<=0.05,t.test(D$gain_M12,mu=0)$p.value,wilcox.test(D$gain_M12,mu=0)$p.value),NA)
      m23 <- ifelse(length(na.omit(D$gain_M23))>2, ifelse(shapiro.test(D$gain_M23)$p.value<=0.05,t.test(D$gain_M23,mu=0)$p.value,wilcox.test(D$gain_M23,mu=0)$p.value),NA)
      m13 <- ifelse(length(na.omit(D$gain_M13))>2, ifelse(shapiro.test(D$gain_M13)$p.value<=0.05,t.test(D$gain_M13,mu=0)$p.value,wilcox.test(D$gain_M13,mu=0)$p.value),NA)
      
      return(data.frame("year" = yr,
                        "gain_M1" = mean(na.omit(D$gain_M1)),
                        "pval_M1" = m1,
                        "stars_M1" = get_stars(m1),
                        "n_M1" = length(na.omit(D$gain_M1)),
                        "gain_M2" = mean(na.omit(D$gain_M2)),
                        "pval_M2" = m2,
                        "stars_M2" = get_stars(m2),
                        "n_M2" = length(na.omit(D$gain_M2)),
                        "gain_M3" = mean(na.omit(D$gain_M3)),
                        "pval_M3" = m3,
                        "stars_M3" = get_stars(m3),
                        "n_M3" = length(na.omit(D$gain_M3)),
             
                         "gain_M12" = mean(na.omit(D$gain_M12)),
                         "pval_M12" = m12,
                         "stars_M12" = get_stars(m12),
                         "n_M12" = length(na.omit(D$gain_M12)),
                         "gain_M23" = mean(na.omit(D$gain_M23)),
                         "pval_M23" = m23,
                         "stars_M23" = get_stars(m23),
                         "n_M23" = length(na.omit(D$gain_M23)),
                         "gain_M13" = mean(na.omit(D$gain_M13)),
                         "pval_M13" = m13,
                         "stars_M13" = get_stars(m13),
                          "n_M13" = length(na.omit(D$gain_M13)))
                        )
    })
   Tab <- do.call(rbind, Tab)
   return(Tab)
  }# end if selection_modalities
  
  if(type == "DS_RS"){
    
  }# end if DS_RS
  
}
