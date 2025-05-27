influentialcases_plots <- function(m,mname2,mname1,size_factor) {
  df<-m[[mname2]]$data
  
  # measure residuals
  ic_distance <-cooks.distance(m[[mname2]], cluster=df$Study)
  
  ########## Not available for robust rma: ###########
  # Studentized Deleted Residual
  ic_rstudent <- rstudent(m[[mname1]], 2, progbar=FALSE, cluster=df$Study, reestimate=TRUE)
  ####################################################
  
  # to chech order of all data
  # cbind(names(ic_distance),as.character(ic_rstudent$cluster$slab))
  
  study_id_map <- setNames(df$StudyID, df$Study)
  study_id_seq <- study_id_map[names(ic_distance)]
  
  threshold_cook <- 4 / length(ic_distance)
  threshold_out <- 2
  
  df_ic <- data.frame(
    Study = names(study_id_seq),
    StudyID = unname(study_id_seq),
    CookDistance = ic_distance,
    StudentizedResiduals = ic_rstudent$cluster$X2,
    Influential = ic_distance > threshold_cook,  # Limite comum de influência
    Outlier = abs(ic_rstudent$cluster$X2) > threshold_out  # Resíduos >2 ou <-2 são suspeitos
  )
  
  p1 <- ggplot(df_ic, aes(x = StudentizedResiduals, y = Study, color = Outlier)) +
    geom_point(size = df_ic$CookDistance * size_factor) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "#A6CEE3")) + 
    geom_vline(xintercept = c(-1, 1)*threshold_out, linetype = "dashed", color = "red") +
    #geom_text(aes(label = ifelse(Outlier, StudyID, "")), hjust = 0.5, size = 3, color = "black", fontface = "bold") +
    labs(title = "Outliers", x = "Stud. Del. Residual", y = "Study") +
    theme_classic() + theme(legend.position = "none")
  
  p2 <- ggplot(df_ic, aes(x = CookDistance, y = Study, color = Influential)) +
    geom_point(size = df_ic$CookDistance * size_factor) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "#A6CEE3")) + 
    geom_vline(xintercept = threshold_cook, linetype = "dashed", color = "red") +
    #geom_text(aes(label = ifelse(Influential, StudyID, "")), hjust = 0.5, size = 3, color = "black", fontface = "bold") +
    labs(title = "Influential Studies", x = "Cook's Distance", y = "Study") +
    theme_classic() + theme(legend.position = "none",axis.title.y = element_blank(), axis.text.y = element_blank())
  
  p3 <- ggplot(df_ic, aes(x = StudentizedResiduals, y = CookDistance, color = Influential)) +
    geom_point(size = df_ic$CookDistance * size_factor) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "#A6CEE3")) + 
    geom_vline(xintercept = c(-1, 1)*threshold_out, linetype = "dashed", color = "red")+
    geom_hline(yintercept = threshold_cook, linetype = "dashed", color = "red")+
    geom_text(aes(label = ifelse(Influential, Study, "")), hjust = 1, vjust = 1, size = 2, color = "black", fontface = "bold") +
    labs(x = "Stud. Del. Residual", y = "Studies Cook's Distance",size=3) +
    theme_classic()+ theme(legend.position = "none")
  
  layout <- "
ABC
AB.
"
  
  figure <- p3 + p1 + p2 + plot_layout(design = layout)
  return(figure)
}