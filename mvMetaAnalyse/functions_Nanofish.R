load_prep_df <- function(file){
  df_original <- read.csv(file,  header = TRUE, stringsAsFactors = TRUE, na.strings = c("NA"))
  
  # check classes
  # str(df_original)
  
  df_original$Outcome = factor(df_original$Outcome,
                               levels = c("Neurochemical outcomes: AChE activity", 
                                          "Neurochemical outcomes: CAT activity",
                                          "Neurochemical outcomes: GPx activity",
                                          "Neurochemical outcomes: GSH content",
                                          "Neurochemical outcomes: GST activity",
                                          "Neurochemical outcomes: Lipid peroxidation",
                                          "Neurochemical outcomes: ROS levels",
                                          "Neurochemical outcomes: SOD activity",
                                          "Motor function: Distance",
                                          "Sensory-motor function: Distance",
                                          "Feeding behaviour: Predatory performance",
                                          "Feeding behaviour: Feeding time",
                                          "Feeding behaviour: Speed during feeding"
                               ))
  
  # choose classes
  df_original$StudyID      <- factor(df_original$StudyID)
  df_original$StudyCluster <- factor(df_original$StudyCluster)
  df_original$EScluster    <- factor(df_original$EScluster)
  df_original$DevelopmentalStage<-factor(sub(".* ", "", 
                                             df_original$DevelopmentalStage), # remove numbers
                                         levels = c("Embrio","Larva","Juvenile","Adult","Unclear"))
  
  # numeric variables transformation (if needed)
  df_original$ExposureDurationDays<-log10(df_original$ExposureDurationDays+1) 
  df_original$Concentration_mg_L  <-log10(df_original$Concentration_mg_L+1) 
  df_original$ParticleSizeMean    <-log10(df_original$ParticleSizeMean+1)
  
  df_original <- select(df_original, c("StudyCluster", "Study", "Year","EScluster", 
                                       "Species","Sex","DevelopmentalStage",
                                       "ParticleMaterial", "ParticleShape", "ParticleSizeMean", "ParticleSizeCat",
                                       "AdministrationRoute", "Concentration_mg_L",
                                       "ExposureDurationCat", "ExposureDurationDays",
                                       "Outcome",
                                       "n_ctrl", "mean_ctrl", "STD_ctrl",
                                       "n_treat", "mean_treat", "STD_treat"))
  return(df_original)
}

set_plot_dimensions <- function(width_choice, height_choice) {
  options(repr.plot.width=width_choice, repr.plot.height=height_choice)
}



### Forest for REmodel and null model
forest_plot <- function(model,outcome,modelID){
  devstage_list = as.character(unique(model$data$DevelopmentalStage))
  devstage_n =    n_distinct(model$data$DevelopmentalStage)
  devstage_rows = vector("list",devstage_n)
  devstage_shade_rows =c()
  k=1
  for (x in devstage_list){
    if (is.na(x))
    {devstage_rows[[k]] = which(rev(is.na(model$data$DevelopmentalStage)))+devstage_n-k
    } else {
      devstage_rows[[k]]= which(rev(model$data$DevelopmentalStage)==x)+devstage_n-k}
    if (k==1|k==3|k==5)
    {devstage_shade_rows = c(devstage_shade_rows,devstage_rows[[k]],max(devstage_rows[[k]])+1)}
    
    k=k+1
  }
  addlines = switch(devstage_n, 4, 3, 2, 1)
  
  pdf(file=paste('Figures/Forest plot_', str_replace(outcome,':','_'), "_",modelID,'.pdf',sep=""),paper="A4")
  
  ### forest plot with extra annotations
  forestplot<-forest.rma(model, 
                         annotate=TRUE, 
                         addpred=TRUE,
                         showweights=FALSE,
                         header=c("Author (s) and Year"),
                         #xlim = c(min(df$yi-df$vi),max(df$yi+df$vi))*c(10,2),
                         #alim, olim, ylim, at, 
                         ylim=c(-1.5,max(unlist(devstage_rows))+devstage_n+addlines),
                         rows=unlist(devstage_rows),
                         steps=5,
                         level=model$level,
                         refline=0,
                         digits=2, 
                         width=(1),
                         xlab = paste("Standartized Mean Difference for",
                                      as.character(outcome)), 
                         slab = Study,
                         mlab = "RE null model",
                         ilab = cbind(
                           #as.character(DevelopmentalStage), 
                           format(10^ExposureDurationDays,digits=1), 
                           as.character(ParticleMaterial),
                           format(round(10^ParticleSizeMean,2),trim=TRUE),
                           format(round(Concentration_mg_L,1), trim=TRUE)
                         ),
                         shade = devstage_shade_rows
                         #colout = df$StudyCluster
  )
  
  text(x=forestplot$ilab.xpos, 
       y=forestplot$ylim[2]-1, 
       c(
         #"Develop.\n stage", 
         "Exp.\n days", 
         "Particle \n material",
         "Particle \n size (µm)",
         "Conc. \n log₁₀[mg/L]"
       ),
       cex=forestplot$cex*0.75,
       font = 2)
  text(x=forestplot$textpos[1],
       y=setdiff(1:1:(max(unlist(devstage_rows))+1),unlist(devstage_rows)),
       rev(devstage_list),
       font=4,
       cex=forestplot$cex*0.8,
       pos=4)
  text(x=forestplot$textpos[1], 
       y=max(unlist(devstage_rows))+2,#forestplot$ylim[1]*2,
       paste("Clusters / Studies / Experiments:",
             paste(n_distinct(model$data$StudyCluster),n_distinct(model$data$Study),n_distinct(model$data$EScluster),sep = " / ")),
       font = 2,
       cex=forestplot$cex*0.75,
       pos=4)
  dev.off() 
  
  return(forestplot)}


### Orchard plot for REmodel
orchardRE_model <- function(model,outcome,I_yposition,heterogeneity){

  switch (length(heterogeneity[[outcome]]),
                 {1},
          
                 {info = paste("italic(I)^{2} ==",heterogeneity[[outcome]][1])
                 I_xposition=0.8},
          
                 {info = paste(list("italic(I)[Total]^{2} ==",
                                    "italic(I)[Study]^{2} == ",
                                    "italic(I)[EScluster]^{2} == "),heterogeneity[[outcome]])
                 I_xposition=c(0.6,0.8,0.7)
                 })
  
  orchard_plot(model, mod = 1,"DevelopmentalStage", xlab = "Standardised mean difference", 
               transfm = "none")+
    annotate(geom="text", x=I_xposition, y=  I_yposition[[outcome]], 
             label= info, parse = TRUE)+
    scale_fill_manual(values="grey") +
    scale_colour_manual(values="grey")}
