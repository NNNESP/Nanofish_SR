# AChE activity 1        CAT activity 2         GPx activity 3*     GSH content 4  
# GST activity 5*         Lipid peroxidation 6    ROS levels 7       SOD activity 8
# Motor function  9                 Sensory-motor function 10 
# Predatory performance 11        [Feeding time 12]*       [Speed during feeding 13]*

# check classes
str(df_original)
# adjust classes
df_original$StudyID      <- factor(df_original$StudyID)
df_original$StudyCluster <- factor(df_original$StudyCluster)
df_original$EScluster    <- factor(df_original$EScluster)
df_original$ParticleShape<-factor(sub(".* ", "",df_original$ParticleShape)) #remove numbers
df_original$DevelopmentalStage<-factor(
  sub(".* ", "", df_original$DevelopmentalStage), # remove numbers
  levels = c("Embrio","Larva","Juvenile","Adult","Unclear"))
df_original$Outcome <- factor(df_original$Outcome, 
                              levels= c("AChE activity", 
                                        "CAT activity",
                                        "GPx activity",
                                        "GSH content",
                                        "GST activity",
                                        "Lipid peroxidation",
                                        "ROS levels",
                                        "SOD activity",
                                        "Motor function",
                                        "Sensory-motor function",
                                        "Predatory performance",
                                        "Feeding time",
                                        "Speed during feeding"))

# numeric variables transformation
# df_original$ExposureDurationDays<-log10(df_original$ExposureDurationDays+1) 
# df_original$Concentration_mg_L  <-log10(df_original$Concentration_mg_L+1) 
# df_original$ParticleSizeMean    <-log10(df_original$ParticleSizeMean+1)

str(df_original)

outcomeslist = sort(unique(df_original$Outcome))
df_original <- metafor::escalc(measure = "SMD", 
                               n1i = n_treat, 
                               n2i = n_ctrl, 
                               m1i = mean_treat, 
                               m2i = mean_ctrl, 
                               sd1i = STD_treat, 
                               sd2i = STD_ctrl, 
                               data = df_original, vtype = "UB") 
write.csv(df_original,file = "mvMetaAnalysis_SMDs.csv")
outcome <- as.character(outcomeslist[[n_outcome]])

mRF <- list() # random effects formula investigation models
mRE <- list() # random effects multilevel models
mMR <- list() # multilevel metaregression models

colormap <-c("Embrio" = "khaki", 
             "Larva" = "lightgreen", 
             "Juvenile" = "lightskyblue", 
             "Adult" = "salmon", 
             "Unclear" = "plum",
             "ExposureDurationDays"= "cornflowerblue",
             "ParticleSizeMean"= "chartreuse4",
             "Concentration_mg_L"= "#fa0079ff",
             "Aminoplast" = "chocolate4",
             "LDPE" = "#2200faff",
             "PE" = "#0089faff",
             "PET" = "#00faedff",
             "PLA" = "#00fa29ff",
             "PS" = "#faf200ff",
             "PVC" = "#fa9a00ff",
             "Publication year"="black",
             "Variance"="navyblue")
df <- df_original[ which(df_original$Outcome == outcome ),]
#df <- df[order(df$DevelopmentalStage,-df$yi),]

print(outcome)
DT::datatable(df)
# choose a hierarchical structure that represents allow representation of the dependencies between data
mRF[["base"]] <- rma.mv(yi, vi, 
                        random = ~ 1 | Study/EScluster/ESid,
                        data=df, method="ML")

nSt <- mRF[["base"]]$s.nlevels[1] # number of studies
nEC <- mRF[["base"]]$s.nlevels[2] # number of ESclusters
nES <- mRF[["base"]]$s.nlevels[3] # number of independent ES (OBS: there is no repeated measurements included)

if (nES == nSt){ 
  RE_formula = formula(~1 | Study) # one ES by Study
} else if (nEC>nSt & nEC<nES){ 
  RE_formula = formula(~1 | Study/EScluster/ESid) # more than one EScluster by Study
} else{
  RE_formula = formula(~1 | Study/ESid) # only one EScluster by Study
}

print(paste0("Considering number of studies: ",nSt,
             ", number of ES clusters: ",nEC,
             " and number of independent observations of ES: ",nES,
             ", the random effects formula for ",outcome,
             " is: ", deparse(RE_formula)))

REterms <- names(mRF[["base"]]$s.names)

## Overfiting correction check
mRF[["base"]] <- rma.mv(yi, vi, 
                        random = RE_formula,
                        data=df, method="ML")

REterms <- names(mRF[["base"]]$s.names)

if (nES == nSt){
  VCVrho <- 0
  VCV <- df$vi
} else {
  VCVrho <- 0.5
  VCV <- vcalc(vi, cluster = interaction(df$EScluster,df$Study), obs = ESid, data = df, rho = VCVrho)
}

mname1 <-"ML"
mname2 <-paste0(mname1,"_RobVar")

mRE[[mname1]] <- rma.mv(yi, VCV,random = RE_formula,
                        data=df,slab = Study, method="REML",test="t", 
                        dfs="contain", level=95)
mRE[[mname1]] <- additional_results(mRE[[mname1]])

mRE[[mname2]] <- robust(mRE[[mname1]], cluster = df$Study, clubSandwich = TRUE)

mRE[[mname2]] <- additional_results(mRE[[mname2]])

I2_ML_RobVar <- i2_ml(mRE[[mname2]])
#I2_ML_RobVar <- i2_ml(mRE[[mname2]], boot = 1000)

extract_report(mRE[[mname2]], outcome, RE_formula)



df_MR <- select(df, c("Year", "StudyID","Study", "EScluster", "ESid","DevelopmentalStage",
                      "ExposureDurationDays","ParticleSizeMean",
                      "Concentration_mg_L","yi", "vi"))

mname1 <- "Conc+"
mname2<-paste0(mname1,"_RobVar")

if (n_outcome!=10){
  mods_formula <- formula(~  0 + ParticleSizeMean
                          + Concentration_mg_L
                          + ExposureDurationDays*DevelopmentalStage)
}else if (n_outcome==10){ #for Sensory-motor function
  #for Sensory-motor function
  mods_formula <- formula(~ 0 + ParticleSizeMean
                          + Concentration_mg_L
                          + ExposureDurationDays)
}

df_MR <- select(df, c("Year", "StudyID","Study", "EScluster", "ESid","DevelopmentalStage",
                      "ExposureDurationDays","ParticleSizeMean",
                      "Concentration_mg_L","yi", "vi"))

#dropped<-"Adult"
#VCV <- VCV[!df_MR$DevelopmentalStage %in% c(dropped),!df_MR$DevelopmentalStage %in% c(dropped) ]
#df_MR <- df_MR[!df_MR$DevelopmentalStage %in% c(dropped), ]

mname1 <- "Conc+"
mname2<-paste0(mname1,"_RobVar")

mMR[[mname1]] <- rma.mv(yi, VCV,random = RE_formula, mods= mods_formula,
                        data=df_MR,slab = Study, method="REML",test="t",
                        dfs="contain", level=95) #it will only keep complete cases

nSt_MR <- mMR[[mname1]]$s.nlevels[1] # number of studies
nEC_MR <- mMR[[mname1]]$s.nlevels[2] # number of ESclusters/independent ES
nES_MR <- mMR[[mname1]]$s.nlevels[3] # number of independent ES

if (length(mMR[[mname1]]$s.nlevels)==3){
  print(paste0("Number of studies: ",nSt_MR,
               ", number of ES clusters: ",nEC_MR,
               " and number of independent observations of ES: ",nES_MR))
}else{
  print(paste0("Number of studies: ",nSt_MR,
               " and number of independent observations of ES: ",nEC_MR))
}

mMR[[mname2]] <- robust(mMR[[mname1]], cluster = df_MR$Study, clubSandwich = TRUE,verbose=TRUE)

mMR[[mname2]][["i2"]] <- i2_ml(mMR[[mname2]])
mMR[[mname2]][["icc"]] <- round(mMR[[mname2]]$sigma2 / sum(mMR[[mname2]]$sigma2), 3)

if (n_outcome!=10){
  anova(mMR[[mname2]], btt=c("DevelopmentalStage","ExposureDurationDays","ParticleSizeMean","Concentration_mg_L"))
}else if (n_outcome==10){ #for Sensory-motor function
  anova(mMR[[mname2]], btt=c("ExposureDurationDays","ParticleSizeMean","Concentration_mg_L"))
}

########

summary(mRE[["ML_RobVar"]])
mRE[["ML_RobVar"]][["icc"]]
I2_ML_RobVar

var.comp(mRE[["ML"]])

########

summary(mMR[["Conc+"]])
summary(mMR[["Conc+_RobVar"]])
m<-mMR[["Conc+_RobVar"]]