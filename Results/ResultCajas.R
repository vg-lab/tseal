totalExperiments <- readRDS("tests/totalExperimentsECG.rds")
infractionExperiments <- min(200,368)
healthExperiments <- min(80,80)
MedicalClasification <- rbind(matrix(1,infractionExperiments,1),matrix(2,healthExperiments,1))


healthExperiments <- healthExperiments + 368

mainExperiments <- totalExperiments[,,1:infractionExperiments]
piedExperiments <- totalExperiments[,,369:healthExperiments]

features = c("Cor","IQR","DM")
maxVars = 20
filter = "la8"

result <- LOOCV(mainExperiments,piedExperiments,method = "linear", f = filter, maxvars = maxVars,features = features,returnClassification = TRUE)
write.csv(result[[2]] == MedicalClasification, "Cor_IQR_DM_20_la8.csv")


features = c("Cor","IQR","PE","DM")
maxVars = 40
filter = "d6"

result <- LOOCV(mainExperiments,piedExperiments,method = "linear", f = filter, maxvars = maxVars,features = features,returnClassification = TRUE)
write.csv(result[[2]] == MedicalClasification, "Cor_IQR_PE_DM_40_d6.csv")

features = c("Var","Cor","IQR","DM")
maxVars = 60
filter = "d6"
result <- LOOCV(mainExperiments,piedExperiments,method = "linear", f = filter, maxvars = maxVars,features = features,returnClassification = TRUE)
write.csv(result[[2]] == MedicalClasification, "Var_Cor_IQR_DM_60_d6.csv")


