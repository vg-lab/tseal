totalExperiments <- readRDS("tests/totalExperiments.rds")
NmainExperiments <- min(200,271)
NpiedExperiments <- min(200,280)
MedicalClasification <- rbind(matrix(1,NmainExperiments,1),matrix(2,NpiedExperiments,1))


NpiedExperiments <- NpiedExperiments + 271
#--------------

Cors <- c(6447,8571,11488,13137,15213,13395,11000,3459,4386,6559,8490,10768,15705,4305,15175,1407,17571,4971,11959,14590)
CorsK <- sapply(Cors,function(x) floor((x-1)/9))
NbK <- combn(1:63,2)
electrodes <- as.vector(sapply(CorsK, function(x) NbK[,x]))
#---------------

electrodesCortex <-c(5,6,18,20,21,22,23,26,27,34,35,40,41,42,43,48,61)


mainExperiments <- totalExperiments[c(5,6,18),,1:NmainExperiments]
piedExperiments <- totalExperiments[c(5,6,18),,272:NpiedExperiments]

MWA <- MultiVaweAnalisys(mainExperiments,piedExperiments,"d6",features = "Cor")
MWA2 <- StepDiscrim(MWA,MedicalClasification,20, features = c("Cor"))

result <- LOOCV(MWA2,MedicalClasification,"linear")
