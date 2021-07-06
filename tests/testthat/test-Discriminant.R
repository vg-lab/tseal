test_that("test trainModel", {
  totalExperiments <- readRDS("../totalExperiments.rds")
  NmainExperiments <- min(10,271)
  NpiedExperiments <- min(10,280)
  MedicalClasification <- rbind(matrix(1,NmainExperiments,1),matrix(2,NpiedExperiments,1))

  electrodes <- 10

  NpiedExperiments <- NpiedExperiments + 271

  mainExperiments <- totalExperiments[1:electrodes,1:20,1:NmainExperiments]
  piedExperiments <- totalExperiments[1:electrodes,1:20,272:NpiedExperiments]

  MWA <- MultiVaweAnalisys(mainExperiments,piedExperiments,"haar")
  MWA2 <- StepDiscrim(MWA,MedicalClasification,20)

  result <- KFCV(MWA2,MedicalClasification,"linear",3)
})
