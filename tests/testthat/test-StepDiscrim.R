test_that("test StepDiscrim", {
  totalExperiments <- readRDS("../totalExperiments.rds")
  NmainExperiments <- min(2,271)
  NpiedExperiments <- min(2,280)
  MedicalClasification <- rbind(matrix(1,NmainExperiments,1),matrix(2,NpiedExperiments,1))

  electrodes <- 10

  NpiedExperiments <- NpiedExperiments + 271

  mainExperiments <- totalExperiments[1:electrodes,1:20,1:NmainExperiments]
  piedExperiments <- totalExperiments[1:electrodes,1:20,272:NpiedExperiments]

  MWA <- MultiVaweAnalisys(mainExperiments,piedExperiments,"haar")
  MWA2 <- StepDiscrim(MWA,MedicalClasification,20)

  m <- read.csv("../Results/StepVar&Cor.csv", header = FALSE)
  expect_equal(MWA2[[1]],as.matrix(m),ignore_attr = TRUE)
})

test_that("test StepDiscrimVar", {
  totalExperiments <- readRDS("../totalExperiments.rds")
  NmainExperiments <- min(2,271)
  NpiedExperiments <- min(2,280)
  MedicalClasification <- rbind(matrix(1,NmainExperiments,1),matrix(2,NpiedExperiments,1))

  electrodes <- 10

  NpiedExperiments <- NpiedExperiments + 271

  mainExperiments <- totalExperiments[1:electrodes,1:20,1:NmainExperiments]
  piedExperiments <- totalExperiments[1:electrodes,1:20,272:NpiedExperiments]

  MWA <- MultiVaweAnalisys(mainExperiments,piedExperiments,"haar")
  MWA2 <- StepDiscrim(MWA,MedicalClasification,20,Cor = FALSE)

  m <- read.csv("../Results/StepVar.csv", header = FALSE)
  expect_equal(MWA2[[1]],as.matrix(m),ignore_attr = TRUE)
})

test_that("test StepDiscrimCor", {
  totalExperiments <- readRDS("../totalExperiments.rds")
  NmainExperiments <- min(2,271)
  NpiedExperiments <- min(2,280)
  MedicalClasification <- rbind(matrix(1,NmainExperiments,1),matrix(2,NpiedExperiments,1))

  electrodes <- 10

  NpiedExperiments <- NpiedExperiments + 271

  mainExperiments <- totalExperiments[1:electrodes,1:20,1:NmainExperiments]
  piedExperiments <- totalExperiments[1:electrodes,1:20,272:NpiedExperiments]

  MWA <- MultiVaweAnalisys(mainExperiments,piedExperiments,"haar")
  MWA2 <- StepDiscrim(MWA,MedicalClasification,20,Var = FALSE)

  m <- read.csv("../Results/StepCor.csv", header = FALSE)
  expect_equal(MWA2[[1]],as.matrix(m),ignore_attr = TRUE)
})
