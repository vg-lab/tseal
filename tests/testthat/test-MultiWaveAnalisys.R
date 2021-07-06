test_that("test MultiWaveAnalisysVarCOR", {
  NVar <- 30
  NCor <- 135

  totalExperiments <- readRDS("../totalExperiments.rds")
  NmainExperiments <- min(2,271)
  NpiedExperiments <- min(2,280)
  MedicalClasification <- rbind(matrix(1,NmainExperiments,1),matrix(2,NpiedExperiments,1))

  electrodes <- 10

  NpiedExperiments <- NpiedExperiments + 271

  mainExperiments <- totalExperiments[1:electrodes,1:20,1:NmainExperiments]
  piedExperiments <- totalExperiments[1:electrodes,1:20,272:NpiedExperiments]

  MWA <- MultiVaweAnalisys(mainExperiments,piedExperiments,"haar")
  m <- read.csv("../Results/MWA.csv", header = FALSE)
  expect_equal(MWA$Values,as.matrix(m),tolerance = 0.01, ignore_attr = TRUE)
  expect_equal(MWA$Vars,NVar)
  expect_equal(MWA$Cors,NCor)
})

test_that("test MultiWaveAnalisysCor", {
  NVar <- 0
  NCor <- 135
  m <- read.csv("../Results/MWA.csv", header = FALSE)
  m <- m[31 : (30+NCor),] # remove the correlations

  totalExperiments <- readRDS("../totalExperiments.rds")
  NmainExperiments <- min(2,271)
  NpiedExperiments <- min(2,280)
  MedicalClasification <- rbind(matrix(1,NmainExperiments,1),matrix(2,NpiedExperiments,1))

  electrodes <- 10

  NpiedExperiments <- NpiedExperiments + 271

  mainExperiments <- totalExperiments[1:electrodes,1:20,1:NmainExperiments]
  piedExperiments <- totalExperiments[1:electrodes,1:20,272:NpiedExperiments]

  MWA <- MultiVaweAnalisys(mainExperiments,piedExperiments,"haar",Var = FALSE)

  expect_equal(MWA$Values,as.matrix(m),tolerance = 0.01,ignore_attr = TRUE)
  expect_equal(MWA$Vars,NVar)
  expect_equal(MWA$Cors,NCor)

})

test_that("test MultiWaveAnalisysVar", {
  NVar <- 30
  NCor <- 0
  m <- read.csv("../Results/MWA.csv", header = FALSE)
  m <- m[1:NVar,] # remove the correlations

  totalExperiments <- readRDS("../totalExperiments.rds")
  NmainExperiments <- min(2,271)
  NpiedExperiments <- min(2,280)
  MedicalClasification <- rbind(matrix(1,NmainExperiments,1),matrix(2,NpiedExperiments,1))

  electrodes <- 10

  NpiedExperiments <- NpiedExperiments + 271

  mainExperiments <- totalExperiments[1:electrodes,1:20,1:NmainExperiments]
  piedExperiments <- totalExperiments[1:electrodes,1:20,272:NpiedExperiments]

  MWA <- MultiVaweAnalisys(mainExperiments,piedExperiments,"haar",Cor = FALSE)

  expect_equal(MWA$Values,as.matrix(m),tolerance = 0.01,ignore_attr = TRUE)
  expect_equal(MWA$Vars,NVar)
  expect_equal(MWA$Cors,NCor)

})





