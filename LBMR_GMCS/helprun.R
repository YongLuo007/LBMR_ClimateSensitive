rm(list = ls())
library(SpaDES); library(raster)
library(magrittr)

inputDir <- file.path(dirname(tempdir()), "LBMR", "inputs") %>% checkPath(create = TRUE)
outputDir <- file.path(dirname(tempdir()), "LBMR", "outputs") 
times <- list(start = 0, end = 10)
parameters <- list(
  .globals = list(verbose = FALSE),
  LBMR = list()
)
# to call the original GM functions
modulesOrig <- list("LBMR", "LBMR_GMOrig")

# to call the climate-sensitivity GM functions
modulesCS <- list("LBMR", "LBMR_GMCS")

firemap <- raster::raster(file.path("~/GitHub/LBMR_ClimateSensitive/LBMR/data/ecoregions.gis"))
firemap[] <- sample(seq(5, 200, by = 5), size = ncell(firemap), replace = TRUE)
CMIAnomalyMap <- firemap
CMIAnomalyMap[] <- rnorm(n = ncell(firemap), mean = 1, sd = 1)

CMINormalMap <- firemap
CMINormalMap[] <- rnorm(n = ncell(firemap), mean  = 10, sd = 2)

objects <- list(successionTimestep = 2,
                rstTimeSinceFire = firemap,
                CMIAnomalyMap = CMIAnomalyMap,
                SpatialDependency = FALSE,
                CMINormalMap = CMINormalMap)
paths <- list(
  cachePath = file.path(outputDir, "cache"),
  modulePath = file.path("~/GitHub/LBMR_ClimateSensitive/"),
  inputPath = inputDir,
  outputPath = outputDir
)

mySim <- simInit(times = times, params = parameters, modules = modulesCS,
                 objects = objects, paths = paths)
mySim <- simInit(times = times, params = parameters, modules = modulesOrig,
                 objects = objects, paths = paths)
mySimOut <- spades(mySim, debug = TRUE)
