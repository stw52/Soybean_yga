#### APSIM Sensitivity Analysis
library(tidyverse)
library(apsimx)
library(sensitivity)

#### Example as followed by Fernando Miguez code

tmp.dir = tempdir()
extd.dir = system.file("extdata", package = "apsimx")

file.copy(file.path(extd.dir, "Wheat.apsimx"), tmp.dir)

pp <- inspect_apsimx("Wheat.apsimx", src.dir = tmp.dir, 
                     node = "Manager", parm = list("SowingFertiliser", 1))

## For simplicity, we create a vector of fertilizer amounts (instead of sampling)
# Code here is for an example of the process, but there is function in package that simplifies this process
ferts <- seq(5, 200, length.out = 7)
col.res <- NULL
for(i in seq_along(ferts)){
  
  edit_apsimx("Wheat.apsimx", src.dir = tmp.dir,
              node = "Other",
              parm.path = pp, parm = "Amount",
              value = ferts[i])
  
  sim <- apsimx("Wheat.apsimx", src.dir = tmp.dir)
  col.res <- rbind(col.res, data.frame(fertilizer.amount = ferts[i], 
                                       wheat.aboveground.wt = mean(sim$Wheat.AboveGround.Wt, na.rm = TRUE)))
}
## Using sens_apsimx

## Identify a parameter of interest
## In this case we want to know the impact of varying the fertilizer amount
## and the plant population
pp1 <- inspect_apsimx("Wheat.apsimx", src.dir = tmp.dir, 
                      node = "Manager", parm = list("SowingFertiliser", 1))
pp1 <- paste0(pp1, ".Amount")

pp2 <- inspect_apsimx("Wheat.apsimx", src.dir = tmp.dir, 
                      node = "Manager", parm = list("SowingRule1", 9))
pp2 <- paste0(pp2, ".Population")

## The names in the grid should (partially) match the parameter path names
grd <- expand.grid(Fertiliser = c(50, 100, 150), Population = c(100, 200, 300))

## This takes 2-3 minutes
sns <- sens_apsimx("Wheat.apsimx", src.dir = tmp.dir,
                   parm.paths = c(pp1, pp2),
                   grid = grd)
summary(sns)

## Simple example: see documentation for other options
X.grid <- parameterSets(par.ranges = list(Fertiliser = c(1, 300), 
                                          Population = c(1, 300)),
                        samples = c(3,3), method = "grid")
X.grid


X.mrrs <- morris(factors = c("Fertiliser", "Population"),
                 r = 3, design = list(type = "oat", levels = 3, grid.jump = 1),
                 binf = c(0, 5), bsup = c(200, 300))
X.mrrs$X



sns2 <- sens_apsimx("Wheat.apsimx", src.dir = tmp.dir,
                    parm.paths = c(pp1, pp2),
                    grid = X.mrrs$X)
## These are the sensitivity results for AboveGround.Wt only
sns.res.ag <- tell(X.mrrs, sns2$grid.sims$Yield)
sns.res.ag

plot(sns.res.ag)


#### Trying this on soils for NY

src_dir <- getwd()

edit_apsimx("soygap_medsoils_only.apsimx",
            node = "Soil",
            soil.child = "Physical",
            root = "NY 2015",
            parm = "Soybean XF",
            value = c(1,1,1,1,1,0.2,0.2))



roots <- c("NY 2015", "NY 2016", "NY 2023", "NY 2024")

for(i in seq_along(roots)) {
  pp1 <- paste0(".Simulations",".",roots[i],".",)
  
  edit_apsimx("soygap_medsoils_only.apsimx",
              node = "Soil",
              root = roots[i], 
              soil.child = "Physical",
              parm.path = 
              parm = ("Soybean XF"),
              value = c(1,1,1,1,1,0.2,0.2),
              edit.tag = paste0("_",roots[i]))


}
medsoil_runs <- apsimx("soygap_medsoils_only.apsimx", src.dir = src_dir)





