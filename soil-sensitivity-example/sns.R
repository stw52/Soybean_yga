## date: 2025-10-10
## apsim version: 2025.3.7681.0

library(apsimx)
library(ggplot2)

## picking a random location in Iowa
random.location <- c(-93.768, 42.019)

## getting an example of a soil profile
random.soil <- get_ssurgo_soil_profile(lonlat = random.location)

## taking a look at the example soil profile
plot(random.soil[[1]])

## running the simulation to make sure things work before
## we tweak them
sim0 <- apsimx('maize_sensitivity.apsimx',
               cleanup = TRUE)

## it seems like things worked :)
head(sim0)

## creating different soil profiles
base.soil <- random.soil[[1]]

## creating a function to make this less painful
## the trick when building a function like this is ensuring functional
## soil profiles. For example, by definition DUL < SAT, so when varying one, 
## it's important to check that it's not crossing any boundaries.
modify_soil_property <- function(soil.profile,
                                 soil.property = c('LL15', 'SAT', 'DUL', 'KS'), ## this can be expanded
                                 multiplier){
  soil.property <- match.arg(soil.property)
  new.soil <- soil.profile
  new.soil$soil[[soil.property]] <- new.soil$soil[[soil.property]] * multiplier 
  return(new.soil)
}

## creating a grid in which I vary DUL. This can be adapted to vary 
## other soil properties too.
sens.grid <- data.frame(DUL.Multiplier = seq(0.8, 1.1, 0.05))
sens.grid$soil.profile <- 1:nrow(sens.grid)

## creating different soil profiles for following the grid
soils.list <- lapply(sens.grid[['DUL.Multiplier']],
                     FUN = function(x){
                       modify_soil_property(base.soil,
                                            soil.property = 'DUL',
                                            multiplier = x)
                     })

## running the simulation
sens0 <- sens_apsimx(file = 'maize_sensitivity.apsimx',
                   grid = sens.grid[, 2, drop = FALSE],
                   parm.paths = 'soil.profile',
                   summary = 'none',
                   soil.profiles = soils.list)

## taking a look at how sensitive the model is to varying DUL
sens1 <- sens0$grid.sims
sens1 <- merge(sens1, sens.grid)

ggplot(sens1)+
  geom_line(aes(x = Date,
                y = Maize.AboveGround.Wt,
                col  = as.character(DUL.Multiplier)))
