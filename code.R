rm(list = ls())
graphics.off()

# Installing INLA

#install.packages("INLA", repos = "https://inla.r-inla-download.org/R/stable", dep = TRUE)

# Loading R packages

library(INLA)
library(fmesher)
library(leaflet)
library(viridis)
library(ggplot2)
library(cowplot)
library(sf)
library(ggspatial)


# Reading data

d <- read.csv("d.csv")
dp <- read.csv("dp.csv")

# Fitting model

coo <- cbind(d$longitude, d$latitude)

#' Make a spatial data frame version of the prevalence data, 
#' solely for making a map
d.sf <- st_transform(st_as_sf(d, coords = c("longitude", "latitude"), crs = 4326), crs = 3857)
head(d.sf$geometry)

mesh <- inla.mesh.2d(loc = coo, max.edge = c(0.1, 5), cutoff = 0.01)

plot(mesh)


spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)
indexs <- inla.spde.make.index("s", spde$n.spde)
A <- inla.spde.make.A(mesh = mesh, loc = coo)

coop <- cbind(dp$longitude,dp$latitude)
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)

stk.e <- inla.stack(tag = "est", data = list(y = d$positive,
                    numtrials = d$examined), A = list(1, A),
                    effects = list(data.frame(b0 = 1, alt = d$alt,
                    temp = d$temp, prec = d$prec, hum = d$hum,
                    pop = d$pop, aqua=d$dist_aqua), s = indexs))

stk.p <- inla.stack(tag = "pred", data = list(y = NA, numtrials = NA),
                    A = list(1, Ap), effects = list(data.frame(b0 = 1,
                    alt = dp$altitude, temp = dp$temp,
                    prec = dp$prec, hum = dp$hum, pop = dp$pop,
                    aqua=dp$dist_aqua), s = indexs))

stk.full <- inla.stack(stk.e, stk.p)

formula <- y ~ 0 + b0 + alt + temp + prec + hum + pop + aqua + f(s, model = spde)

res <- inla(formula, family = "binomial", Ntrials = numtrials,
            control.family = list(link = "logit"),
            control.compute = list(config = TRUE, return.marginals.predictor = TRUE, waic = TRUE),
            data = inla.stack.data(stk.full),
            control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.full)))


# Results

res$waic$waic

summary(res)

#' How much variaton are the fixed effects explaining?
X <- model.matrix(~ alt + temp + prec + hum + pop + dist_aqua, data = d)
B <- res$summary.fixed$mean
var(X %*% B)


#' Calculate the exceedance probabilities from Table 1
covariates <- c("alt", "temp", "prec", "hum", "pop", "aqua")
cbind(sapply(covariates, function(v) 1 - inla.pmarginal(q = 0, marginal = res$marginals.fixed[[v]])))


index <- inla.stack.index(stack = stk.full, tag = "pred")$data

prev_mean <- res$summary.fitted.values[index, "mean"]
prev_ll <- res$summary.fitted.values[index, "0.025quant"]
prev_ul <- res$summary.fitted.values[index, "0.975quant"]

pal <- colorNumeric("viridis", c(0, 1), na.color = "transparent")

# Mean prevalence predicted values
leaflet() %>% addTiles() %>%
  addCircles(lng = coop[, 1], lat = coop[, 2], color = pal(prev_mean)) %>%
  addLegend("bottomright", pal = pal, values = prev_mean, title = "Mean Prev") %>%
  addScaleBar(position = c("bottomleft"))

# 0.025quant prevalence for predicted values
leaflet() %>% addTiles() %>%
  addCircles(lng = coop[, 1], lat = coop[, 2], color = pal(prev_ll)) %>%
  addLegend("bottomright", pal = pal, values = prev_ll, title = "0.025Quant Prev") %>%
  addScaleBar(position = c("bottomleft"))

# 0.975quant prevalence for predicted values
leaflet() %>% addTiles() %>%
  addCircles(lng = coop[, 1], lat = coop[, 2], color = pal(prev_ul)) %>%
  addLegend("bottomright", pal = pal, values = prev_ul, title = "0.975Quant Prev") %>%
  addScaleBar(position = c("bottomleft"))


# Exceedance probabilities

excprob <- sapply(res$marginals.fitted.values[index],
                  FUN = function(marg){1-inla.pmarginal(q = 0.2, marginal = marg)})

pal <- colorNumeric("viridis", c(0, 1), na.color = "transparent")

leaflet() %>% addTiles() %>%
  addCircles(lng = coop[, 1], lat = coop[, 2], color = pal(excprob)) %>%
  addLegend("bottomright", pal = pal, values = excprob, title = "P(prev>0.2)" ) %>%
  addScaleBar(position = c("bottomleft"))

excprob <- sapply(res$marginals.fitted.values[index],
                  FUN = function(marg){1-inla.pmarginal(q = 0.4, marginal = marg)})

leaflet() %>% addTiles() %>%
  addCircles(lng = coop[, 1], lat = coop[, 2], color = pal(excprob)) %>%
  addLegend("bottomright", pal = pal, values = excprob, title = "P(prev>0.4)" ) %>%
  addScaleBar(position = c("bottomleft"))

excprob <- sapply(res$marginals.fitted.values[index],
                  FUN = function(marg){1-inla.pmarginal(q = 0.6, marginal = marg)})

leaflet() %>% addTiles() %>%
  addCircles(lng = coop[, 1], lat = coop[, 2], color = pal(excprob)) %>%
  addLegend("bottomright", pal = pal, values = excprob, title = "P(prev>0.6)" ) %>%
  addScaleBar(position = c("bottomleft"))

excprob <- sapply(res$marginals.fitted.values[index],
                  FUN = function(marg){1-inla.pmarginal(q = 0.8, marginal = marg)})

leaflet() %>% addTiles() %>%
  addCircles(lng = coop[, 1], lat = coop[, 2], color = pal(excprob)) %>%
  addLegend("bottomright", pal = pal, values = excprob, title = "P(prev>0.8)" ) %>%
  addScaleBar(position = c("bottomleft"))


# Plots for spatial field

rang <- apply(mesh$loc[, c(1, 2)], 2, range)
proj <- inla.mesh.projector(mesh, xlim = rang[, 1], ylim = rang[, 2], dims = c(300, 300))
mean_s <- inla.mesh.project(proj, res$summary.random$s$mean)
sd_s <- inla.mesh.project(proj, res$summary.random$s$sd)
df <- expand.grid(x = proj$x, y = proj$y)
df$mean_s <- as.vector(mean_s)
df$sd_s <- as.vector(sd_s)

gmean <- ggplot(df, aes(x = x, y = y, fill = mean_s)) + geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw()

gsd <- ggplot(df, aes(x = x, y = y, fill = sd_s)) + geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw()

plot_grid(gmean, gsd)


#' Quick way to get at the priors 
res.na <- res
res.na$.args$data$y <- rep(NA, length(res.na$.args$data$y))
res.na <- inla.rerun(res.na, plain = TRUE)
summary(res.na)
exp(res.na$summary.hyperpar[,-2])

#' Plot the prior and posterior densities of the hyperparameters
par(mfrow = c(2, 1))
plot(inla.spde.result(res, "s", spde)$marginals.range.nominal$range.nominal.1, type = "l",
     xlim = c(0, 50), , main = "Spatial field range")
lines(inla.spde.result(res.na, "s", spde)$marginals.range.nominal$range.nominal.1, 
      col = 2)
abline(h = 0, col = grey(0.5))
legend("topright", legend = c("Prior", "Posterior"), lty = 1, col = 2:1)

plot(inla.spde.result(res, "s", spde)$marginals.variance.nominal$variance.nominal.1, type = "l",
     xlim = c(0, 10), main = "Spatial field variance")
lines(inla.spde.result(res.na, "s", spde)$marginals.variance.nominal$variance.nominal.1, 
      col = 2)
abline(h = 0, col = grey(0.5))
legend("topright", legend = c("Prior", "Posterior"), lty = 1, col = 2:1)
