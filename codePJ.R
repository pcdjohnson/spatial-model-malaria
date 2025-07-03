#' The repo containing this code was forked from 
#' https://github.com/ncespedesc/spatial-model-malaria
#' which I assume was itself forked from 
#' https://github.com/Paula-Moraga/spatial-model-malaria,
#' which no longer exists. I've made a several changes to the methods
#' which are easy to spot because they're commented with #'
#' instead of #.


#' Clear out objects
rm(list = ls())

# Installing INLA
#install.packages("INLA", repos = "https://inla.r-inla-download.org/R/stable", dep = TRUE)

# Loading R packages
library(INLA)
library(leaflet)
library(viridis)
library(ggplot2)
library(cowplot)
library(sf)
library(ggspatial)

# Reading data
d <- read.csv("d.csv") # prevalence data
dp <- read.csv("dp.csv") # locations we want to predict at

dim(d)
dim(dp)
head(d)
head(dp)

#' Make an observation-level factor for fitting a random effect to
d$obs <- factor(1:nrow(d))
dp$obs <- factor(1:nrow(dp))

#' Altidude has different names between the two data sets, fix this
d$altitude <- d$alt
d$alt <- NULL


#' Altitude in dp has three values <= 0, which can't be log-transformed.
#' Delete them and replace them as half the minimum value of the other 
#' observations
dp$altitude[dp$altitude <= 0] <- NA
dp$altitude[is.na(dp$altitude)] <- min(dp$altitude, na.rm = TRUE)/2

#' Same for dist_aqua in d (1 value) and dp (323 values)
d$dist_aqua[d$dist_aqua <= 0] <- NA
d$dist_aqua[is.na(d$dist_aqua)] <- min(d$dist_aqua, na.rm = TRUE)/2
dp$dist_aqua[dp$dist_aqua <= 0] <- NA
dp$dist_aqua[is.na(dp$dist_aqua)] <- min(dp$dist_aqua, na.rm = TRUE)/2

# coordinates of the surveys
coo <- cbind(d$longitude, d$latitude)

#' Make a spatial data frame version of the prevalence data, 
#' solely for making a map
d.sf <- st_transform(st_as_sf(d, coords = c("longitude", "latitude"), crs = 4326), crs = 3857)
head(d.sf$geometry)

#' Make a map of the survey locations. 
ggplot() +
  annotation_map_tile(type = "cartolight", # available types: rosm::osm.types()
                      zoomin = 0) + 
  geom_sf(data = d.sf, aes(colour = prev, size = sqrt(examined)), inherit.aes = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.25, height = unit(0.1, "cm")) +
  scale_size_continuous(name="sqrt(examined)", range = c(0.2,5)) +
  theme_minimal()

#' Plot scatterplots and histograms of the continuous variables
covariates <- c("altitude", "temp", "prec", "hum", "pop", "dist_aqua")
panel.hist <- function(x, ...) {
  usr <- par("usr")
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, ...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 2)
}
pairs(d[, covariates], upper.panel = panel.cor, diag.panel = panel.hist,
      gap=0, row1attop=FALSE, cex.labels = 1.3)

#' altitude, pop, dist_aqua are all skewed, log them to reduce the influence 
#' of extreme values
for(v in c("altitude", "pop", "dist_aqua")) {
  print(range(d[, v], dp[, v]))
  v.lab <- paste0(v, ".log10")
  d[, v.lab] <- log10(d[, v])
  dp[, v.lab] <- log10(dp[, v])
  covariates[covariates == v] <- v.lab
}; rm(v, v.lab)

#' Much less skewed
pairs(d[, covariates], upper.panel = panel.cor, diag.panel = panel.hist,
      gap=0, row1attop=FALSE, cex.labels = 1.3)


#' Convert all the continuous covariates to SD scores with zero mean. 
#' The makes it easier (kind of) to interpret the effect estimates (betas).
#' It also makes it OK to have a single prior on all the betas. The default
#' fixed effect prior in INLA is N(mean = 0, variance = 1000), and the
#' informativeness of the prior will depend on the scale of the variable.
for(v in covariates) {
  v.lab <- paste0(v, ".sds")
  d[, v.lab] <- scale(d[, v])
  dp[, v.lab] <- scale(dp[, v])
  covariates[covariates == v] <- v.lab
}; rm(v, v.lab)



# Fitting model
mesh <- inla.mesh.2d(loc = coo, max.edge = c(0.1, 5), cutoff = 0.01)
plot(mesh)

# this is quite a dense mesh, but better to err on this side
par(mfrow = c(2, 1))
hist(c(dist(d[, c("longitude", "latitude")])), breaks = 1000)
abline(v = 0.01, col = "red")
hist(c(dist(d[, c("longitude", "latitude")])), breaks = 10000, xlim = c(0, 0.1))
abline(v = 0.01, col = "red")
par(mfrow = c(1, 1))

plot(mesh)
spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)
#' ...to do: check that the default priors make sense. Consider changing 
#' to inla.spde2.pcmatern which is probably easier to set priors for,
#' and gives the range and SD in the summary output.
indexs <- inla.spde.make.index("s", spde$n.spde)
A <- inla.spde.make.A(mesh = mesh, loc = coo)

coop <- cbind(dp$longitude,dp$latitude)
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)

stk.e <- inla.stack(tag = "est", data = list(y = d$positive,
                                             numtrials = d$examined), A = list(1, A, 1),
                    effects = list(data.frame(b0 = 1, altitude.log10.sds = d$altitude.log10.sds,
                                              temp.sds = d$temp.sds, prec.sds = d$prec.sds, hum.sds = d$hum.sds,
                                              pop.log10.sds = d$pop.log10.sds, dist_aqua.log10.sds=d$dist_aqua.log10.sds), 
                                   s = indexs,
                                   obs = d$obs))

stk.p <- inla.stack(tag = "pred", data = list(y = NA, numtrials = NA),
                    A = list(1, Ap, 1), effects = list(data.frame(b0 = 1,
                                                               altitude.log10.sds = dp$altitude.log10.sds, temp.sds = dp$temp.sds,
                                                               prec.sds = dp$prec.sds, hum.sds = dp$hum.sds, pop.log10.sds = dp$pop.log10.sds,
                                                               dist_aqua.log10.sds=dp$dist_aqua.log10.sds), 
                                                       s = indexs,
                                                       obs = dp$obs))

stk.full <- inla.stack(stk.e, stk.p)

#' Prior for the observation-level random effect, which models non-spatial extra-binomial
#' variation
iid.prec.prior <- list(prec = list(prior = "loggamma", param = c(0.5, 0.01)))

formula <- 
  y ~ 0 + b0 + altitude.log10.sds + temp.sds + prec.sds + hum.sds + pop.log10.sds + dist_aqua.log10.sds + 
  f(s, model = spde) + f(obs, model = 'iid', hyper = iid.prec.prior)

res <- inla(formula, family = "binomial", Ntrials = numtrials,
            control.family = list(link = "logit"),
            control.compute = list(config = TRUE, return.marginals.predictor = TRUE, waic = TRUE),            
            data = inla.stack.data(stk.full),
            control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.full)))
#' ...takes a few min

# Results
res$waic$waic
summary(res)
res.summary.fixed <- res$summary.fixed[-1, ]
res.summary.fixed$ID <- factor(rownames(res.summary.fixed), rownames(res.summary.fixed))
ggplot(res.summary.fixed, aes(y = ID, x = exp(mean))) +
  geom_pointrange(aes(xmin = exp(`0.025quant`), xmax = exp(`0.975quant`))) +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_y_discrete(limits = rev(levels(res.summary.fixed$ID))) +
  scale_x_continuous(breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20), 
                     transform = "log") +
  theme_minimal() +
  labs(y = "", x = "Odds ratio estimate",
       title = "Fixed effect estimates \u00B1 95% CI")

#' Calculate the exceedance probabilities from Table 1
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
