
# ---------------- load packages ---------------------------

library(amen)
library(coda)
library(ggplot2)
library(bayesplot)

# ------------------ load data -------------------------------

load(file = "./data/coals_AME_Dec20.RData")
# change the above to point towards the data set

# ------------------ run models ----------------------------------

# this will take a while:

fit_coals1 = ame_rep(Y = Ycoal2, Xdyad = Xdyad2, Xrow = Xnode2,
                      family = "bin", R = 0, symmetric = T,
                      burn = 1000, nscan = 90000, odens = 100, plot = F, print = T)

summary(fit_coals1)
plot(fit_coals1)

fit_coals2 = ame_rep(Y = Ycoal2, Xdyad = Xdyad2, Xrow = Xnode2,
                      family = "bin", R = 0, symmetric = T, 
                      burn = 1000, nscan = 90000, odens = 100, plot = F, print = T, 
                      seed = 2)

summary(fit_coals2)
plot(fit_coals2)

fit_coals3 = ame_rep(Y = Ycoal2, Xdyad = Xdyad2, Xrow = Xnode2,
                       family = "bin", R = 0, symmetric = T, 
                       burn = 1000, nscan = 90000, odens = 100, plot = F, print = T, 
                       seed = 3)

summary(fit_coals3)
plot(fit_coals3)

save.image()

# --------------- check mcmc convergence -------------------

comb_chains = mcmc.list(as.mcmc(fit_coals1$BETA), 
                        as.mcmc(fit_coals2$BETA), 
                        as.mcmc(fit_coals3$BETA))
plot(comb_chains)

(exconv = gelman.diag(comb_chains))

# Values substantially above 1 indicate lack of convergence
# shoot for upper CI values <= 1.03

gelman.plot(comb_chains)

# quick look at results:
mcmc_areas(as.matrix(fit_coals1$BETA), prob = 0.95, pars = vars(-intercept))
mcmc_areas(as.matrix(fit_coals2$BETA), prob = 0.95, pars = vars(-intercept))
mcmc_areas(as.matrix(fit_coals3$BETA), prob = 0.95, pars = vars(-intercept))

# combined posterior distributions from three runs:
coal_post = rbind.data.frame(as.data.frame(fit_coals0b$BETA), 
                            as.data.frame(fit_coals0b2$BETA),
                            as.data.frame(fit_coals0b3$BETA))


# get posterior estimates and credible intervals:
round(Reduce("+", HPDinterval(comb_chains))/3, 2)
sapply(coal_post, function(x) round(mean(x), 2))

# ------------------ plot ---------------------------

theme_jtf = function (base_size = 11, base_family = "") {
  theme_linedraw() %+replace% 
    theme(
      panel.grid.major  = element_line(color = "gray94", size = 0.2),
      panel.grid.minor  = element_line(color = "white")
    )
}


theme_set(theme_jtf())

color_scheme_set("gray")

mcmc_areas(coal_post, 
           pars = c("age_diff.dyad", "time_tog_z.dyad", "grm_rate_tot_z.dyad", "dai_arr_z.dyad"), 
           prob = 0.95) + 
  geom_vline(xintercept = 0, color = "gray30", lty = 2) +
  scale_y_discrete(labels = c("Age difference", "Time together", 
                              "Grooming rate", "Association tie strength")) 
ggtitle("Posterior distributions",
        "with medians and 95% intervals")
