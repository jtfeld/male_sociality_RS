
# still gotta do: bonds w beta male, 
# also change no_alpha_analysis_by_sire2.rds so there's only
# one no_alpha_analysis_by_sire file (2 is the good one)

# ------------- using raw joint arrivals -------------------------

library(dplyr) # for organizing model comparison table
library(lmerTest) # for modeling 
library(MuMIn) # for AICc function
library(sjPlot) # for quick plots of model predictions

full_data = readRDS(file = "./data/full_analysis_by_sire2.rds") # CHANGE

# null model:
m0 = glmer(sire ~ 
             age +
             id1_elo +
             qgr +
             (1|id), family = binomial, 
           data = full_data, na.action = na.fail)

m_n_high_joint_arrs = update(m0, ~ . + n_high_joint_arrs)
m_sum_high_joint_arrs = update(m0, ~ . + sum_high_joint_arrs)

# model comparison table:
AICc(m0, 
     m_n_high_joint_arrs,
     m_sum_high_joint_arrs
) %>% 
  arrange(AICc) %>% 
  mutate(mod = row.names(.),
         delta = AICc(m0) - AICc, 
         weights = MuMIn::Weights(AICc)) %>%
  select(mod, everything()) %>%
  mutate(AICc = round(AICc, 3),
         delta = round(delta, 3),
         weights = round(weights, 3)) 

# results from the best model:
summary(m_n_high_joint_arrs)

# quick plots of model predictions
plot_model(m_n_high_joint_arrs, type = "pred")

# ------------- no alphas in counts of strong assn ties -----------

library(dplyr) # for organizing model comparison table
library(lmerTest) # for modeling 
library(MuMIn) # for AICc function
library(sjPlot) # for quick plots of model predictions

# load no_alphas dataset:
no_alphas = readRDS(file = "./data/no_alpha_analysis_by_sire2.rds") # CHANGE

# null model:
m_0 = glmer(sire ~ 
              age +
              id1_elo +
              qgr +
              (1|id), family = binomial, 
            data = no_alphas, na.action = na.fail) 

#inspect null model:
summary(m_0)

# model comparisons:
m_dai_w_alpha = update(m_0, ~ . + dai_w_alpha)
m_grmtot_w_alpha = update(m_0, ~ . + grmtot_w_alpha)
m_csi_w_alpha = update(m_0, ~ . + csi_w_alpha, nAGQ = 0)
m_csialpha_nhighdai = update(m_csi_w_alpha, ~ . + n_high_dai_arr_noalpha)
m_daialpha_nhighdai = update(m_dai_w_alpha, ~ . + n_high_dai_arr_noalpha, nAGQ = 0)
m_grmalpha_nhighdai = update(m_grmtot_w_alpha, ~ . + n_high_dai_arr_noalpha, nAGQ = 0)
m_n_high_dai_noalpha = update(m_0, ~ . + n_high_dai_arr_noalpha, nAGQ = 0)

# model comparison table:
AICc(m_0, 
     m_csi_w_alpha,
     m_grmtot_w_alpha,
     m_dai_w_alpha,
     m_csialpha_nhighdai,
     m_daialpha_nhighdai,
     m_grmalpha_nhighdai,
     m_n_high_dai_noalpha
) %>% 
  arrange(AICc) %>% 
  mutate(mod = row.names(.),
         delta = AICc(m_0) - AICc, 
         weights = MuMIn::Weights(AICc)) %>%
  select(mod, everything()) %>%
  mutate(AICc = round(AICc, 3),
         delta = round(delta, 3),
         weights = round(weights, 3)) 

# results of best model:
summary(m_csialpha_nhighdai)

# quick plots of model predictions
plot_model(m_csialpha_nhighdai, type = "pred")


# ------------ bond strength with beta males ------------

library(dplyr) # for organizing model comparison table
library(lmerTest) # for modeling 
library(MuMIn) # for AICc function
library(sjPlot) # for quick plots of model predictions

# load no betas dataset:
no_betas = readRDS(file = "./data/no_beta_analysis_by_sire.rds")

# null model:
m_0 = glmer(sire ~ 
              age +
              id1_elo +
              qgr +
              (1|id), family = binomial, 
            data = no_betas, na.action = na.fail) 

summary(m_0)

# update models with bonds w/ beta terms:
m_dai_w_beta = update(m_0, ~ . + dai_w_beta)
m_grmtot_w_beta = update(m_0, ~ . + grmtot_w_beta)
m_csi_w_beta = update(m_0, ~ . + csi_w_beta)

# model comparison tables:
AICc(m_0, 
     m_dai_w_beta,
     m_grmtot_w_beta,
     m_csi_w_beta
) %>% 
  arrange(AICc) %>% 
  mutate(mod = row.names(.),
         delta = AICc(m_0) - AICc, 
         weights = MuMIn::Weights(AICc)) %>%
  select(mod, everything()) %>%
  mutate(AICc = round(AICc, 3),
         delta = round(delta, 3),
         weights = round(weights, 3))

# model summaries:
summary(m_dai_w_beta)
summary(m_grmtot_w_beta)
summary(m_csi_w_beta)

# plots showing no effect of bonds with beta on siring success
plot_model(m_dai_w_beta, type = "pred", terms = "dai_w_beta [all]")
plot_model(m_grmtot_w_beta, type = "pred", terms = "grmtot_w_beta [all]")
plot_model(m_csi_w_beta, type = "pred", terms = "csi_w_beta [all]")

# ------------ Path analysis ----------------------------

# load no_alphas dataset

library(lavaan) # for path analysis
library(semPlot) # for plotting path diagram

no_alphas = readRDS(file = "./data/no_alpha_analysis_by_sire.rds")

head(no_alphas)

# specify model:

mod_path <-'
n_high_dai_arr ~ id1_elo + age
csi_w_alpha ~ id1_elo + age
csi_w_alpha ~ n_high_dai_arr
sire ~ age + id1_elo + qgr + csi_w_alpha + n_high_dai_arr
'

# run model: 
# (note that ordered = "sire" specifies that sire is a binary term)
fit_path <- sem(mod_path, data = no_alphas, ordered = "sire") 

summary(fit_path, fit.measures = TRUE, standardized = T, rsquare = T)

# plot: -----------------------

# initial path plot:
sp = semPaths(fit, 'std', layout = 'tree', 
              edge.label.cex = 1.5, sizeMan = 17, 
              sizeMan2 = 12, curvePivot = TRUE,
              label.cex = 0.75, label.prop = 1,
              nCharNodes = 0, #style = "lisrel",
              fade = T, residuals = F, intercepts = F)

# get default tree layout:
lo = sp$layout

# update layout for clearer plot:
lo[1,2] = 0
lo[2,2] = 0

# get order of node labels:
dput(sp$Arguments$labels)

# updated plot:
semPaths(fit, 'std', layout = lo, 
         edge.label.cex = 1.5, sizeMan = 17, 
         sizeMan2 = 12, curvePivot = TRUE,
         label.cex = 0.75, label.prop = 1,
         nCharNodes = 0, #style = "lisrel",
         fade = T, residuals = F, intercepts = F,
         nodeLabels = c("N strong\nassn. ties", "Bond w\nAlpha", 
                        "Sire\nsuccess", "Male Elo", "Male age", 
                        "Male-Fem.\nRelated."))





