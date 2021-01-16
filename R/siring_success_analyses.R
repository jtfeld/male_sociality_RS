library(dplyr)
library(MuMIn)
library(lmerTest)


# ---------------- Analysis 1 --------------------------------------
full_data = readRDS(file = "./data/full_analysis_by_sire.rds")

m0 = glmer(sire ~ 
             age +
             age2 +
             id1_elo +
             qgr +
             scale(nmales) +
             (1|id), family = binomial, 
           data = full_data, na.action = na.fail) 

summary(m0)

all_mods = dredge(global.model = m0, trace = 2) # , subset = !(time_obs_m & time_obs_atleast1_male)

all_mods = subset(all_mods, !nested(.))

best_mods = get.models(all_mods, subset = delta < 6)

model.sel(best_mods) 

summary(best_mods$`14`)

# specifying best model as null model:
m0 = glmer(sire ~ 
             age +
             id1_elo +
             qgr +
             (1|id), family = binomial, 
           data = full_data, na.action = na.fail)

# ------------------ sociality models: ---------------------------------------

# Grooming models:
m_grm_tottime = update(m0, ~ . + sum_grm_tot) # total grooming time
m_grmtot = update(m0, ~ . + grm_tot_rate) # total grooming time divided by time with males
m_grmmean = update(m0, ~ . + mean_grm_rate) # mean grooming rate across all partners
m_grm_nhightot = update(m0, ~ . + n_high_grm_tot) # number above average grooming rate partners
m_grm_sumhightot = update(m0, ~ . + sum_high_grm_tot) # sum of above average grooming rates
m_grm_nhightime = update(m0, ~ . + n_high_grm_time) # number above average grooming time partners
m_grm_sumhightime = update(m0, ~ . + sum_high_grm_time) # sum of above average grooming times
m_grm_totparners = update(m0, ~ . + tot_grm_partners) # total grooming partners

# SRI (here "DAI") models
m_dai_deg = update(m0, ~ . + dai_deg) # number of dai4 partners
m_dai_str = update(m0, ~ . + dai_str) # sum of dai4 indices
m_n_high_dai = update(m0, ~ . + n_high_dai_arr) # number of above average dai4 partners
m_sum_high_dai = update(m0, ~ . + sum_high_dai_arr) # sum of above average dai4 indices

# CSI models
m_csi = update(m0, ~ . + csi2_top3) # sum of top 3 csis
m_nighcsi = update(m0, ~ . + n_high_csi2) # number of above mean csis
m_sumhighcsi = update(m0, ~ . + sum_high_csi2) # sum of csi values above mean

# AIC table:
AICc(m0, 
     m_grm_tottime,
     m_grmtot, # tot grming rate divided by time with males
     m_grmmean,  # mean grooming rate
     m_grm_nhightot, # number above average grooming rate partners
     m_grm_sumhightot, # sum of above average grooming rates
     m_grm_nhightime, # number above average grooming time partners
     m_grm_sumhightime, # sum of above average grooming times
     m_grm_totparners, # total grooming partners
     # dai-based:
     m_dai_deg, # number of dai4 partners
     m_dai_str, # sum of dai4 indices
     m_n_high_dai, # number of above average dai4 partners
     m_sum_high_dai, # sum of above average dai4 indices
     # csi-based:
     m_csi, # sum of top 3 csis
     m_sumhighcsi, # sum of csi values above mean
     m_nighcsi # number of above mean csis
     
) %>% 
  arrange(AICc) %>% 
  mutate(mod = row.names(.),
         delta = AICc(m0) - AICc, 
         weights = MuMIn::Weights(AICc)) %>%
  select(mod, everything()) %>%
  mutate(AICc = round(AICc, 3),
         delta = round(delta, 3),
         weights = round(weights, 3))

# --------------------------- gregariousness check ---------------

m0 = glmer(sire ~ 
             age +
             id1_elo +
             qgr +
             time_obs_mf_z + # this is total time observed w both males and females
             (1|id), family = binomial, 
           data = full_data, na.action = na.fail)
# re-run the above models and model comparison procedure using the new m0

# ---------------- analysis 2 ----------------------

no_alphas = readRDS(file = "./data/no_alpha_analysis_by_sire.rds")

m_0 = glmer(sire ~ 
              age +
              id1_elo +
              qgr +
              (1|id), family = binomial, #  + grm_given_rate + grm_rec_rate
            data = no_alphas, na.action = na.fail) #, REML = F   nAGQ = 0

summary(m_0)

m_dai_w_alpha = update(m_0, ~ . + dai_w_alpha) # association rate with the alpha
m_grmtot_w_alpha = update(m_0, ~ . + grmtot_w_alpha) # grooming rate with the alpha
m_csi_w_alpha = update(m_0, ~ . + csi_w_alpha, nAGQ = 0) # CSI score with the alpha
m_csialpha_nhighdai = update(m_csi_w_alpha, ~ . + n_high_dai_arr) # csi w/ alpha AND count strong association ties
m_daialpha_nhighdai = update(m_dai_w_alpha, ~ . + n_high_dai_arr, nAGQ = 0) # association rate w/ alpha AND count strong association ties
m_grmalpha_nhighdai = update(m_grmtot_w_alpha, ~ . + n_high_dai_arr, nAGQ = 0) # grooming rate w/ alpha AND count strong association ties
m_n_high_dai_noalpha = update(m_0, ~ . + n_high_dai_arr, nAGQ = 0) # count strong association ties


AICc(m_0, 
     m_csi_w_alpha,
     m_dai_w_alpha,
     m_grmtot_w_alpha,
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
