# TOPOSCORE paper - Analysis script
# v0.2 - JW - Jun 2023

# load required libraries and helper functions
source('tp_helper.R')

## 1. Discovery analysis set ----
log_msg('####### Discovery analysis set #########')
clin_disc <- load_clin(cohort = 'Disc')
met4_disc <- load_microbiome(clin_disc)

### 1.1 CoxPH screen (except Akkermansia) ----
res_surv <- load_or_compute('res_surv2.rds',
                          screen_surv_met4(clin_disc, met4_disc, type = 'OS'))
plot_surv_forest(res_surv, alpha = 0.05)
# select species based on average HR:
res_surv_filt <- res_surv %>% dplyr::filter(HR <= 0.8 | HR >= 1.25)
selected_species <- res_surv_filt$SPECIES
log_msg('%d/%d species selected', length(selected_species), nrow(res_surv))
plot_surv_forest(res_surv_filt %>% dplyr::arrange(HR))
met4_disc_filt <- met4_disc[, c('Sample_id', selected_species)]
hr_annots <- res_surv_filt %>% mutate(HRCAT = ifelse(HR < 1, 'R', 'NR')) %>% 
  dplyr::select(HRCAT, SPECIES)

### 1.2 Correlation screen ----
res <- load_or_compute('res_pairs2.rds', screen_pairs(met4_disc_filt))
# filter based on Fisher Bonferroni-corrected p-values <= 0.05
res_filt <- load_or_compute('res_pairs_filt2.rds', {
  min_p <- bind_rows(list(
    res %>% dplyr::select(VAR = VAR1, FISHER_P),
    res %>% dplyr::select(VAR = VAR2, FISHER_P)
  )) %>% group_by(VAR) %>% summarize(MIN_P = min(FISHER_P)) %>% arrange(-MIN_P)
  sp2 <- min_p[min_p$MIN_P <= 0.05 / nrow(min_p), 'VAR', drop = TRUE]
  res %>% filter(VAR1 %in% sp2 & VAR2 %in% sp2)
})
log_msg('Keeping %d pairs (%d species)', nrow(res_filt), 
        length(unique(c(res_filt$VAR1, res_filt$VAR2))))

### 1.3 Clustering ----
SCORE <- 'fisher_p'
METHOD <- 'ward.D2'
DISTANCE <- 'manhattan'
cc <- cluster_species(res_filt, score = SCORE, method = METHOD, distance = DISTANCE, k = 7) %>% renumber_clusters()
plt_fisher_disc <- plot_score_matrix(res_filt, score = SCORE, method = METHOD, 
                              distance = DISTANCE, annots = list(cc, hr_annots), fontsize = 3)
ggsave(plt_fisher_disc, filename = 'fig_fisher_disc.pdf', width = 30, height = 20, units = "cm")

### 1.4 Definition of SIGB groups ----
cc_names <- unique(cc$CLUSTER)
cc_species <- setNames(lapply(cc_names, function(x) cc[cc$CLUSTER == x, 'SPECIES', drop = TRUE]), cc_names)
SIGB1 <- load_or_compute('sigb12.rds', c(cc_species$C5, cc_species$C6))
SIGB2 <- load_or_compute('sigb22.rds', c(cc_species$C1, cc_species$C2, cc_species$C3))

### 1.5 Toposcoring ----
scores_disc <- compute_toposcore(met4_disc, sigb1 = SIGB1, sigb2 = SIGB2)
pred_disc <- clin_disc %>% left_join(scores_disc, by = 'Sample_id') %>% 
  filter(OS12 != '') %>% mutate(OS12 = factor(OS12, levels = c('NR', 'R')))
roc <- calc_roc(pred_disc$OS12, pred_disc$TOPOB01, verbose = TRUE)
log_msg('ROC AUC = %.2f [%.2f - %.2f]', roc$AUC[1], roc$AUC[2], roc$AUC[3])
youden <- roc$ROC_DF %>% mutate(SENS = TPR, SPEC = 1 - FPR) %>% mutate(J = SENS + SPEC - 1)
ggplot(youden, aes(x = THRESHOLD, y = J)) + geom_point()
ycut_nr <- youden[which.max(youden$J), ] # 0.5351351
ycut_r <- youden[which(youden$THRESHOLD > 0.7 & youden$J > 0.23), ] # 0.7911411
log_msg('Cut-off thresholds = %.4f and %.4f', ycut_nr$THRESHOLD, ycut_r$THRESHOLD)
plt_roc <- plot_roc(roc$ROC_DF) + 
  geom_point(data = ycut_nr, color = 'red') +
  geom_point(data = ycut_r, color = 'green')
ggsave(plt_roc, filename = 'fig_roc.pdf', width = 10, height = 10, units = "cm")
plt_kde_disc <- plot_toposcoreb01_density(scores_disc, clin_disc, 
                                lims = c(ycut_r$THRESHOLD, ycut_nr$THRESHOLD))
ggsave(plt_kde_disc, filename = 'fig_kde_disc.pdf', width = 10, height = 10, units = "cm")


### 1.6 Prediction in discovery cohort (full signature) ----
pred_disc <- assign_prediction(pred_disc, ycut_r$THRESHOLD, ycut_nr$THRESHOLD)
hr_disc <- get_hr(pred_disc, type = 'OS', by = 'PRED')
log_msg('Prediction discovery: HR = %.2f [%.2f-%.2f], p = %.1e', hr_disc[1], hr_disc[2], hr_disc[3], hr_disc[4])
fig_km_disc <- print_plot(plot_mykm(pred_disc, type = 'OS', by = 'PRED'))
ggsave(fig_km_disc2, filename = 'fig_km_disc.pdf', width = 10, height = 10, units = "cm")

library("survival")
library("survminer")
km_fit <- survfit(Surv(OS, Death) ~ PRED, data = pred_disc)
# 绘制Kaplan-Meier曲线
ggsurvplot(km_fit, data = pred_disc, pval = TRUE, conf.int = TRUE,
           xlab = "Time", ylab = "Survival Probability",
           palette = c("#E7B800", "#2E9FDF"))

# ggsurvplot(km_fit, data = pred_disc, 
#            pval = TRUE, conf.int = TRUE,
#            risk.table = TRUE, # Add risk table
#            risk.table.col = "strata", # Change risk table color by groups
#            linetype = "strata", # Change line type by groups
#            surv.median.line = "hv", # Specify median survival
#            ggtheme = theme_bw(), # Change ggplot2 theme
#            palette = c("#E7B800", "#2E9FDF"))

### 1.7 Prediction in discovery cohort (short signature) ----
SIGB1_PCR <- c('Enterocloster_bolteae', 'Clostridium_symbiosum', 'Erysipelatoclostridium_ramosum',
               'Hungatella_hathewayi', 'Veillonella_atypica')
SIGB2_PCR <- c('Anaerostipes_hadrus', 'Blautia_wexlerae', 'Coprococcus_comes', 
               'Dorea_formicigenerans', 'Dorea_longicatena', 'Eubacterium_rectale', 
               'Eubacterium_ventriosum', 'Faecalibacterium_prausnitzii', 
               'Gemmiger_formicilis', 'Phocaeicola_massiliensis', 'Roseburia_hominis', 
               'Roseburia_intestinalis', 'Roseburia_inulinivorans', 
               'Ruminococcus_bicirculans', 'Ruminococcus_lactaris')
scores_disc_short <- compute_toposcore(met4_disc, sigb1 = SIGB1_PCR, sigb2 = SIGB2_PCR)
pred_disc_short <- clin_disc %>% left_join(scores_disc_short, by = 'Sample_id') %>% assign_prediction(ycut_r$THRESHOLD, ycut_nr$THRESHOLD)
hr_disc_short <- get_hr(pred_disc_short, type = 'OS', by = 'PRED')
log_msg('Prediction discovery (short): HR = %.2f [%.2f-%.2f], p = %.1e', hr_disc_short[1], hr_disc_short[2], hr_disc_short[3], hr_disc_short[4])
fig_km_disc_short <- print_plot(plot_mykm(pred_disc_short, type = 'OS', by = 'PRED'))
ggsave(fig_km_disc_short, filename = 'fig_km_disc_short.pdf', width = 10, height = 10, units = "cm")



## 2. Validation analysis set ----
log_msg('####### Validation analysis set #########')
clin_valid <- load_clin(cohort = 'Valid')
met4_valid <- load_microbiome(clin_valid)

### 2.1 Toposcoring ----
scores_valid <- compute_toposcore(met4_valid, sigb1 = SIGB1, sigb2 = SIGB2) 
plt_kde_valid <- plot_toposcoreb01_density(scores_valid, clin_valid, lims = c(ycut_r$THRESHOLD, ycut_nr$THRESHOLD))
ggsave(plt_kde_valid, filename = 'fig_kde_valid.pdf', width = 10, height = 10, units = "cm")


### 2.2 Prediction in validation cohort ----
pred_valid <- clin_valid %>% left_join(scores_valid, by = 'Sample_id') %>% 
  assign_prediction(ycut_r$THRESHOLD, ycut_nr$THRESHOLD)
hr_valid <- get_hr(pred_valid, type = 'OS', by = 'PRED')
log_msg('Prediction validation: HR = %.2f [%.2f-%.2f], p = %.1e', hr_valid[1], hr_valid[2], hr_valid[3], hr_valid[4])
fig_km_valid <- print_plot(plot_mykm(pred_valid, type = 'OS', by = 'PRED'))
ggsave(fig_km_valid, filename = 'fig_km_valid.pdf', width = 10, height = 10, units = "cm")


### 2.3 Prediction in validation cohort (short signature) ----
scores_valid_short <- compute_toposcore(met4_valid, sigb1 = SIGB1_PCR, sigb2 = SIGB2_PCR)
pred_valid_short <- clin_valid %>% left_join(scores_valid_short, by = 'Sample_id') %>% assign_prediction(ycut_r$THRESHOLD, ycut_nr$THRESHOLD)
hr_valid_short <- get_hr(pred_valid_short, type = 'OS', by = 'PRED')
log_msg('Prediction validation (short): HR = %.2f [%.2f-%.2f], p = %.1e', hr_valid_short[1], hr_valid_short[2], hr_valid_short[3], hr_valid_short[4])
fig_km_valid_short <- print_plot(plot_mykm(pred_valid_short, type = 'OS', by = 'PRED'))
ggsave(fig_km_valid_short, filename = 'fig_km_valid_short.pdf', width = 10, height = 10, units = "cm")



## 3. Other analysis sets ----

### 3.1 Healthy donors ----
log_msg('####### HD analysis set #########')
clin_hd <- load_clin(cohort = 'HD') %>% mutate(OS12 = 'Healthy')
met4_hd <- load_microbiome(clin_hd)
scores_hd <- compute_toposcore(met4_hd, sigb1 = SIGB1, sigb2 = SIGB2) 
plot_toposcoreb01_density(scores_hd, clin_hd, lims = c(ycut_r$THRESHOLD, ycut_nr$THRESHOLD))
pred_hd <- clin_hd %>% left_join(scores_hd, by = 'Sample_id') %>% assign_prediction(ycut_r$THRESHOLD, ycut_nr$THRESHOLD)
scores_disc_hd <- bind_rows(scores_disc, scores_hd) 
plt_kde_disc_hd <- plot_toposcoreb01_density(scores_disc_hd, bind_rows(clin_disc, clin_hd), lims = c(ycut_r$THRESHOLD, ycut_nr$THRESHOLD))
ggsave(plt_kde_disc_hd, filename = 'fig_kde_disc_hd.pdf', width = 10, height = 10, units = "cm")


### 3.2 Urothelial analysis set ----
log_msg('####### Urothelial analysis set #########')
clin_uro <- load_clin(cohort = 'Uro')
met4_uro <- load_microbiome(clin_uro)
scores_uro <- compute_toposcore(met4_uro, sigb1 = SIGB1, sigb2 = SIGB2) 
pred_uro <- clin_uro %>% left_join(scores_uro, by = 'Sample_id') %>% assign_prediction(ycut_r$THRESHOLD, ycut_nr$THRESHOLD)
hr_uro <- get_hr(pred_uro, type = 'OS', by = 'PRED')
log_msg('Prediction urothelial: HR = %.2f [%.2f-%.2f], p = %.1e', hr_uro[1], hr_uro[2], hr_uro[3], hr_uro[4])
fig_km_uro <- print_plot(plot_mykm(pred_uro, type = 'OS', by = 'PRED'))
ggsave(fig_km_uro, filename = 'fig_km_uro.pdf', width = 10, height = 10, units = "cm")
scores_uro_short <- compute_toposcore(met4_uro, sigb1 = SIGB1_PCR, sigb2 = SIGB2_PCR)
pred_uro_short <- clin_uro %>% left_join(scores_uro_short, by = 'Sample_id') %>% assign_prediction(ycut_r$THRESHOLD, ycut_nr$THRESHOLD)
hr_uro_short <- get_hr(pred_uro_short, type = 'OS', by = 'PRED')
log_msg('Prediction urothelial (short): HR = %.2f [%.2f-%.2f], p = %.1e', hr_uro_short[1], hr_uro_short[2], hr_uro_short[3], hr_uro_short[4])
fig_km_uro_short <- print_plot(plot_mykm(pred_uro_short, type = 'OS', by = 'PRED'))
ggsave(fig_km_uro_short, filename = 'fig_km_uro_short.pdf', width = 10, height = 10, units = "cm")


### 3.3 RCC analysis set ----
log_msg('####### RCC analysis set #########')
clin_rcc <- load_clin(cohort = 'RCC')
met4_rcc <- load_microbiome(clin_rcc)
scores_rcc <- compute_toposcore(met4_rcc, sigb1 = SIGB1, sigb2 = SIGB2) 
pred_rcc <- clin_rcc %>% left_join(scores_rcc, by = 'Sample_id') %>% assign_prediction(ycut_r$THRESHOLD, ycut_nr$THRESHOLD)
hr_rcc <- get_hr(pred_rcc, type = 'OS', by = 'PRED')
log_msg('Prediction RCC: HR = %.2f [%.2f-%.2f], p = %.1e', hr_rcc[1], hr_rcc[2], hr_rcc[3], hr_rcc[4])
fig_km_rcc <- print_plot(plot_mykm(pred_rcc, type = 'OS', by = 'PRED'))
ggsave(fig_km_rcc, filename = 'fig_km_rcc.pdf', width = 10, height = 10, units = "cm")
scores_rcc_short <- compute_toposcore(met4_rcc, sigb1 = SIGB1_PCR, sigb2 = SIGB2_PCR)
pred_rcc_short <- clin_rcc %>% left_join(scores_rcc_short, by = 'Sample_id') %>% assign_prediction(ycut_r$THRESHOLD, ycut_nr$THRESHOLD)
hr_rcc_short <- get_hr(pred_rcc_short, type = 'OS', by = 'PRED')
log_msg('Prediction RCC (short): HR = %.2f [%.2f-%.2f], p = %.1e', hr_rcc_short[1], hr_rcc_short[2], hr_rcc_short[3], hr_rcc_short[4])
fig_km_rcc_short <- print_plot(plot_mykm(pred_rcc_short, type = 'OS', by = 'PRED'))
ggsave(fig_km_rcc_short, filename = 'fig_km_rcc_short.pdf', width = 10, height = 10, units = "cm")


### 3.4 longitudinal analysis ----
# longitudinal
met4_long <- load_microbiome_longitudinal()
clin_long <- load_clin_longitudinal()
scores_long <- compute_toposcore(met4_long, sigb1 = SIGB1_PCR, sigb2 = SIGB2_PCR) 
clin_long_V0 <- clin_long %>% dplyr::select(Patient_id, Sample_id = Sample_id_V0, AKK_TRICHO = AKK_TRICHO_V0)
pred_long_V0 <- clin_long_V0 %>% left_join(scores_long, by = 'Sample_id') %>% 
  assign_prediction(ycut_r$THRESHOLD, ycut_nr$THRESHOLD) %>%
  dplyr::select(Patient_id, Sample_id_V0 = Sample_id, AKK_TRICHO_V0 = AKK_TRICHO, PRED_V0 = PRED, CAT_V0 = SIGCAT)
clin_long_V3 <- clin_long %>% dplyr::select(Patient_id, Sample_id = Sample_id_V3, AKK_TRICHO = AKK_TRICHO_V3)
pred_long_V3 <- clin_long_V3 %>% left_join(scores_long, by = 'Sample_id') %>% 
  assign_prediction(ycut_r$THRESHOLD, ycut_nr$THRESHOLD) %>%
  dplyr::select(Patient_id, Sample_id_V3 = Sample_id, AKK_TRICHO_V3 = AKK_TRICHO, PRED_V3 = PRED, CAT_V3 = SIGCAT)
pred_long_shifts <- merge(pred_long_V0, pred_long_V3, by = 'Patient_id', all = TRUE)
plt_sankey_pred <- plot_sankey_pred(pred_long_shifts)
ggsave(plt_sankey_pred, filename = 'fig_sankey_pred.pdf', width = 10, height = 10, units = "cm")
desc_sankey_pred(pred_long_shifts)
plt_sankey_sig <- plot_sankey_sigcat(pred_long_shifts)
ggsave(plt_sankey_sig, filename = 'fig_sankey_sig.pdf', width = 10, height = 10, units = "cm")
desc_sankey_sigcat(pred_long_shifts)



# prevalences ----
met4 <- bind_rows(met4_disc, met4_valid)
prevalence_plt <- plot_prevalences(met4)
ggsave(prevalence_plt, filename = 'fig_prevalence.pdf', width = 10, height = 20, units = "cm")
