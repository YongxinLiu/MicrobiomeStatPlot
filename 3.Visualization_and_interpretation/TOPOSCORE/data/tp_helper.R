# R version 4.2.3

################################################################################
# setup
################################################################################

# load required libraries
library(dplyr)    # v1.1.2
library(readxl)   # v1.4.2
library(survival) # v3.5-5
library(ggplot2)  # v3.4.2

# global options
options(dplyr.summarise.inform = FALSE)

# data input directory
#SRCDIR <- '/home/jwojcik/d2t/onedrive/04_Customers/PF04_IGR/PF04.01_TOPOSCORE/standalone_input/'
SRCDIR <- './'
message('Data source directory: ', SRCDIR)

# Expected input files (in SRCDIR):
FILES = list(
  ONCO_CLIN = 'data/DS1_oncology_clinical_data.csv', 
  ONCO_MET4 = 'data/DS2_oncology_microbiome_data.csv',
  HD_CLIN = 'data/DS3_healthy_donor_clinical_data.csv', 
  HD_MET4 = 'data/DS4_healthy_donor_microbiome_data.csv',
  LONGITUDINAL_CLIN = 'data/DS5_longitudinal_clinical_data.csv', 
  LONGITUDINAL_MET4 = 'data/DS6_longitudinal_microbiome_data.csv'
)
# check that they are all in SRCDIR
existence <- vapply(unlist(FILES), function(x) file.exists(file.path(SRCDIR, x)), TRUE)
if (!all(existence)) {
  message('Input files not found: ', paste(unlist(FILES)[which(!existence)], collapse = ', '))
  stop()
}
message('All expected input files found')

# other constants and inits
VARS <- list(
  CENSOR = list(PFS = 'PD', OS = 'Death')
)

# graphical specs
GRAPHICS <- list(
  THEME = ggplot2::theme_light(),
  COLOR = list(
    DEFAULT = 'dark blue',
    SIG = c(ALPHA = 'red', BETA = 'darkorange', GAMMA = 'darkgreen', DELTA = 'chartreuse',
            EPSILON = 'cyan', ETA = 'pink', ZETA = 'brown'),
    OS12 = c(R = 'green', NR = 'red', Healthy = 'blue'),
    OS6bis = c(R = 'green', NR = 'red'),
    OS12bis = c(R = 'green', NR = 'red', `?` = 'gray')
  )
)


  
################################################################################
# util functions
################################################################################


#' Logging util
log_msg <- function(...) {
  msg <- paste(sprintf(...), collapse = '\n')
  message(msg)
}


die_if <- function(expr, ...) {
  if (expr) {
    stop(sprintf(...))
  }
}


#' Read a XLSX file
load_xlsx <- function(filename, dir, sheet, ...) {
  path <- file.path(dir, filename)
  stopifnot(file.exists(path))
  suppressMessages(readxl::read_xlsx(path, sheet = sheet, ...))
}


#' Read a CSV file
read_csv <- function(path) {
  stopifnot(file.exists(path))
  read.csv(path, stringsAsFactors = FALSE)
}


#' Util function to speed-up the re-analysis script
#' If the RDS file exists, read it and return it
#' If not, execute the expression, save to RDS, and return
load_or_compute <- function(rds_file, expr) {
  if (file.exists(rds_file)) {
    log_msg('Results read from existing file %s', rds_file)
    readRDS(rds_file)
  } else {
    out <- expr
    saveRDS(out, rds_file)
    log_msg('Results saved to file %s', rds_file)
    out
  }
}


color_values <- function(var, df) {
  colors <- GRAPHICS$COLOR[[var]]
  if (is.null(colors)) return(topo.colors(length(unique(df[[var]]))))
  colors[as.character(unique(df[[var]]))]
}


shape_values <- function(var, df) {
  shapes <- GRAPHICS$SHAPE[[var]]
  if (is.null(shapes)) return(1:length(unique(df[[var]])))
  shapes[as.character(unique(df[[var]]))]
}


get_met4_presences <- function(met4) {
  data_pres <- as.data.frame(ifelse(met4[, -1] > 0, 1, 0))
  cbind(Sample_id = met4$Sample_id, data_pres)
}



geomean <- function(x) {
  exp(mean(log(x)))
}


inner_legend <- function() {
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1), 
        legend.key.width = unit(1, "lines"), 
        legend.key.height = unit(1, "lines"), 
        plot.margin = unit(c(5, 1, 0.5, 0.5), "lines"))
}


print_plot <- function(plt, apply_theme = TRUE) {
  suppressWarnings(
    if (inherits(plt, 'gList') || inherits(plt, 'gtable') || inherits(plt, 'gTree')) {
      grid::grid.newpage()
      grid::grid.draw(plt)
    } else {
      if (apply_theme && exists('GRAPHICS') && !is.null(GRAPHICS$THEME))
        plt <- plt + GRAPHICS$THEME
      graphics::plot(plt)
    }
  )
}


format_pvalue <- function(x) {
  fo <- if (x < 0.0001) '%.1e' else '%.4f'
  sprintf(fo, x)
}


to_pct <- function(x) {
  sprintf('%.0f%%', 100 * x)
}


################################################################################
# loaders
################################################################################

#' Load clinical data
load_clin <- function(cohort, dir = SRCDIR) {
  clin <- if (cohort == 'HD' || cohort == 'Healthy')
    read_csv(file.path(dir, FILES$HD_CLIN)) 
  else
    read_csv(file.path(dir, FILES$ONCO_CLIN)) %>% filter(Cohort == cohort) 
  log_msg('Reading %s clinical data records from cohort %s', nrow(clin), cohort)
  clin
}


#' Load microbiome data 
load_microbiome <- function(clin, dir = SRCDIR) {
  met4 <- if (unique(clin$Cohort) == 'Healthy') 
    read_csv(file.path(dir, FILES$HD_MET4))
  else
    read_csv(file.path(dir, FILES$ONCO_MET4)) %>% filter(Sample_id %in% clin$Sample_id)
  log_msg('Reading %s microbiome data records', nrow(met4))
  met4
}



load_microbiome_longitudinal <- function(filename = FILES$LONGITUDINAL_MET4, dir = SRCDIR) {
  met4 <- read_csv(file.path(dir, filename))
  log_msg('Reading %s microbiome data records', nrow(met4))
  met4
}


load_clin_longitudinal <- function(filename = FILES$LONGITUDINAL_CLIN, dir = SRCDIR) {
  clin <- read_csv(file.path(dir, filename))
  log_msg('Reading %s clinical data records', nrow(clin))
  clin
}


################################################################################
# survival analysis
################################################################################

get_surv_formula <- function(time_var, censor_var, term_vars, strata_var = NULL) {
  fo_str <- sprintf('Surv(%s, %s) ~ %s', time_var, censor_var, paste(term_vars, collapse = ' + '))
  if (!is.null(strata_var)) fo_str <- sprintf('%s + strata(%s)', fo_str, strata_var)
  as.formula(fo_str)
}



calc_surv_df <- function(data, time_var, censor_var, term_vars = '1') {
  fo <- get_surv_formula(time_var, censor_var, term_vars)
  fit <- survival::survfit(fo, data)
  df <- as.data.frame(broom::tidy(fit))
  if (term_vars != '1' && !'strata' %in% colnames(df))
    df$strata <- paste0(term_vars, '=-')
  df0 <- data.frame(time = 0, estimate = 1, conf.high = 1, conf.low = 1, 
                    strata = unique(df$strata), stringsAsFactors = FALSE)
  df <- bind_rows(df, df0) %>% arrange(strata, time)
  df
}


scale_x_months <- function(data) {
  months <- seq(from = 0, to = max(data$time), by = 12)
  scale_x_continuous(breaks = months)
}


plot_km <- function(data, time_var, censor_var, term_vars = '1',
                    confint = TRUE, censor_ticks = TRUE, med = NULL) {
  
  df <- calc_surv_df(data, time_var, censor_var, term_vars)
  strata <- if ('strata' %in% colnames(df)) 'strata' else NULL
  
  # plot
  plt <- if (is.null(strata))
    ggplot(data = df,  aes(x = time, y = estimate))
  else
    ggplot(data = df,  aes(x = time, y = estimate, color = .data[[strata]], fill = .data[[strata]]))
  
  plt <- plt + geom_step() +
    ylim(0, 1) + ylab(paste('Time to', time_var)) + scale_x_months(df)
  
  if (!is.null(med))
    plt <- plt + 
    scale_color_discrete(labels = sprintf("%s [%.2fM]", med$strata, med$MEDSURV)) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'gray')
  
  # confidence intervals?
  if (confint) plt <- plt +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                show.legend = FALSE, linetype = 'blank', alpha = 0.2)
  
  # censored ticks?
  if (censor_ticks) plt <- plt +
    geom_segment(aes(y = estimate - 0.01, yend = estimate + 0.01, xend = time), 
                 data = df[df$n.censor != 0, ])
  
  plt + GRAPHICS$THEME
  
}


plot_risk_table <- function(km_plt, nb_labels = 10, size = 3) {
  
  data <- km_plt$data
  
  # subset
  .subset_risk_ta <- function(df) {
    time_interval <- max(round(max(df$time) / nb_labels), 1)
    sq <- seq(from = min(df$time), to = max(df$time), by = time_interval)
    i <- vapply(sq, function(x) which(df$time >= x)[1], 1L)
    df[i, ]
  }
  data <- if ('strata' %in% colnames(data)) {
    bind_rows(lapply(unique(data$strata), function(x) .subset_risk_ta(data[data$strata == x, ])))
  } else {
    .subset_risk_ta(data)
  }
  
  # strata?
  plt <- if ('strata' %in% colnames(data)) {
    st <- strsplit(data$strata, split = '=')
    data$strata_light <- vapply(st, function(x) x[2], 'str')
    ggplot(data, aes(x = time, y = strata_light, label = n.risk, color = strata)) +
      theme(axis.text.y = element_text(hjust = 1 )) + ylab('strata')
  } else {
    ggplot(data, aes(x = time, y = 0, label = n.risk)) +
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + ylab('')
  }
  
  plt +
    geom_text(size = size) +
    theme(axis.title.x = element_text(vjust = 1), legend.position = "none") + 
    xlab('Time (months)') + scale_x_months(data) +
    GRAPHICS$THEME
  
}


plot_km_and_risk_table <- function(data, time_var, censor_var, term_vars = '1',
                                   confint = TRUE, censor_ticks = TRUE, hr = NULL, med = NULL,
                                   nb_labels = 10, size = 3) {
  
  th <- theme(
    panel.background = element_rect(fill = NA, color = "gray"),
    panel.grid.major.y =  element_line( color = 'grey', linetype = 3),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x =  element_line( color = 'grey', linetype = 3),
    panel.grid.minor.x =  element_line( color = 'grey', linetype = 3)
  ) 
  
  km_plt <- plot_km(data, time_var, censor_var, term_vars, confint, censor_ticks, med) + th + inner_legend()
  ta_plt <- plot_risk_table(km_plt, nb_labels, size) + th
  
  if (!is.null(hr)) {
    hr_str <- sprintf('HR %.2f [%.2f;%.2f]\np = %s', hr[1], hr[2], hr[3], format_pvalue(hr[4]))
    km_plt <- km_plt +
      guides(color = guide_legend(hr_str), fill = guide_legend(hr_str))
  }
  
  km_grb <- ggplot2::ggplotGrob(km_plt)
  ta_grb <- ggplot2::ggplotGrob(ta_plt)
  
  # build a 3x3 gtable (main panel is in grob layout [7, 5])
  ng <- grid::nullGrob()
  grobs <- list(
    km_grb[7, 1:4],   # KM y-axis
    km_grb[7, 5],     # KM plot
    km_grb[7, 6:ncol(km_grb)],   # KM legend
    ta_grb[7, 1:4],   # risk table y-axis
    ta_grb[7, 5],     # risk table plot
    ng,
    ng,  
    ta_grb[8:9, 5],   # risk table x-axis
    ng
  )
  
  plt <- gtable::gtable_matrix("KM_risk", 
                               grobs = matrix(grobs, ncol = 3, nrow = 3, byrow = TRUE),
                               widths = unit(c(1, 5, 1), "null"), heights = unit(c(5, 1, 0.5), "null")
  )
  
  plt
}


plot_mykm <- function(data, type = 'PFS', by = NULL, max = NULL,...) {
  
  type <- match.arg(type, c('PFS', 'OS'))
  data <- data[!is.na(data[[type]]), ]
  
  # censoring after max?
  if (!is.null(max)) {
    k <- which(data[[type]] >= max)
    if (length(k) > 0) {
      data[k, type] <- max
      data[k, VARS$CENSOR[[type]]] <- 0
    }
  }
  
  if (!is.null(by)) data <- data[!is.na(data[[by]]) & !is.nan(data[[by]]) & data[[by]] != 'NaN', ]
  term_vars <- if (is.null(by)) '1' else by
  
  hr <- if (is.null(by)) NULL else get_hr(data, type, by, term_vars = NULL)
  med <- get_median_surv_time(data, type, by)
  
  plot_km_and_risk_table(data, time_var = type, 
                         censor_var = VARS$CENSOR[[type]], term_vars = term_vars, 
                         hr = hr, med = med, ...)
  
}


categorize_met4_for_surv <- function(clin, met4, by, cutoff = NULL, verbose = FALSE) {
  
  die_if(!by %in% colnames(met4), 'Cannot find %s', by)
  
  if (!is.null(cutoff)) {
    if (verbose) log_msg('Categorizing %s Low/High by cutoff (%f)', by, cutoff)
    met4[[by]] <- ifelse(is.na(met4[[by]]), NA, ifelse(met4[[by]] < cutoff, 'Low', 'High'))
    met4[[by]] <- factor(met4[[by]], levels = c('Low', 'High'))
    met4$MEDIAN <- cutoff
    
  } else {
    met4_median <- median(met4[[by]], na.rm = TRUE)
    if (met4_median != 0) {
      if (verbose) log_msg('Categorizing %s Low/High by median (%f)', by, met4_median)
      met4[[by]] <- ifelse(is.na(met4[[by]]), NA, ifelse(met4[[by]] <= met4_median, 'Low', 'High'))
      met4[[by]] <- factor(met4[[by]], levels = c('Low', 'High'))
      met4$MEDIAN <- met4_median
    } else {
      if (verbose) log_msg('Categorizing %s Absence/Presence', by)
      met4[[by]] <- ifelse(is.na(met4[[by]]), NA, ifelse(met4[[by]] == 0, 'Absence', 'Presence'))
      met4[[by]] <- factor(met4[[by]], levels = c('Absence', 'Presence'))
      met4$MEDIAN <- NA
    }
  }
  
  data <- clin %>% left_join(met4[, c('Sample_id', by, 'MEDIAN')], by = 'Sample_id')
  data[!is.na(data[[by]]), ]
  
}


plot_mykm_met4 <- function(clin, met4, type = 'PFS', by = NULL, cutoff = NULL, ...) {
  data <- categorize_met4_for_surv(clin, met4, by, cutoff, verbose = TRUE)
  plot_mykm(data, type = type, by = by, ...)
}


get_median_surv_time <- function(data, type = 'PFS', by = NULL) {
  
  type <- match.arg(type, c('PFS', 'OS'))
  data <- data[!is.na(data[[type]]), ]
  
  if (!is.null(by)) data <- data[!is.na(data[[by]]) & !is.nan(data[[by]]) & data[[by]] != 'NaN', ]
  term_vars <- if (is.null(by)) '1' else by
  
  df <- calc_surv_df(data, type, censor_var = VARS$CENSOR[[type]], term_vars = term_vars)
  
  .calc_medsurv <- function(df) {
    medsurv <- if (min(df$estimate > 0.5)) {
      NA
    } else {
      df$time[which(df$estimate <= 0.5)[1]]
    }
    data.frame(N = nrow(df), N_EVENTS = sum(df$n.event), MEDSURV = medsurv)
  }
  
  if (is.null(by))
    .calc_medsurv(df)
  else
    df %>% group_by(strata) %>% do(.calc_medsurv(.)) %>% as.data.frame()
  
}


coxph_test <- function(data, type = 'PFS', by, term_vars = NULL) {
  
  type <- match.arg(type, c('PFS', 'OS'))
  data <- data[!is.na(data[[type]]), ]
  
  fo <- get_surv_formula(time_var = type, censor_var = VARS$CENSOR[[type]], 
                         term_vars = c(term_vars, by))
  fit <- survival::coxph(fo, data)
  
  summary(fit)
  
}


get_hr <- function(data, type = 'PFS', by, term_vars = NULL) {
  
  sfit <- coxph_test(data, type, by, term_vars)
  
  hrs <- cbind(
    sfit$conf.int[, c(1, 3, 4), drop = FALSE],
    sfit$coefficients[, 5, drop = FALSE]
  )
  colnames(hrs) <- c('HR', 'HR_LCI', 'HR_UCI', 'P')
  
  hrs[grep(by, rownames(hrs)), , drop = TRUE]
  
}


get_hr_met4 <- function(var, clin, met4, type = 'OS', cutoff = NULL, verbose = TRUE) {
  data <- categorize_met4_for_surv(clin, met4, by = var, cutoff = cutoff, verbose = verbose)
  hr <- get_hr(data, type, var)
  as.data.frame(as.list(hr)) %>% 
    mutate(SPECIES = var, MEDIAN = unique(data$MEDIAN))
}


screen_surv_met4 <- function(clin, met4, type = 'OS', vars = colnames(met4)[-1],
                             verbose = FALSE) {
  
  res_list <- lapply(vars, get_hr_met4, clin = clin, met4 = met4, type = type, verbose = verbose)
  
  bind_rows(res_list) %>% arrange(P)
  
}


plot_surv_forest <- function(res_surv, species = res_surv$SPECIES, alpha = 1,
                             color = 'P') {
  
  data <- res_surv %>% 
    filter(P <= alpha & SPECIES %in% species) %>% 
    mutate(CAT_TYPE = ifelse(is.na(MEDIAN), 'abs/pres', 'low/high'))
  
  if (color == 'P')
    data <- data %>% mutate(P = ifelse(P <= 0.05, 'p<=0.05', 'p>0.05'))
  
  data$SPECIES <- factor(data$SPECIES, levels = data$SPECIES)
  
  ggplot(data, aes(x = SPECIES, y = HR, color = .data[[color]])) +
    geom_segment(aes(y = HR_LCI, yend = HR_UCI, xend = SPECIES)) +
    geom_hline(yintercept = 1, color = 'blue', linetype = 'dotted') +
    scale_color_manual(values = color_values(color, data)) +
    geom_point(aes(shape = CAT_TYPE)) +
    coord_flip() +
    GRAPHICS$THEME + xlab('')
  
}



plot_sankey_pred <- function(data, na.rm = TRUE) {
  
  if (na.rm) data <- data %>% filter(!is.na(PRED_V0) & ! is.na(PRED_V3))
  data <- data %>% 
    dplyr::rename(Baseline = PRED_V0, Month3 = PRED_V3) %>%
    ggsankey::make_long(Baseline, Month3)
  
  ggplot(data, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node))) +
    scale_fill_manual(name = 'Prediction', values = c(`NA` = 'gray', R = 'darkgreen', NR = 'red')) +
    ggsankey::geom_sankey(flow.alpha = 0.5) +
    #ggsankey::geom_sankey_label(aes(label = node)) +
    ggsankey::theme_sankey(base_size = 16) + xlab('')
  
}


desc_sankey_pred <- function(data, na.rm = TRUE) {
  
  if (na.rm) data <- data %>% filter(!is.na(PRED_V0) & ! is.na(PRED_V3))
  
  data %>% 
    group_by(PRED_V0) %>%
    summarize(N_TOT = n(), N_NR_V3 = sum(PRED_V3 == 'NR'), N_R_V3 = sum(PRED_V3 == 'R')) %>%
    ungroup() %>%
    mutate(PC_NR_V3 = to_pct(N_NR_V3 / N_TOT), PC_R_V3 = to_pct(N_R_V3 / N_TOT)) %>%
    as.data.frame()
  
}


desc_sankey_sigcat <- function(data, na.rm = TRUE) {
  
  if (na.rm) data <- data %>% filter(!is.na(CAT_V0) & ! is.na(CAT_V3))
  
  data %>% 
    group_by(CAT_V0) %>%
    summarize(N_TOT = n(), N_SIG1_V3 = sum(CAT_V3 == 'SIG1'), 
              N_Gray_V3 = sum(CAT_V3 == 'Gray'), N_SIG2_V3 = sum(CAT_V3 == 'SIG2')) %>%
    ungroup() %>%
    mutate(PC_SIG1_V3 = to_pct(N_SIG1_V3 / N_TOT), 
           PC_Gray_V3 = to_pct(N_Gray_V3 / N_TOT), PC_SIG2_V3 = to_pct(N_SIG2_V3 / N_TOT)) %>%
    as.data.frame()
  
}


plot_sankey_sigcat <- function(data, na.rm = TRUE) {
  
  if (na.rm) data <- data %>% filter(!is.na(CAT_V0) & ! is.na(CAT_V3))
  data <- data %>% 
    dplyr::rename(Baseline = CAT_V0, Month3 = CAT_V3) %>%
    ggsankey::make_long(Baseline, Month3)
  data$node <- factor(data$node, levels = c('SIG1', 'Gray', 'SIG2'))
  
  ggplot(data, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node))) +
    scale_fill_manual(name = 'Category', values = c(Gray = 'gray', SIG2 = 'darkgreen', SIG1 = 'red')) +
    ggsankey::geom_sankey(flow.alpha = 0.5) +
    #ggsankey::geom_sankey_label(aes(label = node)) +
    ggsankey::theme_sankey(base_size = 16) + xlab('')
  
}


################################################################################
# correlation analysis & clustering
################################################################################

test_pair <- function(pair, met4) {
  
  #log_msg("pair <- c('%s', '%s')", pair[1], pair[2])
  
  met4_filt <- setNames(met4[, pair], c('VAR1', 'VAR2')) %>%
    mutate(CAT1 = sign(VAR1), CAT2 = sign(VAR2)) %>%
    mutate(CAT = paste0(CAT1, CAT2))
  
  n00 <- sum(met4_filt$CAT == '00')
  n01 <- sum(met4_filt$CAT == '01')
  n10 <- sum(met4_filt$CAT == '10')
  n11 <- sum(met4_filt$CAT == '11')
  mat <- matrix(c(n00, n01, n10, n11), nrow = 2)
  ftest <- fisher.test(mat)
  
  met4_filt11 <- met4_filt %>% filter(CAT == '11')
  ptest0 <- if (nrow(met4_filt11) > 2) {
    cor.test(met4_filt11$VAR1, met4_filt11$VAR2, method = 'pearson')
  } else {
    list(estimate = NA, p.value = NA)
  }
  
  data.frame(VAR1 = pair[1], VAR2 = pair[2],
       N00 = n00, N01 = n01, N10 = n10, N11 = n11, FISHER_P = ftest$p.value,
       OR = unname(ftest$estimate), OR_LCI = ftest$conf.int[1], OR_UCI = ftest$conf.int[2],
       PEARSON11_P = ptest0$p.value, PEARSON11_R = ptest0$estimate,
       stringsAsFactors = FALSE)
  
}


screen_pairs <- function(met4, vars = colnames(met4)[-1]) {
  
  grid <- expand.grid(VAR1 = vars, VAR2 = vars, stringsAsFactors = FALSE) %>% 
    filter(VAR1 < VAR2)
  
  res <- grid %>%
    group_by(VAR1, VAR2) %>%
    do(test_pair(pair = c(.$VAR1, .$VAR2), met4)) %>%
    ungroup()
  
  as.data.frame(res)
  
}


get_score_matrix <- function(res, score, vars) {
  
  if (!is.null(vars)) {
    res <- res %>% filter(res$VAR1 %in% vars & res$VAR2 %in% vars)
  }
  
  res_dup <- res
  colnames(res_dup)[1:2] <- colnames(res_dup)[2:1]
  all <- bind_rows(res, res_dup) 
  
  all$SCORE <- switch(score,
                      fisher_p = -log10(all$FISHER_P) * sign(all$OR - 1),
                      pearson11_p = -log10(all$PEARSON11_P) * sign(all$PEARSON11_R),
                      pearson_p = -log10(all$PEARSON_P) * sign(all$PEARSON_R),
                      meta = all$META,
                      pearson_r = all$PEARSON_R
  )
  
  all <- all %>% 
    dplyr::select(SPECIES = VAR1, VAR2, SCORE) %>% 
    tidyr::pivot_wider(names_from = VAR2, values_from = SCORE)
  all <- as.data.frame(all[match(colnames(all)[-1], all$SPECIES), ])
  rownames(all) <- all$SPECIES
  
  all
  
}


cluster_species <- function(res, score = 'fisher_p', vars = NULL, 
                            method = 'complete', distance = 'euclidean', ...) {
  
  # correlation matrix 
  cmat <- get_score_matrix(res, score, vars)[, -1]
  
  # clustering
  cdist <- stats::dist(cmat, method = distance)
  cclus <- stats::hclust(cdist, method = method)
  
  # cut tree and list clusters
  ct <- cutree(cclus, ...)
  
  cluster_df <- data.frame(CLUSTER = paste0('C', ct), SPECIES = names(ct), stringsAsFactors = FALSE)
  rownames(cluster_df) <- NULL
  
  log_msg('Clustering %d markers in %d clusters (%s / %s)',
          nrow(cluster_df), length(unique(cluster_df$CLUSTER)),
          method, distance)
  
  invisible(cluster_df)
  
}


renumber_clusters <- function(cc) {
  
  CLUSTER_MAP <- data.frame(
    CLUSTER = c('C6', 'C2', 'C4', 'C3', 'C5', 'C1', 'C7'),
    NEW_CLUSTER = paste0('C', 1:7), stringsAsFactors = FALSE
  )
  
  cc %>%
    left_join(CLUSTER_MAP, by = 'CLUSTER') %>%
    mutate(CLUSTER = NEW_CLUSTER) %>%
    dplyr::select(-NEW_CLUSTER)
  
}


plot_heatmap <- function(cmat, breaks = NULL, fontsize = 6, annots = get_sigs_df(),
                         method = 'complete', distance = 'euclidean') {
  
  annots <- annots %>% filter(SPECIES %in% rownames(cmat))
  annots <- annots[match(rownames(cmat), annots$SPECIES), ]
  rownames(annots) <- rownames(cmat)
  annots$SPECIES <- NULL
  
  n_colors <- if (is.null(breaks)) 100 else length(breaks)
  color <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(n_colors)
  
  pheatmap::pheatmap(cmat, scale = 'none', revC = TRUE, fontsize = fontsize,
                     color = color, breaks = breaks, annotation_row = annots,
                     clustering_distance_rows = distance, clustering_distance_cols = distance,
                     clustering_method = method, show_colnames = FALSE)
  
}


plot_score_matrix <- function(res, score, vars = NULL, annots = list(get_sigs_df()),
                              method = 'complete', distance = 'euclidean', ...) {
  
  all <- get_score_matrix(res, score, vars)
  
  annots_df <- annots[[1]]
  if (length(annots) > 1)
    for (k in 2:length(annots))
      annots_df <- merge(annots_df, annots[[k]], by = 'SPECIES', all = TRUE)
  
  breaks <- NULL
  if (score == 'pearson_r') breaks <- seq(from = -1, to = 1, by = 0.1)
  
  plot_heatmap(all[, -1], breaks = breaks, annots = annots_df,
               method = method, distance = distance, ...)
  
}


################################################################################
# Toposcore
################################################################################

compute_toposcore <- function(met4, sigb1, sigb2) {
  
  data_pres <- get_met4_presences(met4)
  
  .topo <- function(species, label) {
    species_selected <- intersect(species, colnames(data_pres))
    if (length(species_selected) < length(species)) {
      missing_species <- setdiff(species, species_selected)
      log_msg('%s: %d species not found (%s) -> considered absent', label, length(missing_species),
              paste(missing_species, collapse = ', '))
    }
    
    res <- data.frame(Sample_id = data_pres$Sample_id, N = length(species_selected),
                      COUNT = as.numeric(rowSums(data_pres[, species_selected])),
                      FREQ = as.numeric(rowSums(met4[, species_selected])))
    colnames(res)[2:4] <- c(paste0('N_', label), label, paste0('FREQ_', label))
    res
  }  
  
  res <- merge(.topo(sigb1, 'SIGB1'), .topo(sigb2, 'SIGB2'), by = 'Sample_id')
  
  res <- res %>%
    mutate(SSIGB1 = SIGB1 / length(sigb1), SSIGB2 = SIGB2 / length(sigb2)) %>%
    mutate(TOPOB = SSIGB2 - SSIGB1, TOPOBFREQ = FREQ_SIGB2 - FREQ_SIGB1) %>%
    mutate(TOPOB01 = (TOPOB + 1) / 2)
  
  res
  
}



plot_toposcore_density <- function(scores, clin, var = 'TOPO', color = 'OS12bis', 
                                   xlim = NULL, grey = NULL) {
  
  data <- scores %>% 
    left_join(clin, by = 'Sample_id')
  
  if (color == 'OS12bis')
    data <- data %>% mutate(OS12bis = ifelse(OS12 == '', '?', OS12))
  
  data <- data[!is.na(data[[color]]) & data[[color]] != '', ]
  
  plt <- ggplot(data, aes(x = .data[[var]], color = .data[[color]], fill = .data[[color]])) + 
    geom_density(alpha = 0.25) +
    scale_color_manual(values = color_values(color, data)) +
    scale_fill_manual(values = color_values(color, data)) +
    GRAPHICS$THEME
  
  if (!is.null(grey)) plt <- plt + geom_vline(xintercept = grey, color = 'dark gray', linewidth = 1.5)
  if (!is.null(xlim)) plt <- plt + xlim(xlim)
  
  plt
  
}


plot_toposcoreb01_density <- function(scores, clin, resp = 'OS12', lims = c(0.5, 0.75)) {
  plot_toposcore_density(scores, clin, var = 'TOPOB01', color = resp, xlim = c(0, 1), grey = lims) +
    xlab('Score')
}


assign_prediction <- function(pred, cut_r = 0.75, cut_nr = 0.5, verbose = TRUE) {
  
  if (verbose)
    log_msg('Prediction: NR [%d] < %.2f < gray_zone [%d] < %.2f < R [%d]', 
            length(which(pred$TOPOB01 <= cut_nr)), cut_nr, 
            length(which(pred$TOPOB01 > cut_nr & pred$TOPOB01 < cut_r)),
            cut_r, length(which(pred$TOPOB01 >= cut_r)))
  
  pred <- pred %>%
    mutate(PRED = ifelse(TOPOB01 <= cut_nr, 'NR',
                         ifelse(TOPOB01 >= cut_r, 'R',
                                ifelse(AKK_TRICHO == 'Low', 'R', 'NR')))) %>%
    mutate(SIGCAT = ifelse(TOPOB01 <= cut_nr, 'SIG1',
                           ifelse(TOPOB01 >= cut_r, 'SIG2', 'Gray')))
  
}


get_prediction <- function(pred, resp_var = 'OS12bis') {
  
  tab <- as.data.frame(table(pred[, c(resp_var, 'PRED')])) %>%
    tidyr::pivot_wider(names_from = PRED, values_from = Freq) %>%
    mutate(TOT = NR + R, CORRECT = NA, ACCURACY = NA)
  
  nr <- which(tab[[1]] == 'NR')
  tab$CORRECT[nr] <- tab$NR[nr] / tab$TOT[nr]
  r <- which(tab[[1]] == 'R')
  tab$CORRECT[r] <- tab$R[r] / tab$TOT[r]
  
  tab$ACCURACY[1] <- (tab$NR[nr] + tab$R[r]) / sum(tab$TOT)
  
  as.data.frame(tab)
  
}


calculate_shits <- function(pred_long, baseline_cutoff = 15, month3_cutoffs = c(60, 120)) {
  
  .shift_pred <- function(df) {
    v0 <- df %>% filter(df$DAYS <= baseline_cutoff)
    if (nrow(v0) > 1) v0 <- v0[which.min(abs(v0$DAYS))[1], ]
    if (nrow(v0) == 1)
      res0 <- v0 %>% dplyr::select(SAMPLE_V0 = Sample_id, DAY_V0 = DAYS, SCORE_V0 = TOPOB01, 
                                   CAT_V0 = SIGCAT, AKK_V0 = AKK_TRICHO, PRED_V0 = PRED)
    v3 <- df %>% filter(df$DAYS >= month3_cutoffs[1] & df$DAYS <= month3_cutoffs[2])
    if (nrow(v3) > 1) v3 <- v3[which.min(abs(v3$DAYS - 90))[1], ]
    if (nrow(v3) == 1)
      res3 <- v3 %>% dplyr::select(SAMPLE_V3 = Sample_id, DAY_V3 = DAYS, SCORE_V3 = TOPOB01, 
                                   CAT_V3 = SIGCAT, AKK_V3 = AKK_TRICHO, PRED_V3 = PRED)
    res <- if (nrow(v0) == 1 & nrow(v3) == 1) {
      cbind(res0, res3) %>% mutate(
        CAT_SHIFT = sprintf('%s-%s', CAT_V0, CAT_V3),
        PRED_SHIFT = sprintf('%s-%s', PRED_V0, PRED_V3)
      )
    } else if (nrow(v0) == 1) {
      res0
    } else if (nrow(v3) == 1){
      res3
    } else {
      data.frame(NOTHING = TRUE)
    }
    res
  }
  
  pred_long %>% group_by(Patient_ID) %>% 
    do(.shift_pred(.)) %>% ungroup() %>% 
    filter(is.na(NOTHING)) %>% dplyr::select(-NOTHING) %>%
    as.data.frame()
  
}


calc_roc <- function(true_resp, pred_prob, verbose = FALSE) {
  
  roc_obj <- pROC::roc(as.numeric(true_resp) - 1, pred_prob, 
                       auc = TRUE, ci = TRUE, quiet = !verbose)
  
  roc_df <- data.frame(
    THRESHOLD = rev(roc_obj$thresholds),
    TPR = rev(roc_obj$sensitivities), 
    FPR = rev(1 - roc_obj$specificities)
  )
  
  list(AUC = as.numeric(roc_obj$ci), ROC_DF = roc_df)
}


plot_roc <- function(roc_df) {
  
  ggplot(roc_df, aes(x = FPR, y = TPR)) +
    geom_step() +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'blue') +
    xlim(0, 1) + ylim(0, 1) +
    xlab('1 - Specificity') + ylab('Sensitivity') +
    GRAPHICS$THEME
  
}



plot_prevalences <- function(met4, sigb1 = SIGB1,  
                             sigb2 = SIGB2, pcr = c(SIGB1_PCR, SIGB2_PCR)) {
  
  pres <- get_met4_presences(met4)[, c(sigb1, sigb2, 'Akkermansia_muciniphila')]
  prevalences <- 100 * colSums(pres) / nrow(pres)
  
  data <- data.frame(SPECIES = names(prevalences), PREVALENCE = prevalences,
                     stringsAsFactors = FALSE) %>%
    mutate(TYPE = ifelse(SPECIES == 'Akkermansia_muciniphila', 'Akk',
                         ifelse(SPECIES %in% sigb1, 'SIG1', 'SIG2'))) %>%
    mutate(PCR = ifelse(SPECIES %in% c(pcr, 'Akkermansia_muciniphila'), 'Y', 'N')) %>%
    arrange(PREVALENCE) %>%
    mutate(SPECIES = gsub('_', ' ', SPECIES, fixed = TRUE))
  data$SPECIES <- factor(data$SPECIES, levels = data$SPECIES)
  
  ggplot(data, aes(x = SPECIES, y = PREVALENCE)) +
    geom_bar(aes(fill = TYPE), stat = 'identity') +
    scale_fill_manual(name = '', values = c(SIG1 = 'orange', SIG2 = 'cornflowerblue', Akk = 'pink')) +
    geom_point(y = 0, shape = 18, data = data %>% filter(PCR == 'Y'), color = 'darkgray') +
    ylim(0, 100) + ylab('Prevalence (%)') + xlab('') +
    coord_flip() + expand_limits(x = 0) +
    GRAPHICS$THEME +
    theme(axis.text.y = element_text(size = 6, hjust = 1), panel.grid.major.y = element_blank(),
          panel.border = element_blank(), axis.ticks.y = element_blank(),
          legend.position = 'bottom')
  
}