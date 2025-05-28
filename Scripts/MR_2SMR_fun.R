
#!/usr/bin/env Rscript

# === Two-Sample Mendelian Randomization (MR) Analysis ===
# Input:
#   1. Merged file with shared SNPs between exposure and outcome
#   2. Exposure phenotype name
# Output:
#   MR results, PRESSO results, plots and summary tables

cat(paste("Start Mendelian Randomization analysis at", format(Sys.time(), "%c"), "\n"))  
options(scipen = 999)

# Load libraries
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(devtools)
  library(TwoSampleMR)
  library(na.tools)
})

# Read arguments
args <- commandArgs(trailingOnly = TRUE)
shared_file <- args[1]
pheno_id <- args[2]

# Read shared beta/SE data
Shared_beta_SE <- read.table(shared_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, na.strings = "unknown")

# Extract identifiers
Outcome_id <- strsplit(basename(shared_file), "\\.")[[1]][1]
Pheno_id <- pheno_id
Ph_temp <- strsplit(pheno_id, "\\.")[[1]]
Ph_temp_2 <- Ph_temp[2]

# Adjust rsID if needed (for GIANT exposures)
if (Ph_temp_2 == "giant-ukbb") {
  cat("GIANT phenotype detected\n")
  rsID_file <- paste0("rsID_", pheno_id)
  Shared_rs_exposure <- read.table(rsID_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df <- merge(Shared_beta_SE, Shared_rs_exposure, by.x = "CHR_BP", by.y = "CHR_POS")

  for (i in seq_len(nrow(df))) {
    if (df$effect_allele.exposure[i] != df$effect_allele.outcome[i]) {
      df$beta.exposure[i] <- -df$beta.exposure[i]
      df$effect_allele.exposure[i] <- ifelse(df$effect_allele.exposure[i] == df$Other_Allele[i],
                                             df$effect_allele.outcome[i], df$Other_Allele[i])
    }
  }

  Shared_beta_SE <- df
  Shared_beta_SE$SNP.y <- NULL
  colnames(Shared_beta_SE) <- c("SNP", "effect_allele.outcome", "pval.outcome", "beta.outcome", "se.outcome", 
                                "CHR_BP", "effect_allele.exposure", "beta.exposure", "se.exposure", 
                                "pval.exposure", "Other_allele")
} else {
  cat("Non-GIANT phenotype (assumed NEALE)\n")
}

# Format data for MR
MR_file_exposure <- Shared_beta_SE %>%
  transmute(SNP, beta.exposure, se.exposure, effect_allele.exposure, pval.exposure, id.exposure = Pheno_id)

MR_file_outcome <- Shared_beta_SE %>%
  transmute(SNP, beta.outcome, se.outcome, effect_allele.outcome, pval.outcome, id.outcome = Outcome_id)

PHENO_dat <- format_data(MR_file_exposure, type = "exposure")
GWAS_dat <- format_data(MR_file_outcome, type = "outcome")

# Harmonize and run MR
datos <- harmonise_data(PHENO_dat, GWAS_dat, action = 1)
resultado <- mr(datos)
res_plot <- mr(datos, method_list = c("mr_ivw", "mr_egger_regression"))
cat("MR methods completed\n")

# MR-PRESSO
fitDis <- round(nrow(datos) / 0.05, 0)
presso_res <- run_mr_presso(datos, NbDistribution = fitDis, SignifThreshold = 0.05)
cat("MR-PRESSO completed\n")

# Detect outliers
outliers_alert <- presso_res[[1]][[1]][[6]][2]
if (!is_na(outliers_alert)) {
  cat("Outlier(s) detected\n")
  outliers_list <- unlist(presso_res[[1]][[2]][[3]][1], use.names = FALSE)
  outlier_df <- data.frame(SNP = Shared_beta_SE$SNP[outliers_list])
  write.table(outlier_df, paste0("outliers_", Pheno_id), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
} else {
  cat("No outliers detected\n")
}

# Save session
Rsesion_name <- paste0(Ph_temp[1], "_noOut_Rsession.rds")
RObjects_name <- paste0(Ph_temp[1], "_noOut_Robjects.rds")
save.image(Rsesion_name)

# Recalculate without outliers
if (!is_na(outliers_alert)) {
  out_temp <- sapply(outlier_df$SNP, function(x) grep(x, datos$SNP))
  datos_noOut <- datos[-out_temp, ]
} else {
  datos_noOut <- datos
}
res_noOut <- mr(datos_noOut, method_list = c("mr_ivw", "mr_egger_regression"))

# Sensitivity analyses
res_heterogeneity <- mr_heterogeneity(datos_noOut)
if (res_heterogeneity[2, 8] > 0.05) {
  IVW_fix <- mr_ivw_fe(datos_noOut$beta.exposure, datos_noOut$beta.outcome,
                       datos_noOut$se.exposure, datos_noOut$se.outcome)
  write.table(IVW_fix, paste0(Ph_temp[1], "_IVW_fixed_effects.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

res_pleiotropy <- mr_pleiotropy_test(datos_noOut)
res_single <- mr_singlesnp(datos_noOut)
res_loo <- mr_leaveoneout(datos_noOut)

# Save plots
ggsave(mr_scatter_plot(res_plot, datos_noOut)[[1]], filename = paste0(Pheno_id, "_scatter_plot.png"))
ggsave(mr_forest_plot(res_single)[[1]], filename = paste0(Pheno_id, "_forest_plot.png"))
ggsave(mr_leaveoneout_plot(res_loo)[[1]], filename = paste0(Pheno_id, "_leaveOneOut_plot.png"))
ggsave(mr_funnel_plot(res_single)[[1]], filename = paste0(Pheno_id, "_funnel_plot.png"))

# Save results
OR <- exp(resultado$b)
CIp <- exp(resultado$b + 1.96 * resultado$se)
CIn <- exp(resultado$b - 1.96 * resultado$se)
pval <- resultado$pval

press_vals <- if (!is_na(outliers_alert)) presso_res[[1]][2, ] else presso_res[[1]][1, ]
ORp <- exp(press_vals$Main.MR.results.Causal.Estimate)
CIpp <- exp(press_vals$Main.MR.results.Causal.Estimate + 1.96 * press_vals$Main.MR.results.Sd)
CIpn <- exp(press_vals$Main.MR.results.Causal.Estimate - 1.96 * press_vals$Main.MR.results.Sd)
pvalp <- press_vals$Main.MR.results.P.value

# Summary table
deftable <- data.frame(
  MR_approach = c("MR-Egger", "Random-effect_IVW", "Weight-Median", "MR-PRESSO"),
  OR = round(c(OR[1], OR[3], OR[2], ORp), 4),
  CIp = round(c(CIp[1], CIp[3], CIp[2], CIpp), 4),
  CIn = round(c(CIn[1], CIn[3], CIn[2], CIpn), 4),
  P = round(c(pval[1], pval[3], pval[2], pvalp), 4),
  P_Plei_heter = c(round(res_pleiotropy[1, 7], 4), round(res_heterogeneity[1, 8], 4), NA, NA),
  stringsAsFactors = FALSE
)

# Write outputs
write.table(resultado, paste0("mr_analysis_noOut_", Pheno_id), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(presso_res[[1]][1], paste0("presso_analysis_", Pheno_id), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deftable, paste0("table_MR_noOut_", Pheno_id), sep = "\t", quote = FALSE, row.names = FALSE)
save(datos, resultado, presso_res, file = RObjects_name)
