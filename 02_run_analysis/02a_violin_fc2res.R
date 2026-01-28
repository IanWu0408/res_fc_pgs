

# plot_config.R
# ------------------------------------------------------------
# Config for violin plots only
# ------------------------------------------------------------
PLOT_CFG <- list(
  # Input: your pipeline output summary
  # must contain columns: model, fold, R2_test_full, dR2_test
  fold_metrics_csv = "outputs/RIDGE_GLMNET/RUN_YYYYMMDD_HHMMSS_ROI_v2/summary/all_models_fold_metrics.csv",
  
  # Output folder for plots
  out_dir = "outputs/violin_plots",
  
  # Which models to include (NULL = include all found in csv)
  include_models = NULL,
  
  # Explicit order on x-axis (NULL = use include_models order, else alphabetical)
  model_order = c(
    "cov_only",
    "cov_triple",
    "cov_appraisals",
    "cov_coping",
    "cov_appraisals_coping",
    "cov_triple_appraisals_coping"
    # "cov_allFC"
  ),
  
  # Plot switches
  make_r2_plot  = TRUE,   # violin of R2_test_full across folds
  make_dr2_plot = TRUE,   # violin of dR2_test across folds (exclude cov_only)
  
  # Significance annotations
  do_sig = TRUE,
  sig_method = "t",       # "t" or "wilcox"
  sig_adjust = "holm",    # holm / bonferroni / BH / ...
  
  # Cosmetic (no hardcoded colors)
  jitter_width = 0.12,
  violin_alpha = 0.45,
  point_alpha  = 0.65,
  
  # File names
  r2_png  = "violin_testR2_full_meanSE_sig.png",
  r2_pdf  = "violin_testR2_full_meanSE_sig.pdf",
  dr2_png = "violin_deltaR2_meanSE_sig.png",
  dr2_pdf = "violin_deltaR2_meanSE_sig.pdf",
  
  # Plot titles
  r2_title = "10-fold Test R² distribution by model (Full model)",
  r2_subtitle = "Fold-wise test R² | Mean ± SE shown | Paired tests across folds (Holm-adjusted)",
  dr2_title = "10-fold Test ΔR² distribution by model",
  dr2_subtitle = "ΔR² = R²(full) - R²(cov_only) | Mean ± SE shown | Paired tests across folds (Holm-adjusted)"
)
