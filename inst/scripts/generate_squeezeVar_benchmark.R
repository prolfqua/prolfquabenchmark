#!/usr/bin/env Rscript
# Regenerates inst/extdata/squeezeVar_benchmark.rds, the cached results shown in
# the `SqueezeVar_comparison` vignette.
#
# It compares the two empirical-Bayes variance-moderation backends
#   - prolfqua::squeezeVarRob() : limma classic fitFDist + min_df (the former fork, now deprecated)
#   - limma::squeezeVar()       : the modern estimator prolfqua now uses
# holding everything else fixed (same models, contrasts, downstream stats), on:
#   1. synthetic data with known regulated proteins (30 reps x 2 noise scenarios)
#   2. the real IonStar spike-in dataset (prolfquadata, ECOLI = TP, HUMAN = TN)
# for three moderated facades spanning the df regimes: lm (integer df),
# rlm (fractional df) and lmer_nested (fractional Satterthwaite df).
#
# Runtime is dominated by fitting IonStar lmer_nested on all proteins
# (~10-30 min); that is why the results are cached rather than computed at
# vignette build time. Run from the prolfquabenchmark package root:
#   Rscript inst/scripts/generate_squeezeVar_benchmark.R

suppressMessages({
  library(prolfqua)
  library(prolfquapp)
  library(dplyr)
  library(pROC)
  library(limma)
  library(readxl)
})

OUT <- "inst/extdata/squeezeVar_benchmark.rds"
SEP <- "~"

# ---- shared moderation helper (mirrors prolfqua::moderated_p_limma) ----------
moderate_generic <- function(cd, squeeze_fun) {
  do.call(
    rbind,
    lapply(split(cd, cd$contrast), function(g) {
      g <- g[
        is.finite(g$sigma) &
          is.finite(g$df) &
          is.finite(g$statistic) &
          g$df > 0,
      ]
      if (nrow(g) < 3) {
        return(NULL)
      }
      sv <- squeeze_fun(g$sigma^2, df = g$df, robust = FALSE)
      dfprior <- sv$df.prior
      if (all(is.infinite(dfprior))) {
        dfprior <- mean(g$df) * nrow(g) / 10
      }
      modstat <- g$statistic * g$sigma / sqrt(sv$var.post)
      g$mod.p <- 2 * pt(abs(modstat), df = g$df + dfprior, lower.tail = FALSE)
      g$mod.q <- p.adjust(g$mod.p, method = "BH")
      g$collapsed <- all(is.infinite(sv$df.prior))
      g
    })
  )
}
raw_contrasts <- function(fa) {
  r <- fa$contrast$Contrast$get_contrasts(all = TRUE)
  r[, intersect(
    c("protein_Id", "contrast", "diff", "sigma", "df", "statistic"),
    colnames(r)
  )]
}
# The former prolfqua::squeezeVarRob() was removed in prolfqua 1.6.3. It was
# limma's classic estimator (legacy = TRUE) plus MSqRob's min_df pre-filter;
# reconstructed here (non-robust; production moderation always uses robust = FALSE).
former_fork <- function(v, df, min_df = 1) {
  n <- length(v)
  dfv <- if (length(df) == 1L) rep(df, n) else df
  keep <- is.finite(dfv) & dfv > min_df
  fit <- limma::squeezeVar(v[keep], df = dfv[keep], legacy = TRUE)
  var.post <- if (is.infinite(fit$df.prior)) {
    rep(fit$var.prior, n)
  } else {
    (dfv * v + fit$df.prior * fit$var.prior) / (dfv + fit$df.prior)
  }
  list(df.prior = fit$df.prior, var.prior = fit$var.prior, var.post = var.post)
}
BACKENDS <- list(
  `squeezeVarRob (former)` = function(v, df, robust) former_fork(v, df),
  `limma::squeezeVar (current)` = function(v, df, robust) {
    limma::squeezeVar(v, df = df, robust = robust)
  }
)
CT_SIM <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

# ---- 1. synthetic ------------------------------------------------------------
FC <- list(A = c(D = -12, U = 5, N = 0))
PROP <- list(A = c(D = 10, U = 10))
mkcfg <- function(pep) {
  cfg <- AnalysisConfiguration$new()
  cfg$file_name <- "sample"
  cfg$factors[["group_"]] <- "group"
  cfg$hierarchy[["protein_Id"]] <- c("proteinID", "idtype2")
  if (pep) {
    cfg$hierarchy[["peptide_Id"]] <- "peptideID"
  }
  cfg$set_response("abundance")
  cfg$nr_children <- "nr_peptides"
  cfg
}
sim_data <- function(N, seed, wm, pep) {
  set.seed(seed)
  raw <- sim_lfq_data(Nprot = 300, N = N, PEPTIDE = pep, fc = FC, prop = PROP)
  raw <- raw[raw$group %in% c("Ctrl", "A"), ]
  if (!pep) {
    raw$peptideID <- raw$proteinID
  }
  raw$nr_peptides <- as.numeric(
    !which_missing(raw$abundance, weight_missing = wm)
  )
  if (pep) {
    raw$nr_children <- raw$nr_peptides
  }
  raw <- raw[raw$nr_peptides > 0, ]
  raw$isotopeLabel <- "light"
  raw$qValue <- 0
  truth <- raw |>
    dplyr::filter(.data$group == "A") |>
    dplyr::transmute(
      protein_Id = paste(.data$proteinID, .data$idtype2, sep = SEP),
      is_true = as.integer(mean != 0)
    ) |>
    dplyr::distinct()
  lfq <- LFQData$new(
    setup_analysis(raw, mkcfg(pep)),
    mkcfg(pep)
  )$get_Transformer()$log2()$lfq
  list(lfq = lfq, truth = truth)
}
SIM_METHODS <- list(
  lm = list(pep = FALSE, f = "~ group_"),
  rlm = list(pep = FALSE, f = "~ group_"),
  lmer_nested = list(
    pep = TRUE,
    f = "~ group_ + (1|peptide_Id) + (1|sampleName)"
  )
)
SCEN <- list(moderate = list(N = 6, wm = 0.2), low_df = list(N = 4, wm = 0.5))
R <- 30
run_synthetic_backend <- function(cd, truth, scenario, method, backend) {
  mod <- moderate_generic(cd, function(v, df, robust) {
    BACKENDS[[backend]](v, df, robust)
  })
  if (is.null(mod)) {
    return(NULL)
  }
  d <- dplyr::inner_join(mod, truth, by = "protein_Id")
  d <- d[is.finite(d$mod.p), ]
  if (length(unique(d$is_true)) < 2) {
    return(NULL)
  }
  roc <- pROC::roc(
    d$is_true,
    -log10(pmax(d$mod.p, 1e-300)),
    quiet = TRUE,
    direction = "<"
  )
  called <- d$mod.q <= 0.05
  data.frame(
    scenario = scenario,
    method = method,
    backend = backend,
    AUC = as.numeric(pROC::auc(roc)),
    FDP_at_0.05 = if (sum(called)) mean(d$is_true[called] == 0) else NA,
    power = sum(called & d$is_true == 1) / sum(d$is_true == 1)
  )
}
run_synthetic_method <- function(scenario, params, scenario_index, replicate, method) {
  spec <- SIM_METHODS[[method]]
  simulated <- tryCatch(
    sim_data(params$N, 1000 * scenario_index + replicate, params$wm, spec$pep),
    error = function(e) NULL
  )
  if (is.null(simulated)) {
    return(NULL)
  }
  facade <- tryCatch(
    build_contrast_analysis(simulated$lfq, spec$f, CT_SIM, method = method),
    error = function(e) NULL
  )
  if (is.null(facade)) {
    return(NULL)
  }
  contrasts <- suppressMessages(raw_contrasts(facade))
  metrics <- lapply(
    names(BACKENDS),
    function(backend) {
      run_synthetic_backend(
        contrasts,
        simulated$truth,
        scenario,
        method,
        backend
      )
    }
  )
  list(
    metrics = Filter(Negate(is.null), metrics),
    df = data.frame(
      scenario = scenario,
      method = method,
      frac_df = mean(contrasts$df %% 1 != 0, na.rm = TRUE),
      median_df = median(contrasts$df, na.rm = TRUE)
    )
  )
}
run_synthetic_scenario <- function(scenario, params, scenario_index) {
  results <- unlist(
    lapply(
      seq_len(R),
      function(replicate) {
        lapply(
          names(SIM_METHODS),
          function(method) {
            run_synthetic_method(
              scenario,
              params,
              scenario_index,
              replicate,
              method
            )
          }
        )
      }
    ),
    recursive = FALSE
  )
  results <- Filter(Negate(is.null), results)
  message("synthetic scenario done: ", scenario)
  list(
    metrics = unlist(lapply(results, `[[`, "metrics"), recursive = FALSE),
    df = lapply(results, `[[`, "df")
  )
}
synthetic_results <- Map(
  run_synthetic_scenario,
  names(SCEN),
  SCEN,
  seq_along(SCEN)
)
syn <- unlist(lapply(synthetic_results, `[[`, "metrics"), recursive = FALSE)
syndf <- unlist(lapply(synthetic_results, `[[`, "df"), recursive = FALSE)
syn <- do.call(rbind, syn)
syndf <- do.call(rbind, syndf)
syn_metrics <- syn |>
  group_by(scenario, method, backend) |>
  summarise(
    AUC = mean(AUC),
    FDP_at_0.05 = mean(FDP_at_0.05, na.rm = TRUE),
    power = mean(power),
    n_reps = n() / 1,
    .groups = "drop"
  ) |>
  arrange(scenario, method, desc(backend)) |>
  as.data.frame()
syn_df <- syndf |>
  group_by(scenario, method) |>
  summarise(
    pct_fractional_df = round(100 * mean(frac_df)),
    median_df = round(mean(median_df), 1),
    .groups = "drop"
  ) |>
  as.data.frame()

# ---- 2. real IonStar ---------------------------------------------------------
dd <- file.path(find.package("prolfquadata"), "quantdata")
config <- AnalysisConfiguration$new()
config$file_name <- "raw.file"
config$hierarchy[["protein_Id"]] <- "leading.razor.protein"
config$hierarchy[["peptide_Id"]] <- "sequence"
config$hierarchy_depth <- 1
config$ident_q_value <- "pep"
config$set_response("peptide.intensity")
config$isotope_label <- "isotope"
config$min_peptides_protein <- 2
data_mq <- prolfquapp::tidyMQ_Peptides(file.path(
  dd,
  "MAXQuant_IonStar2018_PXD003881.zip"
))
annotation <- readxl::read_xlsx(file.path(
  dd,
  "annotation_Ionstar2018_PXD003881.xlsx"
))
res <- dplyr::inner_join(data_mq, annotation, by = "raw.file")
config$factors[["dilution."]] <- "sample"
config$factors[["run_Id"]] <- "run_ID"
config$factor_depth <- 1
lfq <- LFQData$new(setup_analysis(res, config), config)
lfq$set_data(
  lfq$data_long() |> dplyr::filter(!grepl("^REV__|^CON__", protein_Id))
)
lfq$filter_proteins_by_peptide_count()
lfq$remove_small_intensities()
tr <- lfq$get_Transformer()
subset_h <- lfq$get_copy()
subset_h$set_data(
  subset_h$data_long() |> dplyr::filter(grepl("HUMAN", protein_Id))
)
subset_h <- subset_h$get_Transformer()$log2()$lfq
lfqNorm <- tr$log2()$robscale_subset(lfqsubset = subset_h)$lfq
agg <- lfqNorm$get_Aggregator()
agg$aggregate()
lfqProt <- agg$lfq_agg
truth <- data.frame(protein_Id = unique(lfqProt$data_long()$protein_Id)) |>
  mutate(
    species = ifelse(
      grepl("ECOLI", protein_Id),
      "ECOLI",
      ifelse(grepl("HUMAN", protein_Id), "HUMAN", "OTHER")
    )
  ) |>
  filter(species != "OTHER") |>
  mutate(is_true = as.integer(species == "ECOLI"))
CT <- c(
  "e_vs_d_1.2" = "dilution.e - dilution.d",
  "d_vs_c_1.25" = "dilution.d - dilution.c",
  "c_vs_b_1.33" = "dilution.c - dilution.b",
  "b_vs_a_1.5" = "dilution.b - dilution.a"
)
ION_METHODS <- list(
  lm = list(src = "prot", f = "~ dilution."),
  rlm = list(src = "prot", f = "~ dilution."),
  lmer_nested = list(
    src = "pep",
    f = "~ dilution. + (1|peptide_Id) + (1|sampleName)"
  )
)
ion <- list()
iondf <- list()
ioncon <- list()
for (m in names(ION_METHODS)) {
  spec <- ION_METHODS[[m]]
  src <- if (spec$src == "prot") lfqProt else lfqNorm
  fa <- tryCatch(
    build_contrast_analysis(src, spec$f, CT, method = m),
    error = function(e) NULL
  )
  if (is.null(fa)) {
    next
  }
  cd <- suppressMessages(raw_contrasts(fa))
  iondf[[m]] <- data.frame(
    method = m,
    pct_fractional_df = round(100 * mean(cd$df %% 1 != 0, na.rm = TRUE)),
    median_df = round(median(cd$df, na.rm = TRUE), 1),
    min_df = round(min(cd$df, na.rm = TRUE), 2)
  )
  mods <- list()
  for (b in names(BACKENDS)) {
    mod <- moderate_generic(cd, function(v, df, robust) {
      BACKENDS[[b]](v, df, robust)
    })
    if (is.null(mod)) {
      next
    }
    mods[[b]] <- mod
    d <- dplyr::inner_join(mod, truth, by = "protein_Id")
    d <- d[is.finite(d$mod.p), ]
    roc <- pROC::roc(
      d$is_true,
      -log10(pmax(d$mod.p, 1e-300)),
      quiet = TRUE,
      direction = "<"
    )
    called <- d$mod.q <= 0.05
    ion[[paste(m, b)]] <- data.frame(
      method = m,
      backend = b,
      AUC = round(as.numeric(pROC::auc(roc)), 3),
      pAUC10 = round(
        as.numeric(pROC::auc(
          roc,
          partial.auc = c(1, 0.9),
          partial.auc.focus = "specificity",
          partial.auc.correct = TRUE
        )),
        3
      ),
      FDP_at_0.05 = round(
        if (sum(called)) mean(d$is_true[called] == 0) else NA,
        3
      ),
      power = round(sum(called & d$is_true == 1) / sum(d$is_true == 1), 3),
      n_called = sum(called)
    )
  }
  if (length(mods) == 2) {
    j <- dplyr::inner_join(
      mods[[1]] |> transmute(protein_Id, contrast, q1 = mod.q),
      mods[[2]] |> transmute(protein_Id, contrast, q2 = mod.q),
      by = c("protein_Id", "contrast")
    )
    ioncon[[m]] <- data.frame(
      method = m,
      n_contrasts = nrow(j),
      pct_calls_concordant = round(
        100 * mean((j$q1 <= 0.05) == (j$q2 <= 0.05)),
        2
      ),
      calls_changed = sum((j$q1 <= 0.05) != (j$q2 <= 0.05)),
      median_abs_dq = round(median(abs(j$q1 - j$q2)), 4)
    )
  }
  message("ionstar method done: ", m)
}

cache <- list(
  synthetic = list(metrics = syn_metrics, df_regime = syn_df, n_reps = R),
  ionstar = list(
    metrics = do.call(rbind, ion),
    concordance = do.call(rbind, ioncon),
    df_regime = do.call(rbind, iondf),
    truth_n = as.integer(table(truth$species)),
    truth_labels = names(table(truth$species))
  ),
  generated = as.character(Sys.Date())
)
saveRDS(cache, OUT)
message("Wrote ", OUT)
