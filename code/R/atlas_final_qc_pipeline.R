# Sepsis.Atlas final cell-QC helpers.
# Input: one preQC Seurat object containing RAW (uncorrected called-cell counts)
#        and RNA (SoupX-corrected counts) assays.
# Order: scDblFinder on RAW -> miQC on RAW singlets -> conservative hard QC.

atlas_final_qc_pipeline_version <- "final-qc-v4-20260803"

atlas_env_flag <- function(name, default = FALSE) {
  x <- tolower(Sys.getenv(name, unset = if (default) "true" else "false"))
  x %in% c("1", "true", "yes", "y", "on")
}

atlas_env_num <- function(name, default, min = -Inf, max = Inf) {
  x <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (!is.finite(x) || x < min || x > max) x <- default
  x
}

atlas_env_int <- function(name, default, min = -Inf, max = Inf) {
  x <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(x) || x < min || x > max) x <- as.integer(default)
  x
}

atlas_require_namespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg, call. = FALSE)
  }
}

atlas_get_assay_counts <- function(obj, assay) {
  x <- tryCatch(
    SeuratObject::GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) SeuratObject::GetAssayData(obj, assay = assay, slot = "counts")
  )
  methods::as(x, "dgCMatrix")
}

atlas_seed_from_string <- function(x, offset = 0L) {
  z <- utf8ToInt(enc2utf8(as.character(x)))
  if (!length(z)) z <- 1L
  seed <- (sum(as.double(z) * seq_along(z)) + as.double(offset)) %% 2147483000
  as.integer(seed + 1L)
}


# flexmix implements logLik as an S4 method (stats4 generic), not the
# stats::logLik S3 generic. Prefer the stored convergence log-likelihood,
# with stats4 dispatch as a fallback.
atlas_flexmix_loglik <- function(model) {
  value <- tryCatch({
    if (methods::is(model, "flexmix") &&
        "logLik" %in% methods::slotNames(model)) {
      as.numeric(methods::slot(model, "logLik"))
    } else {
      as.numeric(stats4::logLik(model))
    }
  }, error = function(e) NA_real_)

  if (length(value) != 1L || !is.finite(value)) NA_real_ else value
}

atlas_pct <- function(counts, feature_mask, totals = NULL) {
  if (is.null(totals)) totals <- Matrix::colSums(counts)
  out <- rep(0, ncol(counts))
  names(out) <- colnames(counts)
  if (any(feature_mask)) {
    numerator <- Matrix::colSums(counts[feature_mask, , drop = FALSE])
    ok <- totals > 0
    out[ok] <- 100 * numerator[ok] / totals[ok]
  }
  out
}

atlas_compute_raw_qc <- function(counts) {
  totals <- Matrix::colSums(counts)
  detected <- Matrix::colSums(counts != 0)
  rn <- rownames(counts)

  mito <- grepl("^MT-", rn, ignore.case = TRUE)
  ribo <- grepl("^(RPS|RPL)[0-9]", rn, ignore.case = TRUE)
  hb <- grepl("^HB(A1|A2|B|D|E1|G1|G2|M|Q1|Z)$", rn, ignore.case = TRUE)

  data.frame(
    cell_barcode = colnames(counts),
    nCount_RAW = as.numeric(totals),
    nFeature_RAW = as.numeric(detected),
    percent.mt_RAW = as.numeric(atlas_pct(counts, mito, totals)),
    percent.ribo_RAW = as.numeric(atlas_pct(counts, ribo, totals)),
    percent.hb_RAW = as.numeric(atlas_pct(counts, hb, totals)),
    stringsAsFactors = FALSE,
    row.names = colnames(counts)
  )
}

atlas_mad_flag <- function(x, side = c("low", "high"), nmads = 5) {
  side <- match.arg(side)
  x <- as.numeric(x)
  ok <- is.finite(x)
  out <- rep(FALSE, length(x))
  if (sum(ok) < 10L) return(out)
  med <- stats::median(x[ok])
  spread <- stats::mad(x[ok], center = med, constant = 1.4826)
  if (!is.finite(spread) || spread <= 0) return(out)
  if (side == "low") out[ok] <- x[ok] < med - nmads * spread
  else out[ok] <- x[ok] > med + nmads * spread
  out
}

atlas_add_general_qc <- function(md) {
  nmads <- atlas_env_num("ATLAS_GENERAL_QC_NMADS", 5, min = 1, max = 20)
  md$qc_flag_low_count_outlier <- atlas_mad_flag(log10(md$nCount_RAW + 1), "low", nmads)
  md$qc_flag_high_count_outlier <- atlas_mad_flag(log10(md$nCount_RAW + 1), "high", nmads)
  md$qc_flag_low_feature_outlier <- atlas_mad_flag(log10(md$nFeature_RAW + 1), "low", nmads)
  md$qc_flag_high_feature_outlier <- atlas_mad_flag(log10(md$nFeature_RAW + 1), "high", nmads)
  md$qc_flag_high_mt_outlier <- atlas_mad_flag(md$percent.mt_RAW, "high", nmads)

  min_counts <- atlas_env_num("ATLAS_GENERAL_MIN_COUNTS", 1, min = 0)
  min_features <- atlas_env_num("ATLAS_GENERAL_MIN_FEATURES", 1, min = 0)
  max_counts <- atlas_env_num("ATLAS_GENERAL_MAX_COUNTS", .Machine$double.xmax, min = 1)
  max_features <- atlas_env_num("ATLAS_GENERAL_MAX_FEATURES", .Machine$double.xmax, min = 1)
  max_mt <- atlas_env_num("ATLAS_GENERAL_MAX_MT", 100, min = 0, max = 100)
  max_ribo <- atlas_env_num("ATLAS_GENERAL_MAX_RIBO", 100, min = 0, max = 100)
  max_hb <- atlas_env_num("ATLAS_GENERAL_MAX_HB", 100, min = 0, max = 100)

  md$general_qc_fail_low_counts <- !is.finite(md$nCount_RAW) | md$nCount_RAW < min_counts
  md$general_qc_fail_low_features <- !is.finite(md$nFeature_RAW) | md$nFeature_RAW < min_features
  md$general_qc_fail_high_counts <- is.finite(md$nCount_RAW) & md$nCount_RAW > max_counts
  md$general_qc_fail_high_features <- is.finite(md$nFeature_RAW) & md$nFeature_RAW > max_features
  md$general_qc_fail_high_mt <- !is.finite(md$percent.mt_RAW) | md$percent.mt_RAW > max_mt
  md$general_qc_fail_high_ribo <- !is.finite(md$percent.ribo_RAW) | md$percent.ribo_RAW > max_ribo
  md$general_qc_fail_high_hb <- !is.finite(md$percent.hb_RAW) | md$percent.hb_RAW > max_hb

  md$general_qc_keep <- !(
    md$general_qc_fail_low_counts |
      md$general_qc_fail_low_features |
      md$general_qc_fail_high_counts |
      md$general_qc_fail_high_features |
      md$general_qc_fail_high_mt |
      md$general_qc_fail_high_ribo |
      md$general_qc_fail_high_hb
  )
  md
}

atlas_run_scdblfinder <- function(counts, library_key) {
  atlas_require_namespace("SingleCellExperiment")
  atlas_require_namespace("scDblFinder")
  atlas_require_namespace("BiocParallel")
  atlas_require_namespace("S4Vectors")

  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  seed <- atlas_seed_from_string(library_key, 1000L)
  set.seed(seed)

  fml <- names(formals(scDblFinder::scDblFinder))
  args <- list(sce)
  if ("clusters" %in% fml) args$clusters <- FALSE
  if ("BPPARAM" %in% fml) args$BPPARAM <- BiocParallel::SerialParam()
  if ("verbose" %in% fml) args$verbose <- FALSE

  explicit_dbr <- suppressWarnings(as.numeric(Sys.getenv("ATLAS_SCDOUBLET_DBR", unset = "")))
  if (is.finite(explicit_dbr) && explicit_dbr > 0 && explicit_dbr < 1 && "dbr" %in% fml) {
    args$dbr <- explicit_dbr
  }
  dbr_per1k <- atlas_env_num("ATLAS_SCDOUBLET_DBR_PER1K", 0.008, min = 0.0001, max = 0.1)
  if ("dbr.per1k" %in% fml) args$dbr.per1k <- dbr_per1k
  dbr_sd <- atlas_env_num("ATLAS_SCDOUBLET_DBR_SD", 0.015, min = 0, max = 1)
  if ("dbr.sd" %in% fml) args$dbr.sd <- dbr_sd

  ans <- tryCatch(
    do.call(scDblFinder::scDblFinder, args),
    error = function(e) e
  )
  if (inherits(ans, "error")) {
    return(list(ok = FALSE, error = conditionMessage(ans), metadata = NULL, package_version = as.character(utils::packageVersion("scDblFinder"))))
  }

  cd <- as.data.frame(SummarizedExperiment::colData(ans), stringsAsFactors = FALSE)
  sc_cols <- grep("^scDblFinder\\.", names(cd), value = TRUE)
  if (!"scDblFinder.class" %in% sc_cols) {
    return(list(ok = FALSE, error = "scDblFinder.class was not produced", metadata = NULL, package_version = as.character(utils::packageVersion("scDblFinder"))))
  }
  cd <- cd[, sc_cols, drop = FALSE]
  rownames(cd) <- colnames(ans)

  list(
    ok = TRUE,
    error = NA_character_,
    metadata = cd,
    package_version = as.character(utils::packageVersion("scDblFinder")),
    seed = seed,
    dbr = if (is.finite(explicit_dbr)) explicit_dbr else NA_real_,
    dbr_per1k = dbr_per1k,
    dbr_sd = dbr_sd,
    stats = S4Vectors::metadata(ans)$scDblFinder.stats
  )
}

atlas_choose_best_miqc_model <- function(sce, library_key) {
  atlas_require_namespace("miQC")
  nstarts <- atlas_env_int("ATLAS_MIQC_NSTARTS", 3L, min = 1L, max = 20L)
  model_type <- Sys.getenv("ATLAS_MIQC_MODEL_TYPE", unset = "linear")
  fits <- vector("list", nstarts)
  errors <- character(nstarts)
  lls <- rep(-Inf, nstarts)

  for (i in seq_len(nstarts)) {
    set.seed(atlas_seed_from_string(library_key, 2000L + i))
    fit <- tryCatch(
      miQC::mixtureModel(
        sce,
        model_type = model_type,
        detected = "detected",
        subsets_mito_percent = "subsets_mito_percent"
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      errors[[i]] <- conditionMessage(fit)
    } else {
      fits[[i]] <- fit
      lls[[i]] <- atlas_flexmix_loglik(fit)
      if (!is.finite(lls[[i]])) lls[[i]] <- -Inf
    }
  }

  if (all(vapply(fits, is.null, logical(1)))) {
    return(list(ok = FALSE, error = paste(unique(errors[nzchar(errors)]), collapse = " | "), model = NULL, nstarts_success = 0L, model_type = model_type))
  }
  finite_ll <- which(is.finite(lls) & !vapply(fits, is.null, logical(1)))
  if (length(finite_ll)) {
    best <- finite_ll[[which.max(lls[finite_ll])]]
  } else {
    best <- which(!vapply(fits, is.null, logical(1)))[1]
  }
  list(
    ok = TRUE,
    error = NA_character_,
    model = fits[[best]],
    logLik = lls[[best]],
    nstarts_success = sum(!vapply(fits, is.null, logical(1))),
    model_type = model_type,
    best_start = best
  )
}

atlas_miqc_decision <- function(sce, model, posterior_cutoff = 0.90,
                                keep_all_below_boundary = TRUE,
                                enforce_left_cutoff = TRUE) {
  atlas_require_namespace("flexmix")
  metrics <- as.data.frame(SummarizedExperiment::colData(sce), stringsAsFactors = FALSE)
  p1 <- flexmix::parameters(model, component = 1)[1]
  p2 <- flexmix::parameters(model, component = 2)[1]
  if (!is.finite(p1) || !is.finite(p2)) stop("miQC model component intercepts are not finite")
  if (p1 > p2) {
    compromised <- 1L
    intact <- 2L
  } else {
    intact <- 1L
    compromised <- 2L
  }

  post <- flexmix::posterior(model)
  if (nrow(post) != nrow(metrics) || ncol(post) < 2L) stop("Unexpected miQC posterior dimensions")
  prob <- as.numeric(post[, compromised])

  keep_posterior <- is.finite(prob) & prob <= posterior_cutoff
  keep_below_boundary <- keep_posterior

  predictions <- tryCatch(stats::fitted(model), error = function(e) NULL)
  if (keep_all_below_boundary && !is.null(predictions) && nrow(predictions) == nrow(metrics) && ncol(predictions) >= intact) {
    intact_prediction <- as.numeric(predictions[, intact])
    rescue <- is.finite(metrics$subsets_mito_percent) &
      is.finite(intact_prediction) &
      metrics$subsets_mito_percent < intact_prediction
    keep_below_boundary[rescue] <- TRUE
  }

  keep_final <- keep_below_boundary
  if (enforce_left_cutoff && any(!keep_final)) {
    min_discard <- min(metrics$subsets_mito_percent[!keep_final], na.rm = TRUE)
    min_index <- which(metrics$subsets_mito_percent == min_discard)[1]
    lib_complexity <- metrics$detected[[min_index]]
    force_discard <- metrics$detected <= lib_complexity & metrics$subsets_mito_percent >= min_discard
    keep_final[force_discard] <- FALSE
  }

  data.frame(
    miQC.prob_compromised = prob,
    miQC.keep_posterior = keep_posterior,
    miQC.keep_below_boundary = keep_below_boundary,
    miQC.keep_model = keep_final,
    stringsAsFactors = FALSE,
    row.names = colnames(sce)
  )
}

atlas_run_miqc <- function(counts, singlet_barcodes, library_key) {
  atlas_require_namespace("SingleCellExperiment")
  atlas_require_namespace("scater")
  atlas_require_namespace("miQC")
  atlas_require_namespace("flexmix")

  min_cells <- atlas_env_int("ATLAS_MIQC_MIN_CELLS", 200L, min = 50L)
  if (length(singlet_barcodes) < min_cells) {
    return(list(ok = FALSE, error = paste0("miQC has only ", length(singlet_barcodes), " singlets; minimum=", min_cells)))
  }
  x <- counts[, singlet_barcodes, drop = FALSE]
  mt_genes <- rownames(x)[grepl("^MT-", rownames(x), ignore.case = TRUE)]
  if (!length(mt_genes)) return(list(ok = FALSE, error = "No mitochondrial genes matched ^MT-"))

  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = x))
  sce <- scater::addPerCellQC(sce, subsets = list(mito = mt_genes))
  cd <- as.data.frame(SummarizedExperiment::colData(sce), stringsAsFactors = FALSE)
  if (!all(c("detected", "subsets_mito_percent") %in% names(cd))) {
    return(list(ok = FALSE, error = "scater did not produce detected and subsets_mito_percent"))
  }
  if (length(unique(cd$detected[is.finite(cd$detected)])) < 3L ||
      length(unique(cd$subsets_mito_percent[is.finite(cd$subsets_mito_percent)])) < 3L) {
    return(list(ok = FALSE, error = "Insufficient variation for a two-component miQC model"))
  }

  fit <- atlas_choose_best_miqc_model(sce, library_key)
  if (!fit$ok) return(fit)

  cutoff <- atlas_env_num("ATLAS_MIQC_POSTERIOR_CUTOFF", 0.90, min = 0.01, max = 0.999)
  decision <- tryCatch(
    atlas_miqc_decision(
      sce,
      fit$model,
      posterior_cutoff = cutoff,
      keep_all_below_boundary = TRUE,
      enforce_left_cutoff = TRUE
    ),
    error = function(e) e
  )
  if (inherits(decision, "error")) {
    return(list(ok = FALSE, error = conditionMessage(decision)))
  }

  posterior_removed_fraction <- 1 - mean(decision$miQC.keep_posterior, na.rm = TRUE)
  boundary_removed_fraction <- 1 - mean(decision$miQC.keep_below_boundary, na.rm = TRUE)
  model_removed_fraction <- 1 - mean(decision$miQC.keep_model, na.rm = TRUE)

  max_auto_removed <- atlas_env_num(
    "ATLAS_MIQC_MAX_AUTO_REMOVED_FRACTION",
    0.50,
    min = 0,
    max = 1
  )
  auto_apply <- is.finite(model_removed_fraction) &&
    model_removed_fraction <= max_auto_removed
  auto_review_reason <- if (auto_apply) {
    "ok"
  } else if (!is.finite(model_removed_fraction)) {
    "miQC_removed_fraction_not_finite"
  } else {
    paste0(
      "miQC_model_removed_fraction_",
      format(round(model_removed_fraction, 4), nsmall = 4),
      "_exceeds_",
      format(round(max_auto_removed, 4), nsmall = 4)
    )
  }

  list(
    ok = TRUE,
    error = NA_character_,
    metadata = decision,
    package_version = as.character(utils::packageVersion("miQC")),
    posterior_cutoff = cutoff,
    model_type = fit$model_type,
    logLik = fit$logLik,
    nstarts_success = fit$nstarts_success,
    best_start = fit$best_start,
    posterior_removed_fraction = posterior_removed_fraction,
    boundary_removed_fraction = boundary_removed_fraction,
    model_removed_fraction = model_removed_fraction,
    max_auto_removed_fraction = max_auto_removed,
    auto_apply = auto_apply,
    auto_review_reason = auto_review_reason
  )
}

atlas_make_exclusion_reasons <- function(md) {
  reasons <- rep("", nrow(md))
  add <- function(mask, label) {
    idx <- which(!is.na(mask) & mask)
    if (!length(idx)) return(invisible(NULL))
    reasons[idx] <<- ifelse(nzchar(reasons[idx]), paste0(reasons[idx], ";", label), label)
  }
  add(!md$scDblFinder.keep, "scDblFinder_doublet")
  add(md$scDblFinder.keep & !md$miQC.keep, "miQC_compromised")
  add(md$general_qc_fail_low_counts, "general_low_counts")
  add(md$general_qc_fail_low_features, "general_low_features")
  add(md$general_qc_fail_high_counts, "general_high_counts")
  add(md$general_qc_fail_high_features, "general_high_features")
  add(md$general_qc_fail_high_mt, "general_high_mt")
  add(md$general_qc_fail_high_ribo, "general_high_ribo")
  add(md$general_qc_fail_high_hb, "general_high_hb")
  reasons[!nzchar(reasons)] <- "kept"
  reasons
}

atlas_write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

atlas_write_csv_gz <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.csv(df, con, row.names = FALSE, na = "")
}
