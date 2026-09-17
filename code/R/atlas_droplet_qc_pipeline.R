# Shared emptyDrops -> valiDrops -> SoupX preprocessing for Sepsis.Atlas.
# The full unfiltered 10x matrix is used for cell calling and ambient RNA estimation.

atlas_droplet_qc_pipeline_version <- "resource-safe-qc-v6-20260802"

atlas_env_flag <- function(name, default = FALSE) {
  x <- tolower(Sys.getenv(name, unset = if (default) "true" else "false"))
  x %in% c("1", "true", "yes", "y", "on")
}

atlas_require_namespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg, call. = FALSE)
  }
}

atlas_get_assay_data <- function(obj, assay = "RNA", layer = "counts") {
  tryCatch(
    SeuratObject::GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) SeuratObject::GetAssayData(obj, assay = assay, slot = layer)
  )
}

atlas_read_10x_gene_expression <- function(raw_dir) {
  x <- Seurat::Read10X(raw_dir, gene.column = 2, unique.features = TRUE)
  if (is.list(x)) {
    if ("Gene Expression" %in% names(x)) x <- x[["Gene Expression"]] else x <- x[[1]]
  }
  x <- methods::as(x, "dgCMatrix")
  totals <- Matrix::colSums(x)
  x[, totals > 0, drop = FALSE]
}


atlas_read_10x_barcode_file <- function(matrix_dir) {
  candidates <- c(
    file.path(matrix_dir, "barcodes.tsv.gz"),
    file.path(matrix_dir, "barcodes.tsv")
  )
  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) return(character())

  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  x <- readLines(con, warn = FALSE)
  x <- sub("\\r$", "", x)
  unique(x[nzchar(x)])
}

atlas_get_cellranger_filtered_barcodes <- function(raw_dir, raw_counts) {
  outs_dir <- dirname(normalizePath(raw_dir, winslash = "/", mustWork = TRUE))
  candidates <- c(
    file.path(outs_dir, "filtered_feature_bc_matrix"),
    file.path(outs_dir, "filtered_gene_bc_matrices")
  )
  filtered_dir <- candidates[dir.exists(candidates)][1]
  if (is.na(filtered_dir)) {
    return(list(
      barcodes = character(),
      source = "not_available",
      filtered_dir = NA_character_,
      n_in_file = 0L,
      n_matched_raw = 0L
    ))
  }

  barcodes <- atlas_read_10x_barcode_file(filtered_dir)
  n_in_file <- length(barcodes)
  barcodes <- colnames(raw_counts)[colnames(raw_counts) %in% barcodes]

  list(
    barcodes = barcodes,
    source = "cellranger_filtered_feature_bc_matrix",
    filtered_dir = normalizePath(filtered_dir, winslash = "/", mustWork = TRUE),
    n_in_file = as.integer(n_in_file),
    n_matched_raw = as.integer(length(barcodes))
  )
}

atlas_get_barcode_rank_fallback <- function(raw_counts, emptydrops, min_candidates = 500L) {
  atlas_require_namespace("DropletUtils")
  atlas_require_namespace("S4Vectors")
  br <- DropletUtils::barcodeRanks(raw_counts)
  br_df <- as.data.frame(br)
  br_df$barcode <- rownames(br_df)
  meta <- S4Vectors::metadata(br)

  thresholds <- suppressWarnings(as.numeric(c(meta$inflection, meta$knee)))
  thresholds <- thresholds[is.finite(thresholds) & thresholds > 0]
  threshold <- if (length(thresholds)) min(thresholds) else NA_real_

  candidates <- character()
  if (is.finite(threshold) && "total" %in% names(br_df)) {
    candidates <- br_df$barcode[is.finite(br_df$total) & br_df$total >= threshold]
  }

  ed_barcodes <- emptydrops$barcode[emptydrops$emptydrops_pass %in% TRUE]
  candidates <- union(candidates, ed_barcodes)

  min_candidates <- suppressWarnings(as.integer(min_candidates))
  if (!is.finite(min_candidates) || min_candidates < 100L) min_candidates <- 500L

  if (length(candidates) < min_candidates) {
    totals <- Matrix::colSums(raw_counts)
    top_n <- min(min_candidates, length(totals))
    top_barcodes <- names(sort(totals, decreasing = TRUE))[seq_len(top_n)]
    candidates <- union(candidates, top_barcodes)
  }

  candidates <- colnames(raw_counts)[colnames(raw_counts) %in% candidates]
  list(
    barcodes = candidates,
    source = "barcodeRanks_plus_emptyDrops",
    filtered_dir = NA_character_,
    n_in_file = NA_integer_,
    n_matched_raw = as.integer(length(candidates)),
    barcode_rank_threshold = threshold
  )
}

atlas_get_validrops_fallback_candidates <- function(raw_dir, raw_counts, emptydrops) {
  min_candidates <- suppressWarnings(as.integer(Sys.getenv(
    "ATLAS_VALIDROPS_FALLBACK_MIN_BARCODES",
    unset = "500"
  )))
  if (!is.finite(min_candidates) || min_candidates < 100L) min_candidates <- 500L

  cr <- atlas_get_cellranger_filtered_barcodes(raw_dir, raw_counts)
  if (length(cr$barcodes) >= min_candidates) {
    cr$barcode_rank_threshold <- NA_real_
    return(cr)
  }

  br <- atlas_get_barcode_rank_fallback(
    raw_counts,
    emptydrops,
    min_candidates = min_candidates
  )
  if (length(cr$barcodes) > 0L) {
    br$barcodes <- union(cr$barcodes, br$barcodes)
    br$barcodes <- colnames(raw_counts)[colnames(raw_counts) %in% br$barcodes]
    br$source <- "cellranger_filtered_too_small_plus_barcodeRanks"
    br$filtered_dir <- cr$filtered_dir
    br$n_in_file <- cr$n_in_file
    br$n_matched_raw <- as.integer(length(br$barcodes))
  }
  br
}

atlas_run_emptydrops <- function(raw_counts, lower = 100L, fdr = 0.01) {
  atlas_require_namespace("DropletUtils")
  ed <- DropletUtils::emptyDrops(raw_counts, lower = as.integer(lower))
  out <- as.data.frame(ed)
  out$barcode <- rownames(out)
  names(out) <- paste0("emptydrops_", tolower(names(out)))
  names(out)[names(out) == "emptydrops_barcode"] <- "barcode"
  out$emptydrops_pass <- !is.na(out$emptydrops_fdr) & out$emptydrops_fdr <= fdr
  out
}

atlas_validrops_bpparam <- function() {
  atlas_require_namespace("BiocParallel")
  workers <- suppressWarnings(as.integer(Sys.getenv("ATLAS_VALIDROPS_WORKERS", unset = "1")))
  if (is.na(workers) || workers < 1L) workers <- 1L
  if (workers == 1L) return(BiocParallel::SerialParam())
  if (.Platform$OS.type == "windows") BiocParallel::SnowParam(workers) else BiocParallel::MulticoreParam(workers)
}

atlas_run_validrops_once <- function(
    counts,
    label_dead,
    rank_barcodes = TRUE,
    alpha = NULL,
    stage_three = TRUE
) {
  args <- list(
    counts = counts,
    label_dead = label_dead,
    rank_barcodes = rank_barcodes,
    stageThree = stage_three,
    bpparam = atlas_validrops_bpparam()
  )
  if (!is.null(alpha)) args$alpha <- alpha

  captured_warnings <- character()
  value <- tryCatch(
    withCallingHandlers(
      do.call(valiDrops::valiDrops, args),
      warning = function(w) {
        captured_warnings <<- c(captured_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  if (inherits(value, "error")) {
    return(list(
      ok = FALSE,
      value = NULL,
      error = conditionMessage(value),
      warnings = unique(captured_warnings)
    ))
  }

  value <- as.data.frame(value, stringsAsFactors = FALSE)
  if (!"barcode" %in% names(value)) {
    return(list(
      ok = FALSE,
      value = NULL,
      error = "valiDrops output has no barcode column",
      warnings = unique(captured_warnings)
    ))
  }

  list(
    ok = TRUE,
    value = value,
    error = NA_character_,
    warnings = unique(captured_warnings)
  )
}

atlas_parse_validrops_retry_alphas <- function() {
  raw <- Sys.getenv(
    "ATLAS_VALIDROPS_RETRY_ALPHAS",
    unset = "0.001,0.01"
  )
  vals <- suppressWarnings(as.numeric(trimws(strsplit(raw, ",", fixed = TRUE)[[1]])))
  vals <- unique(vals[is.finite(vals) & vals > 0 & vals < 0.5])
  if (length(vals) == 0L) vals <- c(0.001, 0.01)
  vals
}

atlas_run_validrops <- function(raw_counts, emptydrops, raw_dir, label_dead = TRUE) {
  atlas_require_namespace("valiDrops")
  atlas_require_namespace("BiocParallel")

  attempts <- list()
  alphas <- atlas_parse_validrops_retry_alphas()

  full_max <- suppressWarnings(as.integer(Sys.getenv(
    "ATLAS_VALIDROPS_FULL_MAX_BARCODES", unset = "1200000"
  )))
  if (!is.finite(full_max) || full_max < 10000L) full_max <- 1200000L
  skip_full <- ncol(raw_counts) > full_max

  if (skip_full) {
    message(
      "  valiDrops full-matrix rank step skipped: ", ncol(raw_counts),
      " non-zero barcodes exceeds ATLAS_VALIDROPS_FULL_MAX_BARCODES=", full_max
    )
    attempts[[1L]] <- data.frame(
      strategy = "rank_barcodes",
      alpha = NA_real_,
      stage_three = TRUE,
      n_input_barcodes = ncol(raw_counts),
      success = FALSE,
      error = "skipped_full_matrix_due_to_barcode_count",
      warnings = NA_character_,
      stringsAsFactors = FALSE
    )
  } else {
    for (alpha in alphas) {
      message("  valiDrops rank_barcodes retry alpha=", format(alpha, scientific = FALSE))
      ans <- atlas_run_validrops_once(
        counts = raw_counts,
        label_dead = label_dead,
        rank_barcodes = TRUE,
        alpha = alpha,
        stage_three = TRUE
      )

      attempts[[length(attempts) + 1L]] <- data.frame(
        strategy = "rank_barcodes",
        alpha = alpha,
        stage_three = TRUE,
        n_input_barcodes = ncol(raw_counts),
        success = ans$ok,
        error = if (ans$ok) NA_character_ else ans$error,
        warnings = paste(ans$warnings, collapse = " | "),
        stringsAsFactors = FALSE
      )

      if (ans$ok) {
        vd <- ans$value
        names(vd)[names(vd) != "barcode"] <- paste0(
          "validrops_",
          names(vd)[names(vd) != "barcode"]
        )
        return(list(
          table = vd,
          status = if (identical(alpha, alphas[[1]])) "ok" else "retry_success",
          strategy = "rank_barcodes",
          alpha = alpha,
          stage_three = TRUE,
          n_input_barcodes = ncol(raw_counts),
          fallback_source = NA_character_,
          fallback_filtered_dir = NA_character_,
          n_cellranger_filtered_barcodes = NA_integer_,
          final_error = NA_character_,
          attempts = do.call(rbind, attempts)
        ))
      }

      message("  valiDrops attempt failed: ", ans$error)
      invisible(gc())
    }
  }

  allow_fallback <- atlas_env_flag(
    "ATLAS_VALIDROPS_ALLOW_EMPTYDROPS_FALLBACK",
    TRUE
  )
  if (!allow_fallback) {
    stop(
      "valiDrops rank_barcodes failed for every configured alpha. ",
      "Set ATLAS_VALIDROPS_ALLOW_EMPTYDROPS_FALLBACK=true to evaluate ",
      "Cell Ranger filtered or barcode-rank candidates with rank_barcodes=FALSE. Last error: ",
      attempts[[length(attempts)]]$error,
      call. = FALSE
    )
  }

  fallback <- atlas_get_validrops_fallback_candidates(
    raw_dir = raw_dir,
    raw_counts = raw_counts,
    emptydrops = emptydrops
  )
  candidate_barcodes <- fallback$barcodes

  min_candidates <- suppressWarnings(as.integer(Sys.getenv(
    "ATLAS_VALIDROPS_FALLBACK_MIN_BARCODES",
    unset = "500"
  )))
  if (!is.finite(min_candidates) || min_candidates < 100L) min_candidates <- 500L
  if (length(candidate_barcodes) < min_candidates) {
    stop(
      "valiDrops rank_barcodes failed and fallback source '", fallback$source,
      "' produced only ", length(candidate_barcodes), " candidate barcodes; ",
      "minimum required is ", min_candidates, ".",
      call. = FALSE
    )
  }

  max_candidates <- suppressWarnings(as.integer(Sys.getenv(
    "ATLAS_VALIDROPS_FALLBACK_MAX_BARCODES",
    unset = "50000"
  )))
  if (!is.finite(max_candidates) || max_candidates < 1000L) max_candidates <- 50000L
  if (length(candidate_barcodes) > max_candidates) {
    stop(
      "valiDrops fallback would evaluate ", length(candidate_barcodes),
      " barcodes from ", fallback$source,
      ", exceeding ATLAS_VALIDROPS_FALLBACK_MAX_BARCODES=",
      max_candidates, ". Review this library manually.",
      call. = FALSE
    )
  }

  message(
    "  valiDrops rank step failed; retrying on ",
    length(candidate_barcodes), " candidates from ", fallback$source,
    " with rank_barcodes=FALSE"
  )
  if (!is.na(fallback$filtered_dir)) {
    message("  Cell Ranger filtered matrix: ", fallback$filtered_dir)
  }
  candidate_counts <- raw_counts[, candidate_barcodes, drop = FALSE]

  ans <- atlas_run_validrops_once(
    counts = candidate_counts,
    label_dead = label_dead,
    rank_barcodes = FALSE,
    alpha = NULL,
    stage_three = TRUE
  )
  attempts[[length(attempts) + 1L]] <- data.frame(
    strategy = paste0(fallback$source, "_rank_barcodes_FALSE"),
    alpha = NA_real_,
    stage_three = TRUE,
    n_input_barcodes = ncol(candidate_counts),
    success = ans$ok,
    error = if (ans$ok) NA_character_ else ans$error,
    warnings = paste(ans$warnings, collapse = " | "),
    stringsAsFactors = FALSE
  )

  stage_three_used <- TRUE
  if (!ans$ok) {
    message(
      "  valiDrops expression stage also failed; final retry with stageThree=FALSE: ",
      ans$error
    )
    ans <- atlas_run_validrops_once(
      counts = candidate_counts,
      label_dead = label_dead,
      rank_barcodes = FALSE,
      alpha = NULL,
      stage_three = FALSE
    )
    stage_three_used <- FALSE
    attempts[[length(attempts) + 1L]] <- data.frame(
      strategy = paste0(fallback$source, "_rank_barcodes_FALSE"),
      alpha = NA_real_,
      stage_three = FALSE,
      n_input_barcodes = ncol(candidate_counts),
      success = ans$ok,
      error = if (ans$ok) NA_character_ else ans$error,
      warnings = paste(ans$warnings, collapse = " | "),
      stringsAsFactors = FALSE
    )
  }

  if (!ans$ok) {
    stop(
      "valiDrops failed after rank retries and filtered/barcode-rank fallback. ",
      "Last error: ", ans$error,
      call. = FALSE
    )
  }

  vd <- ans$value
  names(vd)[names(vd) != "barcode"] <- paste0(
    "validrops_",
    names(vd)[names(vd) != "barcode"]
  )

  list(
    table = vd,
    status = if (stage_three_used) "fallback_success" else "fallback_stage3_skipped",
    strategy = paste0(fallback$source, "_rank_barcodes_FALSE"),
    alpha = NA_real_,
    stage_three = stage_three_used,
    n_input_barcodes = ncol(candidate_counts),
    fallback_source = fallback$source,
    fallback_filtered_dir = fallback$filtered_dir,
    n_cellranger_filtered_barcodes = fallback$n_matched_raw,
    final_error = paste(
      unique(stats::na.omit(vapply(attempts, function(x) x$error[[1]], character(1)))),
      collapse = " | "
    ),
    attempts = do.call(rbind, attempts)
  )
}

atlas_choose_called_barcodes <- function(raw_counts, ed, vd, policy = "validrops", drop_dead = FALSE) {
  vd_pass <- vd$validrops_qc.pass == "pass"
  if (drop_dead && "validrops_label" %in% names(vd)) vd_pass <- vd_pass & vd$validrops_label == "live"
  vd_barcodes <- vd$barcode[!is.na(vd_pass) & vd_pass]
  ed_barcodes <- ed$barcode[ed$emptydrops_pass %in% TRUE]

  called <- switch(
    policy,
    validrops = vd_barcodes,
    intersection = intersect(vd_barcodes, ed_barcodes),
    union = union(vd_barcodes, ed_barcodes),
    stop("Unknown ATLAS_CELL_CALL_POLICY: ", policy)
  )
  called <- colnames(raw_counts)[colnames(raw_counts) %in% called]
  called
}

atlas_precluster_for_soupx <- function(called_counts, sample_id) {
  atlas_require_namespace("Seurat")
  obj <- Seurat::CreateSeuratObject(called_counts, project = sample_id, min.cells = 0, min.features = 0)
  if (ncol(obj) < 50L) return(stats::setNames(rep("all", ncol(obj)), colnames(obj)))

  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  nfeat <- min(2000L, max(200L, nrow(obj) - 1L))
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeat, verbose = FALSE)
  hv <- Seurat::VariableFeatures(obj)
  if (length(hv) < 50L) return(stats::setNames(rep("all", ncol(obj)), colnames(obj)))

  obj <- Seurat::ScaleData(obj, features = hv, verbose = FALSE)
  npcs <- min(30L, length(hv) - 1L, ncol(obj) - 1L)
  if (npcs < 2L) return(stats::setNames(rep("all", ncol(obj)), colnames(obj)))
  obj <- Seurat::RunPCA(obj, features = hv, npcs = npcs, verbose = FALSE)
  dims <- seq_len(min(20L, npcs))
  obj <- Seurat::FindNeighbors(obj, dims = dims, verbose = FALSE)
  obj <- Seurat::FindClusters(
    obj,
    resolution = suppressWarnings(as.numeric(Sys.getenv("ATLAS_SOUX_CLUSTER_RESOLUTION", unset = "0.5"))),
    verbose = FALSE
  )
  stats::setNames(as.character(Seurat::Idents(obj)), colnames(obj))
}

atlas_prepare_soupx_tod <- function(raw_counts, called_counts, sample_id) {
  max_tod <- suppressWarnings(as.integer(Sys.getenv(
    "ATLAS_SOUX_MAX_TOD_BARCODES", unset = "100000"
  )))
  if (!is.finite(max_tod) || max_tod < 10000L) max_tod <- 100000L

  ambient_min <- suppressWarnings(as.integer(Sys.getenv(
    "ATLAS_SOUX_AMBIENT_UMI_MIN", unset = "1"
  )))
  ambient_max <- suppressWarnings(as.integer(Sys.getenv(
    "ATLAS_SOUX_AMBIENT_UMI_MAX", unset = "100"
  )))
  if (!is.finite(ambient_min) || ambient_min < 1L) ambient_min <- 1L
  if (!is.finite(ambient_max) || ambient_max < ambient_min) ambient_max <- 100L

  n_input <- ncol(raw_counts)
  called <- colnames(called_counts)
  target <- max(max_tod, length(called) + 10000L)
  if (n_input <= target) {
    return(list(
      tod = raw_counts,
      n_input = n_input,
      n_used = n_input,
      n_ambient_used = max(0L, n_input - length(called)),
      reduced = FALSE
    ))
  }

  totals <- Matrix::colSums(raw_counts)
  ambient <- names(totals)[
    totals >= ambient_min & totals <= ambient_max &
      !(names(totals) %in% called)
  ]
  n_keep <- max(0L, target - length(called))

  if (length(ambient) > n_keep) {
    code <- utf8ToInt(as.character(sample_id))
    seed <- as.integer((sum(code * seq_along(code)) %% 2147483000L) + 1L)
    old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    ambient <- sample(ambient, n_keep, replace = FALSE)
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }

  keep <- unique(c(called, ambient))
  keep <- colnames(raw_counts)[colnames(raw_counts) %in% keep]
  message(
    sample_id, ": SoupX tod reduced from ", n_input, " to ", length(keep),
    " barcodes (called=", length(called), ", ambient=", length(setdiff(keep, called)), ")"
  )
  list(
    tod = raw_counts[, keep, drop = FALSE],
    n_input = n_input,
    n_used = length(keep),
    n_ambient_used = length(setdiff(keep, called)),
    reduced = TRUE
  )
}

atlas_run_soupx <- function(raw_counts, called_counts, sample_id) {
  atlas_require_namespace("SoupX")
  clusters <- atlas_precluster_for_soupx(called_counts, sample_id)
  tod_info <- atlas_prepare_soupx_tod(raw_counts, called_counts, sample_id)
  sc <- SoupX::SoupChannel(tod = tod_info$tod, toc = called_counts, calcSoupProfile = TRUE)
  sc <- SoupX::setClusters(sc, clusters)

  method <- "autoEstCont"
  status <- "ok"
  auto_error <- NA_character_
  sc <- tryCatch(
    SoupX::autoEstCont(sc, doPlot = FALSE),
    error = function(e) {
      auto_error <<- conditionMessage(e)
      NULL
    }
  )

  if (is.null(sc)) {
    allow <- atlas_env_flag("ATLAS_SOUX_ALLOW_MANUAL_FALLBACK", FALSE)
    if (!allow) {
      stop("SoupX autoEstCont failed for ", sample_id, ": ", auto_error,
           ". Review the sample or set ATLAS_SOUX_ALLOW_MANUAL_FALLBACK=true explicitly.")
    }
    rho <- suppressWarnings(as.numeric(Sys.getenv("ATLAS_SOUX_FALLBACK_RHO", unset = "0.10")))
    if (!is.finite(rho) || rho <= 0 || rho >= 1) stop("Invalid ATLAS_SOUX_FALLBACK_RHO")
    sc <- SoupX::SoupChannel(tod = tod_info$tod, toc = called_counts, calcSoupProfile = TRUE)
    sc <- SoupX::setClusters(sc, clusters)
    sc <- SoupX::setContaminationFraction(sc, rho)
    method <- "manual_fallback"
    status <- paste0("fallback_after_autoEstCont_error: ", auto_error)
  }

  corrected <- SoupX::adjustCounts(sc, roundToInt = TRUE)
  corrected <- methods::as(corrected, "dgCMatrix")
  corrected <- corrected[rownames(called_counts), colnames(called_counts), drop = FALSE]

  rho_by_cell <- rep(NA_real_, ncol(called_counts))
  names(rho_by_cell) <- colnames(called_counts)
  if (!is.null(sc$metaData) && "rho" %in% names(sc$metaData)) {
    z <- sc$metaData$rho
    names(z) <- rownames(sc$metaData)
    rho_by_cell[names(z)] <- z
  }

  list(
    corrected_counts = corrected,
    clusters = clusters,
    rho_by_cell = rho_by_cell,
    method = method,
    status = status,
    auto_error = auto_error,
    tod_input_barcodes = tod_info$n_input,
    tod_used_barcodes = tod_info$n_used,
    ambient_barcodes_used = tod_info$n_ambient_used,
    tod_reduced = tod_info$reduced
  )
}

atlas_safe_write_csv_gz <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(df, path)
}

atlas_new_sample_exclusion <- function(
    sample_id,
    reason,
    n_called_after_droplet_qc,
    minimum_called_cells,
    validrops_strategy = NA_character_,
    validrops_fallback_source = NA_character_,
    qc_summary_path = NA_character_
) {
  message <- paste0(
    sample_id, ": excluded before RDS creation: ", reason,
    " (called=", n_called_after_droplet_qc,
    ", minimum=", minimum_called_cells, ")"
  )
  structure(
    list(
      message = message,
      call = NULL,
      sample_id = sample_id,
      reason = reason,
      n_called_after_droplet_qc = as.integer(n_called_after_droplet_qc),
      minimum_called_cells = as.integer(minimum_called_cells),
      validrops_strategy = validrops_strategy,
      validrops_fallback_source = validrops_fallback_source,
      qc_summary_path = qc_summary_path
    ),
    class = c("atlas_sample_exclusion", "error", "condition")
  )
}

atlas_preprocess_10x_raw <- function(raw_dir, sample_id, qc_output_dir, project_id = NA_character_) {
  for (pkg in c("Matrix", "Seurat", "SeuratObject", "readr", "DropletUtils", "valiDrops", "SoupX", "BiocParallel")) {
    atlas_require_namespace(pkg)
  }

  dir.create(qc_output_dir, recursive = TRUE, showWarnings = FALSE)
  raw_counts <- atlas_read_10x_gene_expression(raw_dir)
  n_raw_barcodes <- ncol(raw_counts)
  lower <- suppressWarnings(as.integer(Sys.getenv("ATLAS_EMPTYDROPS_LOWER", unset = "100")))
  fdr <- suppressWarnings(as.numeric(Sys.getenv("ATLAS_EMPTYDROPS_FDR", unset = "0.01")))
  label_dead <- atlas_env_flag("ATLAS_VALIDROPS_LABEL_DEAD", TRUE)
  drop_dead <- atlas_env_flag("ATLAS_VALIDROPS_DROP_DEAD", FALSE)
  policy <- tolower(Sys.getenv("ATLAS_CELL_CALL_POLICY", unset = "validrops"))

  message(sample_id, ": emptyDrops on ", n_raw_barcodes, " non-zero barcodes")
  ed <- atlas_run_emptydrops(raw_counts, lower = lower, fdr = fdr)
  message(sample_id, ": valiDrops [", atlas_droplet_qc_pipeline_version, "]")
  vd_result <- atlas_run_validrops(
    raw_counts,
    emptydrops = ed,
    raw_dir = raw_dir,
    label_dead = label_dead
  )
  vd <- vd_result$table
  called <- atlas_choose_called_barcodes(raw_counts, ed, vd, policy = policy, drop_dead = drop_dead)
  called_counts <- raw_counts[, called, drop = FALSE]

  prefix <- file.path(qc_output_dir, sample_id)
  atlas_safe_write_csv_gz(ed, paste0(prefix, "__emptyDrops.csv.gz"))
  atlas_safe_write_csv_gz(vd, paste0(prefix, "__valiDrops.csv.gz"))
  atlas_safe_write_csv_gz(
    vd_result$attempts,
    paste0(prefix, "__valiDrops_attempts.csv")
  )
  atlas_safe_write_csv_gz(
    data.frame(barcode = called),
    paste0(prefix, "__called_barcodes.csv.gz")
  )

  minimum_called_cells <- suppressWarnings(as.integer(Sys.getenv(
    "ATLAS_MINIMUM_CALLED_CELLS",
    unset = Sys.getenv("ATLAS_MIN_CALLED_FOR_SOUX", unset = "200")
  )))
  if (!is.finite(minimum_called_cells) || minimum_called_cells < 1L) {
    minimum_called_cells <- 200L
  }

  if (length(called) < minimum_called_cells) {
    exclusion_path <- paste0(prefix, "__droplet_qc_exclusion.csv")
    exclusion_summary <- data.frame(
      project_id = project_id,
      sample_id = sample_id,
      raw_dir = normalizePath(raw_dir, winslash = "/", mustWork = TRUE),
      preprocess_status = "excluded_before_rds",
      exclusion_reason = "called_cells_below_minimum",
      n_raw_barcodes = n_raw_barcodes,
      n_emptydrops_pass = sum(ed$emptydrops_pass %in% TRUE, na.rm = TRUE),
      n_validrops_pass = sum(vd$validrops_qc.pass == "pass", na.rm = TRUE),
      n_called_after_droplet_qc = length(called),
      minimum_called_cells = minimum_called_cells,
      cell_call_policy = policy,
      validrops_status = vd_result$status,
      validrops_strategy = vd_result$strategy,
      validrops_fallback_source = vd_result$fallback_source,
      stringsAsFactors = FALSE
    )
    atlas_safe_write_csv_gz(exclusion_summary, exclusion_path)

    stop(atlas_new_sample_exclusion(
      sample_id = sample_id,
      reason = "called_cells_below_minimum",
      n_called_after_droplet_qc = length(called),
      minimum_called_cells = minimum_called_cells,
      validrops_strategy = vd_result$strategy,
      validrops_fallback_source = vd_result$fallback_source,
      qc_summary_path = exclusion_path
    ))
  }

  message(sample_id, ": SoupX on ", length(called), " called barcodes")
  sx <- atlas_run_soupx(raw_counts, called_counts, sample_id)

  md <- data.frame(barcode = called, row.names = called, stringsAsFactors = FALSE)
  ed_sub <- ed[match(called, ed$barcode), , drop = FALSE]
  vd_sub <- vd[match(called, vd$barcode), , drop = FALSE]
  md <- cbind(md, ed_sub[, setdiff(names(ed_sub), "barcode"), drop = FALSE])
  md <- cbind(md, vd_sub[, setdiff(names(vd_sub), "barcode"), drop = FALSE])
  md$soupx_cluster <- unname(sx$clusters[called])
  md$soupx_contamination_fraction <- unname(sx$rho_by_cell[called])
  md$cell_calling_method <- "valiDrops_with_emptyDrops_audit"
  md$cell_call_policy <- policy
  md$emptydrops_fdr_threshold <- fdr
  md$emptydrops_lower <- lower
  md$validrops_label_dead <- label_dead
  md$validrops_drop_dead <- drop_dead
  md$validrops_status <- vd_result$status
  md$validrops_strategy <- vd_result$strategy
  md$validrops_alpha <- vd_result$alpha
  md$validrops_stage_three <- vd_result$stage_three
  md$validrops_n_input_barcodes <- vd_result$n_input_barcodes
  md$validrops_fallback_source <- vd_result$fallback_source
  md$validrops_fallback_filtered_dir <- vd_result$fallback_filtered_dir
  md$n_cellranger_filtered_barcodes <- vd_result$n_cellranger_filtered_barcodes
  md$validrops_fallback_reason <- vd_result$final_error
  md$soupx_method <- sx$method
  md$soupx_status <- sx$status
  md$soupx_tod_input_barcodes <- sx$tod_input_barcodes
  md$soupx_tod_used_barcodes <- sx$tod_used_barcodes
  md$soupx_ambient_barcodes_used <- sx$ambient_barcodes_used
  md$soupx_tod_reduced <- sx$tod_reduced
  md$n_raw_barcodes <- n_raw_barcodes
  md$n_nonzero_barcodes <- n_raw_barcodes
  md$n_emptydrops_pass <- sum(ed$emptydrops_pass %in% TRUE, na.rm = TRUE)
  md$n_validrops_pass <- sum(vd$validrops_qc.pass == "pass", na.rm = TRUE)
  md$n_called_after_droplet_qc <- length(called)
  md$preqc_stage <- "post_emptyDrops_validrops_SoupX_pre_scDblFinder_miQC"
  md$metadata_embedded_in_loader <- TRUE
  md$input_sample_id <- sample_id

  atlas_safe_write_csv_gz(
    data.frame(
      project_id = project_id,
      sample_id = sample_id,
      raw_dir = normalizePath(raw_dir, winslash = "/", mustWork = TRUE),
      n_raw_barcodes = n_raw_barcodes,
      n_emptydrops_pass = sum(ed$emptydrops_pass %in% TRUE, na.rm = TRUE),
      n_validrops_pass = sum(vd$validrops_qc.pass == "pass", na.rm = TRUE),
      n_called = length(called),
      minimum_called_cells = minimum_called_cells,
      preprocess_status = "completed",
      exclusion_reason = NA_character_,
      cell_call_policy = policy,
      validrops_label_dead = label_dead,
      validrops_drop_dead = drop_dead,
      validrops_status = vd_result$status,
      validrops_strategy = vd_result$strategy,
      validrops_alpha = vd_result$alpha,
      validrops_stage_three = vd_result$stage_three,
      validrops_n_input_barcodes = vd_result$n_input_barcodes,
      validrops_fallback_source = vd_result$fallback_source,
      validrops_fallback_filtered_dir = vd_result$fallback_filtered_dir,
      n_cellranger_filtered_barcodes = vd_result$n_cellranger_filtered_barcodes,
      validrops_fallback_reason = vd_result$final_error,
      soupx_method = sx$method,
      soupx_status = sx$status,
      soupx_tod_input_barcodes = sx$tod_input_barcodes,
      soupx_tod_used_barcodes = sx$tod_used_barcodes,
      soupx_ambient_barcodes_used = sx$ambient_barcodes_used,
      soupx_tod_reduced = sx$tod_reduced,
      stringsAsFactors = FALSE
    ),
    paste0(prefix, "__droplet_qc_summary.csv")
  )

  rm(raw_counts)
  invisible(gc())

  list(
    counts = sx$corrected_counts,
    raw_called_counts = called_counts,
    meta = md,
    read_mode = "Read10X_unfiltered__emptyDrops__valiDrops__SoupX",
    droplet_qc_summary = md[1, c(
      "n_raw_barcodes", "n_emptydrops_pass", "n_validrops_pass",
      "n_called_after_droplet_qc", "cell_call_policy",
      "validrops_status", "validrops_strategy", "validrops_alpha",
      "validrops_stage_three", "validrops_n_input_barcodes",
      "validrops_fallback_source", "n_cellranger_filtered_barcodes",
      "soupx_method", "soupx_status", "soupx_tod_input_barcodes",
      "soupx_tod_used_barcodes", "soupx_ambient_barcodes_used",
      "soupx_tod_reduced"
    ), drop = FALSE]
  )
}

atlas_add_droplet_qc_assays <- function(obj, result) {
  if (!is.null(result$raw_called_counts)) {
    raw_counts <- result$raw_called_counts[rownames(obj), colnames(obj), drop = FALSE]
    obj[["RAW"]] <- SeuratObject::CreateAssayObject(counts = raw_counts)
  }
  Seurat::DefaultAssay(obj) <- "RNA"
  obj
}
