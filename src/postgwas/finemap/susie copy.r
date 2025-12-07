#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(glue)
  library(susieR)
  library(Matrix)
  library(future)
  library(future.apply)
  library(R.utils)
  library(ggplot2)
  library(ggrepel)
  library(parallel)
  library(purrr)
})

## ── 0. Locate and source utilities from same folder ──────────────

get_script_dir <- function() {
  # 1) Rscript
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  fileArg <- grep("^--file=", cmdArgs, value = TRUE)
  if (length(fileArg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", fileArg))))
  }
  # 2) sys.frames (some environments)
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  # 3) Fallback
  getwd()
}

this_dir <- get_script_dir()
source(glue::glue("{this_dir}/utlities.r"))   # assumes this file exists

## ── Helper: safe %||% operator ───────────────────────────────────
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

## ── Helper: message wrapper ──────────────────────────────────────
msg <- function(verbose, ...) {
  if (isTRUE(verbose)) message(...)
}

## ── 1. RAM detection (parallel-safe) ─────────────────────────────
get_free_memory_gb <- function() {
  os <- Sys.info()[["sysname"]]
  
  if (os == "Darwin") {
    out <- suppressWarnings(system("sysctl hw.memsize", intern = TRUE, ignore.stderr = TRUE))
    if (length(out) == 0) return(Inf)
    total <- as.numeric(sub("hw.memsize: ([0-9]+)", "\\1", out))
    if (is.na(total)) return(Inf)
    # Very conservative guess: 25% of total as "free"
    return(0.25 * total / 1024^3)
  }
  
  if (os == "Linux") {
    mem <- tryCatch(readLines("/proc/meminfo"), error = function(e) character(0))
    avail <- grep("^MemAvailable:", mem, value = TRUE)
    if (length(avail) == 0) return(Inf)
    kb <- as.numeric(sub(".* ([0-9]+) kB.*", "\\1", avail))
    return(kb / 1024^2)
  }
  
  # Windows / unknown
  Inf
}

## ── 2. Worker auto-detection ─────────────────────────────────────
auto_detect_workers <- function(
  min_ram_per_worker_gb = 4,
  reserve_cores = 1,
  verbose = TRUE
) {
  phys_cores   <- parallel::detectCores(logical = FALSE) %||% parallel::detectCores()
  max_by_cores <- max(1L, phys_cores - reserve_cores)
  
  free_gb <- tryCatch(get_free_memory_gb(), error = function(e) Inf)
  if (free_gb == Inf || free_gb < 1) {
    if (verbose) message("Free RAM unknown → using ", max_by_cores, " workers (core-based)")
    return(max_by_cores)
  }
  
  max_by_ram <- floor(free_gb / min_ram_per_worker_gb)
  workers    <- max(1L, min(max_by_cores, max_by_ram))
  
  if (verbose) message("Auto workers: ~", round(free_gb, 1), " GB free → ", workers, " workers")
  workers
}

## ── 3. Memory throttle (Linux only) ──────────────────────────────
wait_for_memory <- function(
  threshold_used_pct   = 85,
  check_interval_sec   = 10,
  verbose              = TRUE
) {
  if (Sys.info()[["sysname"]] %in% c("Darwin", "Windows")) {
    return(invisible(TRUE))
  }
  
  repeat {
    mem <- tryCatch(readLines("/proc/meminfo", n = 30), error = function(e) character(0))
    if (length(mem) == 0) { Sys.sleep(check_interval_sec); next }
    
    get_kb <- function(k) {
      as.numeric(sub(".* ([0-9]+) kB.*", "\\1",
                     mem[grepl(paste0("^", k, ":"), mem)])) %||% 0
    }
    
    total     <- get_kb("MemTotal")
    available <- get_kb("MemFree") + get_kb("Buffers") + get_kb("Cached") + get_kb("SReclaimable")
    
    if (total == 0) { Sys.sleep(check_interval_sec); next }
    
    used_pct <- 100 * (1 - available / total)
    if (used_pct < threshold_used_pct) break
    
    if (verbose) {
      message("High memory: ", round(used_pct, 1),
              "% used — waiting ", check_interval_sec, "s...")
    }
    Sys.sleep(check_interval_sec)
  }
  invisible(TRUE)
}

## ── 4. Timeout wrapper ───────────────────────────────────────────
run_with_timeout <- function(expr, timeout, step, tag, verbose = TRUE) {
  if (!is.finite(timeout)) {
    return(tryCatch(expr, error = function(e) NULL))
  }
  
  tryCatch(
    R.utils::withTimeout(expr, timeout = timeout, onTimeout = "error"),
    TimeoutException = function(e) {
      msg(verbose, "Timeout (>", timeout, "s) during ", step, " for ", tag)
      structure(list(timeout = TRUE), class = "timeout")
    },
    error = function(e) {
      msg(verbose, "Error during ", step, " for ", tag, ": ", e$message)
      NULL
    }
  )
}

## ── 5. Fast LD repair ────────────────────────────────────────────
LD_fix_fast <- function(ld.mat, tag = "", verbose = TRUE) {
  ld.mat  <- (ld.mat + t(ld.mat)) / 2
  start_t <- Sys.time()
  msg(verbose, "Fast LD repair started for ", tag)
  
  e <- eigen(ld.mat, symmetric = TRUE)
  e$values[e$values < 1e-8] <- 1e-8
  repaired <- e$vectors %*% diag(e$values) %*% t(e$vectors)
  
  ok <- isTRUE(tryCatch({ chol(repaired); TRUE }, error = function(e) FALSE))
  if (ok) {
    msg(verbose, "Fast LD repair OK (",
        round(difftime(Sys.time(), start_t, units = "secs"), 2), "s)")
    return(repaired)
  }
  
  msg(verbose, "Fast repair failed → falling back to nearPD")
  npd <- Matrix::nearPD(ld.mat, corr = TRUE, keepDiag = TRUE)
  as.matrix(npd$mat)
}

## ── 6. Main locus processor (PASS 1) ─────────────────────────────
process_locus <- function(
  locus,
  df,
  sample_id,
  ld_ref,
  plink,
  analysis_folder,
  lp_threshold,
  L,
  timeout_ld,
  timeout_susie,
  skip_mhc,
  mhc_start,
  mhc_end,
  verbose = TRUE
) {
  
  ## --------------------------------------------------------
  ## LOG REDIRECTION (per-locus)
  ## --------------------------------------------------------

    # Make path absolute to avoid future-globals problems
    analysis_folder <- normalizePath(analysis_folder, mustWork = TRUE)

    tag <- glue("{locus$chr}_{locus$start}_{locus$end}")

    logfile <- file.path(analysis_folder, "logs", glue("{sample_id}_{tag}.log"))

    # Ensure logs folder exists (workers avoid race conditions)
    dir.create(dirname(logfile), recursive = TRUE, showWarnings = FALSE)

    # Cleanup sinks
    while (sink.number() > 0) sink(NULL)
    while (sink.number(type="message") > 0) sink(NULL, type="message")

    # Redirect
    sink(logfile, append = FALSE, type = "output")
    sink(logfile, append = TRUE, type = "message")

    on.exit({
      sink(NULL, type="message")
      sink(NULL)
    }, add = TRUE)

    verbose <- FALSE
    options(warn = 1)
    msg(verbose, "Processing locus ", tag)
  
  ## MHC skip
  if (skip_mhc && locus$chr == "6" &&
      locus$start < mhc_end && locus$end > mhc_start) {
    msg(verbose, "Skipping MHC region")
    return(list(status = "SKIP", reason = "mhc", locus = locus))
  }
  
  wait_for_memory(verbose = verbose)
  
  ## Subset sumstats
  idx <- which(
    df$CHR == as.numeric(locus$chr) &
      df$BP >= locus$start &
      df$BP <= locus$end
  )
  if (length(idx) == 0L) {
    return(list(status = "SKIP", reason = "no_snps", locus = locus))
  }
  
  selected <- data.table::as.data.table(df[idx, ])
  
  if (max(selected$LP, na.rm = TRUE) <= lp_threshold) {
    return(list(status = "SKIP", reason = "low_signal", locus = locus))
  }
  
  ## LD matrix
  ld <- postgwas_ld_matrix(
    variants       = selected$SNP,
    bfile          = ld_ref,
    plink          = plink,
    tag            = gsub(":", "-", tag),
    with_alleles   = FALSE,
    logfile        = glue::glue("{analysis_folder}/ld_matrix_related/{sample_id}_{tag}_ld_matrix_local.log"),
    output_folder  = analysis_folder
  )
  
  selected <- selected[match(rownames(ld), SNP), ]
  data.table::set(selected, j = "variable_index", value = seq_len(nrow(selected)))
  
  if (!identical(selected$SNP, rownames(ld))) {
    stop("SNP order mismatch after LD extraction")
  }
  
  z     <- selected$EZ
  n_eff <- median(selected$NEF, na.rm = TRUE)
  if (!is.finite(n_eff) || n_eff <= 0) n_eff <- 1e5
  
  ## SuSiE
  fitted <- tryCatch(
    susieR::susie_rss(
      z                = z,
      R                = ld,
      n                = n_eff,
      L                = L,
      max_iter         = 100,
      return_correlation = TRUE
    ),
    error = function(e) e
  )
  
  if (inherits(fitted, "error") || is.null(fitted) || !isTRUE(fitted$converged)) {
    return(list(
      status   = "RECOVERY",
      locus    = locus,
      selected = selected,
      ld.mat   = ld,
      n_eff    = n_eff
    ))
  }
  
  ## Annotate + plots + FLAMES
  selected_sus <- annotate_susie(
    fitted   = fitted,
    selected = selected,
    locus    = locus
  )
  
  cred_df <- selected_sus[!is.na(cs) & cs != -1]
  
  plot_susie_loci_ld(
    df            = selected_sus,
    ld            = ld,
    sample_id     = sample_id,
    outdir        = glue::glue("{analysis_folder}/plots/"),
    min_pip_label = 0.01,
    bg_size       = 0.01,
    bg_alpha      = 0.15,
    cs_size       = 1,
    cs_alpha      = 0.2
  )
  
  generate_flames_files(
    fitted  = fitted,
    ld      = ld,
    snp_df  = selected,
    outfile = glue::glue("{analysis_folder}/falmes_output/{sample_id}_{tag}")
  )
  
  ## Save results
  system(glue::glue("mkdir -p {analysis_folder}/locus_files"), ignore.stdout = TRUE)
  system(glue::glue("mkdir -p {analysis_folder}/rds_files"),   ignore.stdout = TRUE)
  
  cred_file <- glue::glue(
    "{analysis_folder}/locus_files/{sample_id}_{tag}_SUSIE_credible_sets.cred1"
  )
  fwrite(cred_df, cred_file, sep = " ", quote = FALSE)
  saveRDS(fitted, glue::glue(
    "{analysis_folder}/rds_files/{sample_id}_{tag}_SUSIE.rds"
  ))
  ## ---- END LOG REDIRECTION ----
  sink(NULL, type = "message")
  sink(NULL)

  return(list(
    status    = "OK",
    vars      = selected_sus,
    cs        = cred_df,
    fitted    = fitted,
    locus     = locus,
    cred_file = cred_file
  ))

}

## ── 7. Recovery (PASS 2, sequential per job) ─────────────────────
recover_locus <- function(
  job,
  sample_id,
  analysis_folder,
  L,
  timeout_ld,
  timeout_susie,
  verbose = TRUE
) {
    locus    <- job$locus
    selected <- job$selected
    ld.mat   <- job$ld.mat
    n_eff    <- job$n_eff
    tag      <- glue::glue("{locus$chr}:{locus$start}-{locus$end}")
    start_rec <- Sys.time()

    ## --------------------------------------------
    ## LOG REDIRECTION (per-locus recovery)
    ## --------------------------------------------
    logfile <- glue::glue("{analysis_folder}/logs/{sample_id}_{tag}_RECOVERY.log")
    dir.create(dirname(logfile), recursive = TRUE, showWarnings = FALSE)

    # Cleanup any previous sinks (paranoia safe)
    while (sink.number() > 0) sink(NULL)
    while (sink.number(type = "message") > 0) sink(NULL, type = "message")

    # Redirect stdout + stderr/messages
    sink(logfile, append = FALSE, type = "output")
    sink(logfile, append = TRUE,  type = "message")

    # Automatically close sinks when function exits
    on.exit({
        sink(NULL, type = "message")
        sink(NULL)
    }, add = TRUE)

    # Disable console output
    verbose <- FALSE
    options(warn = 1)  # warnings go to message sink

    ## Helper: convert SuSiE fit → standard output (+ plotting & FLAMES)
    sus_to_out <- function(fitted, ld_matrix, status, reason = " ", suffix = "") {
      selected_sus <- annotate_susie(
        fitted   = fitted,
        selected = selected,
        locus    = locus
      )
      
      plot_susie_loci_ld(
        df            = selected_sus,
        ld            = ld_matrix,
        sample_id     = sample_id,
        outdir        = glue::glue("{analysis_folder}/plots/"),
        min_pip_label = 0.01,
        bg_size       = 0.01,
        bg_alpha      = 0.15,
        cs_size       = 1,
        cs_alpha      = 0.2
      )
      
      generate_flames_files(
        fitted  = fitted,
        ld      = ld_matrix,
        snp_df  = selected,
        outfile = glue::glue("{analysis_folder}/falmes_output/{sample_id}_{tag}")
      )
      
      cred_df <- selected_sus[!is.na(cs) & cs != -1]
      
      system(glue::glue("mkdir -p {analysis_folder}/locus_files"), ignore.stdout = TRUE)
      system(glue::glue("mkdir -p {analysis_folder}/rds_files"),   ignore.stdout = TRUE)
      
      cred_file <- glue::glue(
        "{analysis_folder}/locus_files/{sample_id}_{tag}_SUSIE_credible_sets{suffix}.cred1"
      )
      fwrite(cred_df, cred_file, sep = " ", quote = FALSE)
      saveRDS(fitted, glue::glue(
        "{analysis_folder}/rds_files/{sample_id}_{tag}_SUSIE{suffix}.rds"
      ))
      
      list(
        status    = status,
        reason    = reason,
        vars      = selected_sus,
        cs        = cred_df,
        fitted    = fitted,
        locus     = locus,
        cred_file = cred_file
      )
    }
    
    ## Step 0 – SuSiE attempt 1 (same LD, higher iter)
    msg(verbose, "   Step 0/4: SuSiE attempt 0 (L=", L, ", max_iter=600)")
    fit0 <- run_with_timeout(
      susieR::susie_rss(
        z        = selected$EZ,
        R        = ld.mat,
        n        = n_eff,
        L        = L,
        max_iter = 600
      ),
      timeout = timeout_susie * 3,
      step    = "SuSiE attempt 1",
      tag     = tag,
      verbose = verbose
    )
    
    if (!is.null(fit0) && !inherits(fit0, "timeout") && isTRUE(fit0$converged)) {
      elapsed <- round(difftime(Sys.time(), start_rec, units = "secs"), 1)
      msg(verbose, "   Recovery SUCCESS (attempt 0) in ", elapsed, "s for ", tag)
      return(
        sus_to_out(
          fitted    = fit0,
          ld_matrix = ld.mat,
          status    = glue("susie_rss success with L {L} and max_iter 600")
        )
      )
    }
    
    ## Step 1 – SuSiE attempt 2 (same LD, lower L)
    msg(verbose, "   Step 1/4: SuSiE attempt 1 (L=5, max_iter=600)")
    fit0b <- run_with_timeout(
      susieR::susie_rss(
        z        = selected$EZ,
        R        = ld.mat,
        n        = n_eff,
        L        = 5,
        max_iter = 600
      ),
      timeout = timeout_susie * 3,
      step    = "SuSiE attempt 2",
      tag     = tag,
      verbose = verbose
    )
    
    if (!is.null(fit0b) && !inherits(fit0b, "timeout") && isTRUE(fit0b$converged)) {
      elapsed <- round(difftime(Sys.time(), start_rec, units = "secs"), 1)
      msg(verbose, "   Recovery SUCCESS (attempt 1) in ", elapsed, "s for ", tag)
      return(
        sus_to_out(
          fitted    = fit0b,
          ld_matrix = ld.mat,
          status    = glue("susie_rss success with L 5 and max_iter 600")
        )
      )
    }
    
    ## Step 2 – LD repair with timeout
    msg(verbose, "   Step 2/4: LD repair (eigen + nearPD fallback)")
    wait_for_memory(verbose = verbose)
    
    ld_repaired <- run_with_timeout(
      expr    = LD_fix_fast(ld.mat, tag = tag, verbose = verbose),
      timeout = timeout_ld,
      step    = "LD repair",
      tag     = tag,
      verbose = verbose
    )
    
    if (is.null(ld_repaired) || inherits(ld_repaired, "timeout")) {
      msg(verbose, "   Recovery FAILED at LD repair for ", tag)
      return(list(
        status = "FAILED",
        locus  = locus,
        reason = "LD_repair_timeout_or_error"
      ))
    }
    
    ## Step 3 – SuSiE attempt 3 (repaired LD, same L)
    msg(verbose, "   Step 3/4: SuSiE attempt 3 (L=", L, ", max_iter=300)")
    fit1 <- run_with_timeout(
      susieR::susie_rss(
        z        = selected$EZ,
        R        = ld_repaired,
        n        = n_eff,
        L        = L,
        max_iter = 300
      ),
      timeout = timeout_susie,
      step    = "SuSiE attempt 3",
      tag     = tag,
      verbose = verbose
    )
    
    if (!is.null(fit1) && !inherits(fit1, "timeout") && isTRUE(fit1$converged)) {
      elapsed <- round(difftime(Sys.time(), start_rec, units = "secs"), 1)
      msg(verbose, "   Recovery SUCCESS (attempt 3) in ", elapsed, "s for ", tag)
      return(
        sus_to_out(
          fitted    = fit1,
          ld_matrix = ld_repaired,
          status    = glue("susie_rss success with repaired LD, L {L}, max_iter 300")
        )
      )
    }
    
    ## Step 4 – SuSiE attempt 4 (repaired LD, reduced L)
    L_new <- max(3L, floor(L / 2))
    msg(verbose, "   Step 4/4: SuSiE attempt 4 (L=", L_new, ", max_iter=600)")
    fit2 <- run_with_timeout(
      susieR::susie_rss(
        z        = selected$EZ,
        R        = ld_repaired,
        n        = n_eff,
        L        = L_new,
        max_iter = 600
      ),
      timeout = timeout_susie * 3,
      step    = "SuSiE attempt 4",
      tag     = tag,
      verbose = verbose
    )
    
    if (!is.null(fit2) && !inherits(fit2, "timeout") && isTRUE(fit2$converged)) {
      elapsed <- round(difftime(Sys.time(), start_rec, units = "secs"), 1)
      msg(verbose, "   Recovery SUCCESS (attempt 4, L=", L_new, ") in ", elapsed, "s for ", tag)
      return(
        sus_to_out(
          fitted    = fit2,
          ld_matrix = ld_repaired,
          status    = glue("susie_rss success with repaired LD, L {L_new}, max_iter 600"),
          suffix    = glue("_L{L_new}")
        )
      )
    }
    
  ## Final failure
  ## ---- END LOG REDIRECTION ----
  sink(NULL, type = "message")
  sink(NULL)

  return(list(
    status = "FAILED",
    locus  = locus,
    reason = "all_recovery_attempts_failed"
  ))

  
}

## ── 8. MAIN FUNCTION (PARALLEL DRIVER) ───────────────────────────
run_susie_finemap_parallel <- function(
  locus_file,
  sumstat_file,
  sample_id,
  ld_ref,
  plink,
  SUSIE_Analysis_folder,
  lp_threshold         = 7.3,
  L                    = 10,
  workers              = "auto",
  min_ram_per_worker_gb = 4,
  verbose              = TRUE,
  timeout_ld_seconds   = 180,
  timeout_susie_seconds = 180,
  skip_mhc             = TRUE,
  mhc_start            = 25e6,
  mhc_end              = 35e6
) {
  dir.create(file.path(SUSIE_Analysis_folder, "logs"), recursive = TRUE, showWarnings = FALSE)
  dir.create(SUSIE_Analysis_folder, recursive = TRUE)
  dir.create(file.path(SUSIE_Analysis_folder, "locus_files"), recursive = TRUE, showWarnings = FALSE)
  
  locus_df <- fread(locus_file)[, .(
    chr   = as.character(CHR),
    start = as.integer(START),
    end   = as.integer(END)
  )]
  
  df <- fread(sumstat_file)
  setnames(
    df,
    old         = c("SNP","REF","ALT","CHR"),
    new         = c("SNP","REF","ALT","CHR"),
    skip_absent = TRUE
  )
  df[, CHR := as.character(CHR)]
  
  options(future.globals.maxSize = +Inf)
  options(future.seed = TRUE)
  
  ## Resolve workers
  if (identical(workers, "auto")) {
    workers <- auto_detect_workers(
      min_ram_per_worker_gb = min_ram_per_worker_gb,
      verbose               = verbose
    )
  } else if (!is.numeric(workers) || workers < 1) {
    stop("workers must be 'auto' or a positive integer")
  }
  workers <- as.integer(workers)
  
  ## Parallel backend (macOS-safe, Linux fast)
  os <- Sys.info()[["sysname"]]
  
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  if (identical(os, "Linux")) {
    msg(verbose, "Linux detected → using multicore backend (fast fork)")
    future::plan(future::multicore, workers = workers, gc = TRUE)
  } else {
    msg(verbose, os, " detected → using multisession backend (PSOCK, no fork)")
    future::plan(future::multisession, workers = workers, gc = TRUE)
  }
  
  msg(
    verbose,
    "Starting PASS 1 with ", workers,
    " workers on ", nrow(locus_df), " loci"
  )
  
  pass1 <- future.apply::future_lapply(
    seq_len(nrow(locus_df)),
    function(i) {
      process_locus(
        locus          = locus_df[i, ],
        df             = df,
        sample_id      = sample_id,
        ld_ref         = ld_ref,
        plink          = plink,
        analysis_folder = SUSIE_Analysis_folder,
        lp_threshold   = lp_threshold,
        L              = L,
        timeout_ld     = timeout_ld_seconds,
        timeout_susie  = timeout_susie_seconds,
        skip_mhc       = skip_mhc,
        mhc_start      = mhc_start,
        mhc_end        = mhc_end,
        verbose        = verbose
      )
    },
    future.seed = TRUE
  )
  
  ## PASS 2: Recovery jobs using fewer workers
  rec_jobs <- purrr::keep(pass1, ~ .x$status == "RECOVERY")
  msg(verbose, length(rec_jobs), " loci need recovery")
  
  rec_results <- list()
  
  if (length(rec_jobs)) {
    rec_workers <- max(1L, floor(workers / 2))
    msg(verbose, "Configuring recovery pass with ", rec_workers, " workers")
    
    future::plan(
      future::multisession,
      workers = rec_workers,
      gc      = TRUE
    )
    
    rec_results <- future.apply::future_lapply(
      rec_jobs,
      function(job) {
        recover_locus(
          job            = job,
          sample_id      = sample_id,
          analysis_folder = SUSIE_Analysis_folder,
          L              = L,
          timeout_ld     = timeout_ld_seconds,
          timeout_susie  = timeout_susie_seconds,
          verbose        = verbose
        )
      },
      future.seed = TRUE
    )
  }
  
  credible_set_df <- aggregate_and_write_results(
    c(pass1, rec_results),
    sample_id,
    SUSIE_Analysis_folder,
    verbose
  )
  
  msg(verbose, "SuSiE fine-mapping finished for ", sample_id)
  invisible(TRUE)
}
