###############################################################################
# Sex sgACC WGCNA (CTRL) — LOOPED: SST, VIP, PyrL2n3, PyrL5n6
# Keeps params consistent: signed, minModuleSize=100, B=50, nPerm=100, MDC perm=1000
# Writes per-celltype outputs to: <results_dir>/pooled/<CELLTYPE>_*
###############################################################################

### ---- 0) Setup / packages ----
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 600)
options(stringsAsFactors = FALSE)

user_lib <- if (.Platform$OS.type == "windows") {
  file.path(Sys.getenv("USERPROFILE"), "AppData/Local/R/win-library",
            paste0(R.version$major, ".", R.version$minor))
} else {
  file.path(Sys.getenv("HOME"), "R",
            paste0(R.version$major, ".", R.version$minor))
}
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c("WGCNA", "fastcluster", "dplyr", "tidyr", "pbapply")
need_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(need_cran)) install.packages(need_cran, dependencies = TRUE)

bioc_pkgs <- c("GO.db", "AnnotationDbi", "impute", "preprocessCore")
need_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(need_bioc)) BiocManager::install(need_bioc, ask = FALSE, update = TRUE)

suppressPackageStartupMessages({
  library(WGCNA)
  library(GO.db)
  library(dplyr)
  library(tidyr)
  library(pbapply)
})

# Force WGCNA cor()
cor <- WGCNA::cor
allowWGCNAThreads()
set.seed(12345)

### ---- 1) Global parameters you want consistent ----
subgroup <- "pooled"
TOP_GENES <- 5000

# WGCNA build params
MIN_MODULE_SIZE <- 100
DEEPSPLIT <- 2
MERGE_CUT_HEIGHT <- 0.25
NETWORK_TYPE <- "signed"
TOM_TYPE <- "signed"
COR_TYPE <- "pearson"

# Subsampling-based robustness
ROB_B <- 50
ROB_PROP_REF <- 0.5
ROB_NPERM <- 100
ROB_SEED0 <- 12345
ROB_ZCUT <- 6

# MDC settings
MDC_PERM_B <- 1000
MDC_SEED <- 12345

### ---- 2) Paths ----
base_dir <- "/Users/arielsnowden/Desktop/sex:sgACC/WGCNA"
results_dir <- file.path(base_dir, "WGCNA_results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

subgroup_dir <- file.path(results_dir, subgroup)
dir.create(subgroup_dir, showWarnings = FALSE, recursive = TRUE)

### ---- 3) Load data (only once) ----
#SST <- read.csv(file.path(base_dir, "log2cpm_filtered_CTRL_SST.csv"), header = TRUE, row.names = 1)
#VIP <- read.csv(file.path(base_dir, "log2cpm_filtered_CTRL_VIP.csv"), header = TRUE, row.names = 1)
#PyrL2n3 <- read.csv(file.path(base_dir, "log2cpm_filtered_CTRL_PyrL2n3.csv"), header = TRUE, row.names = 1)
#PyrL5n6 <- read.csv(file.path(base_dir, "log2cpm_filtered_CTRL_PyrL5n6.csv"), header = TRUE, row.names = 1)

### ---- 3) Load data (PV only) ----

# helpers (define BEFORE use)
as_numeric_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x[!is.finite(x)] <- NA_real_
  x
}
fix_names <- function(x) { colnames(x) <- gsub("\\.", "-", colnames(x)); x }

# expression
PV <- read.csv(file.path(base_dir, "log2cpm_filtered_CTRL_PVALB.csv"),
               header = TRUE, row.names = 1)
PV <- as_numeric_matrix(PV)
PV <- fix_names(PV)

regions <- list(PVALB = PV)

# metadata
meta <- read.csv(file.path(base_dir, "metadata.csv"), header = TRUE)

# Harmonize metadata
names(meta)[names(meta) == "CT"] <- "CellType"
names(meta)[names(meta) == "ID"] <- "SampleID"

meta_ctrl <- meta %>%
  dplyr::filter(DX == "CTRL") %>%
  dplyr::select(ID = SampleID, CellType, Sex, Age, MOD, Intergenic_rate, PMI, Race)

meta_ctrl$Sex <- factor(meta_ctrl$Sex, levels = c("F","M"))
rownames(meta_ctrl) <- meta_ctrl$ID


### ---- 4) Functions ----

prepare_trait_data <- function(expr_data, clinical_data, region_name, subgroup = "pooled") {
  sample_names <- colnames(expr_data)
  clin <- clinical_data[clinical_data$CellType == region_name, , drop = FALSE]
  clin <- clin[match(sample_names, clin$ID), , drop = FALSE]
  
  keep <- !is.na(clin$ID)
  if (subgroup == "male_only")   keep <- keep & clin$Sex == "M"
  if (subgroup == "female_only") keep <- keep & clin$Sex == "F"
  
  sample_names <- sample_names[keep]
  clin <- clin[keep, , drop = FALSE]
  
  trait <- data.frame(row.names = sample_names)
  
  if (subgroup == "pooled") {
    trait$Sex <- ifelse(clin$Sex == "M", 1L, 0L)
  }
  
  add_if <- function(col, scale_it = FALSE) {
    if (col %in% names(clin)) {
      v <- clin[[col]]
      if (scale_it) v <- as.numeric(scale(v))
      trait[[col]] <- v
    }
  }
  add_if("Age", TRUE)
  add_if("PMI", TRUE)
  add_if("Intergenic_rate", TRUE)
  
  nzv <- vapply(trait, function(x) length(unique(na.omit(x))) > 1, logical(1))
  if (!all(nzv)) trait <- trait[, nzv, drop = FALSE]
  
  cc <- complete.cases(trait)
  if (any(!cc)) trait <- trait[cc, , drop = FALSE]
  
  trait
}

prepare_expression_data <- function(expr_data, traitData, top_genes = 5000) {
  keep_samples <- intersect(colnames(expr_data), rownames(traitData))
  if (length(keep_samples) == 0) stop("No overlap between expression and trait IDs.")
  expr_data_clean <- expr_data[, keep_samples, drop = FALSE]
  
  expr_data_clean <- as.matrix(expr_data_clean)
  storage.mode(expr_data_clean) <- "double"
  
  gene_var <- apply(expr_data_clean, 1, var, na.rm = TRUE)
  ok_genes <- is.finite(gene_var) & (gene_var > 0)
  if (sum(ok_genes) == 0) stop("All genes have zero/NA variance.")
  gene_var <- gene_var[ok_genes]
  
  o <- order(gene_var, decreasing = TRUE, na.last = NA)
  n_pick <- min(top_genes, length(o))
  sel_idx <- o[seq_len(n_pick)]
  
  datExpr <- t(expr_data_clean[names(gene_var)[sel_idx], , drop = FALSE])
  storage.mode(datExpr) <- "double"
  datExpr
}

quality_control <- function(datExpr, traitData, out_prefix) {
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) {
    datExpr   <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
    traitData <- traitData[rownames(datExpr), , drop = FALSE]
  }
  
  sampleTree <- hclust(dist(datExpr), method = "average")
  pdf(file.path(subgroup_dir, paste0(out_prefix, "_sample_clustering.pdf")), height = 8, width = 12)
  plot(sampleTree, main = paste("Sample clustering -", out_prefix), xlab = "", sub = "", cex = 0.8)
  dev.off()
  
  list(datExpr = datExpr, traitData = traitData)
}

choose_soft_power <- function(datExpr, out_prefix) {
  powers <- 1:20
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  
  pdf(file.path(subgroup_dir, paste0(out_prefix, "_soft_threshold_selection.pdf")), height = 6, width = 12)
  par(mfrow = c(1,2))
  
  plot(sft$fitIndices[,1],
       -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
       xlab="Soft Threshold (power)", ylab="Scale Free Topology Fit, signed R^2",
       type="n", main = paste("Scale Independence -", out_prefix))
  text(sft$fitIndices[,1],
       -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
       labels=powers, col="red")
  abline(h=0.80, col="red", lty=2)
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
       main = paste("Mean Connectivity -", out_prefix))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
  
  dev.off()
  
  suitable <- sft$fitIndices[sft$fitIndices[,2] > 0.8, 1]
  chosen_power <- if (length(suitable) > 0) suitable[1] else 6
  
  list(chosen_power = chosen_power, sft = sft)
}

construct_network <- function(datExpr, soft_power, out_prefix) {
  net <- blockwiseModules(
    datExpr,
    power              = soft_power,
    networkType        = NETWORK_TYPE,
    TOMType            = TOM_TYPE,
    corType            = COR_TYPE,
    minModuleSize      = MIN_MODULE_SIZE,
    deepSplit          = DEEPSPLIT,
    reassignThreshold  = 0,
    mergeCutHeight     = MERGE_CUT_HEIGHT,
    numericLabels      = FALSE,
    pamRespectsDendro  = FALSE,
    saveTOMs           = FALSE,
    verbose            = 3
  )
  
  moduleColors <- net$colors
  MEs0 <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
  MEs  <- orderMEs(MEs0)
  
  pdf(file.path(subgroup_dir, paste0(out_prefix, "_module_dendrogram.pdf")), height = 8, width = 12)
  for (b in seq_along(net$dendrograms)) {
    plotDendroAndColors(net$dendrograms[[b]],
                        moduleColors[net$blockGenes[[b]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
  }
  dev.off()
  
  assignments <- data.frame(Gene = colnames(datExpr), Module = moduleColors, stringsAsFactors = FALSE)
  write.csv(assignments, file = file.path(subgroup_dir, paste0(out_prefix, "_gene_module_assignments.csv")),
            row.names = FALSE)
  
  write.csv(as.data.frame(MEs),
            file = file.path(subgroup_dir, paste0(out_prefix, "_module_eigengenes.csv")),
            row.names = TRUE)
  
  list(net = net, MEs = MEs, moduleColors = moduleColors)
}

module_trait_analysis <- function(MEs, traitData, out_prefix) {
  moduleTraitCor    <- cor(MEs, traitData, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(MEs))
  moduleTraitPvalue_FDR <- apply(moduleTraitPvalue, 2, function(x) p.adjust(x, method = "BH"))
  
  get_stars <- function(p) {
    out <- matrix("", nrow = nrow(p), ncol = ncol(p))
    out[p < 0.001] <- "***"
    out[p < 0.01 & p >= 0.001] <- "**"
    out[p < 0.05 & p >= 0.01]  <- "*"
    out[p < 0.10 & p >= 0.05] <- "."
    out
  }
  
  draw_heatmap <- function(pmat, tag) {
    stars <- get_stars(pmat)
    textMatrix <- matrix("", nrow = nrow(moduleTraitCor), ncol = ncol(moduleTraitCor))
    for (i in 1:nrow(moduleTraitCor)) for (j in 1:ncol(moduleTraitCor)) {
      textMatrix[i,j] <- paste0(signif(moduleTraitCor[i,j], 2), "\n(",
                                signif(pmat[i,j], 2), ")", stars[i,j])
    }
    pdf(file.path(subgroup_dir, paste0(out_prefix, "_module_trait_", tag, ".pdf")), height = 10, width = 8)
    par(mar = c(8, 10, 4, 2))
    labeledHeatmap(Matrix = moduleTraitCor,
                   xLabels = colnames(traitData),
                   yLabels = names(MEs),
                   colorLabels = TRUE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.8,
                   zlim = c(-1, 1),
                   main = paste("Module–Trait -", out_prefix, "-", tag))
    dev.off()
  }
  
  draw_heatmap(moduleTraitPvalue, "uncorrected")
  draw_heatmap(moduleTraitPvalue_FDR, "FDR_corrected")
  
  tidy <- as.data.frame(as.table(moduleTraitCor))
  names(tidy) <- c("Module","Trait","Correlation")
  tidy$Pvalue_raw <- as.vector(moduleTraitPvalue)
  tidy$Pvalue_FDR <- as.vector(moduleTraitPvalue_FDR)
  tidy$Significant_FDR_0.05 <- tidy$Pvalue_FDR < 0.05
  
  write.csv(tidy, file = file.path(subgroup_dir, paste0(out_prefix, "_module_trait_summary.csv")),
            row.names = FALSE)
  
  invisible(tidy)
}

export_hubs <- function(datExpr, MEs, moduleColors, out_prefix, topN = 25) {
  kME <- signedKME(datExpr, MEs)
  mods <- unique(moduleColors)
  
  hub_list <- lapply(mods, function(m) {
    me_col <- paste0("ME", m)
    if (!me_col %in% colnames(kME)) return(NULL)
    idx <- which(moduleColors == m)
    if (!length(idx)) return(NULL)
    ord <- idx[order(kME[idx, me_col], decreasing = TRUE)]
    top <- head(ord, topN)
    data.frame(Module = m,
               Gene = colnames(datExpr)[top],
               kME = kME[top, me_col],
               stringsAsFactors = FALSE)
  })
  hub_list <- hub_list[!vapply(hub_list, is.null, logical(1))]
  if (length(hub_list)) {
    hub_table <- do.call(rbind, hub_list)
    write.csv(hub_table,
              file = file.path(subgroup_dir, paste0(out_prefix, "_hub_genes_top", topN, "_per_module.csv")),
              row.names = FALSE)
  }
}

# Robustness: Z + medianRank extraction (per subsample)
run_robustness <- function(datExpr, moduleColors, out_prefix) {
  dat_all <- datExpr
  colors_all <- moduleColors
  mod_levels <- sort(unique(colors_all))
  
  Z_mat <- matrix(NA_real_, nrow = length(mod_levels), ncol = ROB_B,
                  dimnames = list(mod_levels, paste0("rep", 1:ROB_B)))
  MR_mat <- matrix(NA_real_, nrow = length(mod_levels), ncol = ROB_B,
                   dimnames = list(mod_levels, paste0("rep", 1:ROB_B)))
  
  all_Z_tables <- vector("list", length = ROB_B)
  names(all_Z_tables) <- paste0("rep", 1:ROB_B)
  
  set.seed(ROB_SEED0)
  pb <- txtProgressBar(min = 0, max = ROB_B, style = 3)
  
  for (b in seq_len(ROB_B)) {
    setTxtProgressBar(pb, b)
    
    idx_ref  <- sample(seq_len(nrow(dat_all)), size = floor(ROB_PROP_REF * nrow(dat_all)))
    idx_test <- setdiff(seq_len(nrow(dat_all)), idx_ref)
    
    dat_ref  <- dat_all[idx_ref, , drop = FALSE]
    dat_test <- dat_all[idx_test, , drop = FALSE]
    
    multiExpr_all <- list(ref = list(data = dat_ref),
                          test = list(data = dat_test))
    colorList_all <- list(ref = colors_all)
    
    presQC <- modulePreservation(
      multiExpr_all,
      colorList_all,
      referenceNetworks = 1,
      nPermutations     = ROB_NPERM,
      networkType       = NETWORK_TYPE,
      corFnc            = "cor",
      corOptions        = "use = 'p'",
      randomSeed        = ROB_SEED0 + b,
      verbose           = 0
    )
    
    tab_Z <- presQC$preservation$Z$ref.ref$inColumnsAlsoPresentIn.test
    tab_Z <- tab_Z[rownames(tab_Z) != "gold", , drop = FALSE]
    all_Z_tables[[b]] <- tab_Z
    
    if ("Zsummary.pres" %in% colnames(tab_Z)) {
      tmp <- tab_Z[, "Zsummary.pres"]
      names(tmp) <- rownames(tab_Z)
      Z_mat[names(tmp), b] <- tmp
    }
    
    # medianRank is in observed
    tab_obs <- presQC$preservation$observed$ref.ref$inColumnsAlsoPresentIn.test
    tab_obs <- tab_obs[rownames(tab_obs) != "gold", , drop = FALSE]
    if ("medianRank.pres" %in% colnames(tab_obs)) {
      tmp2 <- tab_obs[, "medianRank.pres"]
      names(tmp2) <- rownames(tab_obs)
      MR_mat[names(tmp2), b] <- tmp2
    }
  }
  close(pb)
  
  Z_summary <- data.frame(
    Module       = rownames(Z_mat),
    Z_mean       = rowMeans(Z_mat, na.rm = TRUE),
    Z_sd         = apply(Z_mat, 1, sd, na.rm = TRUE),
    n_reps_nonNA = apply(Z_mat, 1, function(x) sum(is.finite(x))),
    stringsAsFactors = FALSE
  )
  
  medianRank_summary <- data.frame(
    Module          = rownames(MR_mat),
    medianRank_mean = rowMeans(MR_mat, na.rm = TRUE),
    medianRank_sd   = apply(MR_mat, 1, sd, na.rm = TRUE),
    n_reps_nonNA    = apply(MR_mat, 1, function(x) sum(is.finite(x))),
    stringsAsFactors = FALSE
  )
  
  # all Z metrics summary
  Z_long <- dplyr::bind_rows(lapply(seq_along(all_Z_tables), function(b) {
    tab <- all_Z_tables[[b]]
    if (is.null(tab)) return(NULL)
    df <- as.data.frame(tab)
    df$Module <- rownames(tab)
    df$Rep <- b
    tidyr::pivot_longer(df,
                        cols = -c(Module, Rep),
                        names_to = "Metric",
                        values_to = "Value")
  }))
  
  Z_metrics_summary <- Z_long %>%
    group_by(Module, Metric) %>%
    summarise(
      Mean         = mean(Value, na.rm = TRUE),
      SD           = sd(Value, na.rm = TRUE),
      n_reps_nonNA = sum(is.finite(Value)),
      .groups = "drop"
    )
  
  # append medianRank into same long-format summary
  MR_long <- data.frame(
    Module       = rownames(MR_mat),
    Metric       = "medianRank.pres",
    Mean         = rowMeans(MR_mat, na.rm = TRUE),
    SD           = apply(MR_mat, 1, sd, na.rm = TRUE),
    n_reps_nonNA = apply(MR_mat, 1, function(x) sum(is.finite(x))),
    stringsAsFactors = FALSE
  )
  Z_metrics_summary <- dplyr::bind_rows(Z_metrics_summary, MR_long)
  
  # preserved modules based on Z_mean threshold (your approach)
  preserved_modules <- Z_summary$Module[Z_summary$Z_mean > ROB_ZCUT]
  
  # save outputs
  write.csv(Z_summary,
            file = file.path(subgroup_dir, paste0(out_prefix, "_modulePreservation_subsampling_Zsummary_reps", ROB_B, ".csv")),
            row.names = FALSE)
  write.csv(medianRank_summary,
            file = file.path(subgroup_dir, paste0(out_prefix, "_modulePreservation_medianRank_reps", ROB_B, ".csv")),
            row.names = FALSE)
  write.csv(Z_metrics_summary,
            file = file.path(subgroup_dir, paste0(out_prefix, "_modulePreservation_Z_allMetrics_plusMedianRank_mean_sd_reps", ROB_B, ".csv")),
            row.names = FALSE)
  
  saveRDS(Z_mat, file = file.path(subgroup_dir, paste0(out_prefix, "_modulePreservation_Zsummary_Zmatrix_reps", ROB_B, ".rds")))
  saveRDS(MR_mat, file = file.path(subgroup_dir, paste0(out_prefix, "_modulePreservation_medianRank_Zmatrix_reps", ROB_B, ".rds")))
  saveRDS(all_Z_tables, file = file.path(subgroup_dir, paste0(out_prefix, "_modulePreservation_Ztables_allReps_reps", ROB_B, ".rds")))
  
  list(Z_summary = Z_summary,
       medianRank_summary = medianRank_summary,
       Z_metrics_summary = Z_metrics_summary,
       preserved_modules = preserved_modules)
}



####MODULE DIFFERENTIAL CONNECTIVITY####
# DC + MDC
#note Ariel: you need to add the sample_shuffled code and rerun 
#maybe you can do this on the data separately w/out rerunning the whole WGCNA
run_connectivity <- function(datExpr, traitData, moduleColors, chosen_power, preserved_modules, out_prefix) {
  stopifnot("Sex" %in% colnames(traitData))
  
  male_ids   <- rownames(traitData)[traitData$Sex == 1]
  female_ids <- rownames(traitData)[traitData$Sex == 0]
  
  datExpr_male   <- datExpr[rownames(datExpr) %in% male_ids, , drop = FALSE]
  datExpr_female <- datExpr[rownames(datExpr) %in% female_ids, , drop = FALSE]
  
  # intramodular connectivity
  kWithin_male <- intramodularConnectivity.fromExpr(
    datExpr     = datExpr_male,
    colors      = moduleColors,
    power       = chosen_power,
    networkType = NETWORK_TYPE,
    corFnc      = "cor",
    corOptions  = "use = 'p'"
  )
  kWithin_female <- intramodularConnectivity.fromExpr(
    datExpr     = datExpr_female,
    colors      = moduleColors,
    power       = chosen_power,
    networkType = NETWORK_TYPE,
    corFnc      = "cor",
    corOptions  = "use = 'p'"
  )
  
  gene_modules <- data.frame(Gene = colnames(datExpr),
                             Module = moduleColors,
                             stringsAsFactors = FALSE)
  
  k_table <- data.frame(
    gene_modules,
    kWithin_male   = kWithin_male$kWithin,
    kWithin_female = kWithin_female$kWithin,
    kWithin_diff   = kWithin_male$kWithin - kWithin_female$kWithin,
    row.names = NULL, check.names = FALSE
  )
  write.csv(k_table, file = file.path(subgroup_dir, paste0(out_prefix, "_kWithin_by_gene_male_female.csv")), row.names = FALSE)
  
  # module-level Welch t-tests on kWithin
  mods <- sort(unique(moduleColors))
  res_list <- lapply(mods, function(m) {
    idx <- which(moduleColors == m)
    if (length(idx) < 5) {
      return(data.frame(Module=m, nGenes=length(idx),
                        mean_kWithin_m=NA_real_, mean_kWithin_f=NA_real_,
                        diff_m_minus_f=NA_real_, t_stat=NA_real_, p_value=NA_real_, cohen_d=NA_real_))
    }
    km <- kWithin_male$kWithin[idx]; km <- km[is.finite(km)]
    kf <- kWithin_female$kWithin[idx]; kf <- kf[is.finite(kf)]
    if (length(km) < 3 || length(kf) < 3) {
      return(data.frame(Module=m, nGenes=length(idx),
                        mean_kWithin_m=mean(km), mean_kWithin_f=mean(kf),
                        diff_m_minus_f=mean(km)-mean(kf),
                        t_stat=NA_real_, p_value=NA_real_, cohen_d=NA_real_))
    }
    tt <- t.test(km, kf, var.equal = FALSE)
    
    m1 <- mean(km); m2 <- mean(kf)
    s1 <- sd(km);   s2 <- sd(kf)
    n1 <- length(km); n2 <- length(kf)
    sp <- sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2))
    d  <- (m1 - m2) / sp
    
    data.frame(Module=m, nGenes=length(idx),
               mean_kWithin_m=m1, mean_kWithin_f=m2,
               diff_m_minus_f=m1 - m2,
               t_stat=unname(tt$statistic),
               p_value=unname(tt$p.value),
               cohen_d=d)
  })
  diffconn_tbl <- do.call(rbind, res_list)
  diffconn_tbl$FDR_BH <- p.adjust(diffconn_tbl$p_value, method = "BH")
  
  write.csv(diffconn_tbl,
            file = file.path(subgroup_dir, paste0(out_prefix, "_diff_connectivity_male_vs_female.csv")),
            row.names = FALSE)
  
  # MDC
  preserved_modules <- preserved_modules[preserved_modules != "grey"]
  if (!length(preserved_modules)) {
    warning("No preserved modules for MDC (after removing grey). Skipping MDC.")
    return(invisible(NULL))
  }
  
  adj_M <- adjacency(datExpr_male,   power = chosen_power, type = "signed")
  adj_F <- adjacency(datExpr_female, power = chosen_power, type = "signed")
  
  compute_MDC <- function(adj1, adj2, genes_idx) {
    sub1 <- adj1[genes_idx, genes_idx, drop = FALSE]
    sub2 <- adj2[genes_idx, genes_idx, drop = FALSE]
    sum(sub1[lower.tri(sub1)]) / sum(sub2[lower.tri(sub2)])
  }
  
  MDC_M_F <- sapply(preserved_modules, function(mod) {
    genes_idx <- which(moduleColors == mod)
    compute_MDC(adj_M, adj_F, genes_idx)
  })
  
  permute_MDC_genes <- function(adj1, adj2, genes_idx, B = 1000) {
    n_genes <- length(genes_idx)
    all_genes <- seq_len(nrow(adj1))
    obs <- compute_MDC(adj1, adj2, genes_idx)
    
    nulls <- replicate(B, {
      random_genes <- sample(all_genes, n_genes)
      compute_MDC(adj1, adj2, random_genes)
    })
    
    (1 + sum(abs(nulls - 1) >= abs(obs - 1))) / (B + 1)
  }
  
  set.seed(MDC_SEED)
  p_gene_M_F <- pbsapply(preserved_modules, function(mod) {
    genes_idx <- which(moduleColors == mod)
    permute_MDC_genes(adj_M, adj_F, genes_idx, B = MDC_PERM_B)
  })
  q_gene_M_F <- p.adjust(p_gene_M_F, method = "BH")
  
  MDC_summary_sex <- data.frame(
    Module          = preserved_modules,
    MDC_Male_Female = MDC_M_F[preserved_modules],
    p_Male_Female   = p_gene_M_F[preserved_modules],
    q_Male_Female   = q_gene_M_F[preserved_modules],
    Scheme          = "Gene-shuffled",
    stringsAsFactors = FALSE
  )
  
  write.csv(MDC_summary_sex,
            file = file.path(subgroup_dir, paste0(out_prefix, "_MDC_Male_vs_Female_geneShuffled.csv")),
            row.names = FALSE)
  
  MDC_summary_sex
}

### ---- 5) Run loop over the 5 cell types ----

### ---- 5) Run loop over the 5 cell types ----
celltypes_to_run <- c("PVALB", "SST", "VIP", "PyrL2n3", "PyrL5n6")

stopifnot(identical(names(regions), "PVALB"))
stopifnot(identical(celltypes_to_run, "PVALB"))

for (ct in celltypes_to_run) {
  cat("\n============================================================\n")
  cat("RUNNING:", ct, "\n")
  cat("============================================================\n")

  expr <- regions[[ct]]
  out_prefix <- ct  # file prefix
  
  # traits (pooled) and keep Sex only like your PV workflow
  trait <- prepare_trait_data(expr, meta_ctrl, region_name = ct, subgroup = subgroup)
  trait <- trait[, "Sex", drop = FALSE]
  
  # expression prep + QC
  datExpr <- prepare_expression_data(expr, trait, top_genes = TOP_GENES)
  qc <- quality_control(datExpr, trait, out_prefix = out_prefix)
  datExpr <- qc$datExpr
  trait <- qc$traitData
  
  # power + modules
  sp <- choose_soft_power(datExpr, out_prefix = out_prefix)
  chosen_power <- sp$chosen_power
  
  net <- construct_network(datExpr, soft_power = chosen_power, out_prefix = out_prefix)
  MEs <- net$MEs
  colors <- net$moduleColors
  
  # module-trait
  module_trait_analysis(MEs, trait, out_prefix = out_prefix)
  
  # hubs
  export_hubs(datExpr, MEs, colors, out_prefix = out_prefix, topN = 25)
  
  # module size summary
  write.csv(as.data.frame(table(colors)),
            file = file.path(subgroup_dir, paste0(out_prefix, "_module_sizes.csv")),
            row.names = FALSE)
  
  # robustness: Z + medianRank
  rob <- run_robustness(datExpr, colors, out_prefix = out_prefix)
  preserved_modules <- rob$preserved_modules
  
  cat("\nPreserved modules (Z_mean >", ROB_ZCUT, "):\n")
  print(preserved_modules)
  
  # DC + MDC
  run_connectivity(datExpr, trait, colors, chosen_power, preserved_modules, out_prefix = out_prefix)
  
  # save a minimal RDS with key objects (handy for Cytoscape scripts later)
  saveRDS(list(
    celltype = ct,
    chosen_power = chosen_power,
    moduleColors = colors,
    MEs = MEs,
    preserved_modules = preserved_modules
  ), file = file.path(subgroup_dir, paste0(out_prefix, "_keyObjects.rds")))
  
  cat("\nDONE:", ct, "\n")
}

cat("\nAll done. Outputs in:\n", subgroup_dir, "\n")

