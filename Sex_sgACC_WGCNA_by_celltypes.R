###############################################################################
# Sex-stratified sgACC WGCNA and module differential connectivity (non-psychiatric subjects)
# Data analyzed were previously generated in Arbabi et al., 2025 - PMID: 39237723
# Cell types: PVALB, SST, VIP, PyrL2n3, PyrL5n6
# Network: signed; minModuleSize = 100; deepSplit = 2
# Robustness/module preservation: 50 split-half repetitions; 100 preservation permutations
# MDC: male/female ratio with 1,000 gene-set permutations
#
###############################################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(WGCNA)
  library(dplyr)
  library(tidyr)
  library(pbapply)
})

allowWGCNAThreads()
set.seed(12345)

### ---- 1) Analysis settings ----

CELL_TYPES <- c("PVALB", "SST", "VIP", "PyrL2n3", "PyrL5n6")
TOP_GENES <- 5000

MIN_MODULE_SIZE <- 100
DEEP_SPLIT <- 2
MERGE_CUT_HEIGHT <- 0.25
NETWORK_TYPE <- "signed"
TOM_TYPE <- "signed"
COR_TYPE <- "pearson"

ROB_B <- 50
ROB_PROP_REF <- 0.5
ROB_NPERM <- 100
ROB_SEED <- 12345
ROB_ZCUT <- 6

MDC_PERM_B <- 1000
MDC_SEED <- 12345

# Set SGACC_WGCNA_DIR to the folder containing metadata.csv and the five
# log2cpm_filtered_CTRL_<CELLTYPE>.csv files. If unset, use the working folder.
BASE_DIR <- Sys.getenv("SGACC_WGCNA_DIR", unset = getwd())
RESULTS_DIR <- file.path(BASE_DIR, "WGCNA_results", "pooled")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

### ---- 2) Input helpers ----

as_numeric_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x[!is.finite(x)] <- NA_real_
  x
}

fix_sample_names <- function(x) {
  colnames(x) <- gsub("\\.", "-", colnames(x))
  x
}

read_expression <- function(cell_type) {
  path <- file.path(
    BASE_DIR,
    paste0("log2cpm_filtered_CTRL_", cell_type, ".csv")
  )
  if (!file.exists(path)) stop("Missing expression file: ", path)
  x <- read.csv(path, header = TRUE, row.names = 1, check.names = FALSE)
  fix_sample_names(as_numeric_matrix(x))
}

metadata_path <- file.path(BASE_DIR, "metadata.csv")
if (!file.exists(metadata_path)) stop("Missing metadata file: ", metadata_path)

metadata <- read.csv(metadata_path, header = TRUE, stringsAsFactors = FALSE)
names(metadata)[names(metadata) == "CT"] <- "CellType"
names(metadata)[names(metadata) == "ID"] <- "SampleID"

required_metadata <- c(
  "SampleID", "CellType", "Sex", "DX", "Age", "PMI", "Intergenic_rate"
)
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata)) {
  stop("Missing metadata columns: ", paste(missing_metadata, collapse = ", "))
}

metadata_ctrl <- metadata %>%
  filter(DX == "CTRL") %>%
  transmute(
    ID = SampleID,
    CellType,
    Sex = factor(Sex, levels = c("F", "M")),
    Age,
    PMI,
    Intergenic_rate
  )

### ---- 3) WGCNA preparation and construction ----

prepare_sample_data <- function(expr_data, clinical_data, cell_type) {
  sample_names <- colnames(expr_data)
  clinical <- clinical_data[clinical_data$CellType == cell_type, , drop = FALSE]
  clinical <- clinical[match(sample_names, clinical$ID), , drop = FALSE]

  keep <- !is.na(clinical$ID)
  clinical <- clinical[keep, , drop = FALSE]
  sample_names <- sample_names[keep]

  sample_data <- data.frame(
    Sex = ifelse(clinical$Sex == "M", 1L, 0L),
    Age = as.numeric(scale(clinical$Age)),
    PMI = as.numeric(scale(clinical$PMI)),
    Intergenic_rate = as.numeric(scale(clinical$Intergenic_rate)),
    row.names = sample_names,
    check.names = FALSE
  )

  # Preserve the complete-case sample set used in the original pipeline.
  sample_data <- sample_data[complete.cases(sample_data), , drop = FALSE]
  if (!nrow(sample_data)) stop("No complete samples for ", cell_type)
  sample_data
}

prepare_expression_data <- function(expr_data, sample_data, top_genes) {
  keep_samples <- intersect(colnames(expr_data), rownames(sample_data))
  if (!length(keep_samples)) stop("No overlapping expression and metadata IDs.")

  x <- expr_data[, keep_samples, drop = FALSE]
  gene_variance <- apply(x, 1, var, na.rm = TRUE)
  gene_variance <- gene_variance[is.finite(gene_variance) & gene_variance > 0]
  if (!length(gene_variance)) stop("All genes have zero or missing variance.")

  selected <- head(names(sort(gene_variance, decreasing = TRUE)), top_genes)
  datExpr <- t(x[selected, , drop = FALSE])
  storage.mode(datExpr) <- "double"
  datExpr
}

quality_control <- function(datExpr, sample_data, prefix) {
  qc <- goodSamplesGenes(datExpr, verbose = 3)
  if (!qc$allOK) {
    datExpr <- datExpr[qc$goodSamples, qc$goodGenes, drop = FALSE]
  }
  sample_data <- sample_data[rownames(datExpr), , drop = FALSE]

  sample_tree <- hclust(dist(datExpr), method = "average")
  pdf(file.path(RESULTS_DIR, paste0(prefix, "_sample_clustering.pdf")),
      height = 8, width = 12)
  plot(sample_tree, main = paste("Sample clustering -", prefix),
       xlab = "", sub = "", cex = 0.8)
  dev.off()

  list(datExpr = datExpr, sample_data = sample_data)
}

choose_soft_power <- function(datExpr, prefix) {
  powers <- 1:20
  sft <- pickSoftThreshold(
    datExpr,
    powerVector = powers,
    networkType = NETWORK_TYPE,
    corFnc = "cor",
    corOptions = list(use = "p"),
    verbose = 5
  )

  pdf(file.path(RESULTS_DIR, paste0(prefix, "_soft_threshold_selection.pdf")),
      height = 6, width = 12)
  par(mfrow = c(1, 2))
  plot(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft threshold (power)",
    ylab = "Scale-free topology fit, signed R^2",
    type = "n",
    main = paste("Scale independence -", prefix)
  )
  text(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers,
    col = "red"
  )
  abline(h = 0.80, col = "red", lty = 2)
  plot(
    sft$fitIndices[, 1], sft$fitIndices[, 5],
    xlab = "Soft threshold (power)", ylab = "Mean connectivity",
    type = "n", main = paste("Mean connectivity -", prefix)
  )
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
  dev.off()

  suitable <- sft$fitIndices[sft$fitIndices[, 2] > 0.8, 1]
  chosen_power <- if (length(suitable)) suitable[1] else 6
  list(chosen_power = chosen_power, diagnostics = sft)
}

construct_network <- function(datExpr, soft_power, prefix) {
  network <- blockwiseModules(
    datExpr,
    power = soft_power,
    networkType = NETWORK_TYPE,
    TOMType = TOM_TYPE,
    corType = COR_TYPE,
    minModuleSize = MIN_MODULE_SIZE,
    deepSplit = DEEP_SPLIT,
    reassignThreshold = 0,
    mergeCutHeight = MERGE_CUT_HEIGHT,
    numericLabels = FALSE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    verbose = 3
  )

  module_colors <- network$colors
  module_eigengenes <- orderMEs(
    moduleEigengenes(datExpr, colors = module_colors)$eigengenes
  )

  pdf(file.path(RESULTS_DIR, paste0(prefix, "_module_dendrogram.pdf")),
      height = 8, width = 12)
  for (block in seq_along(network$dendrograms)) {
    plotDendroAndColors(
      network$dendrograms[[block]],
      module_colors[network$blockGenes[[block]]],
      "Module colors",
      dendroLabels = FALSE,
      hang = 0.03,
      addGuide = TRUE,
      guideHang = 0.05
    )
  }
  dev.off()

  assignments <- data.frame(
    Gene = colnames(datExpr),
    Module = module_colors,
    stringsAsFactors = FALSE
  )
  write.csv(
    assignments,
    file.path(RESULTS_DIR, paste0(prefix, "_gene_module_assignments.csv")),
    row.names = FALSE
  )
  write.csv(
    as.data.frame(module_eigengenes),
    file.path(RESULTS_DIR, paste0(prefix, "_module_eigengenes.csv")),
    row.names = TRUE
  )
  write.csv(
    as.data.frame(table(Module = module_colors)),
    file.path(RESULTS_DIR, paste0(prefix, "_module_sizes.csv")),
    row.names = FALSE
  )

  list(
    network = network,
    module_colors = module_colors,
    module_eigengenes = module_eigengenes
  )
}

### ---- 4) Split-half module robustness ----

run_robustness <- function(datExpr, module_colors, prefix) {
  modules <- sort(unique(module_colors))
  Z_matrix <- matrix(
    NA_real_, nrow = length(modules), ncol = ROB_B,
    dimnames = list(modules, paste0("rep", seq_len(ROB_B)))
  )
  median_rank_matrix <- Z_matrix
  all_Z_tables <- vector("list", ROB_B)
  names(all_Z_tables) <- colnames(Z_matrix)

  set.seed(ROB_SEED)
  progress <- txtProgressBar(min = 0, max = ROB_B, style = 3)
  on.exit(close(progress), add = TRUE)

  for (iteration in seq_len(ROB_B)) {
    setTxtProgressBar(progress, iteration)
    reference_indices <- sample(
      seq_len(nrow(datExpr)),
      size = floor(ROB_PROP_REF * nrow(datExpr))
    )
    test_indices <- setdiff(seq_len(nrow(datExpr)), reference_indices)

    multi_expression <- list(
      reference = list(data = datExpr[reference_indices, , drop = FALSE]),
      test = list(data = datExpr[test_indices, , drop = FALSE])
    )
    color_list <- list(reference = module_colors)

    preservation <- modulePreservation(
      multi_expression,
      color_list,
      referenceNetworks = 1,
      nPermutations = ROB_NPERM,
      networkType = NETWORK_TYPE,
      corFnc = "cor",
      corOptions = "use = 'p'",
      randomSeed = ROB_SEED + iteration,
      verbose = 0
    )

    Z_table <- preservation$preservation$Z$ref.ref$inColumnsAlsoPresentIn.test
    Z_table <- Z_table[rownames(Z_table) != "gold", , drop = FALSE]
    all_Z_tables[[iteration]] <- Z_table
    if ("Zsummary.pres" %in% colnames(Z_table)) {
      values <- Z_table[, "Zsummary.pres"]
      Z_matrix[names(values), iteration] <- values
    }

    observed <- preservation$preservation$observed$ref.ref$inColumnsAlsoPresentIn.test
    observed <- observed[rownames(observed) != "gold", , drop = FALSE]
    if ("medianRank.pres" %in% colnames(observed)) {
      values <- observed[, "medianRank.pres"]
      median_rank_matrix[names(values), iteration] <- values
    }
  }

  Z_summary <- data.frame(
    Module = rownames(Z_matrix),
    Z_mean = rowMeans(Z_matrix, na.rm = TRUE),
    Z_sd = apply(Z_matrix, 1, sd, na.rm = TRUE),
    n_reps_nonNA = rowSums(is.finite(Z_matrix)),
    stringsAsFactors = FALSE
  )
  median_rank_summary <- data.frame(
    Module = rownames(median_rank_matrix),
    medianRank_mean = rowMeans(median_rank_matrix, na.rm = TRUE),
    medianRank_sd = apply(median_rank_matrix, 1, sd, na.rm = TRUE),
    n_reps_nonNA = rowSums(is.finite(median_rank_matrix)),
    stringsAsFactors = FALSE
  )

  Z_long <- bind_rows(lapply(seq_along(all_Z_tables), function(iteration) {
    table_i <- all_Z_tables[[iteration]]
    if (is.null(table_i)) return(NULL)
    out <- as.data.frame(table_i)
    out$Module <- rownames(table_i)
    out$Rep <- iteration
    pivot_longer(out, cols = -c(Module, Rep),
                 names_to = "Metric", values_to = "Value")
  }))
  all_metrics <- Z_long %>%
    group_by(Module, Metric) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      n_reps_nonNA = sum(is.finite(Value)),
      .groups = "drop"
    )
  all_metrics <- bind_rows(
    all_metrics,
    data.frame(
      Module = rownames(median_rank_matrix),
      Metric = "medianRank.pres",
      Mean = rowMeans(median_rank_matrix, na.rm = TRUE),
      SD = apply(median_rank_matrix, 1, sd, na.rm = TRUE),
      n_reps_nonNA = rowSums(is.finite(median_rank_matrix))
    )
  )

  write.csv(
    Z_summary,
    file.path(
      RESULTS_DIR,
      paste0(prefix, "_modulePreservation_subsampling_Zsummary_reps", ROB_B, ".csv")
    ),
    row.names = FALSE
  )
  write.csv(
    median_rank_summary,
    file.path(
      RESULTS_DIR,
      paste0(prefix, "_modulePreservation_medianRank_reps", ROB_B, ".csv")
    ),
    row.names = FALSE
  )
  write.csv(
    all_metrics,
    file.path(
      RESULTS_DIR,
      paste0(prefix, "_modulePreservation_allMetrics_reps", ROB_B, ".csv")
    ),
    row.names = FALSE
  )

  saveRDS(
    Z_matrix,
    file.path(RESULTS_DIR, paste0(prefix, "_modulePreservation_Zmatrix_reps", ROB_B, ".rds"))
  )
  saveRDS(
    median_rank_matrix,
    file.path(RESULTS_DIR, paste0(prefix, "_modulePreservation_medianRank_matrix_reps", ROB_B, ".rds"))
  )

  list(
    Z_summary = Z_summary,
    median_rank_summary = median_rank_summary,
    preserved_modules = Z_summary$Module[Z_summary$Z_mean > ROB_ZCUT]
  )
}

### ---- 5) Differential connectivity and MDC ----

run_connectivity <- function(datExpr, sample_data, module_colors,
                             chosen_power, preserved_modules, prefix) {
  stopifnot("Sex" %in% colnames(sample_data))

  male_ids <- rownames(sample_data)[sample_data$Sex == 1]
  female_ids <- rownames(sample_data)[sample_data$Sex == 0]
  male_expression <- datExpr[male_ids, , drop = FALSE]
  female_expression <- datExpr[female_ids, , drop = FALSE]

  if (nrow(male_expression) < 4 || nrow(female_expression) < 4) {
    stop("At least four samples per sex are required for connectivity analysis.")
  }

  male_connectivity <- intramodularConnectivity.fromExpr(
    datExpr = male_expression,
    colors = module_colors,
    power = chosen_power,
    networkType = NETWORK_TYPE,
    corFnc = "cor",
    corOptions = "use = 'p'"
  )
  female_connectivity <- intramodularConnectivity.fromExpr(
    datExpr = female_expression,
    colors = module_colors,
    power = chosen_power,
    networkType = NETWORK_TYPE,
    corFnc = "cor",
    corOptions = "use = 'p'"
  )

  gene_connectivity <- data.frame(
    Gene = colnames(datExpr),
    Module = module_colors,
    kWithin_male = male_connectivity$kWithin,
    kWithin_female = female_connectivity$kWithin,
    kWithin_diff = male_connectivity$kWithin - female_connectivity$kWithin,
    stringsAsFactors = FALSE
  )
  write.csv(
    gene_connectivity,
    file.path(RESULTS_DIR, paste0(prefix, "_kWithin_by_gene_male_female.csv")),
    row.names = FALSE
  )

  module_statistics <- lapply(sort(unique(module_colors)), function(module) {
    indices <- which(module_colors == module)
    male_k <- male_connectivity$kWithin[indices]
    female_k <- female_connectivity$kWithin[indices]
    male_k <- male_k[is.finite(male_k)]
    female_k <- female_k[is.finite(female_k)]

    if (length(indices) < 5 || length(male_k) < 3 || length(female_k) < 3) {
      return(data.frame(
        Module = module, nGenes = length(indices),
        mean_kWithin_m = mean(male_k), mean_kWithin_f = mean(female_k),
        diff_m_minus_f = mean(male_k) - mean(female_k),
        t_stat = NA_real_, p_value = NA_real_, cohen_d = NA_real_
      ))
    }

    test <- t.test(male_k, female_k, var.equal = FALSE)
    pooled_sd <- sqrt(
      ((length(male_k) - 1) * var(male_k) +
         (length(female_k) - 1) * var(female_k)) /
        (length(male_k) + length(female_k) - 2)
    )
    data.frame(
      Module = module,
      nGenes = length(indices),
      mean_kWithin_m = mean(male_k),
      mean_kWithin_f = mean(female_k),
      diff_m_minus_f = mean(male_k) - mean(female_k),
      t_stat = unname(test$statistic),
      p_value = test$p.value,
      cohen_d = (mean(male_k) - mean(female_k)) / pooled_sd
    )
  })
  module_statistics <- bind_rows(module_statistics)
  module_statistics$FDR_BH <- p.adjust(module_statistics$p_value, method = "BH")
  write.csv(
    module_statistics,
    file.path(RESULTS_DIR, paste0(prefix, "_diff_connectivity_male_vs_female.csv")),
    row.names = FALSE
  )

  preserved_modules <- setdiff(preserved_modules, "grey")
  if (!length(preserved_modules)) {
    warning("No preserved non-grey modules; skipping MDC.")
    return(invisible(NULL))
  }

  male_adjacency <- adjacency(
    male_expression, power = chosen_power, type = NETWORK_TYPE
  )
  female_adjacency <- adjacency(
    female_expression, power = chosen_power, type = NETWORK_TYPE
  )

  compute_mdc <- function(adjacency_1, adjacency_2, indices) {
    network_1 <- adjacency_1[indices, indices, drop = FALSE]
    network_2 <- adjacency_2[indices, indices, drop = FALSE]
    sum(network_1[lower.tri(network_1)]) /
      sum(network_2[lower.tri(network_2)])
  }

  observed_mdc <- vapply(preserved_modules, function(module) {
    compute_mdc(
      male_adjacency,
      female_adjacency,
      which(module_colors == module)
    )
  }, numeric(1))

  permutation_p <- function(module) {
    module_indices <- which(module_colors == module)
    observed <- compute_mdc(male_adjacency, female_adjacency, module_indices)
    null_values <- replicate(MDC_PERM_B, {
      random_indices <- sample(seq_len(ncol(datExpr)), length(module_indices))
      compute_mdc(male_adjacency, female_adjacency, random_indices)
    })
    (1 + sum(abs(null_values - 1) >= abs(observed - 1))) /
      (MDC_PERM_B + 1)
  }

  set.seed(MDC_SEED)
  p_values <- pbsapply(preserved_modules, permutation_p)
  mdc_summary <- data.frame(
    Module = preserved_modules,
    MDC_Male_Female = observed_mdc[preserved_modules],
    p_Male_Female = p_values[preserved_modules],
    q_Male_Female = p.adjust(p_values[preserved_modules], method = "BH"),
    Scheme = "Gene-shuffled",
    stringsAsFactors = FALSE
  )
  write.csv(
    mdc_summary,
    file.path(RESULTS_DIR, paste0(prefix, "_MDC_Male_vs_Female_geneShuffled.csv")),
    row.names = FALSE
  )
  mdc_summary
}

### ---- 6) Run all cell types ----

for (cell_type in CELL_TYPES) {
  message("\nRunning ", cell_type, "...")

  expression_data <- read_expression(cell_type)
  sample_data <- prepare_sample_data(expression_data, metadata_ctrl, cell_type)
  datExpr <- prepare_expression_data(expression_data, sample_data, TOP_GENES)

  qc <- quality_control(datExpr, sample_data, cell_type)
  datExpr <- qc$datExpr
  sample_data <- qc$sample_data

  power_result <- choose_soft_power(datExpr, cell_type)
  chosen_power <- power_result$chosen_power

  network_result <- construct_network(datExpr, chosen_power, cell_type)
  module_colors <- network_result$module_colors
  module_eigengenes <- network_result$module_eigengenes

  robustness <- run_robustness(datExpr, module_colors, cell_type)
  preserved_modules <- robustness$preserved_modules

  message(
    "Preserved modules (mean Z > ", ROB_ZCUT, "): ",
    paste(preserved_modules, collapse = ", ")
  )

  mdc_summary <- run_connectivity(
    datExpr,
    sample_data,
    module_colors,
    chosen_power,
    preserved_modules,
    cell_type
  )

  saveRDS(
    list(
      cell_type = cell_type,
      chosen_power = chosen_power,
      module_colors = module_colors,
      module_eigengenes = module_eigengenes,
      preserved_modules = preserved_modules,
      mdc_summary = mdc_summary
    ),
    file.path(RESULTS_DIR, paste0(cell_type, "_key_objects.rds"))
  )

  message("Completed ", cell_type)
}

message("\nAll outputs written to: ", normalizePath(RESULTS_DIR))
