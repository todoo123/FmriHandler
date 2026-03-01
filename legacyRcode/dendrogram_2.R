library(data.table)
library(Matrix)
library(igraph)
library(htmlwidgets)
library(viridisLite)

# ==============================================================================
# 0. Helper Functions (Dice, Pruning, Plotting)
# ==============================================================================

# 1. Dice Coefficient Calculation
calculate_dice_coefficient <- function(res1, res2, verbose = FALSE) {
  .sum_nC2 <- function(counts) {
    n <- as.numeric(counts)
    sum(n * (n - 1) / 2)
  }
  extract_labels <- function(obj) {
    if (is.list(obj)) {
      if (!is.null(obj$group_labels)) return(as.integer(obj$group_labels))
      if (!is.null(obj$labels)) return(as.integer(obj$labels))
      if (is.vector(obj)) return(as.integer(obj))
      stop("Labels not found.")
    }
    return(as.integer(obj))
  }
  
  L1 <- extract_labels(res1); L2 <- extract_labels(res2)
  if (length(L1) != length(L2)) stop("Length mismatch")
  
  valid_mask_1 <- (L1 > 0L) & (!is.na(L1))
  valid_mask_2 <- (L2 > 0L) & (!is.na(L2))
  
  pairs_A <- .sum_nC2(table(L1[valid_mask_1]))
  pairs_B <- .sum_nC2(table(L2[valid_mask_2]))
  
  common_mask <- valid_mask_1 & valid_mask_2
  if (sum(common_mask) == 0) return(0)
  
  pairs_intersection <- .sum_nC2(table(L1[common_mask], L2[common_mask]))
  denom <- pairs_A + pairs_B
  if (denom == 0) return(0)
  
  return((2 * pairs_intersection) / denom)
}

# 2. Top-K Pruning
prune_topk_by_col_fast <- function(W, k = 200L, clamp_negative = TRUE) {
  stopifnot(inherits(W, "dgCMatrix"))
  n <- ncol(W); p <- W@p; i <- W@i; x <- W@x
  if (clamp_negative) x <- pmax(x, 0)
  keep_i <- vector("list", n); keep_x <- vector("list", n)
  
  for (j in seq_len(n)) {
    lo <- p[j] + 1L; hi <- p[j + 1L]
    if (lo > hi) { keep_i[[j]] <- integer(0); keep_x[[j]] <- numeric(0); next }
    idx <- lo:hi; vals <- x[idx]; rows <- i[idx] + 1L
    nz <- which(is.finite(vals) & vals > 0)
    if (!length(nz)) { keep_i[[j]] <- integer(0); keep_x[[j]] <- numeric(0); next }
    vals <- vals[nz]; rows <- rows[nz]
    if (length(vals) > k) {
      ord <- order(vals, decreasing = TRUE)[seq_len(k)]
      vals <- vals[ord]; rows <- rows[ord]
    }
    keep_i[[j]] <- rows; keep_x[[j]] <- vals
  }
  ii <- unlist(keep_i); jj <- rep.int(seq_len(n), lengths(keep_i)); xx <- unlist(keep_x)
  Wk <- sparseMatrix(i=ii, j=jj, x=xx, dims=dim(W), giveCsparse=TRUE)
  Wk <- Matrix::forceSymmetric(Wk, uplo="U")
  diag(Wk) <- 0
  Matrix::drop0(Wk)
}

# 3. Visualization
plot_voxels_by_cluster <- function(df, cluster_col, i_col=3, j_col=4, k_col=5, 
                                   showlegend=FALSE, point_size=2, max_colors=12, min_cluster_size=50) {
  if (!cluster_col %in% names(df)) stop("Column not found")
  cl <- df[[cluster_col]]
  tab <- table(cl)
  small_cls <- names(tab[tab < min_cluster_size])
  cl2 <- as.character(cl)
  cl2[cl2 %in% small_cls] <- "isolated(0)"
  cl2[cl2 == "0"] <- "isolated(0)"
  
  levels_ord <- names(sort(table(cl2), decreasing=TRUE))
  cl2 <- factor(cl2, levels=levels_ord)
  n_main <- min(sum(levels_ord != "isolated(0)"), max_colors)
  cols <- c(viridisLite::viridis(n_main, option="D"), "grey80")
  names(cols) <- c(levels_ord[levels_ord != "isolated(0)"], "isolated(0)")
  
  .plot_voxels_3d(df, i_col=i_col, j_col=j_col, k_col=k_col, color_fac=cl2, colors=cols, 
                  showlegend=showlegend, point_size=point_size, opacity=1)
}


# ==============================================================================
# 1. Setup & Data Loading
# ==============================================================================
source("/home/hyunseok/code/ROI_detection.R")
source("/home/hyunseok/code/fmri_smoothing.R")
load("/home/hyunseok/data/fmri_list.RData") 
ref_craddock <- read.csv("/home/hyunseok/data/Craddock200_2mm.csv")

# Reference Data 병합 (Dice 계산용 roi 정보 확보)
truncated_list <- fmri_list[[1]]
common_coords <- merge(truncated_list[, c("i","j","k")], ref_craddock[, c("i","j","k")], by=c("i","j","k"))
frontal_ref_craddock <- merge(ref_craddock, common_coords, by=c("i","j","k"))

fmri_comparison <- fmri_list
for(iter in 1:length(fmri_comparison)){
  fmri_comparison[[iter]] <- merge(fmri_comparison[[iter]], frontal_ref_craddock, by=c('i','j','k'))
}

# Parameters
k_fpca <- 9
k_levels <- c(2,4,8,16,32,64,128,256)
# k_levels <- c(10, 50, 100) # Test용
out_dir <- "/home/hyunseok/plot"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_subj <- length(fmri_comparison)
sample_n <- NULL

# ==============================================================================
# 2. Main Loop (Per Subject)
# ==============================================================================
# for (si in seq_len(n_subj)) {
for (si in seq_len(n_subj)) {
  message("\n==============================")
  message("Subject: ", si, " / ", n_subj)
  message("==============================")
  
  # output directory per subject
  subj_dir <- file.path(out_dir, sprintf("subject_%02d", si))
  dir.create(subj_dir, showWarnings = FALSE, recursive = TRUE)

  # --------------------------------------------------------------------------
  # [Shared Step] fPCA (One time per subject)
  # --------------------------------------------------------------------------
  message(">> Step 0: Running fPCA (Shared Object)...")
  fmri_smooth_fpca <- build_fpca_loading_list(
    fmri_comparison[si],
    i_col = 1, j_col = 2, k_col = 3,
    time_cols = 6:483,
    Kpc = k_fpca,
    center = FALSE, scale = FALSE,
    chunk_size = 2000L, verbose = FALSE, use_RS = TRUE
  )
  
  # Craddock 함수 호환용 ijk 컬럼 추가
  fmri_smooth_fpca[[1]]$ijk <- paste(fmri_smooth_fpca[[1]]$i, fmri_smooth_fpca[[1]]$j, fmri_smooth_fpca[[1]]$k, sep = "_")
  
  # Reference ROI 병합 (Dice 및 Plotting용)
  # fPCA 결과는 컬럼 이름이 변경되므로 좌표 기준으로 roi 다시 붙임
  dt_fpca_base <- as.data.table(fmri_smooth_fpca[[1]])
  ref_dt <- as.data.table(fmri_comparison[[si]])[, .(i, j, k, roi)]
  setkey(dt_fpca_base, i, j, k); setkey(ref_dt, i, j, k)
  dt_base <- merge(dt_fpca_base, ref_dt, by = c("i", "j", "k")) # roi 포함됨

  # --------------------------------------------------------------------------
  # [Method A] Fast Greedy -> HTML Save ONLY
  # --------------------------------------------------------------------------
  message(">> Step A: Fast Greedy Clustering (Viz Save Only)...")
  
  # A-1. Build Graph
  W_cos <- build_weighted_adj_all_streaming_parallel(
    df = fmri_smooth_fpca[[1]],
    time_cols = 4:(3+k_fpca), # fPCA loadings
    weight_fun = weight_cosine_ts,
    i_block_size = 10L, r_min = 0, clamp_negative = TRUE,
    chunk_size = 10000L, parallel = "multicore", workers = 20, tame_blas = TRUE, verbose = FALSE
  )
  W <- Matrix::forceSymmetric(W_cos, uplo="U"); diag(W) <- 0

  # 1. 대칭화 (dsCMatrix 생성)
  W_sym <- Matrix::forceSymmetric(W_cos, uplo="U")
  
  # 2. 일반 희소 행렬(dgCMatrix)로 안전하게 변환
  #    Matrix 패키지 최신 버전 호환을 위해 "generalMatrix" -> "CsparseMatrix" 순서로 변환
  W <- as(as(W_sym, "generalMatrix"), "CsparseMatrix")
  
  # 3. 대각 성분 0 처리
  diag(W) <- 0

  Wk <- prune_topk_by_col_fast(W, k = 5000L, clamp_negative = TRUE) # k=5000 pruning
  
  # A-2. Clustering
  g <- graph_from_adjacency_matrix(Wk, mode="undirected", weighted=TRUE, diag=FALSE)
  g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = list(weight="max"))
  fg <- cluster_fast_greedy(g, weights = E(g)$weight)
  
  # A-3. Save HTML per Level
  for (k in k_levels) {
    # Labeling
    memb <- cut_at(fg, no = k)
    
    # Plotting Data 준비
    dt_plot <- copy(dt_base)
    dt_plot[, cluster_label := memb]
    
    # Sampling (Optional)
    if (!is.null(sample_n) && nrow(dt_plot) > sample_n) {
      set.seed(1); dt_plot <- dt_plot[sample(.N, sample_n)]
    }
    
    # Plot
    p <- plot_voxels_by_cluster(
      dt_plot, cluster_col = "cluster_label",
      i_col = 1, j_col = 2, k_col = 3,
      showlegend = FALSE, point_size = 5
    )
    
    # Save
    outfile <- file.path(subj_dir, sprintf("subject%02d_FastGreedy_k%04d.html", si, k))
    saveWidget(widget = p, file = outfile, selfcontained = FALSE, libdir = "libs")
  }
  message("   -> HTML files saved.")

  # --------------------------------------------------------------------------
  # [Method B] Craddock Parcellation -> Dice Metric Save ONLY
  # --------------------------------------------------------------------------
  message(">> Step B: Craddock Parcellation (Dice Metric Save Only)...")
  
  dice_results <- data.frame(subject_id = si, method = "Craddock_L2", k_level = k_levels, dice = NA_real_)
  
  for (idx in seq_along(k_levels)) {
    curr_k <- k_levels[idx]
    
    tryCatch({
      # B-1. Craddock Spectral Clustering
      # * Note: weight_l2_ts as per previous context for Craddock comparison
      res_craddock <- craddock_group_parcellation(
        fmri_list      = fmri_smooth_fpca, 
        K              = curr_k,
        i_col          = 1, j_col = 2, k_col = 3, 
        ijk_col        = k_fpca + 4,        
        time_cols      = 4:(3 + k_fpca),    
        weight_fun     = weight_l2_ts,      # Craddock 방식은 L2 distance + Spectral
        standardize    = FALSE,
        clamp_negative = NULL,
        r_min          = 0.55,
        chunk_size     = 10000L,
        remove_isolated= TRUE,
        kmeans_itermax = 1000,
        seed           = NA,
        verbose        = FALSE
      )
      
      # B-2. Calculate Dice (No Viz)
      # res_craddock$group_labels 순서는 fmri_smooth_fpca[[1]] 순서와 같음
      pred_labels <- res_craddock$group_labels
      
      # dt_base는 이미 roi(참조값)를 가지고 있고 순서도 fPCA 결과와 동일함 (copy로 보존)
      # 단, merge 과정에서 순서가 섞였을 수 있으니 안전하게 좌표로 다시 매칭하거나
      # dt_base 생성 시점을 신뢰해야 함. 여기서는 좌표 기준 병합된 dt_base에 추가.
      
      dt_dice <- copy(dt_base)
      dt_dice[, pred := pred_labels]
      
      # Dice 계산 (0인 라벨 제외는 함수 내부에서 처리)
      d_val <- calculate_dice_coefficient(dt_dice$roi, dt_dice$pred, verbose = FALSE)
      dice_results$dice[idx] <- d_val
      
    }, error = function(e) {
      message(sprintf("   Error at K=%d: %s", curr_k, e$message))
    })
  }
  
  # B-3. Save CSV
  dice_file <- file.path(subj_dir, sprintf("subject%02d_Craddock_Dice.csv", si))
  write.csv(dice_results, dice_file, row.names = FALSE)
  message("   -> Dice CSV saved: ", dice_file)

} # End Subject Loop
