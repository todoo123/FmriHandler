# 시각화 수정
plot_voxels_by_cluster <- function(df,
                                   cluster_col = "cluster_200",
                                   i_col = 3, j_col = 4, k_col = 5,
                                   showlegend = FALSE,
                                   point_size = 2,
                                   max_colors = 12,
                                   min_cluster_size = 50) {
  
  if (!cluster_col %in% names(df))
    stop(sprintf("%s 컬럼이 없습니다.", cluster_col))
  
  cl <- df[[cluster_col]]

  # cluster size 계산
  tab <- table(cl)

  # 작은 클러스터 → isolated 처리
  small_cls <- names(tab[tab < min_cluster_size])
  cl2 <- as.character(cl)
  cl2[cl2 %in% small_cls] <- "isolated(0)"
  cl2[cl2 == "0"] <- "isolated(0)"

  # factor 재정렬: 큰 클러스터부터
  tab2 <- sort(table(cl2), decreasing = TRUE)
  levels_ord <- names(tab2)
  cl2 <- factor(cl2, levels = levels_ord)

  # 색상 팔레트
  library(viridisLite)
  n_main <- min(sum(levels_ord != "isolated(0)"), max_colors)

  main_cols <- viridis(n_main, option = "D")  # perceptually uniform
  cols <- c(main_cols, "grey80")
  names(cols) <- c(levels_ord[levels_ord != "isolated(0)"], "isolated(0)")

  .plot_voxels_3d(
  df,
  i_col = i_col, j_col = j_col, k_col = k_col,
  color_fac = cl2,
  colors = cols,
  showlegend = showlegend,
  point_size = point_size,
  opacity = 1
  )
}




source("/home/hyunseok/code/ROI_detection.R")
source("/home/hyunseok/code/fmri_smoothing.R")
load("/home/hyunseok/data/fmri_list.RData")
ref_craddock<-read.csv("/home/hyunseok/data/Craddock200_2mm.csv")



# 전체 보는 경우
truncated_list <- fmri_list[[1]]
# reference code 에서 내 자료의 좌표값만 살려서 다시 시각화
common_coords <- merge(
  truncated_list[, c("i","j","k")],
  ref_craddock[, c("i","j","k")],
  by = c("i","j","k")
)


frontal_ref_craddock <- merge(ref_craddock, common_coords, by = c("i","j","k"))

fmri_comparison <-fmri_list
# 이미 fmri_list 피험자들끼리는 겹치는 값만을 남겨놓았기때문에 common_coords 변수 그대로 사용가능
for(iter in 1:length(fmri_comparison)){
  fmri_comparison[[iter]] <- merge(
    fmri_comparison[[iter]],
    frontal_ref_craddock,
    by = c('i','j','k')
  )
}

################
# fPCA
################
k <- 9
fmri_smooth_fpca <- build_fpca_loading_list(fmri_comparison[1],
                        i_col = 1, j_col = 2, k_col = 3,
                        time_cols = 6:483,
                        Kpc = k,            # 유지할 주성분(기저) 개수
                        center = FALSE,       # 시간축 평균 함수(열평균) 제거
                        scale = FALSE,       # 필요시 시간축 표준화(보통 FALSE 권장)
                        chunk_size = 2000L,
                        verbose = TRUE,
                        use_RS = TRUE        # RSpectra 사용 여부(대규모 T일 때 빠름)
                        )
fmri_smooth_fpca[[1]][1,]

W_cos <- build_weighted_adj_all_streaming_parallel(
      df              = fmri_smooth_fpca[[1]],
      time_cols       = 4:12,
      weight_fun      = weight_cosine_ts,
      i_block_size    = 10L,
      r_min           = 0,
      clamp_negative  = TRUE,
      chunk_size      = 10000L,
      parallel        = "multicore",     # ← 병렬 모드 전달
      workers         = 20,      # ← 워커 수 전달
      tame_blas       = TRUE,    # ← BLAS 억제
      verbose         = TRUE
    )

# dendrogram

library(Matrix)
library(igraph)

# W: dgCMatrix, symmetric, diag=0 권장, 값은 cosine similarity (대개 0~1)
# 예: W <- W_cos  (이미 만들어 둔 것)
W<-W_cos
# 혹시 비대칭이면 대칭화
if (!isSymmetric(W)) {
  W <- Matrix::forceSymmetric(W, uplo = "U")
  W <- as(W, "dgCMatrix")
}

# 대각 0으로
diag(W) <- 0

library(Matrix)

prune_topk_by_col_fast <- function(W, k = 200L, clamp_negative = TRUE) {
  stopifnot(inherits(W, "dgCMatrix"))
  n <- ncol(W)
  p <- W@p; i <- W@i; x <- W@x
  
  # (선택) 음수 제거
  if (clamp_negative) x <- pmax(x, 0)

  # 저장할 공간을 대략 n*k로 잡고 시작(정확히는 컬럼 nnz가 k보다 작으면 그만큼 덜)
  # R에서 완전 prealloc이 어렵지만, 리스트 축적 대신 벡터 누적을 최소화
  keep_i <- vector("list", n)
  keep_x <- vector("list", n)

  for (j in seq_len(n)) {
    lo <- p[j] + 1L
    hi <- p[j + 1L]
    if (lo > hi) {
      keep_i[[j]] <- integer(0); keep_x[[j]] <- numeric(0)
      next
    }
    idx <- lo:hi
    vals <- x[idx]
    rows <- i[idx] + 1L

    # 0 제거
    nz <- which(is.finite(vals) & vals > 0)
    if (!length(nz)) {
      keep_i[[j]] <- integer(0); keep_x[[j]] <- numeric(0)
      next
    }
    vals <- vals[nz]
    rows <- rows[nz]

    if (length(vals) > k) {
      # 부분정렬(전체정렬보다 빠름)
      # base R엔 argpartition이 없어서 sort가 걸리긴 하지만,
      # k가 작으면 여기서 병목이 생길 수 있어. 그래도 지금은 필요함.
      ord <- order(vals, decreasing = TRUE)[seq_len(k)]
      vals <- vals[ord]
      rows <- rows[ord]
    }
    keep_i[[j]] <- rows
    keep_x[[j]] <- vals
  }

  ii <- unlist(keep_i, use.names = FALSE)
  jj <- rep.int(seq_len(n), times = vapply(keep_i, length, 1L))
  xx <- unlist(keep_x, use.names = FALSE)

  Wk <- sparseMatrix(i = ii, j = jj, x = xx, dims = dim(W), giveCsparse = TRUE)
  Wk <- Matrix::forceSymmetric(Wk, uplo = "U")  # union 대칭화 효과
  Wk <- as(Wk, "dgCMatrix")
  diag(Wk) <- 0
  Matrix::drop0(Wk)
}

cat("Before nnz:", Matrix::nnzero(W), "\n")
cat("avg degree ~", (2*Matrix::nnzero(W))/nrow(W), "\n")
Wk <- prune_topk_by_col_fast(W, k = 5000L, clamp_negative = TRUE)
cat("After  nnz:", Matrix::nnzero(Wk), "\n")
cat("After avg degree ~", 2*Matrix::nnzero(Wk)/nrow(Wk), "\n")


library(igraph)

g <- graph_from_adjacency_matrix(Wk, mode="undirected", weighted=TRUE, diag=FALSE)
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE,
              edge.attr.comb = list(weight="max"))

fg <- cluster_fast_greedy(g, weights = E(g)$weight)

k_levels <- c(2,4,8,16,32,64,128,256,512,1024,2048,4096,8192)
labels_mat <- sapply(k_levels, function(k) cut_at(fg, no = k))
colnames(labels_mat) <- sprintf("cl_k%04d", k_levels)

# 결과 병합
library(data.table)

# 예: subject 1만 먼저
dt_fpca <- as.data.table(fmri_smooth_fpca[[1]])

# labels_mat: rownames가 없다고 가정하면, 반드시 좌표로 붙이는 게 안전
# (labels_mat가 W를 만든 row순서와 완전히 동일하다면 cbind로도 됨)
# 여기서는 "좌표 키로 merge"하는 안전 버전 제공

# labels_mat가 dt_fpca와 같은 row order일 때: (가장 빠름)
# dt_fpca <- cbind(dt_fpca, as.data.table(labels_mat))

# ---- 안전 merge: labels에 좌표를 붙여서 join ----
labels_dt <- cbind(
  dt_fpca[, .(i, j, k)],          # 좌표를 labels에 붙여 key 생성
  as.data.table(labels_mat)
)

setkey(dt_fpca, i, j, k)
setkey(labels_dt, i, j, k)

dt_fpca2 <- labels_dt[dt_fpca]    # dt_fpca의 row 순서 유지하며 라벨 붙임
stopifnot(nrow(dt_fpca2) == nrow(dt_fpca))

# 다시 리스트에 저장
fmri_smooth_fpca[[1]] <- dt_fpca2


sum(is.na(dt_fpca2[[colnames(labels_mat)[1]]]))

# 3D plotting
library(htmlwidgets)

out_dir <- "/home/hyunseok/plot"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 너무 많은 점이면 html이 무거워질 수 있어서 (선택) 샘플링 옵션
# sample_n <- 200000L
sample_n <- NULL

dt_plot <- dt_fpca2
if (!is.null(sample_n) && nrow(dt_plot) > sample_n) {
  set.seed(1)
  dt_plot <- dt_plot[sample(.N, sample_n)]
}

label_cols <- colnames(labels_mat)

for (lc in label_cols) {
  message("Plotting: ", lc)

  p <- plot_voxels_by_cluster(
    dt_plot,
    cluster_col = lc,
    i_col = 1, j_col = 2, k_col = 3,  # dt_plot에서 i,j,k가 1,2,3 컬럼이면
    showlegend = FALSE,
    point_size = 5
  )

  outfile <- file.path(out_dir, sprintf("fpca_voxels_%s.html", lc))
  saveWidget(
    widget = p,
    file = outfile,
    selfcontained = FALSE,
    libdir = file.path(out_dir, "libs")
  )

  message("Saved: ", outfile)
}


table(dt_fpca2$cl_k0002)
table(dt_fpca2$cl_k0004)
table(dt_fpca2$cl_k0008)
table(dt_fpca2$cl_k0016)
table(dt_fpca2$cl_k0032)
table(dt_fpca2$cl_k0064)
table(dt_fpca2$cl_k0128)
table(dt_fpca2$cl_k0256)
table(dt_fpca2$cl_k0512)
table(dt_fpca2$cl_k1024)
table(dt_fpca2$cl_k2048)
table(dt_fpca2$cl_k4096)

mod_vec <- sapply(k_levels, function(k) {
  memb <- cut_at(fg, no = k)           # 길이 = vcount(g)
  modularity(g, membership = memb, weights = E(g)$weight)
})

mod_df <- data.frame(k = k_levels, modularity = mod_vec)

print(mod_df)
