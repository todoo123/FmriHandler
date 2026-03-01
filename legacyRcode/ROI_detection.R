# install.packages("ClusterR")

# ################################################################################
# readme
# author: HyunseokYoon
# ################################################################################
# 해당 파일은 기본적인 correlation ,partial correlation, srvf, l2 distance 기반의 유사도를 이용하여
# raw data 에 대한 ROI detection 을 수행한 코드파일입니다.
# clustering 방법은 Creddock(2011) 방법을 기반으로 하였음.
# 이외에 추가적인 data smoothing 은 사용되지 않았음. 
# gaussian kernel, basis expansion 을 이용한 smoothing 은 각각
# 1. ROI_detection_gaussian_example.R 
# 2. ROI_detection_bspline_example.R 
# 에서 수행하였음.
# 위의 1,2 파일은 smoothing 만을 수행하고, 현재 파일에 있는 clustering function 
# 을 이용하여 clustering 을 수행함.
# ################################################################################




# rm(list = ls())
# ################################################################################
# # data load
# ################################################################################
# setwd('/Users/hyunseokyoon/Desktop/학교/대학원/Project/인문사회융합연구/code')
# setwd('/Users/hyunseokyoon/Desktop/학교/대학원/Project/인문사회융합연구/sample_data')
# load('/Users/hyunseokyoon/Desktop/학교/대학원/Project/인문사회융합연구/code/fmri_list')
# files <- sprintf("%03d.csv", 1:15)
# fmri_list <- list()
# for(i in 1:15){
#   fmri_list[[i]] <- read.csv(files[i])
# }
# 
# ################################################################################
# # data preprocessing
# ################################################################################
# 
# for(i in 1:15){
#   print(dim(fmri_list[[i]]))
# }
# 
# for(i in 1:15){
#   print(length(unique(fmri_list[[i]]$Label)))
# }
# 
# for(i in 1:15){
#   print(length(unique(fmri_list[[i]]$i)))
# }
# 
# for(i in 1:15){
#   print(length(unique(fmri_list[[i]]$j)))
# }
# 
# for(i in 1:15){
#   print(length(unique(fmri_list[[i]]$k)))
# }
# 
# # voxel 중, 어떤 사람에게는 측정되었는데, 어떤 사람에게는 측정되지 않은 것이 있으므로, 모두가 공유하는 voxel 만 filtering
# # 하기로 함
# for (idx in seq_along(fmri_list)) {
#   fmri_list[[idx]]$ijk <- paste(fmri_list[[idx]]$i,
#                                 fmri_list[[idx]]$j,
#                                 fmri_list[[idx]]$k,
#                                 sep = "_")
# }
# 
# common_ijk <- Reduce(intersect, lapply(fmri_list, function(df) df$ijk))
# 
# for (idx in seq_along(fmri_list)) {
#   fmri_list[[idx]] <- fmri_list[[idx]][fmri_list[[idx]]$ijk %in% common_ijk, ]
# }
# 
# for(i in 1:15){
#   print(dim(fmri_list[[i]]))
# }
# 
# ################################################################################
# # data visualization
# ################################################################################
# # 개인 별 voxel fmri signal 이 어떤지 간략하게 visualization
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# dim(fmri_list[[1]])
# 
# # voxel index 선택 (예: 10번째 voxel)
# v_idx <- 30000
# # 리스트 안에서 모든 사람의 해당 voxel 시그널 꺼내기
# signal_list <- list()
# 
# for (id in seq_along(fmri_list)) {
#   df <- fmri_list[[id]]
#   
#   # 해당 voxel (v_idx 행), 시간별 값(6~485 열) 추출
#   signal <- as.numeric(as.character(unlist(df[v_idx, 6:483])))
#   
#   # 데이터프레임 생성
#   signal_df <- data.frame(
#     time   = 1:length(signal),
#     signal = signal,
#     person = paste0("P", id)
#   )
#   
#   # 리스트에 저장
#   signal_list[[id]] <- signal_df
# }
# 
# signal_df <- dplyr::bind_rows(signal_list)
# 
# ggplot(signal_df, aes(x = time, y = signal, color = person, group = person)) +
#   geom_line(alpha = 0.5) +
#   theme_minimal() + 
#   labs(
#     title = paste0('voxel: ', v_idx),
#     x = "Time",
#     y = "Signal"
#   )
# 
# # 결론: white noise 에 가까워 보인다. 
# # 그래서 질문은, 혹시 이런 신호를 처리하기 전에 해당 신호를 좀 평활화하거나 차분하거나 하는 절차가 필요한지?

################################################################################
# ROI detection: using parcellation
################################################################################

################################################################################
## 1: making weighted adjacency matrix
################################################################################

# ---- Packages ----
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(matrixStats)
  library(fdasrvf)
  library(future.apply)
  library(ClusterR)
})


# ------------------------------------------------------------------------------
# 함수명: make_26_offsets
# 역할  : (0,0,0) 기준으로 26 neighborhood 3d 좌표를 반환
# 입력  : 
#   -
# 출력  :
#   - data.table with columns (dx, dy, dz)
# 주의  : 
#   -
# ------------------------------------------------------------------------------
make_26_offsets <- function() {
  offs <- as.data.table(expand.grid(dx = -1:1, dy = -1:1, dz = -1:1))
  offs <- offs[!(dx == 0 & dy == 0 & dz == 0)]
  # 26 rows, columns: dx, dy, dz
  offs[]
}


# ------------------------------------------------------------------------------
# 함수명: build_edges_26
# 역할  : voxel 좌표를 이용해 26-이웃 edge 리스트 생성
# 입력  :
#   - df    : data.frame (i,j,k 좌표 포함)
#   - i_col : i 좌표 열 인덱스 (기본값=3)
#   - j_col : j 좌표 열 인덱스 (기본값=4)
#   - k_col : k 좌표 열 인덱스 (기본값=5)
# 출력  :
#   - data.table with columns (i_idx, j_idx), i_idx < j_idx
# 주의  : 중복 엣지 제거, 1-based 인덱스 반환
# ------------------------------------------------------------------------------
build_edges_26 <- function(df, i_col = 3, j_col = 4, k_col = 5) {
  # data check
  stopifnot(ncol(df) >= max(i_col, j_col, k_col))
  dt <- as.data.table(df)
  # Keep only coordinate + row index
  coords <- dt[, .(i = as.integer(.SD[[1]]),
                   j = as.integer(.SD[[2]]),
                   k = as.integer(.SD[[3]])),
               .SDcols = c(i_col, j_col, k_col)]
  coords[, idx := .I]
  setkey(coords, i, j, k)
  
  offs <- make_26_offsets()
  
  # For each offset, shift coords and join back to original coords to find neighbors
  edge_list <- rbindlist(lapply(seq_len(nrow(offs)), function(r) {
    dx <- offs$dx[r]; dy <- offs$dy[r]; dz <- offs$dz[r]
    
    # 원래 좌표를 (dx,dy,dz) 만큼 이동
    left <- coords[, .(i = i + dx, j = j + dy, k = k + dz, idxL = idx)]
    
    # 이동된 좌표(left)와 실제 좌표(coords)를 (i,j,k) 기준으로 매칭
    m <- left[coords, on = .(i, j, k), nomatch = 0L]
    
    # 매칭된 결과에서 node ID 쌍 만들기
    m <- m[, .(i_idx = pmin(idxL, idx), j_idx = pmax(idxL, idx))]
    
    # 중복 제거, 자기자신 제거
    unique(m[i_idx < j_idx])
  }), use.names = TRUE)
  
  # Deduplicate across all offsets
  # output: 26 neighbor 에 해당하는 node index
  return(unique(edge_list, by = c("i_idx", "j_idx")))
}

# ------------------------------------------------------------------------------
# 함수명: row_standardize
# 역할  : 각 voxel(time series row)을 평균 0, 표준편차 1로 표준화
# 입력  :
#   - X : numeric matrix (N x T), 행 = voxel, 열 = time point
# 출력  :
#   - Z : 표준화된 행렬 (row-wise zero mean / unit variance)
# 주의  :
#   - NA 값은 무시하고 계산 (na.rm = TRUE)
#   - 분산이 0이거나 비유한 값은 1e-8로 대체하여 나눗셈 안정화
# ------------------------------------------------------------------------------
row_standardize <- function(X) {
  mu <- rowMeans(X, na.rm = TRUE)
  sd <- matrixStats::rowSds(X, na.rm = TRUE)
  sd[!is.finite(sd) | sd == 0] <- 1e-8
  Z <- sweep(X, 1, mu, FUN = "-")
  Z <- sweep(Z, 1, sd, FUN = "/")
  return(Z)
}

# ------------------------------------------------------------------------------
# 함수명: weight_corr_by_rows
# 역할  : voxel time series 행렬과 edge list를 이용해 엣지 가중치 계산
#          (기본 가중치 = Pearson correlation; 대규모 데이터에 대비해 chunk 처리)
# 입력  :
#   - X               : numeric matrix (N x T), 행 = voxel, 열 = time point
#   - edges           : data.table with columns (i_idx, j_idx)
#   - prestandardized : FALSE면 row-wise 표준화(row_standardize) 수행, TRUE면 입력을 그대로 사용
#   - chunk_size      : 한 번에 처리할 edge batch 크기 (기본=50,000)
# 출력  :
#   - w : numeric vector, 각 edge에 대응하는 Pearson correlation 값
# 주의  :
#   - X는 반드시 행 = voxel, 열 = time 구조여야 함
#   - row-wise 표준화 후 correlation = 두 행 벡터의 내적 / (T-1)
#   - 메모리 절약과 속도 개선을 위해 edge들을 chunk 단위로 분할 계산
# ------------------------------------------------------------------------------
weight_corr_by_rows <- function(X, edges, prestandardized = FALSE, chunk_size = 50000L) {
  # vectorized correlation coefficient computation function
  
  stopifnot(is.matrix(X))
  Tlen <- ncol(X)
  if (!prestandardized) {
    Z <- row_standardize(X)
  } else {
    Z <- X
  }
  nE <- nrow(edges)
  w <- numeric(nE)
  # correlation of two standardized row vectors = dot / (T-1)
  denom <- (Tlen - 1L)
  idx_seq <- seq.int(1L, nE, by = chunk_size)
  for (start in idx_seq) {
    end <- min(start + chunk_size - 1L, nE)
    ids <- start:end
    A <- Z[edges$i_idx[ids], , drop = FALSE]
    B <- Z[edges$j_idx[ids], , drop = FALSE]
    w[ids] <- rowSums(A * B) / denom
  }
  return(w)
}

# ------------------------------------------------------------------------------
# 함수명: dist_to_weight
# 역할  : 거리 벡터 d를 엣지 가중치(weight)로 변환
#          (거리 클수록 weight 작아짐; kernel 또는 변환 적용)
# 입력  :
#   - d      : numeric vector of distances
#   - method : 변환 방식 선택 ("rbf" = Gaussian kernel, "inv" = inverse transform, "none" = 그대로 사용)
#   - sigma  : rbf 변환의 scale 파라미터 (NULL 시 pilot 표본으로 추정)
#   - scale  : inv 변환의 scale 파라미터 (NULL 시 pilot 표본으로 추정)
#   - pilot  : NULL이면 d 일부를 샘플링하여 sigma/scale 추정에 사용
# 출력  :
#   - w : numeric vector, 변환된 weight
# 주의  :
#   - sigma/scale이 NULL이면 pilot 표본에서 양수값의 중앙값(median)으로 추정
#   - 추정값이 유효하지 않으면 기본값 1.0 사용
#   - method="none"은 distance를 그대로 반환하므로 일반적인 similarity 용도에는 비권장
# ------------------------------------------------------------------------------
dist_to_weight <- function(d, method = c("rbf","inv","none"),
                           sigma = NULL, scale = NULL,
                           pilot = NULL) {
  # distance 값을 weight 로 사용할 수 있도록 해줌 - distance 크면 weight 는 작게
  method <- match.arg(method)
  if (method == "rbf") {
    # gaussian kernel
    if (is.null(sigma)) {
      # 파일럿 표본의 중앙값으로 스케일 추정
      if (is.null(pilot)) pilot <- d[sample.int(length(d), min(10000L, length(d)))]
      sigma <- stats::median(pilot[pilot > 0], na.rm = TRUE)
      if (!is.finite(sigma) || sigma <= 0) sigma <- 1.0
    }
    w <- exp(-(d^2) / (2 * sigma^2))
  } else if (method == "inv") {
    # inverse transform
    if (is.null(scale)) {
      if (is.null(pilot)) pilot <- d[sample.int(length(d), min(10000L, length(d)))]
      scale <- stats::median(pilot[pilot > 0], na.rm = TRUE)
      if (!is.finite(scale) || scale <= 0) scale <- 1.0
    }
    w <- 1 / (1 + d / scale)
  } else {
    w <- d  # 비권장: 유사도가 아니라 거리 그대로
  }
  return(w)
}


# ------------------------------------------------------------------------------
# 함수명: weight_pcorr_by_rows
# 역할  : 행=voxel, 열=time 행렬 X와 edge list로 "리지 부분상관" 가중치 계산
#          - Woodbury(우드버리) 항등식으로 T×T 선형계만 풀어 고차원(N>>T)에서도 안전
#          - edges를 chunk로 나누어 메모리 사용 제어
# 입력  :
#   - X               : numeric matrix (N x T), 행=voxel(변수), 열=time
#   - edges           : data.table/data.frame with (i_idx, j_idx) [1-based]
#   - prestandardized : FALSE면 행 표준화(평균0, 표준편차1) 수행
#   - lambda          : ridge 강도(>0). 기본=1.0 (행 표준화 후 대각≈1 기준)
#   - chunk_size      : 한 번에 계산할 edge 개수(기본=50,000)
#   - mem_saving      : TRUE면 각 chunk에서 등장하는 node만 부분해 풀어 메모리 절약
#                       FALSE면 전체 S=solve(B,t(U))를 한 번에 만들어 더 빠르게 계산
#   - eps             : 수치 안정성용 작은 값
# 출력  :
#   - w : numeric vector, 각 edge(i,j)의 리지 부분상관 추정치 ρ̃_{ij} ∈ [-1,1]
#
# 수식 요약:
#   Z = row-standardize(X), U = Z / sqrt(T-1)  (N x T)
#   B = I_T + (1/λ) * (U^T U)   (T x T, PD)
#   S = solve(B, t(U))          (T x N)  [mem_saving=FALSE일 때 한 번에]
#   a_i = u_i^T (solve(B) u_i),  b_ij = u_i^T (solve(B) u_j)
#   ρ̃_{ij} = [ (1/λ) * b_ij ] / sqrt( (1 - (1/λ) a_i) (1 - (1/λ) a_j) )
# ------------------------------------------------------------------------------
weight_pcorr_by_rows <- function(X, edges,
                                 prestandardized = TRUE,
                                 lambda = 1.0,
                                 chunk_size = 50000L,
                                 mem_saving = FALSE,
                                 eps = 1e-12) {
  stopifnot(is.matrix(X), lambda > 0, is.data.frame(edges),
            all(c("i_idx","j_idx") %in% names(edges)))
  
  N <- nrow(X); Tlen <- ncol(X)
  if (Tlen < 2L) stop("시간축 길이(T)는 2 이상이어야 합니다.")
  if (nrow(edges) == 0L) return(numeric(0L))
  
  # -- 행 표준화 (평균 0, 표준편차 1)
  .row_standardize <- function(A, eps = 1e-12) {
    mu <- rowMeans(A)
    A0 <- A - mu
    sd <- sqrt(pmax(rowSums(A0 * A0) / (ncol(A) - 1L), eps))
    A0 / sd
  }
  Z <- if (!prestandardized) .row_standardize(X, eps) else X
  
  # -- U, Woodbury 준비: B = I_T + (1/λ) U^T U  (T x T)
  U <- Z / sqrt(Tlen - 1L)             # N x T
  L <- 1.0 / lambda
  B <- diag(Tlen) + L * crossprod(U)   # crossprod(U) = t(U) %*% U (T x T)
  cholB <- chol(B)                      # 안정적 촐레스키
  
  nE <- nrow(edges)
  w <- numeric(nE)
  w[] <- NA_real_
  
  if (!mem_saving) {
    # ---------- 전체 S = solve(B, t(U))를 한 번에 구해 캐시 (빠름, 메모리↑) ----------
    # S_full: T x N, ST: N x T
    S_full <- backsolve(cholB, forwardsolve(t(cholB), t(U)))
    ST <- t(S_full)
    # a_vec: 길이 N, a_i = u_i^T solve(B) u_i
    a_vec <- rowSums(U * ST)
    one_minus_La <- pmax(1 - L * a_vec, eps)
    
    idx_seq <- seq.int(1L, nE, by = chunk_size)
    for (start in idx_seq) {
      end <- min(start + chunk_size - 1L, nE)
      ids <- start:end
      ii <- edges$i_idx[ids]
      jj <- edges$j_idx[ids]
      
      # 동일 노드 방지(있다면 0 처리)
      same <- which(ii == jj)
      if (length(same)) {
        w[ids][same] <- 0
        keep <- setdiff(seq_along(ii), same)
      } else {
        keep <- seq_along(ii)
      }
      if (length(keep) == 0L) next
      
      ii_k <- ii[keep]; jj_k <- jj[keep]
      
      # b_ij = sum_k U[i,k] * ST[j,k]
      bij <- rowSums(U[ii_k, , drop = FALSE] * ST[jj_k, , drop = FALSE])
      denom <- sqrt(pmax(one_minus_La[ii_k], eps) * pmax(one_minus_La[jj_k], eps))
      w_chunk <- (L * bij) / denom
      w_chunk <- pmin(pmax(w_chunk, -1), 1)
      w[ids][keep] <- w_chunk
    }
    
  } else {
    # ---------- 메모리 절약 모드: 각 chunk의 unique node만 부분해 풀기 (메모리↓) ----------
    idx_seq <- seq.int(1L, nE, by = chunk_size)
    for (start in idx_seq) {
      end <- min(start + chunk_size - 1L, nE)
      ids <- start:end
      ii <- edges$i_idx[ids]
      jj <- edges$j_idx[ids]
      
      # 동일 노드 방지
      same <- which(ii == jj)
      if (length(same)) {
        w[ids][same] <- 0
        keep <- setdiff(seq_along(ii), same)
      } else {
        keep <- seq_along(ii)
      }
      if (length(keep) == 0L) next
      
      ii_k <- ii[keep]; jj_k <- jj[keep]
      nodes <- unique(c(ii_k, jj_k))
      
      # U_sub: k x T, solve(B, t(U_sub)): T x k  ->  ST_sub: k x T
      U_sub <- U[nodes, , drop = FALSE]
      S_sub <- backsolve(cholB, forwardsolve(t(cholB), t(U_sub)))   # T x k
      ST_sub <- t(S_sub)                                            # k x T
      
      # a_sub: 길이 k
      a_sub <- rowSums(U_sub * ST_sub)
      one_minus_La_sub <- pmax(1 - L * a_sub, eps)
      
      # 원래 인덱스 → sub 인덱스 매핑
      pos_i <- match(ii_k, nodes)
      pos_j <- match(jj_k, nodes)
      
      bij <- rowSums(U_sub[pos_i, , drop = FALSE] * ST_sub[pos_j, , drop = FALSE])
      denom <- sqrt(pmax(one_minus_La_sub[pos_i], eps) * pmax(one_minus_La_sub[pos_j], eps))
      w_chunk <- (L * bij) / denom
      w_chunk <- pmin(pmax(w_chunk, -1), 1)
      w[ids][keep] <- w_chunk
    }
  }
  
  # 안정화
  w[!is.finite(w)] <- NA_real_
  return(w)
}

# ------------------------------------------------------------------------------
# 함수명: weight_cosine_ts
# 역할  : voxel time series 간 코사인 유사도 기반 엣지 가중치 계산
#          (옵션에 따라 표준화 적용 후 cosine 계산, (cos+1)/2 변환으로 [0,1] 매핑)
# 입력  :
#   - X           : numeric matrix (N x T), 행 = voxel, 열 = time point
#   - edges       : data.table/data.frame with columns (i_idx, j_idx)
#   - standardize : TRUE면 각 row(time series)를 평균0/표준편차1로 표준화 후 코사인 계산
#   - chunk_size  : 한 번에 처리할 edge batch 크기 (기본=50,000)
# 출력  :
#   - numeric vector: 각 edge에 대응하는 weight 값
# 주의  :
#   - 코사인 유사도 cos(x,y) = (x·y) / (||x|| ||y||)
#   - (cos+1)/2 로 [0,1]로 매핑 → 0(완전 반대) ~ 1(완전 같은 방향)
#   - 표준화 옵션이 TRUE일 경우, 분산 0 또는 NaN/Inf는 1e-8로 보정
#   - 수치 안정화를 위해 분모에 eps를 더하고, 결과를 [0,1]로 클램핑
#   - 이 함수는 "엣지별 벡터"를 반환합니다. 희소 행렬/가중치 행렬 구성은 별도 단계에서 처리하세요.
# ------------------------------------------------------------------------------
weight_cosine_ts <- function(X, edges,
                             standardize = FALSE,
                             chunk_size = 50000L) {
  # --- 표준화(행 단위: 각 시계열 평균0/표준편차1) ---
  if (standardize) {
    mu <- rowMeans(X, na.rm = TRUE)
    sd <- matrixStats::rowSds(X, na.rm = TRUE); sd[!is.finite(sd) | sd == 0] <- 1e-8
    X <- sweep(sweep(X, 1, mu, "-"), 1, sd, "/")
  }
  
  nE <- nrow(edges)
  w  <- numeric(nE)
  eps <- 1e-12  # 분모 0 방지
  
  # --- 청크 단위로 처리 (메모리 절약) ---
  for (s in seq.int(1L, nE, by = chunk_size)) {
    e <- min(s + chunk_size - 1L, nE)
    A <- X[edges$i_idx[s:e], , drop = FALSE]
    B <- X[edges$j_idx[s:e], , drop = FALSE]
    
    # 내적 및 노름
    dot   <- rowSums(A * B)                          # x·y
    nA2   <- rowSums(A * A)                          # ||x||^2
    nB2   <- rowSums(B * B)                          # ||y||^2
    denom <- sqrt(pmax(nA2, 0)) * sqrt(pmax(nB2, 0)) # ||x||·||y||
    denom <- denom + eps
    
    cosv  <- dot / denom                             # [-1,1] 이론상
    # 수치 안전 클램핑
    cosv  <- pmin(pmax(cosv, -1.0), 1.0)
    
    # (cos + 1)/2 → [0,1]
    w[s:e] <- 0.5 * (cosv + 1.0)
  }
  
  # 최종 안전 클램핑
  w <- pmin(pmax(w, 0.0), 1.0)
  return(w)
}

# ------------------------------------------------------------------------------
# 함수명: weight_l2_ts
# 역할  : voxel time series 간 L2 거리 기반 엣지 가중치 계산
#          (옵션에 따라 표준화 적용 후 L2 계산, kernel 변환으로 weight 변환)
# 입력  :
#   - X           : numeric matrix (N x T), 행 = voxel, 열 = time point
#   - edges       : data.table with columns (i_idx, j_idx)
#   - standardize : TRUE면 각 row(time series)를 평균0/표준편차1로 표준화 후 L2 계산
#   - kernel      : 거리→가중치 변환 방식 ("rbf", "inv", "none")
#   - sigma       : rbf kernel의 scale 파라미터 (NULL 시 pilot 표본으로 추정)
#   - scale       : inv 변환의 scale 파라미터 (NULL 시 pilot 표본으로 추정)
#   - chunk_size  : 한 번에 처리할 edge batch 크기 (기본=50,000)
# 출력  :
#   - numeric vector: 각 edge에 대응하는 weight 값
# 주의  :
#   - L2 distance = ||xi - xj||₂ / √T 로 계산 (시계열 길이에 따른 스케일 보정)
#   - 표준화 옵션이 TRUE일 경우, 분산 0 또는 NaN/Inf는 1e-8로 보정
#   - 큰 네트워크에서도 메모리 효율을 위해 chunk 단위로 edge를 나눠서 처리
# ------------------------------------------------------------------------------
weight_l2_ts <- function(X, edges,
                         standardize = TRUE,
                         kernel = c("rbf","inv","none"),
                         sigma = NULL, scale = NULL,
                         chunk_size = 50000L) {
  if (standardize) {
    mu <- rowMeans(X, na.rm = TRUE)
    sd <- matrixStats::rowSds(X, na.rm = TRUE); sd[!is.finite(sd) | sd == 0] <- 1e-8
    X <- sweep(sweep(X, 1, mu, "-"), 1, sd, "/")
  }
  nE <- nrow(edges)
  d <- numeric(nE)
  # L2 거리 = ||xi - xj||_2 / sqrt(T)  (스케일 공정화를 위해 √T로 나눔)
  # 따라서 엄밀한 의미로는 L2 distance 가 아니고, 
  # RMS 거리 = L2 거리 / sqrt(T) (길이 T에 따른 스케일 보정) 라고 불러야 한다.
  # 이는 시계열의 길이에 따른 편차를 보정하기 위해서인데, 
  # fMRI data 는 해당 처리를 했으므로 표준화의 영향은 무시해도 된다.
  Tlen <- ncol(X); denom <- sqrt(Tlen)
  for (s in seq.int(1L, nE, by = chunk_size)) {
    e <- min(s + chunk_size - 1L, nE)
    A <- X[edges$i_idx[s:e], , drop = FALSE]
    B <- X[edges$j_idx[s:e], , drop = FALSE]
    d[s:e] <- sqrt(rowSums((A - B)^2)) / denom
  }
  return(dist_to_weight(d, method = match.arg(kernel), sigma = sigma, scale = scale))
}

# ------------------------------------------------------------------------------
# 함수명: weight_srvf_ts
# 역할  : SRVF(Square Root Velocity Function) 기반 시계열 유사도/가중치 계산
#          - warp="none"    : SRVF 간 단순 RMS L2 거리 (빠름)
#          - warp="dtw"     : DTW 정합 기반 RMS 거리 (근사, 매우 느림)
#          - warp="elastic" : fdasrvf 기반 SRVF-elastic 정합 (정석, 매우 느림)
# 입력  :
#   - X          : numeric matrix (N x T), 행 = voxel, 열 = time point
#   - edges      : data.table with columns (i_idx, j_idx)
#   - dt         : time step (기본=1.0), 도함수 계산 시 사용
#   - smooth     : TRUE면 3-point 이동평균 사전 스무딩
#   - warp       : warping 방식 선택 ("none","dtw","elastic")
#   - kernel     : 거리→가중치 변환 ("rbf","inv","none")
#   - sigma/scale: kernel 파라미터 (NULL이면 pilot 표본으로 추정)
#   - eps        : SRVF 안정화 상수 (기본=1e-8)
#   warp = "elastic" 병렬화 parameter
#   - chunk_size : warp="none" 시 배치 크기 (기본=20000)
#   - parallel   : warp="elastic" 시 병렬 옵션 ("none","multisession","multicore")
#   - workers    : 병렬 워커 수 (기본=가용 코어-1)
#   - edge_chunk : warp="elastic" 시 edge를 나눌 청크 크기 (기본=200)
# 출력  :
#   - numeric vector: edge별 weight
# 주의  :
#   - warp="elastic"은 SRVF 이론에 충실하나 매우 느림 (에지마다 SRVF 정렬 수행)
#   - warp="dtw"는 근사 정합으로 SRVF 재매개변환 규칙(√γ̇)을 정확히 반영하지 않음
#   - 대규모 네트워크에서는 템플릿 정렬/샘플링 전략 병행 권장
# ------------------------------------------------------------------------------

weight_srvf_ts <- function(X, edges,
                           dt = 1.0,
                           smooth = FALSE,
                           warp = c("none","dtw","elastic"),
                           kernel = c("rbf","inv","none"),
                           sigma = NULL, scale = NULL,
                           eps = 1e-8,
                           chunk_size = 20000L,          # warp="none" 용
                           # ---- parallel options (elastic용) ----
                           parallel = c("none","multisession","multicore"),
                           workers = max(1L, parallel::detectCores() - 1L),
                           edge_chunk = 200L,            # elastic용: 에지 청크 크기
                           # ---- elastic.distance 옵션 ----
                           elastic_lambda = 0,
                           elastic_pen    = "roughness"
){
  warp    <- match.arg(warp)
  kernel  <- match.arg(kernel)
  parallel <- match.arg(parallel)
  
  # (선택) 간단 이동평균 스무딩
  if (smooth) {
    ma1 <- function(v) stats::filter(v, rep(1/3,3), sides = 2, circular = FALSE)
    X <- t(apply(X, 1, function(r) { z <- as.numeric(ma1(r)); z[is.na(z)] <- r[is.na(z)]; z }))
  }
  
  nE <- nrow(edges)
  d  <- numeric(nE)
  
  # --- warp 별 계산 ---
  if (warp %in% c("none","dtw")) {
    # 도함수(중심 차분) → SRVF
    diff_center <- function(v) {
      Tlen <- length(v)
      d <- c(v[2] - v[1],
             (v[3:Tlen] - v[1:(Tlen-2)])/2,
             v[Tlen] - v[Tlen-1])
      d / dt
    }
    D <- t(apply(X, 1, diff_center))             # N x (T-1)
    Q <- D / sqrt(abs(D) + eps)                  # SRVF
    
    if (warp == "none") {
      # SRVF 공간 RMS L2 거리
      Tq <- ncol(Q); denom <- sqrt(Tq)
      for (s in seq.int(1L, nE, by = chunk_size)) {
        e <- min(s + chunk_size - 1L, nE)
        Ai <- Q[edges$i_idx[s:e], , drop = FALSE]
        Aj <- Q[edges$j_idx[s:e], , drop = FALSE]
        d[s:e] <- sqrt(rowSums((Ai - Aj)^2)) / denom
      }
    } else {  # warp == "dtw"
      if (!requireNamespace("dtw", quietly = TRUE)) {
        stop("warp='dtw'는 'dtw' 패키지가 필요합니다. install.packages('dtw')")
      }
      for (s in seq_len(nE)) {
        qi <- Q[edges$i_idx[s], ]; qj <- Q[edges$j_idx[s], ]
        al <- dtw::dtw(qi, qj, distance.only = FALSE)
        ii <- al$index1; jj <- al$index2
        d[s] <- sqrt(mean((qi[ii] - qj[jj])^2))
      }
    }
    
  } else {  # ---- warp == "elastic": fdasrvf::elastic.distance 사용 ----
    if (!requireNamespace("fdasrvf", quietly = TRUE)) {
      stop("warp='elastic'은 'fdasrvf' 패키지가 필요합니다. install.packages('fdasrvf')")
    }
    
    # RStudio/Windows에서 multicore 안정성 보정
    if (parallel != "none") {
      if (parallel == "multicore" && .Platform$OS.type == "windows") {
        warning("Windows에서는 multicore 미지원 → multisession으로 대체합니다.")
        parallel <- "multisession"
      }
      if (parallel == "multicore" && Sys.getenv("RSTUDIO") == "1") {
        warning("RStudio 환경에서는 multicore가 불안정합니다. multisession으로 전환합니다.")
        parallel <- "multisession"
      }
      strategy <- if (parallel == "none") "sequential" else parallel
      future::plan(strategy, workers = workers)
      # 중첩 스레딩 억제
      old_omp <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
      old_mkl <- Sys.getenv("MKL_NUM_THREADS", unset = NA)
      Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
      on.exit({
        if (!is.na(old_omp)) Sys.setenv(OMP_NUM_THREADS = old_omp) else Sys.unsetenv("OMP_NUM_THREADS")
        if (!is.na(old_mkl)) Sys.setenv(MKL_NUM_THREADS = old_mkl) else Sys.unsetenv("MKL_NUM_THREADS")
        future::plan("sequential")
      }, add = TRUE)
    }
    
    # 시간축(원시 X 기준)
    tvec <- seq(0, 1, length.out = ncol(X))
    # 에지 인덱스를 청크로 분할 → 워커 오버헤드 감소
    idx <- seq_len(nE)
    chunk_ids <- split(idx, ceiling(idx / edge_chunk))
    
    # 큰 객체를 명시 인자로 전달하고, 글로벌 캡처 차단
    dist_chunks <- if (parallel == "none") {
      lapply(chunk_ids, function(ids) {
        out <- numeric(length(ids))
        for (k in seq_along(ids)) {
          s <- ids[k]
          fi <- as.numeric(X[edges$i_idx[s], ])
          fj <- as.numeric(X[edges$j_idx[s], ])
          # SRVF-elastic 거리 (스칼라)
          out[k] <- fdasrvf::elastic.distance(fi, fj, time = tvec,
                                              lambda = elastic_lambda,
                                              pen    = elastic_pen)
        }
        out
      })
    } else {
      future.apply::future_lapply(
        X = chunk_ids,
        FUN = function(ids, X_, edges_, tvec_, elastic_lambda_, elastic_pen_) {
          out <- numeric(length(ids))
          for (k in seq_along(ids)) {
            s <- ids[k]
            fi <- as.numeric(X_[edges_$i_idx[s], ])
            fj <- as.numeric(X_[edges_$j_idx[s], ])
            out[k] <- fdasrvf::elastic.distance(fi, fj, time = tvec_,
                                                lambda = elastic_lambda_,
                                                pen    = elastic_pen_)
          }
          out
        },
        X_ = X,
        edges_ = edges,
        tvec_ = tvec,
        elastic_lambda_ = elastic_lambda,
        elastic_pen_    = elastic_pen,
        future.packages = "fdasrvf",
        future.globals  = FALSE,   # ★ 대형 글로벌 전송 방지
        future.seed     = TRUE
      )
    }
    
    # 원래 순서로 합치기
    d[unlist(chunk_ids, use.names = FALSE)] <- unlist(dist_chunks, use.names = FALSE)
  }
  
  # 거리 → 가중치 변환
  return(dist_to_weight(d, method = kernel, sigma = sigma, scale = scale))
}
# ------------------------------------------------------------------------------
# 함수명: postprocess_weights
# 역할  : edge weight 벡터에 후처리를 적용
#          - 음수값 보정 (0으로 클램핑)
#          - 임계값(r_min) 미만은 0으로 설정
# 입력  :
#   - w              : numeric vector (edge weights)
#   - clamp_negative : TRUE면 음수 weight를 0으로 보정
#   - r_min          : threshold 값 (NULL이면 skip)
# 출력  :
#   - 후처리된 numeric vector (edge weights)
# ------------------------------------------------------------------------------
postprocess_weights <- function(w, clamp_negative = TRUE, r_min = NULL) {
  # 음수 weight를 0으로 보정
  if (clamp_negative) 
    w <- pmax(w, 0)
  
  # 임계값 r_min 미만은 0으로 설정
  if (!is.null(r_min)) 
    w[w < r_min] <- 0
  
  return(w)
}


# ------------------------------------------------------------------------------
# 함수명: assemble_sparse_adj
# 역할  : edge list + weight 로부터 sparse weighted adjacency matrix 생성
#          - 유효한 weight만 남기고
#          - 상삼각(i<j)만 사용하여 대칭 행렬 구성
# 입력  :
#   - N     : 정점(노드, voxel) 개수
#   - edges : data.table (i_idx, j_idx)  # edge endpoint indices
#   - w     : numeric vector (edge weights; length = nrow(edges))
# 출력  :
#   - generalMatrix 클래스 sparse weighted adjacency matrix (N x N, symmetric)
# ------------------------------------------------------------------------------
assemble_sparse_adj <- function(N, edges, w) {
  # 1) 유효한 weight만 필터링 (유한값 & 0이 아님)
  keep <- which(is.finite(w) & w != 0)
  if (length(keep) == 0L) {
    # 모두 0 또는 NA이면 영행렬 반환
    return(Matrix::Matrix(0, nrow = N, ncol = N, sparse = TRUE))
  }
  
  # edge endpoint index와 weight 추출
  i <- edges$i_idx[keep]
  j <- edges$j_idx[keep]
  x <- w[keep]
  
  # 2) 안전장치: 상삼각(i<j) edge만 유지
  keep_up <- which(i < j)
  i <- i[keep_up]; j <- j[keep_up]; x <- x[keep_up]
  
  # 3) sparse symmetric matrix 구성
  #    - 상삼각 입력만 주면 대칭으로 확장됨
  W <- Matrix::sparseMatrix(
    i = i, j = j, x = x,
    dims = c(N, N), symmetric = TRUE
  )
  
  # generalMatrix 형식으로 강제 변환
  W <- methods::as(W, "generalMatrix")
  
  return(W)
}


# ------------------------------------------------------------------------------
# 함수명: build_weighted_adj_26
# 역할  : voxel 시계열 데이터로부터 26-이웃 기반 weighted adjacency matrix 생성
#          - voxel 좌표에서 edge list를 만들고, 사용자 지정 weight 함수로 가중치 계산
# 입력  :
#   - df             : data.frame (voxel 데이터; 좌표 + 시계열 포함)
#   - time_cols      : 시계열 데이터 열 인덱스 (기본=6:(ncol(df)-2))
#   - i_col, j_col, k_col : voxel 좌표 열 인덱스 (기본=3,4,5)
#   - weight_fun     : 가중치 계산 함수 (기본=weight_l2_ts; 예: weight_pcorr_by_rows, weight_corr_by_rows, weight_srvf_ts 등)
#   - clamp_negative : TRUE면 음수 weight를 0으로 보정 (유사도 척도일 경우 불필요)
#   - r_min          : threshold; 이 값보다 작은 weight는 0으로 설정
#   - chunk_size     : weight_fun에 전달할 청크 크기 (대규모 edge 처리용)
#   - verbose        : TRUE면 진행 상황 메시지 출력
#   - ...            : weight_fun에 전달할 추가 인자 (예: kernel, warp 등)
# 출력  :
#   - sparse matrix (dgCMatrix): N × N weighted adjacency matrix
# 주의  :
#   - 26-neighborhood edge list는 build_edges_26()에 의해 생성됨
#   - weight_fun은 반드시 (X, edges, ...) 서명을 지원해야 함
#   - 후처리로 r_min threshold 적용; 필요시 clamp_negative도 가능
# ------------------------------------------------------------------------------
build_weighted_adj_26 <- function(df,
                                  time_cols = 6:(ncol(df) - 2),
                                  i_col = 3, j_col = 4, k_col = 5,
                                  weight_fun = weight_corr_ts,   # 기본: L2 기반
                                  clamp_negative = FALSE,      # 거리→유사도면 음수 없음
                                  r_min = NULL,
                                  chunk_size = 50000L,
                                  verbose = TRUE,
                                  ...) {                       # ← 가중치 함수로 추가 인자 전달
  if (verbose) message("Step 1/3: Building 26-neighborhood edge list...")
  edges <- build_edges_26(df, i_col, j_col, k_col)
  N <- nrow(df); if (verbose) message(sprintf("  Voxels: %d, edges: %d", N, nrow(edges)))
  
  if (verbose) message("Step 2/3: Extracting time-series matrix...")
  X <- as.matrix(as.data.frame(df)[, time_cols, drop = FALSE]); storage.mode(X) <- "double"
  
  if (verbose) message("Step 3/3: Computing weights via custom function ...")
  # 가중치 함수 호출: 필요 파라미터는 ...로 전달
  w_raw <- weight_fun(X, edges, chunk_size = chunk_size, ...)
  
  # (선택) 후처리: threshold
  w <- w_raw
  if (!is.null(r_min)) w[w < r_min] <- 0
  
  assemble_sparse_adj(N, edges, w)
}

# ------------------------------------------------------------------------------
# 헬퍼: (i_lo .. i_hi) 구간의 상삼각 엣지(i<j) 생성
#  - N^2 크기의 행렬을 만들지 않고 인덱스 수열만 생성
# warning: 현재 데이터에 대한 전체 수열 생성 시에 memory 과부하 오므로 서버 사용할 것
# ------------------------------------------------------------------------------
build_edges_block <- function(i_lo, i_hi, N) {
  stopifnot(1L <= i_lo, i_lo <= i_hi, i_hi < N)
  # i는 i_lo..i_hi, 각 i에 대해 (N-i)개의 j
  i_idx <- rep.int(i_lo:i_hi, times = (N - i_lo):(N - i_hi))
  # j는 i+1..N 구간을 차례대로 이어붙인 수열
  # 각 i마다 길이가 (N - i)
  len_per_i <- (N - i_lo):(N - i_hi)
  j_start    <- (i_lo + 1L):(i_hi + 1L)
  j_idx <- unlist(Map(function(js, len) seq.int(js, js + len - 1L),
                      js = j_start, len = len_per_i), use.names = FALSE)
  data.table::data.table(i_idx = i_idx, j_idx = j_idx)
}


# ------------------------------------------------------------------------------
# 함수명: build_weighted_adj_all_streaming
# 역할  : 모든 voxel 쌍을 'i 인덱스 블록' 단위로 스트리밍 처리하여
#         사용자 정의 weight_fun으로 가중치를 계산하고 희소 인접행렬 생성
# 인터페이스(중요):
#   - weight_fun(X, edges, chunk_size=..., ...) 형태여야 함
#     * X: (N x T) 시계열 행렬 (double)
#     * edges: data.table(i_idx, j_idx)  (1-based, i<j)
#     * chunk_size: (옵션) 내부에서 사용할 청크 크기
# 입력:
#   - df             : voxel 데이터프레임 (행=voxel)
#   - time_cols      : 시계열 열 인덱스
#   - weight_fun     : 가중치 계산 함수 (예: weight_corr_ts, weight_corr_by_rows 등)
#   - i_block_size   : i 인덱스를 이 크기만큼 묶어 스트리밍 처리
#   - r_min          : 임계값 미만은 0으로 컷
#   - clamp_negative : TRUE면 음수 가중치를 0으로 클램프
#   - chunk_size     : weight_fun에 그대로 전달
#   - verbose        : 진행 출력
#   - ...            : weight_fun에 전달할 추가 인자
# 반환:
#   - dgCMatrix (대칭 희소 인접행렬, 상삼각으로부터 대칭화)
# warning: 미완성
# ------------------------------------------------------------------------------
build_weighted_adj_all_streaming <- function(df,
                                             time_cols = 6:(ncol(df) - 2),
                                             weight_fun = weight_corr_ts,
                                             i_block_size = 2000L,   # 메모리/시간 트레이드오프
                                             r_min = NULL,
                                             clamp_negative = FALSE,
                                             chunk_size = 50000L,
                                             verbose = TRUE,
                                             ...) {
  N <- nrow(df)
  if (N < 2L) {
    if (verbose) message("N < 2: 빈 그래프 반환")
    return(Matrix::Matrix(0, nrow = N, ncol = N, sparse = TRUE))
  }
  
  # Step 1. 시계열 행렬 준비
  if (verbose) message("Step 1/3: Extracting time-series matrix...")
  X <- as.matrix(as.data.frame(df)[, time_cols, drop = FALSE])
  storage.mode(X) <- "double"
  
  # 누적 컨테이너(트리플릿)
  Ii_list <- list()
  Jj_list <- list()
  Xx_list <- list()
  chunk_id <- 0L
  total_kept <- 0L
  
  # Step 2. i-블록 스트리밍
  if (verbose) {
    message("Step 2/3: Streaming over i-blocks...")
    message(sprintf("  Voxels: %d, i_block_size: %d", N, i_block_size))
  }
  
  i_starts <- seq.int(1L, N - 1L, by = i_block_size)
  n_blocks <- length(i_starts)
  for (b in seq_along(i_starts)) {
    i_lo <- i_starts[b]
    i_hi <- min(i_lo + i_block_size - 1L, N - 1L)
    
    if (verbose) {
      message(sprintf("  Block %d/%d: i in [%d, %d]", b, n_blocks, i_lo, i_hi))
    }
    
    # 엣지 생성 (상삼각, i<j)
    edges_block <- build_edges_block(i_lo, i_hi, N)
    
    # 가중치 계산
    w_raw <- weight_fun(X, edges_block, chunk_size = chunk_size, ...)
    
    # 후처리(임계/음수클램프)
    w <- w_raw
    if (!is.null(r_min)) w[w < r_min] <- 0
    if (isTRUE(clamp_negative)) w <- pmax(w, 0)
    
    # 유효 가중치만 보관
    keep <- which(is.finite(w) & w != 0)
    if (length(keep)) {
      Ii_list[[length(Ii_list) + 1L]] <- edges_block$i_idx[keep]
      Jj_list[[length(Jj_list) + 1L]] <- edges_block$j_idx[keep]
      Xx_list[[length(Xx_list) + 1L]] <- w[keep]
      total_kept <- total_kept + length(keep)
    }
    
    # 메모리 완화를 위해 간헐적으로 리스트 컴팩트
    chunk_id <- chunk_id + 1L
    if (chunk_id %% 10L == 0L) {
      gc()
    }
  }
  
  # Step 3. 희소행렬 조립
  if (verbose) {
    message("Step 3/3: Assembling sparse adjacency...")
    message(sprintf("  Nonzero weights kept (upper triangle): %d", total_kept))
  }
  
  if (total_kept == 0L) {
    if (verbose) message("  No edges survived thresholding; returning zero matrix.")
    return(Matrix::Matrix(0, nrow = N, ncol = N, sparse = TRUE))
  }
  
  i_all <- do.call(c, Ii_list)
  j_all <- do.call(c, Jj_list)
  x_all <- do.call(c, Xx_list)
  
  # 상삼각으로 생성 후 대칭화
  A_upper <- Matrix::sparseMatrix(
    i = i_all, j = j_all, x = x_all,
    dims = c(N, N), symmetric = FALSE, giveCsparse = TRUE
  )
  A_sym <- Matrix::forceSymmetric(A_upper, uplo = "U")
  A_sym <- as(A_sym, "dgCMatrix")  # 표준 dgCMatrix로 캐스팅
  
  # (선택) 0차수 노드가 많다면 여기서 후처리 가능
  # d <- Matrix::rowSums(A_sym); ...
  
  if (verbose) {
    nnz <- length(A_sym@x)
    message(sprintf("  Done. nnz (symmetric): %d", nnz))
  }
  A_sym
}

build_weighted_adj_all_streaming_parallel <- function(
    df,
    time_cols     = 6:(ncol(df) - 2),
    weight_fun    = weight_corr_by_rows,
    i_block_size  = 2000L,
    r_min         = NULL,
    clamp_negative= FALSE,
    chunk_size    = 50000L,
    verbose       = TRUE,
    # --- 병렬 옵션 ---
    parallel      = c("none","multisession","multicore"),
    workers       = max(1L, parallel::detectCores() - 1L),
    # --- BLAS/OMP 억제 (권장) ---
    tame_blas     = TRUE,
    ...
){
  parallel <- match.arg(parallel)
  N <- nrow(df)
  if (N < 2L) return(Matrix::Matrix(0, nrow = N, ncol = N, sparse = TRUE))

  if (verbose) message("[W_allpairs|par] Extract X ...")
  X <- as.matrix(as.data.frame(df)[, time_cols, drop = FALSE])
  storage.mode(X) <- "double"

  i_starts <- seq.int(1L, N - 1L, by = i_block_size)
  blocks <- lapply(i_starts, function(i_lo) {
    i_hi <- min(i_lo + i_block_size - 1L, N - 1L)
    list(i_lo = i_lo, i_hi = i_hi)
  })

  if (verbose) message(sprintf("[W_allpairs|par] Blocks: %d (i_block_size=%d)", length(blocks), i_block_size))
  
  # <<< 추가: ...를 리스트로 캡처해서 워커에 명시 전달 >>>
  dots <- list(...)

  # 3) 병렬 세팅 (+ cgroups 우회 폴백 유지)
  if (parallel != "none") {
    if (tame_blas) {
      old_omp <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
      old_mkl <- Sys.getenv("MKL_NUM_THREADS", unset = NA)
      Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
      on.exit({
        if (!is.na(old_omp)) Sys.setenv(OMP_NUM_THREADS = old_omp) else Sys.unsetenv("OMP_NUM_THREADS")
        if (!is.na(old_mkl)) Sys.setenv(MKL_NUM_THREADS = old_mkl) else Sys.unsetenv("MKL_NUM_THREADS")
      }, add = TRUE)
    }
    
    plan_req <- if (parallel == "multicore" && .Platform$OS.type == "windows") "multisession" else parallel
    ok <- TRUE; err <- NULL
    tryCatch({
      future::plan(plan_req, workers = workers)
    }, error = function(e) { ok <<- FALSE; err <<- e })
    
    if (!ok) {
      if (grepl("CGroups|cgroups|cpuset", conditionMessage(err), ignore.case = TRUE)) {
        options(parallelly.cgroups.enabled = FALSE)
        Sys.setenv(R_PARALLELLY_CGROUPS = "FALSE")
      }
      ok <- TRUE; err <- NULL
      tryCatch({
        future::plan("multisession", workers = workers)
        print("multisession activated")
        print(workers)
      }, error = function(e2) { ok <<- FALSE; err <<- e2 })
    }
    
    if (!ok) {
      warning(sprintf("[parallel] plan() failed -> falling back to serial: %s",
                      conditionMessage(err)))
      parallel <- "none"
    } else {
      on.exit(future::plan("sequential"), add = TRUE)
    }
  }
  
  # 4) 블록 워커 (여기에 안전호출 로직을 인라인)
  block_worker <- function(bl, X_, weight_fun_, chunk_size_, r_min_, clamp_negative_, dots_) {
    # --- 로컬 안전 호출 ---
    safe_weight_call_local <- function(fun, X, edges, chunk_size, dots) {
      fn_formals <- names(formals(fun))
      pass <- dots[intersect(names(dots), fn_formals)]
      do.call(fun, c(list(X = X, edges = edges, chunk_size = chunk_size), pass))
    }
    
    i_lo <- bl$i_lo; i_hi <- bl$i_hi
    edges_block <- build_edges_block(i_lo, i_hi, nrow(X_))
    
    w_raw <- safe_weight_call_local(weight_fun_, X_, edges_block, chunk_size_, dots_)
    w <- w_raw
    if (!is.null(r_min_)) w[w < r_min_] <- 0
    if (isTRUE(clamp_negative_)) w <- pmax(w, 0)
    
    keep <- which(is.finite(w) & w != 0)
    if (!length(keep)) {
      list(i = integer(0), j = integer(0), x = numeric(0))
    } else {
      list(i = edges_block$i_idx[keep], j = edges_block$j_idx[keep], x = w[keep])
    }
  }

  # 5) 실행
  if (parallel == "none") {
    parts <- lapply(blocks, block_worker,
                    X_ = X, weight_fun_ = weight_fun,
                    chunk_size_ = chunk_size,
                    r_min_ = r_min, clamp_negative_ = clamp_negative,
                    dots_ = dots)

  } else {
    parts <- future.apply::future_lapply(
      blocks, block_worker,
      X_ = X, weight_fun_ = weight_fun,
      chunk_size_ = chunk_size,
      r_min_ = r_min, clamp_negative_ = clamp_negative,
      dots_ = dots,
      future.globals = c("build_edges_block"), # 워커로 내보낼 전역 함수 지정
      future.seed    = TRUE
    )
  }
  
  # 6) 희소행렬 조립
  i_all <- unlist(lapply(parts, `[[`, "i"), use.names = FALSE)
  j_all <- unlist(lapply(parts, `[[`, "j"), use.names = FALSE)
  x_all <- unlist(lapply(parts, `[[`, "x"), use.names = FALSE)
  
  if (!length(i_all)) {
    if (verbose) message("[W_allpairs|par] no nonzeros -> zero matrix")
    return(Matrix::Matrix(0, nrow = N, ncol = N, sparse = TRUE))
  }
  
  if (verbose) message(sprintf("[W_allpairs|par] nnz(upper): %d", length(x_all)))
  A_upper <- Matrix::sparseMatrix(i = i_all, j = j_all, x = x_all,
                                  dims = c(N, N), symmetric = FALSE, giveCsparse = TRUE)
  A_sym <- Matrix::forceSymmetric(A_upper, uplo = "U")
  as(A_sym, "dgCMatrix")
}

# 코드 사용 예시
# Assuming your data.frame is named `df` and has 45590 rows x 485 cols as described:
#   col1: label (char)
#   col2: PFC_flag (char, all "O")
#   col3: i (int)
#   col4: j (int)
#   col5: k (int)
#   col6:col483 time series (float)
#   col484: subject ID (int)
#   col485: ijk string (char)

# W_corr <- build_weighted_adj_26(
#        fmri_list[[1]],
#        time_cols = 6:483,       # or 6:(ncol(df)-2)
#        i_col = 3, j_col = 4, k_col = 5,
#        clamp_negative = TRUE,   # set negative correlations to 0
#        r_min = 0.0,             # or e.g. 0.2 to sparsify weak edges
#        chunk_size = 50000L,
#        verbose = TRUE
#      )
#      
# W_pcorr <- build_weighted_adj_26(
#        fmri_list[[1]],
#        time_cols = 6:483,       # or 6:(ncol(df)-2)
#        i_col = 3, j_col = 4, k_col = 5,
#        weight_fun = weight_l2_ts,
#        clamp_negative = TRUE,   # set negative correlations to 0
#        r_min = 0.0,             # or e.g. 0.2 to sparsify weak edges
#        chunk_size = 50000L,
#        verbose = TRUE
#      )
# 
# W_l2 <- build_weighted_adj_26(
#   fmri_list[[1]],
#   time_cols = 6:483,
#   i_col = 3, j_col = 4, k_col = 5,
#   weight_fun = weight_l2_ts,
#   kernel = "rbf",          # "rbf" | "inv" | "none"
#   sigma = NULL,            # NULL이면 파일럿에서 자동 추정
#   chunk_size = 50000L,
#   verbose = TRUE
# )
# 
# # 졸라 느림 - 절대 안돼;;;; 
# options(future.globals.maxSize = 3 * 1024^3)
# W_srvf <- build_weighted_adj_26(
#   fmri_list[[1]],
#   time_cols = 6:483,
#   i_col = 3, j_col = 4, k_col = 5,
#   weight_fun = weight_srvf_ts,
#   warp = "elastic",
#   parallel = "multisession",
#   workers = 7,
#   edge_chunk = 200,
#   elastic_lambda = 0,
#   elastic_pen = "roughness",
#   kernel = "rbf"
# )

# str(W)  # dgCMatrix, symmetric weighted adjacency
#
# ----- Custom weight (for future FDA) -----
# Define your own function with the same signature as `weight_corr_by_rows`:
# my_weight_fun <- function(X, edges, ...) {
#   # X: N x T matrix of time-series (or coefficients if you first project onto basis)
#   # edges: data.table with i_idx, j_idx
#   # return numeric vector of length nrow(edges)
# }
# Then pass `weight_fun = my_weight_fun` to build_weighted_adj_26().


################################################################################
## 2. spectral clustering
################################################################################
suppressPackageStartupMessages({
  library(Matrix)
  library(RSpectra)
})
# 병목
# ------------------------------------------------------------------------------
# 함수명: spectral_clustering_ncut
# 역할  : 정규화 인접행렬 S = D^{-1/2} W D^{-1/2}의 상위 K 고유벡터(Ng–Jordan–Weiss)
#         를 이용해 그래프의 스펙트럴 클러스터링(NCut)을 수행
# 입력  :
#   - W              : 대칭 희소 인접행렬 (dgCMatrix/dsCMatrix/sparseMatrix, N x N)
#   - K              : 원하는 클러스터 개수
#   - remove_isolated: TRUE면 차수 0(고립) 노드를 클러스터링에서 제외 후 라벨 0으로 복원
#   - kmeans_nstart  : k-means nstart (여러 초기값 시도)
#   - kmeans_itermax : k-means 최대 반복 횟수
#   - seed           : 난수 시드(고유분해/군집 결과 재현성)
#   - verbose        : 진행 메시지 출력 여부
# 출력  :
#   - list(
#       labels       : 길이 N의 정수 벡터 (1..K = 클러스터, 0 = 고립/미할당)
#       centers      : k-means 중심(K x K)
#       eigvec       : 활성 부분그래프의 상위 K 고유벡터 (n_active x K)
#       eigval       : 상위 K 고유값
#       active_index : 활성 노드 인덱스(차수>0)
#       S_diag       : S의 대각 일부(간단 검사용)
#       info         : 리스트(N, n_active, K)
#     )
# 주의  :
#   - 입력 W가 비대칭이면 상삼각 기준으로 forceSymmetric 적용
#   - S = D^{-1/2} W D^{-1/2}는 대칭/희소이며 최대 고유값은 1
#   - 연결 성분 수만큼 고유값=1이 존재 (K는 활성 노드 수 이상이어야 함)
#   - RSpectra::eigs_sym 호환을 위해 S를 dgCMatrix로 캐스팅
#   - Ng–Jordan–Weiss: 고유벡터 행 정규화 후 k-means 수행
# ------------------------------------------------------------------------------
spectral_clustering_ncut <- function(
    W,                      # symmetric dgCMatrix (N x N)
    K,                      # 원하는 클러스터 개수
    remove_isolated = TRUE, # 차수 0 노드 제거/표식
    kmeans_nstart = 10,
    kmeans_itermax = 100,
    seed = 42,
    verbose = TRUE,
    screeplot = FALSE,
    screeplot_name = NULL
){
  stopifnot(inherits(W, "dsCMatrix") || inherits(W, "dgCMatrix") || inherits(W, "sparseMatrix"))
  # 대칭 보장
  if (!isSymmetric(W)) {
    if (verbose) message("Input not symmetric; forceSymmetric(U) 적용")
    W <- forceSymmetric(W, uplo = "U")
  }
  W <- as(W, "dgCMatrix")
  
  N <- nrow(W)
  d <- Matrix::rowSums(W)
  idx_active <- which(d > 0)
  
  if (length(idx_active) < K) {
    stop(sprintf("활성 노드 수(%d)가 K(%d)보다 작습니다. K를 줄이거나 그래프 연결성을 늘려주세요.", length(idx_active), K))
  }
  
  # 고립 노드 처리
  if (remove_isolated && length(idx_active) < N) {
    if (verbose) message(sprintf("차수 0(고립) 노드 %d개 제외 후 진행", N - length(idx_active)))
  }
  
  # 활성 부분그래프
  Wsub <- W[idx_active, idx_active, drop = FALSE]
  dsub <- Matrix::rowSums(Wsub)
  
  # 정규화 인접행렬 S = D^{-1/2} W D^{-1/2}
  invsqrt_d <- 1 / sqrt(dsub)
  invsqrt_d[!is.finite(invsqrt_d)] <- 0
  Dmhalf <- Diagonal(x = invsqrt_d)
  # 희소 곱셈: S는 여전히 희소 & 대칭
  S <- Dmhalf %*% Wsub %*% Dmhalf
  S <- forceSymmetric(S, uplo = "U")
  
  # >>> 이 줄 추가! (dsCMatrix -> dgCMatrix로 변환)
  S <- as(S, "dgCMatrix")
  
  # 상위 고유벡터 K개 (가장 큰 고유값 기준). S의 최대 고유값은 1이며,
  # 연결 성분 수만큼 고유값=1이 존재합니다.
  if (verbose) message("RSpectra 고유분해 진행 (which='LA') ...")
  # set.seed(seed)
  eig <- RSpectra::eigs_sym(S, k = K, which = "LA")  # Largest Algebraic
  U   <- eig$vectors  # n_active x K
  
  # test
  if(screeplot == TRUE){
    source("/home/hyunseok/code/scree_plot.R")
    scree_plot_rspectra(
    eig,
    title = "Scree (% of positive-sum)",
    normalize = "possum",
    save_path = paste0("/home/hyunseok/plot/scree_possum_",screeplot_name,".png"),
    width = 1600, height = 1000, dpi = 200
    )
  }
  # test
  

  # 행별 정규화 (Ng–Jordan–Weiss)
  row_norms <- sqrt(rowSums(U^2))
  row_norms[row_norms == 0 | !is.finite(row_norms)] <- 1
  Y <- U / row_norms
  # k-means
  if (verbose) message("k-means 클러스터링 ...")
  # set.seed(seed)
  # km <- stats::kmeans(Y, centers = K, nstart = kmeans_nstart, algorithm = 'Hartigan-Wong', iter.max = kmeans_itermax)
  km <- ClusterR::MiniBatchKmeans(
    data        = Y,
    clusters    = K,
    batch_size  = 5000,
    max_iters   = kmeans_itermax,
    num_init    = kmeans_nstart,
    initializer = "kmeans++",
    verbose     = FALSE
  )
  pred <- predict_MBatchKMeans(
    data      = Y,
    CENTROIDS = km$centroids,
    fuzzy     = FALSE
  )
  # 전체 노드로 라벨 복원
  labels <- integer(N)
  labels[] <- 0L  # 0 = 고립 노드(또는 미할당)
  # labels[idx_active] <- km$cluster
  labels[idx_active] <- pred
  
  # 반환
  list(
    labels = labels,       # 길이 N, 0은 고립/미할당
    centers = km$centroids,  # k-means 중심 (K x p)
    eigvec = U,            # 활성 부분의 고유벡터 (n_active x K)
    eigval = eig$values,   # 고유값 K개
    active_index = idx_active,
    S_diag = diag(S)[1:min(5, nrow(S))], # 간단한 체크용
    info = list(
      N = N,
      n_active = length(idx_active),
      K = K
    )
  )
}


# 사용 예시
# W: 앞서 만든 weighted adjacency (dgCMatrix, symmetric)
# 권장 K: 150~200(해석 용이), 더 미세한 분할은 600~1000도 가능

# res <- spectral_clustering_ncut(
#   W_pcorr,
#   K = 100,
#   remove_isolated = TRUE,
#   kmeans_nstart = 10,
#   kmeans_itermax = 300,
#   seed = 2025,
#   verbose = TRUE
# )

# res <- spectral_clustering_ncut(
#   W_srvf,
#   K = 100,
#   remove_isolated = TRUE,
#   kmeans_nstart = 10,
#   kmeans_itermax = 100,
#   seed = 2025,
#   verbose = TRUE
# )
# 
# table(res$labels)         # 각 클러스터 크기(0은 고립 노드)
# head(res$eigval)          # 상위 고유값 확인 (1 근처 값이 성분 수만큼 나올 수 있음)
# 
# # 라벨을 원래 데이터프레임에 붙이기
# fmri_list[[1]]$cluster_group100_pcorr <- res$labels



################################################################################
## 3. visualization
################################################################################
# ------------------------------------------
# 3D voxel scatter by cluster / structural label
# 시각화 부분은 함수 주석을 따로 추가하지 않았습니다.
# ------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(plotly)
  library(grDevices)  # hcl.colors
})


# 팔레트: 범주 n개에 대해 구분 가능한 색 생성
# n<=12면 Set3 비슷한 톤, 그 이상이면 HCL 기반의 "Dynamic" 팔레트
.discrete_palette <- function(n) {
  if (n <= 12) {
    # Set3 느낌으로 hcl 기반 톤 매칭
    cols <- hcl.colors(max(12, n), palette = "Pastel1")
    cols[seq_len(n)]
  } else {
    hcl.colors(n, palette = "Dynamic")  # 많은 범주에 적당
  }
}

# 공통: plot_ly에서 빠르게 찍기 위한 유틸
# df: data.frame/data.table
# i_col, j_col, k_col: 좌표 열 이름(문자) 또는 인덱스
# color_fac: 색상에 쓸 팩터(길이 nrow(df))
# showlegend: 논리값(클러스터가 많으면 FALSE 권장)
# point_size: 점 크기
# .plot_voxels_3d <- function(df, i_col = 3, j_col = 4, k_col = 5,
#                             color_fac, showlegend = FALSE, point_size = 2, colors = NULL) {
#   DT <- as.data.table(df)
#   # 좌표 벡터
#   xi <- DT[[i_col]]; yj <- DT[[j_col]]; zk <- DT[[k_col]]
  
#   # 팩터 및 팔레트
#   f <- as.factor(color_fac)
#   levs <- levels(f)
#   nlev <- length(levs)
#   pal <- .discrete_palette(nlev)
  
#   # plotly는 범주마다 trace를 하나씩 만듭니다(성능상 많은 범주면 legend OFF 권장)
#   p <- plot_ly(type = "scatter3d", mode = "markers")
#   for (ii in seq_len(nlev)) {
#     idx <- which(f == levs[ii])
#     if (length(idx) == 0) next
#     p <- add_markers(
#       p,
#       x = xi[idx], y = yj[idx], z = zk[idx],
#       name = levs[ii],
#       marker = list(size = point_size, opacity = 1, color = pal[ii]),
#       hoverinfo = "text",
#       text = paste0("(", xi[idx], ", ", yj[idx], ", ", zk[idx], ")\n", "group: ", levs[ii]),
#       showlegend = showlegend
#     )
#   }
  
#   p <- layout(
#     p,
#     scene = list(
#       xaxis = list(title = "i"),
#       yaxis = list(title = "j"),
#       zaxis = list(title = "k"),
#       aspectmode = "data"  # 좌표 비율 유지
#     ),
#     legend = list(orientation = "v")
#   )
#   p
# }
.plot_voxels_3d <- function(df, i_col = 3, j_col = 4, k_col = 5,
                            color_fac, showlegend = FALSE, point_size = 2,
                            colors = NULL, opacity = 1) {
  DT <- data.table::as.data.table(df)

  xi <- DT[[i_col]]; yj <- DT[[j_col]]; zk <- DT[[k_col]]

  f <- as.factor(color_fac)
  levs <- levels(f)
  nlev <- length(levs)

  # --- 팔레트 결정: colors가 있으면 그걸 쓰고, 없으면 기본 팔레트 ---
  pal <- NULL
  if (!is.null(colors)) {
    # named vector면 level 이름으로 매칭
    if (!is.null(names(colors))) {
      pal <- unname(colors[levs])
    } else {
      # unnamed vector면 길이 체크
      stopifnot(length(colors) >= nlev)
      pal <- colors[seq_len(nlev)]
    }
  } else {
    pal <- .discrete_palette(nlev)
  }

  p <- plotly::plot_ly(type = "scatter3d", mode = "markers")
  for (ii in seq_len(nlev)) {
    idx <- which(f == levs[ii])
    if (!length(idx)) next

    p <- plotly::add_markers(
      p,
      x = xi[idx], y = yj[idx], z = zk[idx],
      name = levs[ii],
      marker = list(size = point_size, opacity = opacity, color = pal[ii]),
      hoverinfo = "text",
      text = paste0("(", xi[idx], ", ", yj[idx], ", ", zk[idx], ")\n", "group: ", levs[ii]),
      showlegend = showlegend
    )
  }

  plotly::layout(
    p,
    scene = list(
      xaxis = list(title = "i"),
      yaxis = list(title = "j"),
      zaxis = list(title = "k"),
      aspectmode = "data"
    ),
    legend = list(orientation = "v")
  )
}

# ----------------------------
# 1) 스펙트럴 클러스터 색상 버전
# ----------------------------
plot_voxels_by_cluster <- function(df,
                                   cluster_col = "cluster_200",
                                   i_col = 3, j_col = 4, k_col = 5,
                                   showlegend = FALSE,
                                   point_size = 2) {
  if (!cluster_col %in% names(df)) stop(sprintf("%s 컬럼이 없습니다.", cluster_col))
  colv <- df[[cluster_col]]
  
  # 0(고립/미할당) 처리: 별도 색(회색)으로 묶고, 나머지는 팩터
  iso_idx <- which(colv == 0)
  grp <- as.character(colv)
  grp[iso_idx] <- "isolated(0)"
  
  .plot_voxels_3d(
    df,
    i_col = i_col, j_col = j_col, k_col = k_col,
    color_fac = grp,
    showlegend = showlegend,
    point_size = point_size
  )
}

# ----------------------------
# 2) 구조적 레이블(데이터 1열) 색상 버전
# ----------------------------
plot_voxels_by_struct_label <- function(df,
                                        label_col = NULL,
                                        i_col = 3, j_col = 4, k_col = 5,
                                        showlegend = TRUE,
                                        point_size = 2) {
  # label 컬럼 이름 자동 추론 (없으면 첫 번째 열)
  if (is.null(label_col)) {
    if ("label" %in% names(df)) label_col <- "label" else label_col <- names(df)[1]
  }
  if (!label_col %in% names(df)) stop(sprintf("'%s' 컬럼이 없습니다.", label_col))
  
  grp <- as.factor(df[[label_col]])
  .plot_voxels_3d(
    df,
    i_col = i_col, j_col = j_col, k_col = k_col,
    color_fac = grp,
    showlegend = showlegend,     # 구조 라벨 수가 너무 많으면 FALSE로 바꾸세요
    point_size = point_size
  )
}


# 사용 예시
# # df <- fmri_list[[1]]
# # 1) 스펙트럴 클러스터 색상
# p1 <- plot_voxels_by_cluster(
#   fmri_list[[1]],
#   cluster_col = "cluster_200",
#   i_col = 3, j_col = 4, k_col = 5,
#   showlegend = FALSE,   # 200개 클러스터면 범례 OFF 권장
#   point_size = 5
# )
# 
# # 2) 구조적 레이블 색상 (데이터 1열이 label이라면 자동 인식)
# p2 <- plot_voxels_by_struct_label(
#   fmri_list[[1]],
#   label_col = "Label",  # 없으면 생략 가능; 자동으로 첫 열 사용
#   i_col = 3, j_col = 4, k_col = 5,
#   showlegend = TRUE,    # 레이블 종류 많으면 FALSE로
#   point_size = 5
# )
# 
# 
# p3 <- plot_voxels_by_cluster(
#   fmri_list[[1]],
#   cluster_col = "cluster_24",
#   i_col = 3, j_col = 4, k_col = 5,
#   showlegend = TRUE,   # 200개 클러스터면 범례 OFF 권장
#   point_size = 5
# )
# 
# p1
# p2
# p3

################################################################################
## 4. inter person clustering
################################################################################
# 이후 코드는 주석 추가 따로 하지 않은 상태입니다.
# 최종단 clustering
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(RSpectra)
})

## ---------------------------------------------------------
# (옵션) 기존 build_weighted_adj_26 대신, 미리 만든 edges를 재사용 (일반화 버전)
# ---------------------------------------------------------
build_weighted_adj_with_edges <- function(
    df, edges,
    time_cols = 6:(ncol(df)-2),
    weight_fun = weight_corr_by_rows,   # ← 원하는 가중치 함수로 교체 가능
    clamp_negative = NULL,              # ← NULL이면 자동 결정: 음수 있으면 TRUE, 아니면 FALSE
    r_min = NULL,
    chunk_size = 50000L,
    verbose = FALSE,
    ...
){
  # 시간행렬
  X <- as.matrix(as.data.frame(df)[, time_cols, drop = FALSE])
  storage.mode(X) <- "double"
  
  # 가중치 계산: 어떤 weight_fun이든 동일 시그니처로 호출
  # (예: weight_corr_by_rows / weight_l2_ts / weight_srvf_ts)
  w_raw <- weight_fun(X, edges, chunk_size = chunk_size, ...)
  
  # clamp_negative 자동 결정 (기본: 음수가 있으면 TRUE, 아니면 FALSE)
  if (is.null(clamp_negative)) {
    clamp_negative <- any(w_raw < 0, na.rm = TRUE)
  }
  
  # 후처리(음수 clamp / 임계치)
  w <- postprocess_weights(w_raw, clamp_negative = clamp_negative, r_min = r_min)
  
  # 희소 대칭 행렬
  assemble_sparse_adj(nrow(df), edges, w)
}


# ---------------------------------------------------------
# 모든 subject가 동일 voxel 집합/순서를 갖도록 정렬
# (col 485: "i_j_k" 키 기준)
# ---------------------------------------------------------
align_by_ijk <- function(fmri_list, ijk_col = 485) {
  ref <- as.character(fmri_list[[1]][[ijk_col]])
  for (s in seq_along(fmri_list)) {
    cur <- as.character(fmri_list[[s]][[ijk_col]])
    if (!identical(sort(ref), sort(cur))) {
      stop(sprintf("fmri_list[[%d]]: ijk 키 집합이 기준과 다릅니다.", s))
    }
    # ref 순서로 재정렬
    ord <- match(ref, cur)
    if (anyNA(ord)) stop(sprintf("fmri_list[[%d]]: 매칭 실패한 ijk가 있습니다.", s))
    fmri_list[[s]] <- fmri_list[[s]][ord, , drop = FALSE]
  }
  fmri_list
}

# ---------------------------------------------------------
# edge 단위 동일-클러스터 여부(0/1) 벡터 만들기
#  - labels: 길이 N (0은 고립/미할당)
#  - both>0 인 edge만 집계에 포함 (count 증가), 그 외는 포함 X
# ---------------------------------------------------------
edge_agreement_01 <- function(labels, edges) {
  li <- labels[edges$i_idx]; lj <- labels[edges$j_idx]
  use <- (li > 0L & lj > 0L)
  agree <- integer(nrow(edges))
  agree[use] <- as.integer(li[use] == lj[use])
  list(agree = agree, use = use)
}

# ---------------------------------------------------------
# 메인: Craddock 2단계 그룹 파이프라인 (가중치 함수 주입형)
# ---------------------------------------------------------
craddock_group_parcellation <- function(
    fmri_list,
    K = 200,
    # 열 지정
    i_col = 3, j_col = 4, k_col = 5, ijk_col = 485,
    time_cols = 6:483,
    # 개인 가중치 옵션
    weight_fun = weight_corr_by_rows,    # ← 어떤 가중치 함수든 교체 가능
    clamp_negative = NULL,               # ← NULL=자동(음수 있으면 clamp), TRUE/FALSE로 강제 가능
    r_min = NULL,
    chunk_size = 50000L,
    # NCut 옵션
    remove_isolated = TRUE, kmeans_nstart = 5, kmeans_itermax = 300, seed = 2025,
    verbose = TRUE,
    ...                                   # ← weight_fun에 전달할 추가 인자들 (kernel, sigma, dt, smooth, warp 등)
){
  S <- length(fmri_list)
  if (verbose) message(sprintf("Subjects: %d", S))
  # 0) ijk 정렬 보장
  fmri_list <- align_by_ijk(fmri_list, ijk_col = ijk_col)
  
  # 1) 26-이웃 edge를 한 번만 구성 (기준: 첫 번째 subject)
  if (verbose) message("Build 26-neighborhood edges (once)...")
  edges <- build_edges_26(fmri_list[[1]], i_col = i_col, j_col = j_col, k_col = k_col)
  N <- nrow(fmri_list[[1]])
  E <- nrow(edges)
  if (verbose) message(sprintf("  Voxels: %d, Edges: %d", N, E))
  
  # 2) 개인별 라벨 + edge agreement 누적
  labels_list <- vector("list", S)
  sum_agree   <- numeric(E)   # \sum_s 1{same cluster}
  sum_count   <- numeric(E)   # \sum_s 1{both>0}
  
  for (s in seq_len(S)) {
    if (verbose) message(sprintf("Subject %d/%d: build W & NCut ...", s, S))
    df <- fmri_list[[s]]
    
    # (a) 개인 가중 그래프 — weight_fun & ... 그대로 전달
    W_s <- build_weighted_adj_with_edges(
      df, edges,
      time_cols       = time_cols,
      weight_fun      = weight_fun,
      clamp_negative  = clamp_negative,   # NULL이면 자동 판단(상관=clamp, L2/SRVF=그대로)
      r_min           = r_min,
      chunk_size      = chunk_size,
      verbose         = FALSE,
      ...
    )
    # (b) 개인 NCut → 라벨
    res_s <- spectral_clustering_ncut(
      W_s, K = K,
      remove_isolated = remove_isolated,
      kmeans_nstart   = kmeans_nstart,
      kmeans_itermax  = kmeans_itermax,
      seed            = seed, verbose = FALSE
    )
    
    labels_list[[s]] <- res_s$labels
    
    # (c) edge 단위 동일-클러스터 여부 집계
    ag <- edge_agreement_01(res_s$labels, edges)
    sum_agree <- sum_agree + ag$agree
    sum_count <- sum_count + as.numeric(ag$use)
  }
  
  # 3) 그룹 coincidence 가중치 (edge별 평균 동의율)
  if (verbose) message("Compute group coincidence weights ...")
  w_group <- numeric(E)
  nz <- which(sum_count > 0)
  w_group[nz] <- sum_agree[nz] / sum_count[nz]
  
  # 4) 그룹 그래프 W_group 생성 후 최종 NCut
  if (verbose) message("Assemble group graph & run final NCut ...")
  W_group <- assemble_sparse_adj(N, edges, w_group)
  W_group <- as(forceSymmetric(W_group, uplo = "U"), "dgCMatrix")
  
  res_group <- spectral_clustering_ncut(
    W_group, K = K,
    remove_isolated = remove_isolated,
    kmeans_nstart   = kmeans_nstart,
    kmeans_itermax  = kmeans_itermax,
    seed            = seed, verbose = TRUE
  )
  
  list(
    group_labels   = res_group$labels,    # 길이 N (0=고립)
    group_graph    = W_group,             # dgCMatrix (coincidence 가중 그래프)
    edges          = edges,               # 26-이웃 edge 목록
    edge_weights   = w_group,             # edge별 동의율 [0,1]
    subject_labels = labels_list,         # 개인별 라벨 리스트
    info = list(N = N, E = E, S = S, K = K)
  )
}


# -------------------------------------------------------------------------
# 메인: Craddock 2단계 그룹 파이프라인 (All-pairs, i-블록 스트리밍, 병렬 옵션 지원)
#  - build_weighted_adj_all_streaming_parallel() 를 내부에서 호출
#  - .safe_weight_call 는 build_* 쪽에서 이미 사용하도록 패치되어 있다고 가정
# -------------------------------------------------------------------------
craddock_group_parcellation_allpairs <- function(
    fmri_list,
    K = 200,
    # 공통 열 지정
    ijk_col       = 485,
    time_cols     = 6:483,
    # 개인 그래프 옵션 (모든 노드쌍 스트리밍 빌드)
    weight_fun    = weight_corr_by_rows,   # weight_l2_ts, weight_pcorr_by_rows, weight_srvf_ts 등 교체 가능
    clamp_negative= NULL,                  # NULL=자동(음수 있으면 0 클램프), TRUE/FALSE 강제
    r_min         = NULL,                  # 임계치 미만 컷(가중치 희소화)
    i_block_size  = 2000L,                 # i-블록 크기 (메모리/시간 트레이드오프)
    chunk_size    = 50000L,                # weight_fun 내부 에지 배치 크기
    # NCut 옵션
    remove_isolated = TRUE,
    kmeans_nstart   = 5,
    kmeans_itermax  = 300,
    seed            = 2025,
    verbose         = TRUE,
    # ---- 병렬 옵션 (개인 그래프 단계) ----
    parallel        = c("none","multisession","multicore"),
    workers         = max(1L, parallel::detectCores() - 1L),
    tame_blas       = TRUE,                # 병렬-안의 BLAS/MKL 다중스레딩 억제,
    # ----- plot -------
    screeplot      = FALSE,
    screeplot_name = NULL,
    ...
){
  parallel <- match.arg(parallel)
  
  stopifnot(is.list(fmri_list), length(fmri_list) >= 1L)
  S <- length(fmri_list)
  if (verbose) message(sprintf("[AllPairs] Subjects: %d", S))
  
  # -- 0) ijk 정렬 보장 (모든 subject가 동일 순서/집합)
  fmri_list <- align_by_ijk(fmri_list, ijk_col = ijk_col)
  
  # 공통 N 확인
  N <- nrow(fmri_list[[1]])
  if (N < 2L) stop("N < 2: 군집 불가")
  if (verbose) message(sprintf("[AllPairs] Voxels: %d (full upper-triangle will be streamed)", N))
  
  # -- 1) 개인별 그래프 W_s (전쌍 스트리밍) 만들고 NCut → 라벨 수집
  if (verbose) message("[AllPairs] (1/3) Subject-level W_s build (all-pairs; streamed) & NCut...")
  labels_list <- vector("list", S)
  
  for (s in seq_len(S)) {
    if (verbose) message(sprintf("  - Subject %d/%d: build W_s via streaming%s ...",
                                 s, S, if (parallel=="none") "" else paste0(" [", parallel, ", workers=", workers, "]")))
    df <- fmri_list[[s]]
    
    # === 변경된 부분: 병렬 스트리밍 빌더 사용 ===
    W_s <- build_weighted_adj_all_streaming_parallel(
      df              = df,
      time_cols       = time_cols,
      weight_fun      = weight_fun,
      i_block_size    = i_block_size,
      r_min           = r_min,
      clamp_negative  = clamp_negative,
      chunk_size      = chunk_size,
      parallel        = parallel,     # ← 병렬 모드 전달
      workers         = workers,      # ← 워커 수 전달
      tame_blas       = tame_blas,    # ← BLAS 억제
      verbose         = FALSE,
      ...
    )
    
    # 개인 NCut → 라벨
    res_s <- spectral_clustering_ncut(
      W_s, K = K,
      remove_isolated = remove_isolated,
      kmeans_nstart   = kmeans_nstart,
      kmeans_itermax  = kmeans_itermax,
      seed            = seed, verbose = FALSE,
      screeplot = TRUE,
      screeplot_name = screeplot_name
    )
    labels_list[[s]] <- as.integer(res_s$labels)
  }
  
  # -- 2) 그룹 동의율(coincidence) 그래프 W_group을 “라벨만”으로 i-블록 스트리밍 조립
  if (verbose) message("[AllPairs] (2/3) Assemble group coincidence graph (streaming over i-blocks)...")
  
  Ii_list <- list(); Jj_list <- list(); Xx_list <- list()
  total_kept <- 0L
  
  i_starts <- seq.int(1L, N - 1L, by = i_block_size)
  n_blocks <- length(i_starts)
  
  for (b in seq_along(i_starts)) {
    i_lo <- i_starts[b]
    i_hi <- min(i_lo + i_block_size - 1L, N - 1L)
    if (verbose) message(sprintf("  - Block %d/%d: i in [%d, %d]", b, n_blocks, i_lo, i_hi))
    
    # 상삼각 블록 에지 인덱스 생성 (i<j)
    edges_block <- build_edges_block(i_lo, i_hi, N)  # (i_idx, j_idx)
    nb <- nrow(edges_block)
    if (nb == 0L) next
    
    # 라벨로부터 에지 동의/카운트 누적
    sum_agree <- numeric(nb)
    sum_count <- numeric(nb)
    
    for (s in seq_len(S)) {
      L  <- labels_list[[s]]
      li <- L[edges_block$i_idx]
      lj <- L[edges_block$j_idx]
      use <- (li > 0L & lj > 0L)
      sum_agree[use] <- sum_agree[use] + as.numeric(li[use] == lj[use])
      sum_count[use] <- sum_count[use] + 1.0
    }
    
    # 블록 가중치 = 평균 동의율 (분모 0은 0)
    w_block <- numeric(nb)
    nz <- which(sum_count > 0)
    if (length(nz)) w_block[nz] <- sum_agree[nz] / sum_count[nz]
    
    # 0이 아닌 가중치만 저장
    keep <- which(is.finite(w_block) & w_block != 0)
    if (length(keep)) {
      Ii_list[[length(Ii_list) + 1L]] <- edges_block$i_idx[keep]
      Jj_list[[length(Jj_list) + 1L]] <- edges_block$j_idx[keep]
      Xx_list[[length(Xx_list) + 1L]] <- w_block[keep]
      total_kept <- total_kept + length(keep)
    }
    
    if ((b %% 10L) == 0L) gc()
  }
  
  if (verbose) message(sprintf("  * kept upper-tri nonzeros: %d", total_kept))
  
  # 상삼각 트리플릿으로 그룹 그래프 조립 후 대칭화
  if (total_kept == 0L) {
    if (verbose) message("  * No edges survived; returning zero matrix and empty labels.")
    W_group <- Matrix::Matrix(0, nrow = N, ncol = N, sparse = TRUE)
  } else {
    i_all <- do.call(c, Ii_list)
    j_all <- do.call(c, Jj_list)
    x_all <- do.call(c, Xx_list)
    
    A_upper <- Matrix::sparseMatrix(
      i = i_all, j = j_all, x = x_all,
      dims = c(N, N), symmetric = FALSE, giveCsparse = TRUE
    )
    W_group <- Matrix::forceSymmetric(A_upper, uplo = "U")
    W_group <- methods::as(W_group, "dgCMatrix")
  }
  
  # -- 3) 그룹 NCut
  if (verbose) message("[AllPairs] (3/3) Final NCut on group coincidence graph ...")
  res_group <- spectral_clustering_ncut(
    W_group, K = K,
    remove_isolated = remove_isolated,
    kmeans_nstart   = kmeans_nstart,
    kmeans_itermax  = kmeans_itermax,
    seed            = seed, verbose = TRUE, 
    screeplot = FALSE
  )
  
  list(
    group_labels   = as.integer(res_group$labels),  # 길이 N (0=고립/미할당)
    group_graph    = W_group,                       # dgCMatrix (coincidence weights)
    subject_labels = labels_list,                   # 각 개인 라벨
    info = list(
      N = N, S = S, K = K, i_block_size = i_block_size,
      mode = "allpairs_streamed", parallel = parallel, workers = workers
    )
  )
}

################################################################################
## 5. clustering result estimation
################################################################################

## =========================================================
##  A. Edge-기반 보조 함수들 (빠르고 메모리 안전)
## =========================================================

# 라벨 -> 에지 단위 동일클러스터(0/1) 벡터
labels_to_edge01 <- function(labels, edges) {
  li <- labels[edges$i_idx]; lj <- labels[edges$j_idx]
  as.integer(li > 0L & lj > 0L & li == lj)
}

# Dice(A,B) on edges (0/1 벡터)
dice_edge01 <- function(a, b) {
  sa <- sum(a); sb <- sum(b); inter <- sum(a & b)
  if ((sa + sb) == 0) return(NA_real_)
  2 * inter / (sa + sb)
}

## =========================================================
##  B. Dice-LOOCV (Two-level 방식의 효율적 LOO 재구성)
##     [수정됨] out$edges 유무에 따라 희소(26-이웃) 모드와
##             전체 쌍(all-pairs) 스트리밍 모드를 자동 전환
## =========================================================
loocv_dice_from_out <- function(
    out,
    K               = NULL,
    remove_isolated = TRUE,
    kmeans_nstart   = 10,
    kmeans_itermax  = 1000,
    seed            = 2025,
    verbose         = TRUE,
    # --- All-pairs 모드용 스트리밍 파라미터 ---
    i_block_size    = 2000L
){
  stopifnot(!is.null(out$subject_labels)) # 공통 필수 요소
  
  # --- 모드 자동 감지 ---
  # out$edges가 있으면 Sparse(26-neighbor) 모드
  # out$edges가 없으면 All-pairs(streaming) 모드
  is_sparse_mode <- !is.null(out$edges)
  
  labels_list <- out$subject_labels
  if (is.null(K)) K <- out$info$K
  N <- out$info$N
  S <- length(labels_list)
  
  if (verbose) {
    message(sprintf("LOOCV: subjects=%d, K=%d, N=%d", S, K, N))
    message(sprintf("Mode: %s", 
                    if(is_sparse_mode) "Sparse (26-neighbor, using out$edges)" 
                    else "All-pairs (streaming, i_block_size=%d)", i_block_size))
  }
  
  dices <- numeric(S)
  
  # --- 로직 분기 ---
  
  if (is_sparse_mode) {
    # ==========================================================
    # MODE 1: 희소(Sparse) 모드 (craddock_group_parcellation)
    # ==========================================================
    edges <- out$edges
    E <- nrow(edges)
    
    # 각 피험자의 edge 동의(agree/use) 미리 계산 (빠름)
    ag_list <- lapply(labels_list, function(lbl) edge_agreement_01(lbl, edges))
    sum_agree <- Reduce(`+`, lapply(ag_list, `[[`, "agree"))
    sum_count <- Reduce(`+`, lapply(ag_list, function(x) as.numeric(x$use)))
    
    if (verbose) message(sprintf("Pre-calculated agreement for %d sparse edges.", E))
    
    for (m in seq_len(S)) {
      if (verbose) message(sprintf("  LOO (Sparse) %d/%d ...", m, S))
      
      # 1. LOO 그룹 그래프 (Wm) 생성 (엣지 기반)
      agree_m <- ag_list[[m]]$agree
      use_m   <- as.numeric(ag_list[[m]]$use)
      
      num <- sum_agree - agree_m
      den <- pmax(0, sum_count - use_m)
      
      w <- numeric(E)
      nz <- which(den > 0)
      w[nz] <- num[nz] / den[nz]
      
      # [호환성] 'assemble_sparse_adj' (오래된/비병렬 함수) 사용
      Wm <- assemble_sparse_adj(N, edges, w) 
      
      # 2. LOO 그룹 클러스터링 (Wm)
      # [호환성] 'spectral_clustering_ncut' (오래된/비병렬 함수) 사용
      res_m <- spectral_clustering_ncut(
        Wm, K = K,
        remove_isolated = remove_isolated,
        kmeans_nstart   = kmeans_nstart,
        kmeans_itermax  = kmeans_itermax,
        seed            = seed,
        verbose         = FALSE
      )
      
      # 3. Dice 계산 (엣지 기반)
      a <- labels_to_edge01(res_m$labels, edges)     # 그룹(LOO) 인접
      b <- labels_to_edge01(labels_list[[m]], edges) # 당사자 개인 인접
      dices[m] <- dice_edge01(a, b)
    }
    
  } else {
    # ==========================================================
    # MODE 2: 전체 쌍(All-pairs) 모드 (craddock_group_parcellation_allpairs)
    # ==========================================================
    if (verbose) message("Using i-block streaming to build LOO graphs. This may take time.")
    
    for (m in seq_len(S)) {
      if (verbose) message(sprintf("  LOO (All-pairs) %d/%d ...", m, S))
      
      # 1. LOO 그룹 그래프 (Wm) 생성 (스트리밍 기반)
      #    (craddock_group_parcellation_allpairs의 2단계 로직을 여기서 재현)
      
      Ii_list <- list(); Jj_list <- list(); Xx_list <- list()
      total_kept <- 0L
      
      i_starts <- seq.int(1L, N - 1L, by = i_block_size)
      n_blocks <- length(i_starts)
      
      # m번 피험자를 제외한 S-1명의 인덱스
      loo_indices <- setdiff(seq_len(S), m)
      
      for (b in seq_along(i_starts)) {
        i_lo <- i_starts[b]
        i_hi <- min(i_lo + i_block_size - 1L, N - 1L)
        
        # (이 함수는 스크립트 전역에 이미 존재해야 함)
        edges_block <- build_edges_block(i_lo, i_hi, N) 
        nb <- nrow(edges_block)
        if (nb == 0L) next
        
        sum_agree <- numeric(nb)
        sum_count <- numeric(nb)
        
        # S-1 명의 피험자에 대해서만 동의율 계산
        for (s in loo_indices) { 
          L  <- labels_list[[s]]
          li <- L[edges_block$i_idx]
          lj <- L[edges_block$j_idx]
          use <- (li > 0L & lj > 0L)
          sum_agree[use] <- sum_agree[use] + as.numeric(li[use] == lj[use])
          sum_count[use] <- sum_count[use] + 1.0 # 유효한 피험자 수 카운트
        }
        
        w_block <- numeric(nb)
        nz <- which(sum_count > 0) # 분모가 0이 아닌 (유효한) 엣지
        if (length(nz)) w_block[nz] <- sum_agree[nz] / sum_count[nz]
        
        keep <- which(is.finite(w_block) & w_block != 0)
        if (length(keep)) {
          Ii_list[[length(Ii_list) + 1L]] <- edges_block$i_idx[keep]
          Jj_list[[length(Jj_list) + 1L]] <- edges_block$j_idx[keep]
          Xx_list[[length(Xx_list) + 1L]] <- w_block[keep]
          total_kept <- total_kept + length(keep)
        }
      } # end (block b)
      
      # Wm 조립 (All-pairs의 그룹 그래프 생성 방식과 동일)
      if (total_kept == 0L) {
        Wm <- Matrix::Matrix(0, nrow = N, ncol = N, sparse = TRUE)
      } else {
        i_all <- do.call(c, Ii_list); Jj_list_c <- do.call(c, Jj_list); x_all <- do.call(c, Xx_list)
        A_upper <- Matrix::sparseMatrix(
          i = i_all, j = Jj_list_c, x = x_all,
          dims = c(N, N), symmetric = FALSE, giveCsparse = TRUE
        )
        Wm <- Matrix::forceSymmetric(A_upper, uplo = "U")
        Wm <- methods::as(Wm, "dgCMatrix")
      }
      
      # 2. LOO 그룹 클러스터링 (Wm)
      # [호환성] 'spectral_clustering_ncut' (공통 사용 함수)
      res_m <- spectral_clustering_ncut(
        Wm, K = K,
        remove_isolated = remove_isolated,
        kmeans_nstart   = kmeans_nstart,
        kmeans_itermax  = kmeans_itermax,
        seed            = seed,
        verbose         = FALSE
      )
      
      # 3. Dice 계산 (레이블 기반 All-pairs)
      # (이 함수는 스크립트 전역에 이미 존재해야 함)
      dices[m] <- overlap_from_labels(res_m$labels, labels_list[[m]], dice = TRUE)
      
    } # end (subject m)
  } # end (else: all-pairs mode)
  
  list(
    per_subject = dices,
    mean = mean(dices, na.rm = TRUE),
    sd   = stats::sd(dices,  na.rm = TRUE)
  )
}

## =========================================================
##  C. 실루엣(silhouette) 계산 (그래프-제한 유사도 버전 - edge 가 있는 애들끼리만 비교한것)
##     - W: 희소 유사도 그래프(26-이웃, 비음수 권장)
##     - labels: 평가할 라벨(그룹 라벨 고정)
## =========================================================

# faster silhouette score function
# 모든 외부 군집에 대해서 평균내는 식으로 계산 - 원래 score 보다 과대추정되지만
# 조금 더 빠름

# silhouette_from_sparseW_labels <- function(W, labels) {
#   # 안전장치 및 삼중표기 변환
#   W <- forceSymmetric(as(W, "dgCMatrix"), uplo = "U")
#   W <- as(W, "dgTMatrix")
#   
#   # sparse matrix 에서 값이 존재하는 i,j,value 만 추출
#   i <- W@i + 1L; j <- W@j + 1L; x <- W@x
#   L <- as.integer(labels)
#   N <- length(L)
#   
#   # 라벨 0(고립/미할당) 관련 에지는 제외
#   ok   <- (L[i] > 0L & L[j] > 0L)
#   i <- i[ok]; j <- j[ok]; x <- x[ok]
#   
#   # 연결된 node i 와 node j 에 할당된 cluster 가 같은 경우
#   # 대부분 연결되어 있으면 같은 cluster 에 할당되었다:
#   # sum(same) / length(same) = 0.83...
#   # 이런 양상 때문에 LOOCV 는 일반적으로 높게 나오고, siluoette 은 낮게 나왔을
#   # 수 있겠다. 
#   # 군집 할당에 있어서 가중치가 아니라 이미 정의된 연결성이 영향을 많이 미치게 되면
#   # 실제 가중치의 유사성으로 판단하는 실루엣 점수가 낮게 나올 수 있겠다.
#   # 그리고 해당 siluoette은 가장 가까운 군집과의 거리 가 아니라, 속하지 않은
#   # 모든 node 와의 거리의 평균을 구한 것이기 때문에 오히려 과대추정되어있을 
#   # 가능성이 높다.
#   
#   same <- (L[i] == L[j])
#   
#   sum_in  <- numeric(N); cnt_in  <- integer(N)
#   sum_out <- numeric(N); cnt_out <- integer(N)
#   
#   # 에지마다 양 끝 노드에 모두 누적
#   for (idx in seq_along(x)) {
#     ii <- i[idx]; jj <- j[idx]; val <- x[idx]
#     if (same[idx]) {
#       sum_in[ii]  <- sum_in[ii]  + val; cnt_in[ii]  <- cnt_in[ii]  + 1L
#       sum_in[jj]  <- sum_in[jj]  + val; cnt_in[jj]  <- cnt_in[jj]  + 1L
#     } else {
#       sum_out[ii] <- sum_out[ii] + val; cnt_out[ii] <- cnt_out[ii] + 1L
#       sum_out[jj] <- sum_out[jj] + val; cnt_out[jj] <- cnt_out[jj] + 1L
#     }
#   }
#   
#   Ai <- ifelse(cnt_in  > 0L, sum_in  / cnt_in,  NA_real_)
#   Bi <- ifelse(cnt_out > 0L, sum_out / cnt_out, NA_real_)
#   
#   denom <- pmax(Ai, Bi)
#   si <- ifelse(is.finite(denom) & denom > 0,
#                (Ai - Bi) / denom,
#                NA_real_)
#   # 라벨 0 제외
#   keep <- which(L > 0L & is.finite(si))
#   overall <- if (length(keep)) mean(si[keep]) else NA_real_
#   
#   by_cluster <- tapply(si[keep], L[keep], mean, na.rm = TRUE)
#   
#   list(overall = overall, by_cluster = by_cluster, per_node = si)
# }


# real silhuoette score
silhouette_from_sparseW_labels <- function(W, labels) {
  # 연결된 node i 와 node j 에 할당된 cluster 가 같은 경우
  # 대부분 연결되어 있으면 같은 cluster 에 할당되었다:
  # sum(same) / length(same) = 0.83...
  # 이런 양상 때문에 LOOCV 는 일반적으로 높게 나오고, siluoette 은 낮게 나왔을
  # 수 있겠다. 
  # 군집 할당에 있어서 가중치가 아니라 이미 정의된 연결성이 영향을 많이 미치게 되면
  # 실제 가중치의 유사성으로 판단하는 실루엣 점수가 낮게 나올 수 있겠다.
  # 그리고 해당 siluoette은 '가장 가까운 군집과의 거리' 가 아니라, 속하지 않은
  # 모든 node 와의 거리의 평균을 구한 것이기 때문에 오히려 과대추정되어있을 
  # 가능성이 높다.
  
  # --- 0) 안전장치 & 희소행렬 준비 ---
  W <- forceSymmetric(as(W, "dgCMatrix"), uplo = "U")
  W <- as(W, "dgTMatrix")
  
  i <- W@i + 1L; j <- W@j + 1L; x <- W@x
  L <- as.integer(labels)
  N <- length(L)
  
  # 라벨 0(미할당/고립)과 연결된 간선 제거
  ok <- (L[i] > 0L & L[j] > 0L)
  i <- i[ok]; j <- j[ok]; x <- x[ok]
  if (length(x) == 0L) {
    return(list(overall = NA_real_, by_cluster = numeric(), per_node = rep(NA_real_, N)))
  }
  
  # 중복 카운트 방지: 상삼각만 사용
  keep_upper <- i < j
  i <- i[keep_upper]; j <- j[keep_upper]; x <- x[keep_upper]
  if (length(x) == 0L) {
    return(list(overall = NA_real_, by_cluster = numeric(), per_node = rep(NA_real_, N)))
  }
  
  # --- 1) per-node, per-cluster 평균 유사도 계산 (핵심) ---
  # 이 블록이 Bi = "타 클러스터 평균 중 최대값" 계산을 가능하게 함.
  # 각 간선을 양방향으로 펼쳐서 (node, neighbor_cluster) 단위로 집계
  # 클러스터 라벨을 1..K로 리맵(간단/안전)
  uniq <- sort(unique(L[L > 0L]))
  mapL <- integer(max(L))
  mapL[uniq] <- seq_along(uniq)
  C <- integer(N); C[L > 0L] <- mapL[L[L > 0L]]  # 0은 그대로 0
  ci <- C[i]; cj <- C[j]
  
  nodes <- c(i, j)
  neighC <- c(cj, ci)
  wts   <- c(x, x)
  
  # (node, neighC) 쌍을 하나의 정수 key로 인코딩해서 집계
  key <- nodes + (neighC - 1L) * N  # 1..N + (cluster-1)*N
  
  # 합/개수 집계
  sum_by_key <- tapply(wts, key, sum)
  cnt_by_key <- tapply(rep.int(1L, length(key)), key, sum)
  
  keys <- as.integer(names(sum_by_key))
  nds  <- ((keys - 1L) %% N) + 1L         # node id 복원
  cls  <- ((keys - 1L) %/% N) + 1L         # 이웃 클러스터 id 복원
  mns  <- as.numeric(sum_by_key) / as.numeric(cnt_by_key)  # (node, cluster) 평균 유사도
  
  # 노드별로 (자기 클러스터 평균 Ai, 타 클러스터 평균 중 최대 Bi) 계산
  Ai <- rep(NA_real_, N)
  Bi <- rep(NA_real_, N)
  
  # nds를 기준으로 인덱스 묶기
  idx_by_node <- split(seq_along(nds), nds)
  for (s in names(idx_by_node)) {
    u <- as.integer(s)
    if (C[u] == 0L) next  # 미할당/고립 노드 skip
    idxs <- idx_by_node[[s]]
    cls_u <- cls[idxs]
    mns_u <- mns[idxs]
    
    # 자기 클러스터 평균 유사도
    in_pos <- which(cls_u == C[u])
    if (length(in_pos)) Ai[u] <- mns_u[in_pos[1L]]  # (node, own_cluster)는 보통 1개
    
    # 타 클러스터 평균 유사도 중 최대값
    out_pos <- which(cls_u != C[u])
    if (length(out_pos)) Bi[u] <- max(mns_u[out_pos]) else Bi[u] <- NA_real_
  }
  
  # --- 2) 실루엣 계산 (유사도 버전) ---
  # si = (Ai - Bi) / max(Ai, Bi)
  denom <- pmax(Ai, Bi)
  si <- ifelse(is.finite(denom) & denom > 0, (Ai - Bi) / denom, NA_real_)
  
  # --- 3) 요약 통계 ---
  keep <- which(L > 0L & is.finite(si))
  overall <- if (length(keep)) mean(si[keep]) else NA_real_
  by_cluster <- if (length(keep)) tapply(si[keep], L[keep], mean, na.rm = TRUE) else numeric()
  
  list(overall = overall, by_cluster = by_cluster, per_node = si)
}

silhouette_cluster_strict <- function(W, labels) {
  # W: voxel-by-voxel similarity matrix (대칭, 실수; 희소 가능)
  # labels: 정수형 벡터 (길이 N). 0은 비할당으로 간주하며 K는 labels>0의 고유값 개수.
  # 수식:
  #   a_k = [ z_k^T S z_k ] / [ n_k (n_k - 1) ]
  #   b_k = [ z_k^T S 1  - z_k^T S z_k ] / [ n_k (N - n_k) ]
  #   s_k = (a_k - b_k) / max(a_k, b_k)
  #   si(C) = (1/K) * sum_k s_k
  #
  # 구현 주의:
  # - "엣지만" 쓰지 않습니다. 유사도 0인 쌍(관측 안 된 쌍)도 분모에 포함됩니다.
  # - S는 대칭으로 강제, 대각은 0으로 맞춥니다(수식의 i != j 반영).
  
  # --- 준비: 행렬 형태 표준화 & 대칭/대각 정리 ---
  if (inherits(W, "matrix")) {
    # 밀집인 경우
    S <- Matrix::Matrix(W, sparse = FALSE)
  } else {
    # 희소/기타는 dgC로 캐스팅
    S <- methods::as(W, "dgCMatrix")
  }
  # 대칭 강제(상삼각 신뢰), 대각 0
  S <- Matrix::forceSymmetric(S, uplo = "U")
  Matrix::diag(S) <- 0
  
  N <- nrow(S)
  L <- as.integer(labels)
  if (length(L) != N) stop("length(labels) must equal nrow(W).")
  
  # 유효 클러스터(>0)
  uniq <- sort(unique(L[L > 0L]))
  K <- length(uniq)
  if (K == 0) {
    return(list(
      overall = NA_real_, per_cluster = numeric(),
      a = numeric(), b = numeric(), n_k = integer()
    ))
  }
  
  # 공통 항: S %*% 1 (클러스터별 b_k 분자에 필요)
  one <- rep.int(1, N)
  S_one <- as.numeric(S %*% one)
  
  a <- numeric(K)
  b <- numeric(K)
  s <- numeric(K)
  nk <- integer(K)
  names(a) <- names(b) <- names(s) <- names(nk) <- as.character(uniq)
  
  for (kk in seq_len(K)) {
    k <- uniq[kk]
    z <- as.numeric(L == k)          # 0/1 지시 벡터
    n_k <- sum(z)
    nk[kk] <- n_k
    
    # 분모가 0 되는 경우(클러스터 크기 0/1 또는 전부 등) 방지
    if (n_k <= 1L || n_k >= N) {
      a[kk] <- NA_real_
      b[kk] <- NA_real_
      s[kk] <- NA_real_
      next
    }
    
    # 분자 계산
    Sz   <- as.numeric(S %*% z)        # S z
    numA <- sum(z * Sz)                # z^T S z  (대각 0이므로 i!=j 합과 동일)
    numB <- sum(z * S_one) - numA      # z^T S 1 - z^T S z
    
    # 분모
    denA <- n_k * (n_k - 1)
    denB <- n_k * (N - n_k)
    
    # a_k, b_k
    a_k <- numA / denA
    b_k <- numB / denB
    a[kk] <- a_k
    b[kk] <- b_k
    
    # s_k = (a_k - b_k)/max(a_k, b_k)  (정확히 수식대로)
    denom <- max(a_k, b_k)
    s[kk] <- if (is.finite(denom) && denom != 0) (a_k - b_k) / denom else NA_real_
  }
  
  overall <- mean(s, na.rm = TRUE)  # si(C) = (1/K) * sum_k s_k  (NA 제외 평균)
  
  list(
    overall = overall,
    per_cluster = s,   # s_k
    a = a,             # a_k
    b = b,             # b_k
    n_k = nk           # 각 클러스터 크기
  )
}

## =========================================================
##  D. 피험자별 SI: 그룹 라벨 고정, 각 피험자 유사도로 평가
##     - weight_fun 은 (X, edges, ...) -> weight 벡터 형태
## =========================================================
silhouette_per_subjects <- function(
    fmri_list, edges, labels_group,
    weight_fun,
    time_cols       = 6:483,
    clamp_negative  = NULL,  # NULL=자동(clamp 필요시), TRUE/FALSE 강제
    r_min           = NULL,  # 상관식 유사도면 0.5 같은 threshold도 가능
    chunk_size      = 50000L,
    verbose         = TRUE,
    ...
){
  S <- length(fmri_list)
  sils <- numeric(S)
  
  for (s in seq_len(S)) {
    if (verbose) message(sprintf("Silhouette: subject %d/%d ...", s, S))
    W_s <- build_weighted_adj_with_edges(
      df              = fmri_list[[s]],
      edges           = edges,
      time_cols       = time_cols,
      weight_fun      = weight_fun,
      clamp_negative  = clamp_negative,
      r_min           = r_min,
      chunk_size      = chunk_size,
      verbose         = FALSE,
      ...
    )
    ## 2) 기존 함수로 실루엣 계산
    # si <- silhouette_from_sparseW_labels(W_s, labels_group)
    # creddock 방법 그대로 따르기 위함
    si <- silhouette_cluster_strict(W_s, labels_group)
    sils[s] <- si$overall
  }
  
  list(
    per_subject = sils,
    mean = mean(sils, na.rm = TRUE),
    sd   = stats::sd(sils,  na.rm = TRUE)
  )
}

## =========================================================
##  D'. 피험자별 SI (All-pairs): 그룹 라벨 고정, 각 피험자 유사도로 평가
##      - edges 불필요, 전쌍 그래프를 스트리밍/병렬로 만들어 사용
## =========================================================
silhouette_per_subjects_allpairs <- function(
    fmri_list,                 # subject별 data.frame list
    labels_group,              # 그룹 라벨 (길이 N; 0=고립/미할당 허용)
    weight_fun,                # weight_* 함수
    time_cols       = 6:483,
    clamp_negative  = NULL,    # NULL=자동, TRUE/FALSE 강제
    r_min           = NULL,    # 임계 미만 컷 (상관식이면 0~0.2 등)
    i_block_size    = 2000L,
    chunk_size      = 50000L,
    verbose         = TRUE,
    # 병렬 옵션
    parallel        = c("none","multisession","multicore"),
    workers         = max(1L, parallel::detectCores() - 1L),
    tame_blas       = TRUE,
    ...
){
  parallel <- match.arg(parallel)
  S <- length(fmri_list)
  sils <- numeric(S)
  
  for (s in seq_len(S)) {
    if (verbose) message(sprintf("Silhouette(allpairs): subject %d/%d ...", s, S))
    
    # 1) 전쌍 그래프 생성 (스트리밍/병렬)
    W_s <- build_weighted_adj_all_streaming_parallel(
      df              = fmri_list[[s]],
      time_cols       = time_cols,
      weight_fun      = weight_fun,
      i_block_size    = i_block_size,
      r_min           = r_min,
      clamp_negative  = clamp_negative,
      chunk_size      = chunk_size,
      parallel        = parallel,
      workers         = workers,
      tame_blas       = tame_blas,
      verbose         = FALSE,
      ...
    )
    
    ## 2) 기존 함수로 실루엣 계산
    # si <- silhouette_from_sparseW_labels(W_s, labels_group)
    # creddock 방법 그대로 따르기 위함
    si <- silhouette_cluster_strict(W_s, labels_group)
    sils[s] <- si$overall
  }
  
  list(
    per_subject = sils,
    mean = mean(sils, na.rm = TRUE),
    sd   = stats::sd(sils,  na.rm = TRUE)
  )
}

## =========================================================
##  E. Creddock reference file 를 이용해 DICE metric 계산
## =========================================================
# reference data 를 이용해 DICE metric 계산
# =========================================================
# Dice via cluster-based adjacency (ref vs pred)
# 1) (i,j,k) 공통 좌표로 정렬
# 2) 각 라벨로 "같은 클러스터면 1"인 binary adjacency (sparse) 생성
# 3) 두 adjacency의 Dice 계산 (상삼각 기준, 라벨 불변)
# ---------------------------------------------------------
# Options:
#  - treat_zero_as_bg: 0/NA 라벨을 배경으로 취급(연결 없음)
#  - return_adj: TRUE면 adjacency(상삼각) 희소행렬을 반환
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
})

dice_from_cluster_adjacency <- function(
    ref_df, pred_df,
    ijk_cols        = c("i","j","k"),
    ref_label_col   = "roi",
    pred_label_col  = "cluster",
    treat_zero_as_bg = TRUE,
    return_adj       = FALSE
){
  stopifnot(all(ijk_cols %in% names(ref_df)),
            all(ijk_cols %in% names(pred_df)),
            ref_label_col %in% names(ref_df),
            pred_label_col %in% names(pred_df))
  
  # -- (1) 공통 (i,j,k)로 정렬 (순서는 안 맞을 수 있으니 정렬해서 맞춤)
  common <- inner_join(
    ref_df  %>% select(all_of(ijk_cols)),
    pred_df %>% select(all_of(ijk_cols)),
    by = ijk_cols
  ) %>% distinct() %>%
    arrange(across(all_of(ijk_cols)))
  
  if (nrow(common) < 2L) stop("공통 좌표가 너무 적습니다 (>=2 필요).")
  
  ref_aligned <- common %>% left_join(ref_df  %>% select(all_of(ijk_cols), !!ref_label_col),  by = ijk_cols)
  pred_aligned<- common %>% left_join(pred_df %>% select(all_of(ijk_cols), !!pred_label_col), by = ijk_cols)
  
  ref_labels  <- ref_aligned[[ref_label_col]]
  pred_labels <- pred_aligned[[pred_label_col]]
  N <- nrow(common)
  
  # -- (2) 라벨 -> adjacency (sparse, symmetric, binary)
  make_adj_from_labels <- function(labels) {
    # 배경(0/NA) 제외 여부
    if (treat_zero_as_bg) {
      valid <- !is.na(labels) & labels != 0
      lab_use <- labels[valid]
      idx_use <- which(valid)
    } else {
      valid <- !is.na(labels)
      lab_use <- labels[valid]
      idx_use <- which(valid)
    }
    if (length(idx_use) == 0L) {
      # 모두 배경이면 0행렬 반환
      return(Matrix(0, N, N, sparse = TRUE))
    }
    
    # 각 노드-클러스터 멤버십 Z (N x K, one-hot, sparse)
    # factor로 클러스터 인덱스 압축 (라벨 값 자체와 무관 → 라벨 불변)
    f <- as.integer(factor(lab_use))
    Z <- sparseMatrix(i = idx_use, j = f, x = 1, dims = c(N, max(f)))
    
    # 같은 클러스터면 A_ij > 0 (대각 0으로)
    A <- tcrossprod(Z)                    # N x N (정수)
    diag(A) <- 0
    (A > 0) + 0                          # binary sparse (0/1)
  }
  
  A_ref  <- make_adj_from_labels(ref_labels)
  A_pred <- make_adj_from_labels(pred_labels)
  
  # -- (3) Dice 계산: 상삼각만 사용(무방향 에지 1회 카운트)
  U_ref  <- triu(A_ref,  k = 1)
  U_pred <- triu(A_pred, k = 1)
  
  # 교집합: 둘 다 1인 위치 수
  # (희소 논리 연산을 위해 >0로 변환한 뒤 곱 대신 논리 AND)
  inter <- sum((U_ref  > 0) & (U_pred > 0))
  s_ref <- sum(U_ref  > 0)
  s_pred<- sum(U_pred > 0)
  
  dice <- if ((s_ref + s_pred) == 0) NA_real_ else (2 * inter) / (s_ref + s_pred)
  
  out <- list(
    n_nodes   = N,
    n_edges_ref  = s_ref,
    n_edges_pred = s_pred,
    tp_edges  = inter,
    dice_edge = dice
  )
  if (return_adj) {
    out$adj_ref_upper  <- U_ref
    out$adj_pred_upper <- U_pred
  }
  out
}

# # ---------------------------
# # 사용 예
# # ---------------------------
# # reference data 를 이용하여 DICE metric 을 계산
# # 0.45... 15개 데이터만 사용해서 그런가?
# 
# res <- dice_from_cluster_adjacency(
#   frontal_ref_craddock, fmri_comparison[[1]],
#   ijk_cols = c("i","j","k"),
#   ref_label_col  = "roi",
#   pred_label_col = "cluster_group121_corr",
#   treat_zero_as_bg = TRUE,
#   return_adj = FALSE
# )

## =========================================================
##  F. similarity metric 별 clustering 결과 비교
## =========================================================

# 벡터 n에 대해 sum_{k} C(n_k,2)
.choose2_sum <- function(n) {
  n <- as.numeric(n)
  sum(n * (n - 1) / 2)
}

# 두 라벨 벡터(정수; 0은 미할당)를 직접 비교해
# |A∩B|/(|A|+|B|) 또는 Dice (=2*...)를 계산
overlap_from_labels <- function(L1, L2, dice = FALSE) {
  L1 <- as.integer(L1); L2 <- as.integer(L2)
  
  # |A|, |B| : 각 라벨링의 "같은 클러스터" 엣지 수
  nA <- .choose2_sum(table(L1[L1 > 0L]))
  nB <- .choose2_sum(table(L2[L2 > 0L]))
  
  # |A∩B| : 두 라벨링 모두에서 같은 클러스터인 쌍의 수
  idx <- (L1 > 0L & L2 > 0L)
  if (!any(idx)) return(NA_real_)
  tab <- table(L1[idx], L2[idx])   # n_{cd}
  nI  <- .choose2_sum(tab)         # sum_{c,d} C(n_cd,2)
  
  denom <- nA + nB
  if (denom == 0) return(NA_real_)
  if (dice) 2 * nI / denom else nI / denom
}

# out_* 객체에서 group-level 라벨만 뽑아 계산
get_group_labels <- function(out) {
  if (!is.null(out$labels)) return(as.integer(out$labels))
  if (!is.null(out$group_labels)) return(as.integer(out$group_labels))
  if (!is.null(out$consensus) && !is.null(out$consensus$labels))
    return(as.integer(out$consensus$labels))
  stop("group-level labels가 필요합니다.")
}

overlap_between_outs <- function(outA, outB, dice = FALSE) {
  LA <- get_group_labels(outA)
  LB <- get_group_labels(outB)
  overlap_from_labels(LA, LB, dice = dice)
}

# 편의: 여러 결과물을 한 번에 pairwise 비교
pairwise_overlap_outs <- function(named_out_list, dice = FALSE) {
  stopifnot(is.list(named_out_list), length(named_out_list) >= 1L)
  nm <- names(named_out_list); m <- length(named_out_list)
  Ls <- lapply(named_out_list, get_group_labels)
  
  M <- matrix(NA_real_, m, m, dimnames = list(nm, nm))
  diag(M) <- 1
  
  if (m >= 2L) {
    for (i in seq_len(m - 1L)) {
      for (j in seq.int(i + 1L, m)) {  # 안전한 시퀀스
        M[i, j] <- M[j, i] <- overlap_from_labels(Ls[[i]], Ls[[j]], dice = dice)
      }
    }
  }
  M
}


# ## =========================================================
# ##  G. 예시 실행: corr / L2 / SRVF 각각에 대해 LOOCV & SI
# ## =========================================================
# # fmri_list: subject별 데이터프레임 리스트 (모두 동일 voxel 집합/좌표/ijk)
# 
# # correlation 기반
# out_corr_smooth <- craddock_group_parcellation(
#   fmri_list_smooth, # gaussian kernel 로 smoothing 함
#   K = 100,
#   i_col = 3, j_col = 4, k_col = 5, ijk_col = 485,
#   time_cols = 6:483,
#   weight_fun = weight_corr_by_rows,
#   # prestandardized=FALSE는 weight_fun 인자로 ...에 전달
#   prestandardized = TRUE,
#   clamp_negative = NULL,   # 자동: 음수 있으면 clamp
#   r_min = 0.0,
#   chunk_size = 10000L,
#   remove_isolated = TRUE,
#   kmeans_itermax = 1000,
#   seed = 2025,
#   verbose = TRUE
# )
# 
# # partial correlation 기반
# out_pcorr_smooth <- craddock_group_parcellation(
#   fmri_list_smooth,
#   K = 100,
#   i_col = 3, j_col = 4, k_col = 5, ijk_col = 485,
#   time_cols = 6:483,
#   weight_fun = weight_pcorr_by_rows,
#   # prestandardized=FALSE는 weight_fun 인자로 ...에 전달
#   prestandardized = TRUE,
#   clamp_negative = NULL,   # 자동: 음수 있으면 clamp
#   r_min = 0.0,
#   chunk_size = 10000L,
#   remove_isolated = TRUE,
#   kmeans_itermax = 1000,
#   seed = 2025,
#   verbose = TRUE
# )
# 
# 
# # l2 distance with rbf
# out_l2_smooth <- craddock_group_parcellation(
#   fmri_list_smooth,
#   K = 100,
#   i_col = 3, j_col = 4, k_col = 5, ijk_col = 485,
#   time_cols = 6:483,
#   weight_fun = weight_l2_ts,
#   standardize = TRUE,      # weight_l2_ts 인자
#   kernel = "rbf",          # "rbf" | "inv" | "none"(비권장: 거리 그대로)
#   sigma = NULL,            # NULL이면 파일럿에서 자동 추정
#   clamp_negative = NULL,   # 자동(보통 비음수라 clamp 안 함)
#   r_min = 0.0,
#   chunk_size = 10000L,
#   remove_isolated = TRUE,
#   kmeans_itermax = 1000,
#   seed = 2025,
#   verbose = TRUE
# )
# 
# # srvf distance with rbf
# out_srvf_smooth <- craddock_group_parcellation(
#   fmri_list_smooth,
#   K = 100,
#   i_col = 3, j_col = 4, k_col = 5, ijk_col = 485,
#   time_cols = 6:483,
#   weight_fun = weight_srvf_ts,
#   warp = "none",           # "dtw"는 매우 느림
#   smooth = FALSE,
#   kernel = "rbf",
#   sigma = NULL,
#   dt = 1.0,
#   clamp_negative = NULL,   # 자동
#   r_min = 0.0,
#   chunk_size = 10000L,     # q 길이 = T-1
#   remove_isolated = TRUE,
#   kmeans_itermax = 1000,
#   seed = 2025,
#   verbose = TRUE
# )
# 
# 
# # 최종 그룹 라벨을 모든 데이터프레임에 부여 (원하면 하나에만 붙여도 OK)
# for (i in seq_along(fmri_list)) {
#   fmri_list_smooth[[i]]$cluster_group100_corr <- out_corr_smooth$group_labels
# }
# for (i in seq_along(fmri_list)) {
#   fmri_list_smooth[[i]]$cluster_group100_l2 <- out_l2_smooth$group_labels
# }
# for (i in seq_along(fmri_list)) {
#   fmri_list_smooth[[i]]$cluster_group100_srvf <- out_srvf_smooth$group_labels
# }
# for (i in seq_along(fmri_list)) {
#   fmri_list_smooth[[i]]$cluster_group100_pcorr <- out_pcorr_smooth$group_labels
# }
# 
# # 확인
# table(fmri_list[[1]]$group_labels)
# table(fmri_list_smooth[[1]]$cluster_group100_corr)
# 
# p_corr <- plot_voxels_by_cluster(
#   fmri_list_smooth[[1]],
#   cluster_col = "cluster_group100_corr",
#   i_col = 3, j_col = 4, k_col = 5,
#   showlegend = TRUE,   # 200개 클러스터면 범례 OFF 권장
#   point_size = 5
# )
# 
# p_pcorr <- plot_voxels_by_cluster(
#   fmri_list[[1]],
#   cluster_col = "cluster_group200_pcorr",
#   i_col = 3, j_col = 4, k_col = 5,
#   showlegend = TRUE,   # 200개 클러스터면 범례 OFF 권장
#   point_size = 5
# )
# 
# p_l2 <- plot_voxels_by_cluster(
#   fmri_list[[1]],
#   cluster_col = "cluster_group200_l2",
#   i_col = 3, j_col = 4, k_col = 5,
#   showlegend = TRUE,   # 200개 클러스터면 범례 OFF 권장
#   point_size = 5
# )
# 
# p_srvf  <- plot_voxels_by_cluster(
#   fmri_list[[1]],
#   cluster_col = "cluster_group200_srvf",
#   i_col = 3, j_col = 4, k_col = 5,
#   showlegend = TRUE,   # 200개 클러스터면 범례 OFF 권장
#   point_size = 5
# )
# 
# 
# ## 1) Dice-LOOCV
# loo_corr <- loocv_dice_from_out(out_corr_smooth, verbose = TRUE)
# loo_pcorr <- loocv_dice_from_out(out_pcorr_smooth, verbose = TRUE)
# loo_l2   <- loocv_dice_from_out(out_l2_smooth,   verbose = TRUE)
# loo_srvf <- loocv_dice_from_out(out_srvf_smooth, verbose = TRUE)
# 
# ## 2) 실루엣(그룹 라벨 고정, 피험자별 유사도에서 평가)
# si_corr <- silhouette_per_subjects(
#   fmri_list_smooth, out_corr_smooth$edges, out_corr_smooth$group_labels,
#   weight_fun     = weight_corr_by_rows,
#   time_cols      = 6:483,
#   clamp_negative = TRUE,    # 상관 유사도는 보통 음수 clamp
#   r_min          = 0.0,
#   chunk_size     = 50000L,
#   prestandardized = TRUE,  # weight_corr_by_rows 인자 전달
#   verbose        = TRUE
# )
# 
# si_pcorr <- silhouette_per_subjects(
#   fmri_list_smooth, out_pcorr_smooth$edges, out_pcorr_smooth$group_labels,
#   weight_fun     = weight_pcorr_by_rows,
#   time_cols      = 6:483,
#   clamp_negative = TRUE,    # 상관 유사도는 보통 음수 clamp
#   r_min          = 0.0,
#   chunk_size     = 50000L,
#   prestandardized = TRUE,  # weight_corr_by_rows 인자 전달
#   verbose        = TRUE
# )
# 
# si_l2 <- silhouette_per_subjects(
#   fmri_list_smooth, out_l2_smooth$edges, out_l2_smooth$group_labels,
#   weight_fun     = weight_l2_ts,
#   time_cols      = 6:483,
#   standardize    = TRUE,
#   kernel         = "rbf",
#   sigma          = NULL,    # 파일럿 자동
#   clamp_negative = FALSE,   # 거리→유사도는 비음수
#   r_min          = 0.0,
#   chunk_size     = 50000L,
#   verbose        = TRUE
# )
# 
# si_srvf <- silhouette_per_subjects(
#   fmri_list_smooth, out_srvf_smooth$edges, out_srvf_smooth$group_labels,
#   weight_fun     = weight_srvf_ts,
#   time_cols      = 6:483,
#   warp           = "none",
#   smooth         = FALSE,
#   kernel         = "rbf",
#   sigma          = NULL,
#   dt             = 1.0,
#   clamp_negative = FALSE,
#   r_min          = 0.0,
#   chunk_size     = 20000L,
#   verbose        = TRUE
# )
# 
# ## 3) 요약 테이블
# metrics_summary <- data.frame(
#   method   = c("corr(rt)", "pcorr(rs)", "L2", "SRVF"),
#   dice_mean = c(loo_corr$mean, loo_pcorr$mean, loo_l2$mean, loo_srvf$mean),
#   dice_sd   = c(loo_corr$sd, loo_pcorr$sd, loo_l2$sd, loo_srvf$sd),
#   si_mean   = c(si_corr$mean, si_pcorr$mean, si_l2$mean, si_srvf$mean),
#   si_sd     = c(si_corr$sd, si_pcorr$sd, si_l2$sd, si_srvf$sd)
# )
# 
# print(metrics_summary)
# 
# ## =========================================================
# ##  G. 예시 실행: corr / pcorr / L2 / SRVF 이 얼마나 비슷하게 
# ##    clustering 했는가? Dice metric 사용하여 확인
# ## =========================================================
# 
# # 한 번에 행렬로
# M_dice    <- pairwise_overlap_outs(list(corr = out_corr_smooth, pcorr = out_pcorr_smooth, l2 = out_l2_smooth, srvf = out_srvf_smooth), dice = TRUE)
# M_dice
