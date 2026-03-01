# # basis expansion with fourier basis
# setwd('/Users/hyunseokyoon/Desktop/학교/대학원/Project/인문사회융합연구/code')
# setwd('/Users/hyunseokyoon/Desktop/학교/대학원/Project/인문사회융합연구/sample_data')
# load('/Users/todoo/Desktop/학교/대학원/Project/인문사회융합연구/code/fmri_list')
# 
# 
# matplot(as.matrix(fmri_list[[1]][1:100, 6:483]), type = 'l')
# 
# for(i in 1:100){
#   ts.plot(as.numeric(fmri_list[[1]][i, 6:483]), ylim = c(-200, 200))
# }
# 
# 
# library(tseries)
# output <- 0
# for(i in 1:length(fmri_list)){
#   for(j in 1:dim(fmri_list[[1]])[1]){
#     res<-adf.test(as.numeric(fmri_list[[i]][j, 6:483]))
#     if(is.na(res$p.value)){
#       next
#     }
#     if(res$p.value < 0.05){
#       output <- output + 1
#     }
#   }
#   print(i)
# }
# output # 123660
# 15*dim(fmri_list[[1]])[1] # 533850
# # 꽤 많은 signal 이 정상성 가정을 충족하지 않음.
# # wavelet basis 를 사용하는 것이 좋겠음

# 원래 사용하던 smoothing 방법
# gaussian kernel
# FWHM -> sigma (in samples)
################################################################################
# gaussian kernel smoothing
################################################################################
fwhm_to_sigma <- function(fwhm_samples) fwhm_samples / 2.354820045

# 1D Gaussian kernel (sigma in samples), truncated at ±(truncate * sigma)
gaussian_kernel_1d <- function(sigma_samples, truncate = 3) {
  if (sigma_samples <= 0) return(1)  # no smoothing
  r <- ceiling(truncate * sigma_samples)
  x <- -r:r
  k <- exp(-0.5 * (x / sigma_samples)^2)
  k / sum(k)
}

# 1D convolution with reflect padding
conv1d_reflect <- function(x, k) {
  if (length(k) == 1L) return(x)
  n <- length(x)
  r <- (length(k) - 1L) %/% 2L
  left  <- x[seq(r, 1, by = -1)]
  right <- x[seq(n - r + 1, n)]
  padded <- c(left, x, right)
  y <- stats::filter(padded, k, sides = 2, circular = FALSE)
  y[(r + 1):(r + n)]
}

# 스무딩(가우시안; 입력: FWHM in samples)
smooth_time_by_fwhm_samples <- function(ts, fwhm_samples, truncate = 3) {
  sigma <- fwhm_to_sigma(fwhm_samples)
  k <- gaussian_kernel_1d(sigma, truncate = truncate)
  conv1d_reflect(as.numeric(ts), k)
}


# 한 개의 data.frame(한 피험자)에 대해: 행(=voxel)별로 시계열 스무딩
smooth_df_timeseries <- function(df,
                                 time_cols = 6:483,
                                 fwhm_samples = 3,
                                 truncate = 3,
                                 chunk_size = 20000L,
                                 verbose = TRUE) {
  
  stopifnot(all(time_cols %in% seq_len(ncol(df))))
  X <- as.matrix(df[, time_cols, drop = FALSE])
  n_rows <- nrow(X)
  Tlen   <- ncol(X)
  
  if (verbose) {
    message(sprintf("Smoothing %d voxels × %d timepoints (FWHM=%g samples)", 
                    n_rows, Tlen, fwhm_samples))
  }
  
  # 미리 커널 한 번만 생성
  sigma <- fwhm_to_sigma(fwhm_samples)
  k <- gaussian_kernel_1d(sigma, truncate = truncate)
  
  # 출력 버퍼
  X_out <- matrix(NA_real_, nrow = n_rows, ncol = Tlen)
  
  # chunk 단위 처리
  idx_starts <- seq(1L, n_rows, by = chunk_size)
  for (s in seq_along(idx_starts)) {
    i1 <- idx_starts[s]
    i2 <- min(i1 + chunk_size - 1L, n_rows)
    if (verbose) message(sprintf("  chunk %d/%d: rows %d..%d", s, length(idx_starts), i1, i2))
    
    # 각 행에 스무딩 적용
    for (i in i1:i2) {
      X_out[i, ] <- conv1d_reflect(X[i, ], k)
    }
  }
  
  # 결과를 원래 df에 반영
  df[, time_cols] <- X_out
  df
}

# fmri_list 전체(모든 피험자)에 적용
smooth_all_subjects <- function(fmri_list,
                                time_cols = 6:483,
                                fwhm_samples = 3,
                                truncate = 3,
                                chunk_size = 20000L,
                                verbose = TRUE) {
  out_list <- vector("list", length(fmri_list))
  for (s in seq_along(fmri_list)) {
    if (verbose) message(sprintf("\nSubject %d / %d", s, length(fmri_list)))
    out_list[[s]] <- smooth_df_timeseries(fmri_list[[s]],
                                          time_cols = time_cols,
                                          fwhm_samples = fwhm_samples,
                                          truncate = truncate,
                                          chunk_size = chunk_size,
                                          verbose = verbose)
  }
  out_list
}

# -------- 사용 예시 --------
# 1) 한 피험자만
# fmri_list[[1]] <- smooth_df_timeseries(fmri_list[[1]], time_cols = 6:483, fwhm_samples = 3)

# 2) 전 피험자 일괄
# fmri_list_smooth <- smooth_all_subjects(fmri_list, time_cols = 6:483, fwhm_samples = 7)
# index <- fmri_list_smooth[[1]]$cluster_group200_corr == 1
# matplot(t(as.matrix(fmri_list_smooth[[1]][index,6:483])), type = 'l')
# matplot(t(as.matrix(fmri_list[[1]][index,6:483])), type = 'l')
# 
# 
# matplot(t(as.matrix(fmri_list_smooth[[1]][1:10,6:483])), type = 'l')
# matplot(t(as.matrix(fmri_list[[1]][1:10,6:483])), type = 'l')



################################################################################
# bspline basis expansion smoothing
################################################################################
# load('/Users/hyunseokyoon/Desktop/학교/대학원/Project/인문사회융합연구/code/fmri_list')
# =============================================================================
# B-spline 기반 basis expansion 계수표 생성기 (chunk 처리)
# - 입력: fmri_list (리스트; 각 원소는 data.frame)
# - 출력: fmri_basis_coeff_list (리스트; 각 원소는 data.frame, 앞 5열 메타+계수열)
# =============================================================================
library(fda)
library(MASS)

build_bspline_coeff_list <- function(
    fmri_list,
    i_col = 1, j_col = 2, k_col = 3,
    time_cols = 6:483,
    Kint = 60, norder = 4, lambda = 1e-8, Lfd = 2,
    chunk_size = 2000L, verbose = TRUE
){
  out_list <- vector("list", length(fmri_list))
  for (s in seq_along(fmri_list)) {
    df <- fmri_list[[s]]
    meta <- df[, c(i_col, j_col, k_col)]
    X <- as.matrix(df[, time_cols]); storage.mode(X) <- "double"
    N <- nrow(X); Tlen <- ncol(X); pts <- seq(0, 1, length.out = Tlen)
    
    breaks <- seq(0, 1, length.out = Kint + 2)
    basis  <- create.bspline.basis(c(0,1), norder = norder, breaks = breaks)
    fp     <- fdPar(basis, Lfdobj = Lfd, lambda = lambda)
    nbasis <- basis$nbasis
    if (isTRUE(verbose)) message(sprintf("[Subj %d] N=%d T=%d nbasis=%d", s, N, Tlen, nbasis))
    
    coef_blocks <- vector("list", ceiling(N / chunk_size))
    start <- 1L; b <- 1L
    while (start <= N) {
      end <- min(N, start + chunk_size - 1L)
      fit <- smooth.basis(pts, t(X[start:end, , drop = FALSE]), fp)
      coef_blocks[[b]] <- t(fit$fd$coefs)     # (chunk x nbasis)
      start <- end + 1L; b <- b + 1L
    }
    coef_mat <- do.call(rbind, coef_blocks)
    colnames(coef_mat) <- sprintf("coef_%02d", seq_len(nbasis))
    
    out_df <- cbind(meta, as.data.frame(coef_mat))
    out_df$ijk <- with(out_df, paste(i, j, k, sep = "_"))  # ← 마지막에 붙임
    
    # 재구성에 필요할 정보 저장
    attr(out_df, "basis") <- basis
    attr(out_df, "pts")   <- pts
    out_list[[s]] <- out_df
  }
  names(out_list) <- names(fmri_list)
  out_list
}
# 
# 
# # =========================
# # 사용
# # =========================
# fmri_basis_coeff_list <- build_bspline_coeff_list(
#   fmri_list,
#   label_col = 1, pfc_col = 2, i_col = 3, j_col = 4, k_col = 5,
#   time_cols = 6:483,
#   Kint = 80, norder = 4, lambda = 1e-8, Lfd = 2,
#   chunk_size = 2000L, verbose = TRUE
# )
# 
# # =========================
# # test
# # =========================
# # ---------------------------
# # 1) 재구성 (예: 피험자 1)
# # ---------------------------
# s <- 1
# coef_cols <- grepl("^coef_", names(fmri_basis_coeff_list[[s]]))
# C     <- as.matrix(fmri_basis_coeff_list[[s]][, coef_cols])    # (N x nbasis)
# basis <- attr(fmri_basis_coeff_list[[s]], "basis")
# pts   <- attr(fmri_basis_coeff_list[[s]], "pts")
# 
# fdobj <- fd(t(C), basis)
# Xhat  <- t(eval.fd(pts, fdobj))         # 근사 시계열
# Xorig <- as.matrix(fmri_list[[s]][, 6:483])
# Res   <- Xorig - Xhat                   # 잔차
# 
# # ---------------------------
# # 2) 샘플 몇 개 골라서 plot + residual
# # ---------------------------
# set.seed(1)
# idx <- sample(seq_len(nrow(Xorig)), 6)
# 
# par(mfrow=c(3,2), mar=c(3,3,2,1))   # 각 곡선: 왼쪽=원본vs근사, 오른쪽=잔차
# for (i in idx) {
#   # (a) 원본 vs 근사
#   ylim <- range(Xorig[i,], Xhat[i,], finite=TRUE)
#   plot(pts, Xorig[i,], type="l", lwd=1,
#        ylim=ylim, main=paste("Curve", i),
#        xlab="t", ylab="signal")
#   lines(pts, Xhat[i,], lwd=2, lty=2, col="red")
#   legend("topright", c("original","fit"), lty=c(1,2), col=c("black","red"), bty="n", cex=0.8)
#   
#   # (b) Residual
#   plot(pts, Res[i,], type="l", lwd=1,
#        main=paste("Residual", i),
#        xlab="t", ylab="residual")
#   abline(h=0, lty=3)
# }
# 
# 
# 
# # =========================
# # clustering
# # =========================
# 
# ## =========================================================
# ##  F. 예시 실행: corr / L2 / SRVF 각각에 대해 LOOCV & SI
# ## =========================================================
# # fmri_list: subject별 데이터프레임 리스트 (모두 동일 voxel 집합/좌표/ijk)
# 
# # correlation 기반
# out_corr_smooth <- craddock_group_parcellation(
#   fmri_basis_coeff_list, # bspline 으로 smoothing 함
#   K = 100,
#   i_col = 3, j_col = 4, k_col = 5, ijk_col = 90,
#   time_cols = 6:89,
#   weight_fun = weight_corr_by_rows,
#   # prestandardized=FALSE는 weight_fun 인자로 ...에 전달
#   prestandardized = TRUE,
#   clamp_negative = NULL,   # 자동: 음수 있으면 clamp
#   r_min = 0.0,
#   chunk_size = 10000L,
#   remove_isolated = TRUE,
#   kmeans_itermax = 1000,
#   verbose = TRUE
# )
# 
# # partial correlation 기반
# out_pcorr_smooth <- craddock_group_parcellation(
#   fmri_basis_coeff_list,
#   K = 100,
#   i_col = 3, j_col = 4, k_col = 5, ijk_col = 90,
#   time_cols = 6:89,
#   weight_fun = weight_pcorr_by_rows,
#   # prestandardized=FALSE는 weight_fun 인자로 ...에 전달
#   prestandardized = TRUE,
#   clamp_negative = NULL,   # 자동: 음수 있으면 clamp
#   r_min = 0.0,
#   chunk_size = 10000L,
#   remove_isolated = TRUE,
#   kmeans_itermax = 1000,
#   verbose = TRUE
# )
# 
# 
# # 최종 그룹 라벨을 모든 데이터프레임에 부여 (원하면 하나에만 붙여도 OK)
# for (i in seq_along(fmri_basis_coeff_list)) {
#   fmri_basis_coeff_list[[i]]$cluster_group100_corr <- out_corr_smooth$group_labels
# }
# for (i in seq_along(fmri_basis_coeff_list)) {
#   fmri_basis_coeff_list[[i]]$cluster_group100_pcorr <- out_pcorr_smooth$group_labels
# }
# 
# 
# # 확인
# table(fmri_basis_coeff_list[[1]]$cluster_group100_corr)
# 
# p_corr <- plot_voxels_by_cluster(
#   fmri_basis_coeff_list[[1]],
#   cluster_col = "cluster_group100_corr",
#   i_col = 3, j_col = 4, k_col = 5,
#   showlegend = TRUE,   # 200개 클러스터면 범례 OFF 권장
#   point_size = 5
# )
# 
# p_pcorr <- plot_voxels_by_cluster(
#   fmri_basis_coeff_list[[1]],
#   cluster_col = "cluster_group100_pcorr",
#   i_col = 3, j_col = 4, k_col = 5,
#   showlegend = TRUE,   # 200개 클러스터면 범례 OFF 권장
#   point_size = 5
# )
# 
# ## 1) Dice-LOOCV
# loo_corr <- loocv_dice_from_out(out_corr_smooth, verbose = TRUE)
# loo_pcorr <- loocv_dice_from_out(out_pcorr_smooth, verbose = TRUE)
# 
# ## 2) 실루엣(그룹 라벨 고정, 피험자별 유사도에서 평가)
# si_corr <- silhouette_per_subjects(
#   fmri_basis_coeff_list, out_corr_smooth$edges, out_corr_smooth$group_labels,
#   weight_fun     = weight_corr_by_rows,
#   time_cols      = 6:89,
#   clamp_negative = TRUE,    # 상관 유사도는 보통 음수 clamp
#   r_min          = 0.0,
#   chunk_size     = 50000L,
#   prestandardized = TRUE,  # weight_corr_by_rows 인자 전달
#   verbose        = TRUE
# )
# 
# si_pcorr <- silhouette_per_subjects(
#   fmri_basis_coeff_list, out_pcorr_smooth$edges, out_pcorr_smooth$group_labels,
#   weight_fun     = weight_pcorr_by_rows,
#   time_cols      = 6:89,
#   clamp_negative = TRUE,    # 상관 유사도는 보통 음수 clamp
#   r_min          = 0.0,
#   chunk_size     = 50000L,
#   prestandardized = TRUE,  # weight_corr_by_rows 인자 전달
#   verbose        = TRUE
# )
# 
# 
# ## 3) 요약 테이블
# metrics_summary <- data.frame(
#   method   = c("corr(rt)", "pcorr(rs)"),
#   dice_mean = c(loo_corr$mean, loo_pcorr$mean),
#   dice_sd   = c(loo_corr$sd, loo_pcorr$sd),
#   si_mean   = c(si_corr$mean, si_pcorr$mean),
#   si_sd     = c(si_corr$sd, si_pcorr$sd)
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
# M_dice    <- pairwise_overlap_outs(list(corr = out_corr_smooth, pcorr = out_pcorr_smooth), dice = TRUE)
# M_dice
# 
# 
# # ------------------------------------------------------------
# # 클러스터별 재구성 곡선 시각화 & 파일 저장
# # ------------------------------------------------------------
# s      <- 1   # 피험자 번호
# dat    <- fmri_basis_coeff_list[[s]]
# basis  <- attr(dat, "basis")
# pts    <- attr(dat, "pts")
# 
# coef_cols <- grepl("^coef_", names(dat))   # 계수 열
# cluster_col <- "cluster_group100_pcorr"      # 원하는 cluster 변수명
# 
# outdir <- "bspline_pcorr_cluster100_1-100"                  # 저장 디렉토리
# if (!dir.exists(outdir)) dir.create(outdir)
# 
# for (k in 1:100) {
#   idx <- dat[[cluster_col]] == k
#   if (sum(idx) == 0) next   # 빈 클러스터면 skip
#   
#   C <- as.matrix(dat[idx, coef_cols])
#   fdobj <- fd(t(C), basis)
#   Xhat  <- t(eval.fd(pts, fdobj))    # (n_voxel x T)
#   
#   png(file.path(outdir, sprintf("cluster_%02d.png", k)),
#       width = 1200, height = 800)
#   matplot(pts, t(Xhat), type = "l", lty = 1,
#           xlab = "t", ylab = "signal",
#           main = sprintf("Subject %d - Cluster %d", s, k))
#   dev.off()
# }


################################################################################
# fPCA basis expansion smoothing
################################################################################
# library(fda)
# library(MASS)

build_fpca_loading_list <- function(
    fmri_list,
    i_col = 1, j_col = 2, k_col = 3,
    time_cols = 6:483,
    Kpc = 60,            # 유지할 주성분(기저) 개수
    center = FALSE,       # 시간축 평균 함수(열평균) 제거
    scale = FALSE,       # 필요시 시간축 표준화(보통 FALSE 권장)
    chunk_size = 2000L,
    verbose = TRUE,
    use_RS = FALSE        # RSpectra 사용 여부(대규모 T일 때 빠름)
){
  stopifnot(is.list(fmri_list))
  out_list <- vector("list", length(fmri_list))
  
  for (s in seq_along(fmri_list)) {
    df <- fmri_list[[s]]
    meta <- df[, c(i_col, j_col, k_col)]
    X <- as.matrix(df[, time_cols, drop = FALSE]); storage.mode(X) <- "double"
    N <- nrow(X); Tlen <- ncol(X)
    pts <- seq(0, 1, length.out = Tlen)
    
    # (1) 시간축 센터링/스케일링
    mu_t <- if (center) colMeans(X) else rep(0, Tlen)
    Xc   <- sweep(X, 2L, mu_t, FUN = "-")
    if (scale) {
      sd_t <- apply(Xc, 2L, sd)
      sd_t[sd_t == 0 | !is.finite(sd_t)] <- 1
      Xc <- sweep(Xc, 2L, sd_t, FUN = "/")
    }
    
    # (2) 시간–시간 공분산(= T×T)과 고유분해
    #     cov(Xc)는 행=관측(N), 열=변수(T) 기준이라 T×T 결과가 나옴
    #     큰 N에서도 T가 수백이라면 직접 S = crossprod(Xc)/(N-1)이 효율적
    S <- crossprod(Xc) / (N - 1)   # T × T
    
    Kkeep <- min(Kpc, Tlen)
    if (isTRUE(verbose)) {
      message(sprintf("[Subj %d] N=%d T=%d  -> keep %d PCs", s, N, Tlen, Kkeep))
    }
    
    eigvals <- NULL
    E <- NULL
    if (use_RS) {
      # 빠른 상위 K개 고유쌍
      EIG <- tryCatch(
        RSpectra::eigs_sym(S, k = Kkeep, which = "LA"),
        error = function(e) NULL
      )
      if (!is.null(EIG)) {
        eigvals <- as.numeric(EIG$values)
        E <- as.matrix(EIG$vectors)  # T × Kkeep, 열 직교정규
      }
    }
    if (is.null(E)) {
      # fallback: base::eigen (전체 분해)
      EIG <- eigen(S, symmetric = TRUE)
      eigvals <- as.numeric(EIG$values[seq_len(Kkeep)])
      E <- as.matrix(EIG$vectors[, seq_len(Kkeep), drop = FALSE])
    }
    
    # (3) 로딩(score) 계산: (X - mu_t) %*% E  (청크 처리)
    loading_blocks <- vector("list", ceiling(N / chunk_size))
    start <- 1L; b <- 1L
    while (start <= N) {
      end <- min(N, start + chunk_size - 1L)
      X_chunk <- X[start:end, , drop = FALSE]
      Xc_chunk <- sweep(X_chunk, 2L, mu_t, FUN = "-")
      if (scale) Xc_chunk <- sweep(Xc_chunk, 2L, sd_t, FUN = "/")
      # 로딩(=score): (chunk × T) %*% (T × K) = (chunk × K)
      loading_blocks[[b]] <- Xc_chunk %*% E
      start <- end + 1L; b <- b + 1L
    }
    loading_mat <- do.call(rbind, loading_blocks)
    colnames(loading_mat) <- sprintf("loading_%02d", seq_len(ncol(loading_mat)))
    
    # 출력 df (메타 + 로딩)
    out_df <- cbind(meta, as.data.frame(loading_mat))
    # ijk 키 덧붙이기
    out_df$ijk <- with(out_df, paste(i, j, k, sep = "_"))
    
    # 재구성 및 진단에 필요한 정보 저장
    # f_hat = mu_t + (loading_row %*% t(E))  (scale=TRUE면 역스케일 필요)
    attr(out_df, "eigenvectors") <- E            # T × K
    attr(out_df, "eigenvalues")  <- eigvals      # 길이 K
    attr(out_df, "mu_t")         <- mu_t         # 길이 T
    if (scale) attr(out_df, "sd_t") <- sd_t
    attr(out_df, "pts")          <- pts
    attr(out_df, "var_explained") <- eigvals / sum(diag(S))  # 대략적 비율
    out_list[[s]] <- out_df
  }
  names(out_list) <- names(fmri_list)
  out_list
}

# helper
reconstruct_fpca_curves <- function(fpca_df, original_df, start, end) {
  eigenfun   <- attr(fpca_df, "eigenvectors")
  eigenval   <- attr(fpca_df, "eigenvalues")
  loading <- as.matrix(original_df[,start:end]) %*% eigenfun
  recon_fun <- eigenfun %*% t(loading)
  return(recon_fun)
}


# usage
# fpca_list <- build_fpca_loading_list(
#   fmri_list,
#   time_cols = 6:483,
#   Kpc = 70,                # 예: 상위 70개만 유지
#   center = TRUE, scale = FALSE,
#   use_RS = FALSE,
#   chunk_size = 2000L, verbose = TRUE
# )
# 
# 
# eigenfun   <- attr(fpca_list[[1]], "eigenvectors")
# eigenval   <- attr(fpca_list[[1]], "eigenvalues")
# 
# loading <- as.matrix(fmri_list[[1]][,6:483]) %*% eigenfun
# sum(loading < 0)
# 
# recon_fun <- eigenfun %*% t(loading)
# 
# # 복원
# s <- 3
# Xhat <- reconstruct_fpca_curves(fpca_list[[s]], fmri_list[[s]], start = 6, end = 483)
# dim(Xhat)
# # 시각화(한 곡선)
# x_orig <- as.matrix(fmri_list[[s]][1, 6:483])
# x_hat  <- Xhat[,1]
# ts.plot(cbind(t(x_orig), x_hat), lty = c(1,1), col = c("black", "red"),
#         xlab = "Time", ylab = "Value", main = "Original (black) vs fPCA-smoothed (red)")
# 
# 
