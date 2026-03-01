source("/home/hyunseok/code/ROI_detection.R")
source("/home/hyunseok/code/fmri_smoothing.R")
load("/home/hyunseok/data/fmri_list.RData")
ref_craddock<-read.csv("/home/hyunseok/data/Craddock200_2mm.csv")

# 10000 개 sampling 해서 보는 경우 - 계산상 부담 때문
# truncated_list <- fmri_list[[1]][sample(c(1:dim(fmri_list[[1]])[1]),10000),]
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
fmri_smooth_fpca <- build_fpca_loading_list(fmri_comparison,
                        i_col = 1, j_col = 2, k_col = 3,
                        time_cols = 6:483,
                        Kpc = k,            # 유지할 주성분(기저) 개수
                        center = FALSE,       # 시간축 평균 함수(열평균) 제거
                        scale = FALSE,       # 필요시 시간축 표준화(보통 FALSE 권장)
                        chunk_size = 2000L,
                        verbose = TRUE,
                        use_RS = TRUE        # RSpectra 사용 여부(대규모 T일 때 빠름)
                        )

test_list_fpca <- fmri_smooth_fpca

K <- 20
# fPCA - L2
## restricted
restricted_fPCA_l2 <- craddock_group_parcellation(
  test_list_fpca, # gaussian kernel 로 smoothing 함
  K = K,
  i_col = 1, j_col = 2, k_col = 3, ijk_col = k + 4,
  time_cols = 4:(k + 3),
  weight_fun = weight_l2_ts,
  # prestandardized=FALSE는 weight_fun 인자로 ...에 전달
  standardize = FALSE,
  clamp_negative = NULL,   # 자동: 음수 있으면 clamp
  r_min = 0.55, # creddock 이 제시한 제한 사항
  chunk_size = 10000L,
  remove_isolated = TRUE,
  kmeans_itermax = 1000,
  seed = 2025,
  verbose = TRUE
)

test_list_fpca[[1]]$cluster9_l2 <- restricted_fPCA_l2$group_labels

p <- plot_voxels_by_cluster(
  test_list_fpca[[1]],
  cluster_col = "cluster9_l2",
  i_col = 1, j_col = 2, k_col = 3,
  showlegend = FALSE,   # 200개 클러스터면 범례 OFF 권장
  point_size = 5
)

out_dir <- "/home/hyunseok/plot"
outfile <- file.path(out_dir, "fPCA_restricted_l2_voxels_cluster20.html")
library(htmlwidgets)
saveWidget(
  widget = p,
  file   = outfile,
  selfcontained = FALSE,              # ★ pandoc 없어도 됨
  libdir = file.path(out_dir, "libs") # ★ 필요한 js/css 여기로 빠짐
)

loo_result <- loocv_dice_from_out(restricted_fPCA_l2, i_block_size = 40L)
print("restricted_fPCA_l2")
print(loo_result)

si_score <- silhouette_per_subjects_allpairs(
  test_list_fpca,
  restricted_fPCA_l2$group_labels,
  weight_fun     = weight_l2_ts,
  time_cols      = 4:(k + 3),
  clamp_negative = TRUE,    # 상관 유사도는 보통 음수 clamp
  # creddock setting
  r_min          = 0.55,
  i_block_size   = 40L,
  chunk_size     = 10000L,
  parallel       = "multicore",
  workers        = 20,
  standardize    = FALSE,  # weihght_corr_by_rows 인자 전달
  verbose        = FALSE
)

print(paste0("restricted_fPCA_l2 with cluster number: ", k))
print(si_score$per_subject)
print(si_score$mean)


## allpairs
allpairs_fPCA_l2 <- craddock_group_parcellation_allpairs(
  test_list_fpca,
  K               = K,
  ijk_col         = (k + 4),
  time_cols       = 4:(k + 3),
  weight_fun      = weight_l2_ts,
  standardized    = TRUE,
  r_min           = 0.55,
  i_block_size    = 10L,
  chunk_size      = 10000L,
  parallel        = "multicore",
  workers         = 20,
  tame_blas       = TRUE,
  seed            = 2025,
  verbose         = FALSE,
  screeplot = TRUE, screeplot_name = "fPCA - l2"
)

test_list_fpca[[1]]$cluster9_l2 <- allpairs_fPCA_l2$group_labels

p <- plot_voxels_by_cluster(
  test_list_fpca[[1]],
  cluster_col = "cluster9_l2",
  i_col = 1, j_col = 2, k_col = 3,
  showlegend = FALSE,   # 200개 클러스터면 범례 OFF 권장
  point_size = 5
)

out_dir <- "/home/hyunseok/plot"
outfile <- file.path(out_dir, "fPCA_allpairs_l2_voxels_cluster20.html")
library(htmlwidgets)
saveWidget(
  widget = p,
  file   = outfile,
  selfcontained = FALSE,              # ★ pandoc 없어도 됨
  libdir = file.path(out_dir, "libs") # ★ 필요한 js/css 여기로 빠짐
)


loo_result <- loocv_dice_from_out(allpairs_fPCA_l2, i_block_size = 40L)
print("allpairs_fPCA_l2")
print(loo_result)

si_score <- silhouette_per_subjects_allpairs(
  test_list_fpca, allpairs_fPCA_l2$group_labels,
  weight_fun     = weight_l2_ts,
  time_cols      = 4:(k + 3),
  clamp_negative = TRUE,    # 상관 유사도는 보통 음수 clamp
  r_min          = 0.55,
  i_block_size    = 10L,
  chunk_size     = 10000L,
  parallel = "multicore",
  workers = 20,
  standardized   = TRUE,  # weihght_corr_by_rows 인자 전달
  verbose        = FALSE
)
print(paste0("allpairs_fPCA_l2 with cluster number: ", k))
print(si_score$per_subject)
print(si_score$mean)

# fPCA - cosine
## restricted
restricted_fPCA_cos <- craddock_group_parcellation(
  test_list_fpca, # gaussian kernel 로 smoothing 함
  K = K,
  i_col = 1, j_col = 2, k_col = 3, ijk_col = (k + 4),
  time_cols = 4:(k + 3),
  weight_fun = weight_l2_ts,
  # prestandardized=FALSE는 weight_fun 인자로 ...에 전달
  standardize = FALSE,
  clamp_negative = NULL,   # 자동: 음수 있으면 clamp
  r_min = 0.55, # creddock 이 제시한 제한 사항
  chunk_size = 10000L,
  remove_isolated = TRUE,
  kmeans_itermax = 1000,
  seed = 2025,
  verbose = TRUE
)

test_list_fpca[[1]]$cluster9_cos <- restricted_fPCA_cos$group_labels

p <- plot_voxels_by_cluster(
  test_list_fpca[[1]],
  cluster_col = "cluster9_cos",
  i_col = 1, j_col = 2, k_col = 3,
  showlegend = FALSE,   # 200개 클러스터면 범례 OFF 권장
  point_size = 5
)

out_dir <- "/home/hyunseok/plot"
outfile <- file.path(out_dir, "fPCA_restricted_cos_voxels_cluster20.html")
library(htmlwidgets)
saveWidget(
  widget = p,
  file   = outfile,
  selfcontained = FALSE,              # ★ pandoc 없어도 됨
  libdir = file.path(out_dir, "libs") # ★ 필요한 js/css 여기로 빠짐
)

loo_result <- loocv_dice_from_out(restricted_fPCA_cos, i_block_size = 40L)
print("restricted_fPCA_cos")
print(loo_result)

si_score <- silhouette_per_subjects_allpairs(
  test_list_fpca,
  restricted_fPCA_cos$group_labels,
  weight_fun     = weight_l2_ts,
  time_cols      = 4:(k + 3),
  clamp_negative = TRUE,    # 상관 유사도는 보통 음수 clamp
  # creddock setting
  r_min          = 0.55,
  i_block_size   = 40L,
  chunk_size     = 10000L,
  parallel       = "multicore",
  workers        = 20,
  standardize    = FALSE,  # weihght_corr_by_rows 인자 전달
  verbose        = FALSE
)

print(paste0("restricted_fPCA_cos with cluster number: ", k))
print(si_score$per_subject)
print(si_score$mean)

## allpairs
allpairs_fPCA_cos <- craddock_group_parcellation_allpairs(
  test_list_fpca,
  K               = K,
  ijk_col         = k + 4,
  time_cols       = 4:(k + 3),
  weight_fun      = weight_cosine_ts,
  standardize     = FALSE,
  r_min           = 0.55,
  i_block_size    = 10L,
  chunk_size      = 10000L,
  parallel        = "multicore",
  workers         = 20,
  tame_blas       = TRUE,
  seed            = 2025,
  verbose         = FALSE,
  screeplot = TRUE, screeplot_name = "fPCA - cosine"
)

test_list_fpca[[1]]$cluster9_cos <- allpairs_fPCA_cos$group_labels

p <- plot_voxels_by_cluster(
  test_list_fpca[[1]],
  cluster_col = "cluster9_cos",
  i_col = 1, j_col = 2, k_col = 3,
  showlegend = FALSE,   # 200개 클러스터면 범례 OFF 권장
  point_size = 5
)

out_dir <- "/home/hyunseok/plot"
outfile <- file.path(out_dir, "fPCA_allpairs_cos_voxels_cluster20.html")
library(htmlwidgets)
saveWidget(
  widget = p,
  file   = outfile,
  selfcontained = FALSE,              # ★ pandoc 없어도 됨
  libdir = file.path(out_dir, "libs") # ★ 필요한 js/css 여기로 빠짐
)

loo_result <- loocv_dice_from_out(allpairs_fPCA_cos, i_block_size = 40L)
print("allpairs_fPCA_cos")
print(loo_result)


si_score <- silhouette_per_subjects_allpairs(
  test_list_fpca, allpairs_fPCA_cos$group_labels,
  weight_fun     = weight_cosine_ts,
  time_cols      = 4:(k + 3),
  clamp_negative = TRUE,    # 상관 유사도는 보통 음수 clamp
  r_min          = 0.55,
  i_block_size    = 10L,
  chunk_size     = 10000L,
  parallel = "multicore",
  workers = 20,
  standardized   = TRUE,  # weihght_corr_by_rows 인자 전달
  verbose        = FALSE
)
print(paste0("allpairs_fPCA_cos with cluster number: ", k))
print(si_score$per_subject)
print(si_score$mean)

# rm(out_allpairs_corr)