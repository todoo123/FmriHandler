scree_plot_rspectra <- function(
  eig,
  title      = "Scree plot (top-K eigenvalues)",
  show_cum   = TRUE,
  normalize  = c("none", "possum", "absum"),
  save_path  = NULL,          # 예: "plots/scree_k50.png"
  width      = 1600,          # px for png/jpeg, inches for pdf
  height     = 1000,
  dpi        = 150,           # png/jpeg 에만 적용
  bg         = "white"
) {
  stopifnot(is.list(eig), !is.null(eig$values))
  normalize <- match.arg(normalize)

  # 1) eigenvalues (ensure decreasing, real parts)
  vals <- as.numeric(Re(eig$values))
  ord  <- order(vals, decreasing = TRUE)
  vals <- vals[ord]

  # 2) normalize 옵션 처리
  if (normalize == "none") {
    y <- vals
    ylab <- "Eigenvalue"
    cum <- NULL
  } else {
    denom <- switch(normalize,
                    possum = sum(pmax(vals, 0)),
                    absum  = sum(abs(vals)))
    if (is.na(denom) || denom <= 0) {
      warning("Normalization denominator <= 0; falling back to raw eigenvalues.")
      y <- vals; ylab <- "Eigenvalue"; cum <- NULL
    } else {
      prop <- vals / denom
      y    <- 100 * prop
      ylab <- "Percentage (%)"
      cum  <- 100 * cumsum(prop)
    }
  }

  k <- seq_along(vals)

  # --- helper: open/close graphics device if save_path given ---
  dev_opened <- FALSE
  close_dev  <- function() if (dev_opened) grDevices::dev.off()

  if (!is.null(save_path)) {
    # 디렉토리 없으면 생성
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ext <- tolower(tools::file_ext(save_path))
    if (ext %in% c("png", "jpg", "jpeg")) {
      # width/height in px
      if (ext == "png") {
        grDevices::png(filename = save_path, width = width, height = height,
                       units = "px", res = dpi, bg = bg, type = "cairo")
      } else {
        grDevices::jpeg(filename = save_path, width = width, height = height,
                        units = "px", res = dpi, quality = 95, bg = bg, type = "cairo")
      }
      dev_opened <- TRUE
    } else if (ext == "pdf") {
      # width/height in inches for pdf
      grDevices::pdf(file = save_path, width = width/96, height = height/96, onefile = FALSE, paper = "special", bg = bg, useDingbats = FALSE)
      dev_opened <- TRUE
    } else {
      warning("Unknown extension: use .png, .jpg/.jpeg, or .pdf. No file will be saved.")
    }
    on.exit(close_dev(), add = TRUE)
  }

  # 3) draw plot
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  plot(k, y, type = "b", pch = 19, xlab = "Component (ranked)",
       ylab = ylab, main = title)
  grid()

  if (show_cum && !is.null(cum)) {
    par(new = TRUE)
    plot(k, cum, type = "s", lwd = 2, axes = FALSE, xlab = "", ylab = "", ylim = c(0, 100))
    axis(side = 4); mtext("Cumulative (%)", side = 4, line = 3)
    legend("topright",
           legend = c(if (normalize == "none") "Eigenvalue" else "Percentage", "Cumulative %"),
           lty = c(1,1), pch = c(19, NA), lwd = c(1,2), bty = "n")
  }

  invisible(list(values = vals, y = y, cum = cum, saved = !is.null(save_path), path = save_path))
}
