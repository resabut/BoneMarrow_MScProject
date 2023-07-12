
mod_plotXDEHm <- function(testobj, showRowName = FALSE, cellWidthTotal = 250, cellHeightTotal = 400, showCluster = FALSE, colann = NULL, rowann = NULL, annotation_colors = NULL, subsampleCell = TRUE, numSubsampleCell = 1000, sep = NA, break.0 = TRUE) {
  testvar = testobj$testvar
  fit <- testobj$populationFit
  if ("XDEType" %in% names(testobj)) {
    XDEType <- testobj$XDEType
  } else {
    XDEType <- getXDEType(testobj)
  }
  XDEType <- XDEType[rownames(testobj$covariateGroupDiff)]
  if (subsampleCell) {
    id <- round(seq(1, ncol(fit[[1]]), length.out = numSubsampleCell))
    for (i in 1 : length(fit)) {
      fit[[i]] <- fit[[i]][, id]
    }
    if (sum(XDEType == "meanSig", na.rm = T) > 0) {
      meanid <- which(XDEType == "meanSig")
      max <- apply(abs(testobj$covariateGroupDiff[-meanid, id, drop = F]), 1, max)
      FitDiff.scale1 <- testobj$covariateGroupDiff[-meanid, id, drop = F] / max
      FitDiff.scale <- rbind(FitDiff.scale1, testobj$covariateGroupDiff[meanid, id, drop = F])
      FitDiff.scale <- FitDiff.scale[rownames(testobj$covariateGroupDiff), , drop = F]
    } else {
      max <- apply(abs(testobj$covariateGroupDiff[, id, drop = F]), 1, max)
      FitDiff.scale <- testobj$covariateGroupDiff[, id, drop = F] / max
    }
    FitDiff.sd <- scalematrix(testobj$covariateGroupDiff[, id, drop = F])
    colnames(FitDiff.scale) <- paste0("FitDiff:cell", seq(1, ncol(FitDiff.scale)))
    testobj$pseudotime <- sort(sample(testobj$pseudotime, numSubsampleCell))
    print("subsample done!")
  } else {
    if (sum(XDEType == "meanSig", na.rm = T) > 0) {
      meanid <- which(XDEType == "meanSig")
      max <- apply(abs(testobj$covariateGroupDiff[-meanid, , drop = F]), 1, max)
      FitDiff.scale1 <- testobj$covariateGroupDiff[-meanid, , drop = F] / max
      FitDiff.scale <- rbind(FitDiff.scale1, testobj$covariateGroupDiff[meanid, , drop = F])
      FitDiff.scale <- FitDiff.scale[rownames(testobj$covariateGroupDiff), , drop = F]
    } else {
      max <- apply(abs(testobj$covariateGroupDiff), 1, max)
      FitDiff.scale <- testobj$covariateGroupDiff / max
    }
    colnames(FitDiff.scale) <- paste0("FitDiff:cell", seq(1, ncol(FitDiff.scale)))
    FitDiff.sd <- scalematrix(testobj$covariateGroupDiff)
  }
  max <- apply(abs(testobj$covariateGroupDiff), 1, max)
  alluniformdiff <- testobj$covariateGroupDiff / max
  oridata <- testobj$covariateGroupDiff
  fit.bak = fit
  clu <- testobj$cluster
  rownames(testobj$cellanno) <- testobj$cellanno[, 1]
  testobj$cellanno <- testobj$cellanno[names(testobj$pseudotime), ]
  if ("expr.ori" %in% names(testobj)) {
    testobj$expr <- testobj$expr.ori[, names(testobj$pseudotime)]
  } else {
    testobj$expr <- testobj$expr[, names(testobj$pseudotime)]
  }
  fit.scale <- do.call(cbind, fit)
  fit.scale <- fit.scale[names(testobj$cluster), ]
  fit.scale <- scalematrix(fit.scale)
  colnames(fit.scale) <- paste0(rep(names(fit), each = ncol(fit.scale) / length(fit)), ";cell", seq(1, ncol(fit.scale)))
  changepoint <- sapply(names(clu), function(i) {
    ap <- which(FitDiff.sd[i, -ncol(FitDiff.sd)] * FitDiff.sd[i, -1] < 0)
    ap[which.min(abs(ap - ncol(FitDiff.sd) / 2))]
  })
  res <- data.frame(clu = clu, cor = sapply(names(clu), function(i) cor(FitDiff.scale[i, seq(1, ncol(FitDiff.scale))], seq(1, ncol(FitDiff.scale)))), changepoint = changepoint, XDEType = XDEType[names(clu)])
  res1 <- res[res$XDEType == "trendSig", ]
  res2 <- res[res$XDEType == "bothSig", ]
  res3 <- res[res$XDEType == "other", ]
  res4 <- res[res$XDEType == "meanSig", ]
  o1 <- rownames(res1)[order(res1$clu, res1$cor > 0, res1$changepoint)]
  pn <- rowMeans(alluniformdiff[rownames(res2), , drop = FALSE])
  o2 <- rownames(res2)[order(pn > 0, res2$clu, res2$cor > 0, res2$changepoint)]
  o3 <- rownames(res3)[order(res3$clu, res3$cor > 0, res3$changepoint)]
  o4 <- rownames(res4)[order(res4$clu)]
  res <- res[c(o1, o2, o4), ]
  rle <- rle(paste0(res$clu, res$XDEType))$lengths
  clu[rownames(res)] <- rep(1 : length(rle), rle)
  fit.scale <- fit.scale[rownames(res), ]
  FitDiff.scale <- FitDiff.scale[rownames(res), ]
  cellanno <- testobj$cellanno
  expr = testobj$expr
  expr <- expr[, names(testobj$pseudotime)]
  tmp <- lapply(names(fit), function(i) {
    b <- cellanno[, 2] %in% rownames(testobj$design)[testobj$design[, testvar] == sub(".*_", "", i)]
    a <- cellanno[b, 1]
    expr[rownames(fit.scale), colnames(expr) %in% a]
  })
  expr.scale <- do.call(cbind, tmp)
  expr.scale <- scalematrix(expr.scale)
  expr.scale <- expr.scale[rownames(fit.scale), ]
  if (is.null(colann)) {
    colann <- data.frame(pseudotime = testobj$pseudotime[colnames(expr.scale)], group = as.character(testobj$design[cellanno[match(colnames(expr.scale), cellanno[, 1]), 2], testvar]), expression = "Original", stringsAsFactors = F)
    col.group = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann$group)) + 1)
    names(col.group) = c("NA", unique(colann$group))
  }
  rownames(colann) = colnames(expr.scale)
  col.expression = brewer.pal(n = 8, name = "Pastel2")[1 : 3]
  names(col.expression) = c("Original", "ModelFitted", "ModeledGroupDiff")
  col.pseudotime = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann$pseudotime)))
  names(col.pseudotime) = unique(colann$pseudotime)
  if (is.null(rowann)) {
    rowann = data.frame(cluster = factor(as.character(clu), levels = as.character(1 : max(clu))), XDEType = factor(as.character(XDEType[names(clu)]), levels = as.character(unique(res[, 4]))), stringsAsFactors = F)
    rownames(rowann) = names(clu)
  }
  rowann <- rowann[rownames(fit.scale), , drop = F]
  rowann[, "XDEType"] <- factor(as.character(rowann[, "XDEType"]), levels = c("trendSig", "meanSig", "bothSig", "nonXDE", "other"))
  if (length(unique(clu)) < 8) {
    col.clu = brewer.pal(8, "Set1")[1 : length(unique(clu))]
  } else {
    col.clu = colorRampPalette(brewer.pal(8, "Set1"))(length(unique(clu)))
  }
  set.seed(12345)
  col.clu <- sample(col.clu)
  names(col.clu) = levels(rowann$clu)
  if (is.null(colann) | is.null(annotation_colors)) {
    col.meanDiff = c("blue", "red")
    names(col.meanDiff) <- c("Positive", "Negative")
    col.XDEType = brewer.pal(8, "Set3")[1 : 3]
    names(col.XDEType) = c("trendSig", "bothSig", "meanSig")
    annotation_colors = list(pseudotime = col.pseudotime, group = col.group, expression = col.expression, cluster = col.clu, XDEType = col.XDEType)
  }
  col.gs <- c("pink", "skyblue")
  names(col.gs) <- c("No", "Yes")
  col.limmaPb <- c("pink", "skyblue")
  names(col.limmaPb) <- c("nonDiff", "Diff")
  annotation_colors[["gs"]] <- col.gs
  annotation_colors[["limmaPb"]] <- col.limmaPb
  col.signalType <- brewer.pal(8, "Set3")[1 : 3]
  names(col.signalType) <- c("trend only", "mean only", "both")
  annotation_colors[["signalType"]] <- col.signalType
  cpl = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  if (break.0) {
    cpl <- c(cpl[1 : 40], cpl[60 : 100])
  }
  plist <- list()
  if (!is.na(sep)) {
    rownames(expr.scale) <- sub(sep, "", rownames(expr.scale))
    rownames(rowann) <- sub(sep, ":.*", rownames(rowann))
    rownames(oridata) <- sub(sep, "", rownames(oridata))
  }
  p1data <- expr.scale
  p1data[p1data > quantile(as.vector(p1data), 0.95, na.rm = T)] <- quantile(as.vector(p1data), 0.95, na.rm = T)
  p1data[p1data < quantile(as.vector(p1data), 0.05, na.rm = T)] <- quantile(as.vector(p1data), 0.05, na.rm = T)
  col_fun = circlize::colorRamp2(seq(-max(abs(p1data)), max(abs(p1data)), length.out = 50), colorRampPalette(c("blue3", "skyblue", "white", "pink", "red3"))(50))
  pt_col_fun = circlize::colorRamp2(seq(1, ceiling(max(colann$pseudotime))), colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(max(colann$pseudotime)))
  annotation_colors$pseudotime <- pt_col_fun
  ht1 <- ComplexHeatmap::Heatmap(p1data, cluster_rows = F, cluster_columns = FALSE, show_row_names = showRowName, show_column_names = FALSE, col = col_fun, heatmap_legend_param = list(legend_direction = "horizontal"), top_annotation = HeatmapAnnotation(df = colann, col = annotation_colors, show_annotation_name = F), left_annotation = rowAnnotation(df = rowann, col = annotation_colors), width = 2)
  colann.fit1 <- data.frame(pseudotime = rep(1 : ncol(fit[[1]]), length(fit)), group = gsub(sub("_.*", "_", names(fit)[1]), "", sub(";.*", "", colnames(fit.scale))), expression = "ModelFitted", stringsAsFactors = F)
  colann.fit2 <- data.frame(pseudotime = seq(1, ncol(FitDiff.scale)), group = "NA", expression = "ModeledGroupDiff", stringsAsFactors = F)
  colann.fit <- rbind(colann.fit1, colann.fit2)
  col.pseudotime = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann.fit$pseudotime)))
  names(col.pseudotime) = unique(colann.fit$pseudotime)
  annotation_colors$pseudotime <- col.pseudotime
  fit.scale <- cbind(fit.scale, FitDiff.scale)
  rownames(colann.fit) = colnames(fit.scale)
  if (!is.na(sep)) {
    rownames(fit.scale) <- sub(sep, "", rownames(fit.scale))
  }
  p2data <- fit.scale[, rownames(colann.fit)[colann.fit[, "expression"] != "ModeledGroupDiff"]]
  p2data[p2data > quantile(as.vector(p2data), 0.99, na.rm = T)] <- quantile(as.vector(p2data), 0.99, na.rm = T)
  p2data[p2data < quantile(as.vector(p2data), 0.01, na.rm = T)] <- quantile(as.vector(p2data), 0.01, na.rm = T)
  col_fun = circlize::colorRamp2(seq(-max(abs(p2data)), max(abs(p2data)), length.out = 50), colorRampPalette(c("blue3", "skyblue", "white", "pink", "red3"))(50))
  pt_col_fun = circlize::colorRamp2(seq(1, length(unique(colann.fit$pseudotime))), colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann.fit$pseudotime))))
  annotation_colors$pseudotime <- pt_col_fun
  ht2 <- ComplexHeatmap::Heatmap(p2data, cluster_rows = F, cluster_columns = FALSE, show_row_names = showRowName, show_column_names = FALSE, heatmap_legend_param = list(legend_direction = "horizontal"), col = col_fun, top_annotation = HeatmapAnnotation(df = colann.fit[colnames(p2data), ], col = annotation_colors, show_annotation_name = F), width = 2)
  p3data <- fit.scale[, rownames(colann.fit)[colann.fit[, "expression"] == "ModeledGroupDiff"]]
  p3data[p3data > quantile(as.vector(p3data), 0.99, na.rm = T)] <- quantile(as.vector(p3data), 0.99, na.rm = T)
  p3data[p3data < quantile(as.vector(p3data), 0.01, na.rm = T)] <- quantile(as.vector(p3data), 0.01, na.rm = T)
  col_fun = circlize::colorRamp2(seq(-max(abs(p3data)), max(abs(p3data)), length.out = 50), colorRampPalette(c("blue3", "skyblue", "white", "pink", "red3"))(50))
  ht3 <- ComplexHeatmap::Heatmap(p3data, cluster_rows = F, cluster_columns = FALSE, show_row_names = showRowName, show_column_names = FALSE, heatmap_legend_param = list(legend_direction = "horizontal"), col = col_fun, top_annotation = HeatmapAnnotation(df = colann.fit[colnames(p3data), ], col = annotation_colors, show_annotation_name = F), width = 1)
  p4data <- fit.scale[, rownames(colann.fit)[colann.fit[, "expression"] == "ModeledGroupDiff"]]
  p4data <- (p4data - rowMeans(p4data)) / apply(p4data, 1, sd)
  p4data[p4data > quantile(as.vector(p4data), 0.99, na.rm = T)] <- quantile(as.vector(p4data), 0.99, na.rm = T)
  p4data[p4data < quantile(as.vector(p4data), 0.01, na.rm = T)] <- quantile(as.vector(p4data), 0.01, na.rm = T)
  p4data[rownames(p4data) %in% sub(":.*", "", rownames(rowann)[rowann$XDEType == "meanSig"]), ] <- 0
  col_fun = circlize::colorRamp2(seq(-max(abs(p4data)), max(abs(p4data)), length.out = 50), colorRampPalette(c("blue3", "skyblue", "white", "pink", "red3"))(50))
  ht4 <- ComplexHeatmap::Heatmap(p4data, cluster_rows = F, cluster_columns = FALSE, show_row_names = showRowName, show_column_names = FALSE, heatmap_legend_param = list(legend_direction = "horizontal"), col = col_fun, top_annotation = HeatmapAnnotation(df = colann.fit[colnames(p4data), ], col = annotation_colors, show_annotation_name = F), width = 1)
  rownames(alluniformdiff) <- sub(":.*", "", rownames(alluniformdiff))
  p5data <- rowMeans(alluniformdiff[rownames(fit.scale), , drop = FALSE]) %*% matrix(1, nrow = 1, ncol = ncol(p4data))
  rownames(p5data) <- rownames(p4data)
  p5data[rownames(p5data) %in% sub(":.*", "", rownames(rowann)[rowann$XDEType == "trendSig"]), ] <- 0
  col_fun = circlize::colorRamp2(seq(-max(abs(p5data)), max(abs(p5data)), length.out = 50), colorRampPalette(c("blue3", "skyblue", "white", "pink", "red3"))(50))
  ht5 <- ComplexHeatmap::Heatmap(p5data, cluster_rows = F, cluster_columns = FALSE, show_row_names = showRowName, show_column_names = FALSE, heatmap_legend_param = list(legend_direction = "horizontal"), col = col_fun, width = 1)
  htlist <- ht1 + ht2 + ht4 + ht5
  draw(htlist, merge_legend = FALSE, annotation_legend_side = "right", heatmap_legend_side = "bottom")
}
