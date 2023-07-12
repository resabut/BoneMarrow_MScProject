plot_lineage <- function(lineage_to_plot) {
    cell_colors_clust <- cell_pal(lineage_to_plot$ct, hue_pal())
    plot(reducedDims(lineage_to_plot)$UMAP,
         col = cell_colors_clust,
         pch=16,
         asp = 1,
         xlab = 'UMAP 1',
         ylab = 'UMAP 2')
    lines(SlingshotDataSet(lineage_to_plot), lwd=2, type = 'lineages', col = 'black')
    lines(SlingshotDataSet(lineage_to_plot)@curves$Lineage1$s, lwd=2, col = 'red')
    legend('bottomright', legend = levels(factor(lineage_to_plot$ct)), col = hue_pal()(length(levels(factor
                                                                                                     (lineage_to_plot$ct)))), pch = 16)
}