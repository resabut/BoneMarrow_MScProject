lamian_downstream <- function(Res) {
    stat <- Res$statistics
    stat <- stat[order(stat[, 1], -stat[, 3]),]
    # genes with FDR.overall < 0.05 cutoff
    diffgene <- rownames(stat[stat[, grep ('^fdr.*overall$', colnames(stat))]< 0.05, ])
    # population level stats
    Res$populationFit <- getPopulationFit(Res, gene = diffgene, type = 'variable')
    # group difference for a given covariate
    Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = diffgene)
    # cluster genes
    Res$cluster <- clusterGene(Res, gene = diffgene, type = 'variable', k = 5)
    return(Res)
}