VisualizeUmap2 <- function(genes, type, SeuratObj, ncol = NULL, pt.size = 0.05, min.cutoff = -5, max.cutoff = 5, order = TRUE, slot = "data", alpha = 0.5){
  if(type == "ATAC"){
    genes = intersect(rownames(SeuratObj$ACTIVITY),unique(genes))
    gene_expr = GetAssayData(object = SeuratObj$ACTIVITY)[genes, ]
  }
  if(type == "RNA"){
    genes = intersect(rownames(SeuratObj),unique(genes))
    gene_expr = GetAssayData(object = SeuratObj, slot = "scale.data")[genes, ]
  }
  umap_df = SeuratObj@reductions$umap@cell.embeddings
  umapplots = lapply(genes, function(x){
    umap_expr = merge(umap_df, as.data.frame(gene_expr[x, ]), by.x = 0, by.y = 0)
    row.names(umap_expr) = umap_expr[,1]
    umap_expr = umap_expr[,-1]
    colnames(umap_expr)[3] = "Gene"
    # sort
    if (order == TRUE) {
        umap_expr = umap_expr[order(umap_expr[,3]),]
    }
    # set min max cutoff
    umap_expr[umap_expr[,3]>max.cutoff,] <- max.cutoff
    umap_expr[umap_expr[,3]<min.cutoff,] <- min.cutoff


    if(type == "RNA"){
      legendtitle = "Expression level Z-score"}
    if(type == "ATAC"){
      legendtitle = "RP score"}
    p = ggplot(umap_expr, aes(x = UMAP_1, y = UMAP_2, color = Gene)) +
      geom_point(aes(colour = Gene), size = pt.size, alpha = alpha, na.rm = TRUE) +
      scale_color_gradient2(name = legendtitle, low = muted("blue"), mid = "white", high = muted("red"), limits=c(min.cutoff,max.cutoff)) +
      labs(title = x) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"))
    return(p)
  })
  combinedplot = CombinePlots(umapplots, ncol = ncol)
  return(combinedplot)
}

# example of usage
#p <- VisualizeUmap2(tfs, "RNA", rna.aftercc$RNA, ncol = 2, order = TRUE, slot = slot, min.cutoff = min.cutoff, max.cutoff = max.cutoff, pt.size = pt.size, alpha = alpha)
#ggsave(filename="selected.tfs.png",plot=p, height=6, width=10)
