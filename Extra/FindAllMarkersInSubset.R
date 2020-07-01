MAESTROFindAllMarkersInSubset <- function(seuratObject, subset.clusters, min.pct = 0.1, logfc.threshold = 0.25, p.val.adj.cutoff = 1e-5, test.use = "presto", only.pos = FALSE ) {
    subset.object <- SubsetData(object=seuratObject, subset.name = "orig.ident", accept.value = subset.orig.ident, ident.use = subset.clusters)
    res <- data.frame()
    for ( c1 in subset.clusters ) {
        c.1 <- c(c1)
        c.2 <- setdiff(subset.clusters, c.1)
        cluster.genes <- FindMarkersMAESTRO(object = subset.object, min.pct = min.pct, logfc.threshold = logfc.threshold, test.use = test.use, only.pos = only.pos, ident.1=c.1, ident.2=c.2)
	cluster.genes <- cluster.genes[cluster.genes$p_val_adj<p.val.adj.cutoff, ]
	cluster.genes$cluster <- c1
	cluster.genes$gene <- row.names(cluster.genes)
	res.cluster.genes <- as_tibble(cluster.genes)
        res <- rbind(res, res.cluster.genes)
    }
    return(list(RNA=subset.object, genes=res))
}

# example of usage

#res.0w <- MAESTROFindAllMarkersInSubset(rna.aftercc$RNA, min.pct = genes.pct, logfc.threshold = genes.logfc, p.val.adj.cutoff = genes.cutoff, test.use = genes.test.use, only.pos = only.pos, subset.clusters = cluster.0w)
#top10.0w <- res.0w$genes %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


