# Usage: used to search for each transcription start site (TSS) and identify every possible link between a peak and its target gene

library(SnapATAC)
library(GenomicRanges)

# Read data and set pvalue cutoff to be 1e-3
atac_snap <- readRDS('../output/EB1/atac_snap_EB1_pmat.rds')
df_TSS <- read.table('../data/database/hg38_chr_TSS.txt', header = T, sep = '\t', quote = '', stringsAsFactors = F)
pvalue_cutoff <- 3

chr1 <- c()
start1 <- c()
end1 <- c()
chr2 <- c()
start2 <- c()
end2 <- c()
gene_name <- c()
pvalue <- c()

# Iterate for each TSS and identify peaks related to this TSS
for (i in 1:nrow(df_TSS)) {
  if (i %% 1000 == 0) {
    cat(df_TSS[i, 1], ': ', i, '/', nrow(df_TSS), '\n', sep = '')
  }
  if (df_TSS[i, 1] %in% colnames(atac_snap@gmat)) {
    TSS.loci = GRanges(df_TSS[i, 2], IRanges(df_TSS[i, 3], df_TSS[i, 3] + 1))
    
    gene.loci = resize(TSS.loci, width=500000, fix="center")
    gene.val = atac_snap@gmat[,df_TSS[i, 1]]
    data.use = atac_snap@pmat
    peak.use = atac_snap@peak
    ncell = nrow(atac_snap@pmat)
    idy = queryHits(findOverlaps(peak.use, gene.loci)) # Find peaks within a certain distance to the TSS
    if ((x=length(idy)) != 0L & length(idy) > 1){
      data.use = data.use[,idy]
      peak.use = peak.use[idy]
      idy = which(Matrix::colSums(data.use) > 0) 
      if ((x=length(idy)) != 0L & length(idy) > 1) {
        data.use = data.use[,idy]
        peak.use = peak.use[idy]
        peaks.id = seq(ncol(data.use))
        models = lapply(peaks.id, function(t){
          if ('x' %in% row.names(summary(stats::glm(y ~ x, data = data.frame(y=data.use[,t], x=gene.val), family = binomial(link='logit')))[["coefficients"]])) {
            summary(stats::glm(y ~ x, data = data.frame(y=data.use[,t], x=gene.val), family = binomial(link='logit')))[["coefficients"]]["x",]
          } else {
            as.data.frame(t(data.frame(c(0, 0, 0, 1))), row.names = 'x')
          }
        }) # Use the snapATAC test for associating genes and peaks
        models <- do.call(rbind, models)
        colnames(models) <- c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')
        models[models[,"z value"] < 0,"Pr(>|z|)"] = 1
        peak.use$beta = models[,"Estimate"]
        peak.use$zvalue = models[,"z value"]
        peak.use$stde = models[,"Std. Error"]
        peak.use$Pval = models[,"Pr(>|z|)"]
        peak.use$logPval = -log10(peak.use$Pval)
        
        pairs.df = as.data.frame(peak.use)
        rm(peak.use, data.use)
        pairs.df = data.frame(
          chr1=pairs.df[,"seqnames"],
          start1=pairs.df[,"start"],
          end1=pairs.df[,"end"],
          pval=pairs.df[,"logPval"]
        )
        pairs.df.filtered <- pairs.df[which(pairs.df$pval > pvalue_cutoff), ] # Select peaks that satisfy the cutoff for p-values
        if (nrow(pairs.df.filtered) > 0) {
          chr1 <- append(chr1, pairs.df.filtered$chr1)
          start1 <- append(start1, pairs.df.filtered$start1)
          end1 <- append(end1, pairs.df.filtered$end1)
          chr2 <- append(chr2, rep(df_TSS[i, 2], nrow(pairs.df.filtered)))
          start2 <- append(start2, rep(df_TSS[i, 3], nrow(pairs.df.filtered)))
          end2 <- append(end2, rep(df_TSS[i, 3] + 1, nrow(pairs.df.filtered)))
          gene_name <- append(gene_name, rep(df_TSS[i, 1], nrow(pairs.df.filtered)))
          pvalue <- append(pvalue, pairs.df.filtered$pval)
        }
      }
    }
  }
}

df_total <- data.frame(chr1, start1, end1, chr2, start2, end2, gene_name, pvalue)
write.table(df_total, '../output/EB1/EB1_gene_peak.txt', sep = '\t', quote = F, row.names = F, col.names = F)

