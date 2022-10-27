# Function to carry out differential expression analyisis with DESeq
# Source: http://joey711.github.io/phyloseq-extensions/DESeq2.html

require(phyloseq)
require(DESeq2)
require(ggplot2)
phyloseq_dea <- function(pseq,
                         test = "Wald",
                         fit_type = "parametric",
                         alpha = 0.05) {

    dds <- phyloseq_to_deseq2(pseq, ~Event)
    dds <- DESeq(dds, test = test, fitType = fit_type)

    res <- results(dds, cooksCutoff = FALSE)
    sigtab <- res[which(res$padj < alpha), ]
    sigtab <- cbind(as(sigtab, "data.frame"),
                    as(tax_table(pseq)[rownames(sigtab), ], "matrix"))
    return(sigtab)
}

plot_dea <- function(sigtab) {
    # Phylum order
    x <- tapply(sigtab$log2FoldChange,
                sigtab$Phylum,
                function(x) max(x))
    x <- sort(x, TRUE)
    sigtab$Phylum <- factor(as.character(sigtab$Phylum), levels = names(x))

    # Genus order
    x <- tapply(sigtab$log2FoldChange,
                sigtab$Genus,
                function(x) max(x))
    x <- sort(x, TRUE)
    sigtab$Genus <- factor(as.character(sigtab$Genus), levels = names(x))

    ggplot(sigtab,
          aes(x = Genus,
              y = log2FoldChange,
              color = Phylum)) +
            geom_point(size = 6) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = -90,
                                             hjust = 0,
                                             vjust = 0.5))

}