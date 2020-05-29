library(scater)

## Palette taken from `scater`
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

## Load the SCE
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda",
     verbose=T)
# sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo
rm(chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all)
# First drop "ambig.lowNtrxts" (93 nuclei)
sce.nac.all <- sce.nac.all[ ,sce.nac.all$cellType.final != "ambig.lowNtrxts"]
sce.nac.all.orig <- sce.nac.all
sce.nac.all <- sce.nac.all[ ,sce.nac.all$cellType.final == c("MSN.D1.1","MSN.D1.2","MSN.D1.3","MSN.D1.4","MSN.D2.1","MSN.D2.2")]
#sce.nac.all <- sce.nac.all[ ,sce.nac.all$cellType.final == c("Inhib.1","Inhib.2","Inhib.3","Inhib.4")]

sce.nac.all$cellType.final <- droplevels(sce.nac.all$cellType.final)
# Create list of genes to plot
genes2plot = list(
  'experiment1' = c('DRD1', 'TAC1','RXFP2','GABRQ'),
  'experiment2' = c('DRD1', 'TAC1','CRHR2','RXFP1'),
  'experiment3' = c('DRD1', 'DRD2','TAC1','PENK'),
  'experiment4' = c('DRD1', 'DRD2','CRHR2','HTR7'),
  'experiment5' = c('DRD1', 'DRD2','RELN','TAC1','RXFP2'),
  'experiment6' = c('PVALB', 'GAD1','PTHLH','KIT'),
  'msn_markers' = c('PPP1R1B','BCL11B','PDYN','PENK','TAC1'),
  'inhib_markers' = c('GAD1','GAD2','PPP1R1B','PDYN','PENK','TAC1')
  # etc.
)

pdf("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/bam_pseudobulk_abby/RNAscopeAll.pdf", height=13, width=3)
for(i in 1:length(genes2plot)){
  print(
    plotExpression(sce.nac.all, exprs_values = "logcounts", features=genes2plot[[i]],
    x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7, ncol=1,
    add_legend=F) 
     + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                 geom = "crossbar", width = 0.3,
                                 colour=rep(tableau20[3:8], length(genes2plot[[i]]))) 
    +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
    ggtitle(label=paste0(names(genes2plot)[i], " probes"))
  )
}

dev.off()

