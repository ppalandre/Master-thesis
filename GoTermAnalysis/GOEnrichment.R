setwd("../GoTermAnalysis")

if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

BiocManager::install(c("GO.db",  "clusterProfiler", "AnnotationForge", "enrichplot"), force = TRUE)
install.packages('tidyr')


# Install all your R dependencies

library(AnnotationForge)
library(tidyr)


createAnnotation <- function(goTable, symbolTable, outDir, genus="Synechocystis", species="sp.strainPCC6803"){
    mapping_data <- read.table(goTable, sep="\t", header=TRUE, quote="\"")

    if(!("EVIDENCE" %in% colnames(mapping_data))){
        mapping_data$EVIDENCE <- "experimental"
    }
    colnames(mapping_data) <- c("GID", "GO", "EVIDENCE")
    symbol_data <- read.table(symbolTable, sep="\t", header=TRUE, quote="\"")
    colnames(symbol_data) <- c("GID", "SYMBOL")
    
    dir.create(outDir, recursive=TRUE)
    rv <- makeOrgPackage(
        gene_info=symbol_data,
        go=mapping_data,
        version="0.1",
        maintainer="DR <schdruzzi@gmail.com>",
        author="DR <schdruzzi@gmail.com>",
        #outputDir = snakemake@output[["annotation_db"]], ###
        outputDir = "annotation_db",
        tax_id = "42",
        genus = genus,
        species=species,
        goTable = "go"
    )
    return(rv)
}

package <- createAnnotation(
    goTable="./INPUT2_GO-mapping/output_GOALL_terms_PHILLIPP_NoDupes.txt",      # pfad zur tabelle mit den mappings locus_tag zu Go Term
    symbolTable="./INPUT2_GO-mapping/output_symbols_PHILLIPP.txt",   # pfad zur tabelle mapping ProteinID zu locus_tag
    outDir="annotation_db",        # hier wird die database hininstalliert
    #outDir=getwd(),        # hier wird die database hininstalliert
    genus="Synechocystis",              # völlig egal (string)
    species="sp.strainPCC6803"       # völlig egal (string)
)

install.packages(package, repos=NULL, type="source")
library(basename(package), character.only = TRUE)

library(clusterProfiler)

install.packages("viridis")
library(viridis)

######################################################################################################################

log2cutoff <- 2.5
padjcutoff <- 0.05

minGSSize <- 30
maxGSSize <- 500


# Hier noch schauen dass die Tabelle richtig geladen wird

defile <- "./INPUT1_DESeq2_Results/msms/msms_cold.txt" # Pfad zu deiner DESeq Tabelle
detable <- read.table(defile,sep="\t",header=TRUE, quote="\"")
detable <- na.omit(detable)


up <- detable$log2FoldChange >= log2cutoff & detable$padj < padjcutoff
up <- detable[up, ]
down <- detable$log2FoldChange <= -log2cutoff & detable$padj < padjcutoff
down <- detable[down, ]

## eigentlich Funktion für die GO analysis, up-regulated genes
egoBPup <- enrichGO(gene = as.character(up$"UniProt"),
                    universe = as.character(detable$"UniProt"),
                    OrgDb = basename(package),
                    keyType = "GID",
                    ont = "ALL",
                    pAdjustMethod = "BH", # Benjamini Hochberg, padjuted wegen multiple testing problem
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    minGSSize = minGSSize,
                    maxGSSize = maxGSSize,
                    readable = FALSE)
summary <- data.frame(egoBPup )

# Hier wird einfach nur eine leere Tabelle gebaut falls keine GO terms enriched sind
if (dim(summary)[1] == 0){
  df <- data.frame(matrix(ncol = 10, nrow = 0))
  x <- c("ONTOLOGY", "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
  colnames(df) <- x
  summary <- df
}
write.table(summary , file = "./OUTPUT1_GO-enrichement-RAW/msms/GOALLoutput_msms-cold_lfc0&padj0-05", row.names=FALSE, sep="\t") # Pfad  zum output file für hochregulierte gene

## eigentlich Funktion für die GO analysis, down-regulated genes
egoBPdown <- enrichGO(gene = as.character(down$"Gene.Name"),
                    universe = as.character(detable$"Gene.Name"),
                    OrgDb = basename(package),
                    keyType = "GID",
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    minGSSize = minGSSize,
                    maxGSSize = maxGSSize,
                    readable = TRUE)

summary <- data.frame(egoBPdown )

# Hier wird einfach nur eine leere Tabelle gebaut falls keine GO terms enriched sind
if (dim(summary)[1] == 0){
  df <- data.frame(matrix(ncol = 10, nrow = 0))
  x <- c("ONTOLOGY", "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
  colnames(df) <- x
  summary <- df
}
write.table(summary , file = "./OUTPUT1_GO-enrichement-RAW/GOALLoutput_COLD_DOWN_CDS", row.names=FALSE, sep="\t") # Pfad  zum output file für runterregulierte gene


######################################################################################################################


svg("./OUTPUT2_GO-enrichement-PLOTS/msms/GOALLbarplot_msms-cold_lfc0&padj0-05.svg")
#jpeg("./OUTPUT2_GO-enrichement-PLOTS/msms/GOALLbarplot_msms-cold_lfc0&padj0-05.jpeg", quality = 75)
barplot(summary, showCategory=50,font.size=8, title="") + viridis::scale_fill_viridis()


dev.off() 

