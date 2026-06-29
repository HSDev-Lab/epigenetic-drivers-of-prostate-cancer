library(tidyverse)
library(annotatr)

DMRS_grouped <- read_csv("/Users/tobymilne-clark/Documents/Data/Atlas/Prostate_Atlas/DMRs_grouped_mean_with_coords_Prostate_Atlas.csv")

annotate_dmrs <- function(DMRS, annotation = "hg38"){
  colnames(DMRS)[1:3] <- c("chr","start","end")
  if(annotation == "hg19"){
    annots = c('hg19_basicgenes',"hg19_lncrna_gencode","hg19_enhancers_fantom")
    annotations = annotatr::build_annotations(genome = 'hg19', annotations = annots)
  } else if(annotation == "hg38"){
    annots = c('hg38_basicgenes',"hg38_lncrna_gencode","hg38_enhancers_fantom")
    annotations = annotatr::build_annotations(genome = 'hg38', annotations = annots)
    
    print("Downloading EnsDb.Hsapiens.v101")
    
    ah <- AnnotationHub::AnnotationHub()
    AnnotationHub::query(ah, "EnsDb.Hsapiens.v101")
    edb <- ah[["AH83216"]]
    txs <- ensembldb::transcripts(edb, columns = c("tx_id", "tx_biotype", "tx_id_version", "gc_content", "gene_name", "gene_id"))
    ENST <- as.data.frame(txs)
  }
  DMR.granges <- DMRS %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE,starts.in.df.are.0based = TRUE)
  dm_annotated = annotatr::annotate_regions(
    regions = DMR.granges,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
  DMRs_annotated <- as.data.frame(dm_annotated) %>% dplyr::distinct(annot.symbol, .keep_all = TRUE)
  if(annotation == "hg38"){
    DMRs_annotated <- DMRs_annotated %>% dplyr::left_join(ENST, by = c("annot.tx_id" = "tx_id_version")) %>%
      dplyr::distinct(gene_id, .keep_all = TRUE)
  }
  DMRs_annotated
}

extract_proteins <- function(annotated_DMRs, save_dir = getwd(), seperate_networks = "TRUE", sep_col = "celltypes_expanded"){
  
  ah <- AnnotationHub::AnnotationHub()
  AnnotationHub::query(ah, "EnsDb.Hsapiens.v101")
  edb <- ah[["AH83216"]]
  ensembldb::mapIds(edb, keys=annotated_DMRs$gene_name, column="PROTEINID", keytype="SYMBOL", multiVals="CharacterList")
  
  if(seperate_networks == "TRUE"){
    for(i in unique(annotated_DMRs[[sep_col]])){
      gene <- annotated_DMRs %>% filter(sep_col == i) %>% dplyr::select(gene_name)
      protein <- unname(unlist(mapIds(edb, keys=gene$gene_name, column="PROTEINID", keytype="SYMBOL", multiVals="list"))) 
      print(paste0("Proteins_",i,"_cells.csv"))
      gene %>% write_csv(paste0(save_dir,"/Genes_",i,"_cells.csv"),col_names = FALSE)
      protein %>% na.omit() %>% as.data.frame() %>% write_csv(paste0(save_dir,"Proteins_",i,"_cells.csv"),col_names = FALSE)
    }
  } else {
    gene <- annotated_DMRs %>% dplyr::select(gene_name)
    protein <- unname(unlist(mapIds(edb, keys=gene$gene_name, column="PROTEINID", keytype="SYMBOL", multiVals="list"))) 
    gene %>% write_csv(paste0(save_dir,"/Genes.csv"),col_names = FALSE)
    protein %>% na.omit() %>% write_csv(paste0(save_dir,"/Proteins.csv"),col_names = FALSE)
  }
}
  
  
  
  
  
  










