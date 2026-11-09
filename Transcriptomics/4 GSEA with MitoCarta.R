pacman::p_load(
  tidyverse, janitor, patchwork, ggrepel, openxlsx, knitr, readxl,
  DESeq2, clusterProfiler, enrichplot, ComplexHeatmap
)

this.path::here() %>% setwd()

rm(list = ls())



# Overlap with our DE for mVAT --------------------------------------------
mito_genes <- read_excel("../Transcriptomics/data/Mouse.MitoCarta3.0.xls", sheet = "A Mouse MitoCarta3.0") %>%
  pull(Symbol) %>%
  unique()

length(mito_genes) # 1,140







# Prepare MitoCarta metabolism gene set ----------------------------------
pathways <- read_excel("data/Mouse.MitoCarta3.0.xls", sheet = "C MitoPathways") %>% 
  clean_names()

mito_genes <- pathways |> 
  separate_longer_delim(genes, delim = ", ") |> 
  pull(genes) |> 
  unique()




# Terminal leaf pathways: any hierarchy string that isn't a prefix of another
hierarchy <- pathways %>% 
  select(mito_pathway, mito_pathway_hierarchy) %>% 
  deframe()

is_leaf <- function(h, all_h) {
  prefix <- paste0(h, " > ")
  !any(startsWith(all_h[all_h != h], prefix))
}

pathways$is_leaf <- vapply(hierarchy, is_leaf, logical(1), all_h = hierarchy)

term2gene_leaf <- pathways %>%
  filter(is_leaf) %>%
  separate_rows(genes, sep = ",\\s*") %>%
  transmute(
    term = mito_pathway,
    gene = genes
  ) %>%
  distinct()

n_distinct(term2gene_leaf$gene) # 977

term2gene_leaf |> filter(term == "CII subunits")  # Sdha, Sdhb, Sdhc, Sdhd

ngenes <- sapply( split(term2gene_leaf$gene, term2gene_leaf$term), length )
table(ngenes)

MASS::truehist(ngenes)
length(ngenes)    # 122 testable terms
sum(ngenes >=  4) # 103 testable terms including CII subunits
sum(ngenes >=  5) #  93 testable terms
sum(ngenes >=  9) #  58 testable terms
sum(ngenes >= 10) #  50 testable terms

rm(pathways, is_leaf, ngenes)


# Load differential-expression results ----------------------------------
load("results/DEG.rda")



# Run GSEA with leaf nodes of mitocarta 3.0 -----------------------------
gsea_list <- list()

for(depot in names(big.res)){
  
  set.seed(123)
  print(depot)
  
  ## Diet effect
  ranked <- big.res[[ depot ]] %>% 
    rename_with( ~ str_replace(., paste0(depot, "_"), "")) |> 
    filter(SYMBOL %in% mito_genes) %>%                         # Restrict to just the mitochondrial genes
    mutate(score = -1*log10(Diet.p)*sign(Diet.LFC) ) %>% 
    dplyr::filter(!is.na(score)) %>% 
    arrange(desc(score)) %>% 
    dplyr::select(SYMBOL, score) %>% 
    drop_na() %>% 
    deframe()
  
  gsea_list[[ paste0(depot, "_diet") ]] <- GSEA(ranked, TERM2GENE = term2gene_leaf, minGSSize = 4, eps = 0, 
                 pvalueCutoff = 1, sampleSize = 1001, nPermSimple = 1e7, pAdjustMethod = "fdr")
  
  rm(ranked)
  
  ## CL effect
  ranked <- big.res[[ depot ]] |> 
    rename_with( ~ str_replace(., paste0(depot, "_"), "")) |> 
    filter(SYMBOL %in% mito_genes) %>%                         # Restrict to just the mitochondrial genes
    mutate(score = -1*log10(CL.p)*sign(CL.LFC) ) |> 
    dplyr::filter(!is.na(score)) |> 
    arrange(desc(score)) |> 
    dplyr::select(SYMBOL, score) |> 
    drop_na() |> 
    deframe()
  
  gsea_list[[ paste0(depot, "_CL") ]] <- GSEA(ranked, TERM2GENE = term2gene_leaf, minGSSize = 4, eps = 0, 
                                              pvalueCutoff = 1, sampleSize = 1001, nPermSimple = 1e7, pAdjustMethod = "fdr")
  
  rm(ranked, depot)
}

save(gsea_list, file = "results/GSEA_MitoCarta.rda")

openxlsx::write.xlsx(gsea_list, file = "results/GSEA_MitoCarta_tmp.xlsx")


# GSEA in wide format with key columns ------------------------------------
gsea_list_light <- list()

for(combo in names(gsea_list)){
  
  obj <- gsea_list[[combo]]
  
  obj <- obj@result |> 
    mutate(ncore = lengths(strsplit(core_enrichment, "/"))) |> 
    select(ID, Description, setSize, NES, FDR = p.adjust) |> 
    remove_rownames()
  
  colnames(obj)[-c(1,2)] <- paste0(combo, ".", colnames(obj)[-c(1,2)])
  
  gsea_list_light[[ combo ]] <- obj
  rm(obj, combo)
}

gsea_list_light <- Reduce("full_join", gsea_list_light)

openxlsx::write.xlsx(gsea_list_light, file = "results/GSEA_MitoCarta_wide_tmp.xlsx")
