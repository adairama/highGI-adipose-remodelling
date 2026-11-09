# SETUP --------------------------------------------------------------------
pacman::p_load(tidyverse, janitor, patchwork, ggrepel, DESeq2, openxlsx, knitr, readxl,
               clusterProfiler)
this.path::here() |> setwd()
rm(list=ls())


# load the DEG results ----------------------------------------------------
load("results/DEG.rda")
gtf_map <- gtf_map |> select(ENTREZID, SYMBOL) |> remove_rownames()


# GSEA with KEGG ----------------------------------------------------------
gsea_list <- list()

for(depot in names(big.res)){
  
  print(depot)
  
  ## Diet effect
  ranked <- big.res[[ depot ]] |> 
    mutate(score = -1*log10(Diet.p)*sign(Diet.LFC) ) |> 
    dplyr::filter(!is.na(score)) |> 
    arrange(desc(score)) |> 
    left_join(gtf_map) |> 
    dplyr::select(ENTREZID, score) |> 
    drop_na() |> 
    deframe()

  gsea_list[[ paste0(depot, "_diet") ]] <- gseKEGG(ranked, organism = "rno", pvalueCutoff = 1.01) |> 
    setReadable_gtf(id2symbol = deframe(gtf_map))
  
  rm(ranked)
  
  ## CL effect
  ranked <- big.res[[ depot ]] |> 
    mutate(score = -1*log10(CL.p)*sign(CL.LFC) ) |> 
    dplyr::filter(!is.na(score)) |> 
    arrange(desc(score)) |> 
    left_join(gtf_map) |> 
    dplyr::select(ENTREZID, score) |> 
    drop_na() |> 
    deframe()
  
  gsea_list[[ paste0(depot, "_CL") ]] <- gseKEGG(ranked, organism = "rno", pvalueCutoff = 1.01) |> 
    setReadable_gtf(id2symbol = deframe(gtf_map))
  
  rm(ranked, depot)
}

save(gsea_list, file = "results/GSEA_KEGG.rda")

openxlsx::write.xlsx(gsea_list, file = "results/GSEA_KEGG.xlsx")


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

openxlsx::write.xlsx(gsea_list_light, file = "results/GSEA_KEGG_wide_tmp.xlsx")



# UpSet plot for Diet effect ----------------------------------------------

require(ComplexHeatmap)

diet <- list(
  mVAT_up   = wide |> filter(mVAT_diet_FDR <= 0.05, mVAT_diet_NES > 0) |> pull(Description),
  SAT_up    = wide |> filter(SAT_diet_FDR  <= 0.05, SAT_diet_NES  > 0) |> pull(Description),
  iBAT_up   = wide |> filter(iBAT_diet_FDR <= 0.05, iBAT_diet_NES > 0) |> pull(Description),
  
  mVAT_down = wide |> filter(mVAT_diet_FDR <= 0.05, mVAT_diet_NES < 0) |> pull(Description),
  SAT_down  = wide |> filter(SAT_diet_FDR  <= 0.05, SAT_diet_NES  < 0) |> pull(Description),
  iBAT_down = wide |> filter(iBAT_diet_FDR <= 0.05, iBAT_diet_NES < 0) |> pull(Description)
)

m <- make_comb_mat(diet)

UpSet(m, 
      top_annotation   = upset_top_annotation(m, add_numbers = TRUE),
      right_annotation = upset_right_annotation(m, add_numbers = TRUE),
      set_order = c("mVAT_up", "SAT_up", "iBAT_up", "mVAT_down", "SAT_down", "iBAT_down")
      )

ggplot(wide, aes(x = mVAT_diet_NES, y = SAT_diet_NES)) +
  geom_point()


# UpSet plot for CL effect ----------------------------------------------

CL <- list(
  mVAT_up   = wide |> filter(mVAT_CL_FDR <= 0.05, mVAT_CL_NES > 0) |> pull(Description),
  SAT_up    = wide |> filter(SAT_CL_FDR  <= 0.05, SAT_CL_NES  > 0) |> pull(Description),
  iBAT_up   = wide |> filter(iBAT_CL_FDR <= 0.05, iBAT_CL_NES > 0) |> pull(Description),
  
  mVAT_down = wide |> filter(mVAT_CL_FDR <= 0.05, mVAT_CL_NES < 0) |> pull(Description),
  SAT_down  = wide |> filter(SAT_CL_FDR  <= 0.05, SAT_CL_NES  < 0) |> pull(Description),
  iBAT_down = wide |> filter(iBAT_CL_FDR <= 0.05, iBAT_CL_NES < 0) |> pull(Description)
)

m <- make_comb_mat(CL)

UpSet(m, 
      top_annotation   = upset_top_annotation(m, add_numbers = TRUE),
      right_annotation = upset_right_annotation(m, add_numbers = TRUE),
      set_order = c("mVAT_up", "SAT_up", "iBAT_up", "mVAT_down", "SAT_down", "iBAT_down")
)
