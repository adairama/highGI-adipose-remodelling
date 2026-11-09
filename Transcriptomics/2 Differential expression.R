# Setup -------------------------------------------------------------------
pacman::p_load(tidyverse, janitor, patchwork, ggrepel, DESeq2, openxlsx, knitr, readxl)
this.path::here() |> setwd()
rm(list=ls())


# Gene filtering and running differential expression ---------------------
load("results/dds27.rda")

# For each tissue depot, we restrict to genes that are expressed (at least 10 counts per sample) in 3 or more samples.

dds_list <- list()

for(depot in c("mVAT", "SAT", "iBAT")){
  
  dds2 <- dds[ , which(dds$tissue == depot)]
  keep <- rowSums(counts(dds2) >= 10) >= 3
  dds2 <- dds2[keep, ]
  
  design(dds2) <- ~ -1 + group
  dds2 <- DESeq(dds2)
  
  dds_list[[ depot ]] <- dds2
  rm(depot, dds2, keep)
}


# Number of genes tested in each depot. 
#  14,472 genes in the intersect of all three depots. 
#  17,455 genes in the universe.

ggvenn::ggvenn( lapply(dds_list, rownames), show_percentage = FALSE )

all_genes <- Reduce("union", lapply(dds_list, rownames))
length(all_genes)


## Extract DEG results -----------------------------------------------------
extract_results <- function(obj, name = "group", 
                            num, denom, prefix,
                            FC_cut = log2(1.5), FDR_cut = 0.05){
  
  out <- results(obj, contrast = c(name, num, denom)) %>% 
    data.frame() %>%
    rownames_to_column("SYMBOL") %>%
    select(SYMBOL, baseMean, LFC=log2FoldChange, p=pvalue, FDR=padj) %>%
    arrange(p) %>% 
    mutate( sig = sign(LFC) * (abs(LFC) > FC_cut & FDR < FDR_cut),
            sig = factor(sig, levels = c(-1, 0, 1))) 
  
  colnames(out)[-c(1,2)] <- paste(prefix, colnames(out)[-c(1,2)], sep = ".")
  return(out)
}

big.res  <- lapply(dds_list, function(obj){
  
  res <- NULL
  res[[1]] <- extract_results(obj, num = "HGI.S",  denom = "LGI.S", prefix = "Diet") 
  res[[2]] <- extract_results(obj, num = "HGI.CL", denom = "HGI.S", prefix = "CL")
  
  res <- Reduce("full_join", res)
  return(res)
})

gtf_map <- rowData(dds) |> data.frame()


## Number of significant genes for Diet effect
sapply(big.res, function(x) table(x$Diet.sig))
#     mVAT   SAT  iBAT
# -1  2219  1630   220
#  1  3217  2026   162


## Number of significant genes for CL effect
sapply(big.res, function(x) table(x$CL.sig))
#     mVAT   SAT iBAT
# -1  3040  1596 2473
#  1  2043   949 3403


## Format and write results to Excel (wide format)
for(depot in names(big.res)){
  df <- big.res[[ depot ]]
  colnames(df)[-1] <- paste0(depot, "_", colnames(df)[-1])
  big.res[[ depot ]] <- df
}

wide <- Reduce("full_join", big.res)
wide <- left_join(wide, gtf_map)

wide |> filter( grepl("^mt-", SYMBOL) )


# Overlap with MitoCarta 3.0 ----------------------------------------------
mc <- read_excel("data/Mouse.MitoCarta3.0.xls", sheet = "C MitoPathways") |> 
  separate_longer_delim(Genes, delim = ", ") |> 
  distinct(SYMBOL = Genes) |> 
  mutate(mitocarta = 1)

wide <- left_join(wide, mc) |> 
  mutate(mitocarta = replace_na(mitocarta, 0))

table(wide$mitocarta)

save(big.res, wide, gtf_map, file = "results/DEG.rda")
write.xlsx(wide, "results/DEG_results_wide_tmp.xlsx")



# Reproduce Fig 3 of https://www.science.org/doi/10.1126/sciadv.ae --------

dim(mc)  # Total mitochondrial gene: 1037

for(x in paste0(levels(dds$tissue), "_Diet.sig")){
  cat("\n\n", x, "\n")
  table(wide$mitocarta, wide[ , x], exclude = NULL) |> print()
}

for(x in paste0(levels(dds$tissue), "_CL.sig")){
  cat("\n\n", x, "\n")
  table(wide$mitocarta, wide[ , x], exclude = NULL) |> print()
}



# Upset plot for Diet effect ----------------------------------------------
require(ComplexHeatmap)

diet <- list(
  
  mVAT_up   = wide |> filter(mVAT_Diet.sig ==  1) |> pull(SYMBOL),
  mVAT_down = wide |> filter(mVAT_Diet.sig == -1) |> pull(SYMBOL),

  SAT_up    = wide |> filter(SAT_Diet.sig ==  1) |> pull(SYMBOL),
  SAT_down  = wide |> filter(SAT_Diet.sig == -1) |> pull(SYMBOL),
  
  iBAT_up    = wide |> filter(iBAT_Diet.sig ==  1) |> pull(SYMBOL),
  iBAT_down  = wide |> filter(iBAT_Diet.sig == -1) |> pull(SYMBOL)
)
  
m <- make_comb_mat(diet)

UpSet(m, 
      top_annotation   = upset_top_annotation(m, add_numbers = TRUE),
      right_annotation = upset_right_annotation(m, add_numbers = TRUE),
      set_order = c("mVAT_up", "SAT_up", "iBAT_up", "mVAT_down", "SAT_down", "iBAT_down")
)


# Upset plot for CL effect ----------------------------------------------

CL <- list(
  
  mVAT_up   = wide |> filter(mVAT_CL.sig ==  1) |> pull(SYMBOL),
  mVAT_down = wide |> filter(mVAT_CL.sig == -1) |> pull(SYMBOL),
  
  SAT_up    = wide |> filter(SAT_CL.sig ==  1) |> pull(SYMBOL),
  SAT_down  = wide |> filter(SAT_CL.sig == -1) |> pull(SYMBOL),
  
  iBAT_up    = wide |> filter(iBAT_CL.sig ==  1) |> pull(SYMBOL),
  iBAT_down  = wide |> filter(iBAT_CL.sig == -1) |> pull(SYMBOL)
)

m <- make_comb_mat(CL)

UpSet(m, 
      top_annotation   = upset_top_annotation(m, add_numbers = TRUE),
      right_annotation = upset_right_annotation(m, add_numbers = TRUE),
      set_order = c("mVAT_up", "SAT_up", "iBAT_up", "mVAT_down", "SAT_down", "iBAT_down")
)
