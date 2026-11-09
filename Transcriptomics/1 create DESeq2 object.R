# Setup --------------------------------------------------------------------------
pacman::p_load(tidyverse, janitor, patchwork, ggrepel, DESeq2, openxlsx, knitr, readxl)
this.path::here() |> setwd()
rm(list=ls())


# Read in counts ----------------------------------------------------------
counts <- read_tsv("data/nfcore/star_salmon/salmon.merged.gene_counts.tsv")

identical( counts$gene_id, counts$gene_name )

counts <- counts |> 
  column_to_rownames("gene_id") |> 
  select(-gene_name) |> 
  as.matrix()

dim(counts) 
# [1] 37440    27



# mitochondrial encoded genes ------------------------------------------------------

mt_map <- c(
  # rRNA
  "mt-Rnr1" = "mt-Rnr1",
  "mt-Rnr2" = "mt-Rnr2",
  
  # tRNA
  "Trnf"  = "mt-Tf",
  "Trnv"  = "mt-Tv",
  "Trnl1" = "mt-Tl1",
  "Trni"  = "mt-Ti",
  "Trnq"  = "mt-Tq",
  "Trnm"  = "mt-Tm",
  "Trnw"  = "mt-Tw",
  "Trna"  = "mt-Ta",
  "Trnn"  = "mt-Tn",
  "Trnc"  = "mt-Tc",
  "Trny"  = "mt-Ty",
  "Trns1" = "mt-Ts1",
  "Trnd"  = "mt-Td",
  "Trnk"  = "mt-Tk",
  "Trng"  = "mt-Tg",
  "Trnr"  = "mt-Tr",
  "Trnh"  = "mt-Th",
  "Trns2" = "mt-Ts2",
  "Trnl2" = "mt-Tl2",
  "Trne"  = "mt-Te",
  "Trnt"  = "mt-Tt",
  "Trnp"  = "mt-Tp",
  
  # protein coding
  "ND1"  = "mt-Nd1",
  "ND2"  = "mt-Nd2",
  "COX1" = "mt-Co1",
  "COX2" = "mt-Co2",
  "ATP8" = "mt-Atp8",
  "ATP6" = "mt-Atp6",
  "COX3" = "mt-Co3",
  "ND3"  = "mt-Nd3",
  "ND4L" = "mt-Nd4l",
  "ND4"  = "mt-Nd4",
  "ND5"  = "mt-Nd5",
  "ND6"  = "mt-Nd6",
  "CYTB" = "mt-Cytb"
)

common <- intersect(names(mt_map), rownames(counts))   # all 37 found
length(common)

rowSums(counts)[ common ] |> sort(decreasing = TRUE)
#     COX1     COX3     ATP6     CYTB      ND4      ND2     COX2      ND1      ND5      ND3     ND4L      ND6     ATP8 
# 24200874 10627669  7797068  6692578  6460052  5778440  4872517  4265232  4142338  1654288   832723   270579   226341 
# mt-Rnr2  mt-Rnr1     Trnm    Trns2     Trnp     Trnf     Trnv    Trnl1     Trni     Trnq     Trnw     Trna     Trnn 
#    8793     8157        6        3        1        0        0        0        0        0        0        0        0 
# Trnc     Trny    Trns1     Trnd     Trnk     Trng     Trnr     Trnh    Trnl2     Trne     Trnt 
#    0        0        0        0        0        0        0        0        0        0        0 

counts <- counts[ rowSums(counts) > 0, ]  
dim(counts)
# [1] 26766    27


# Mapping for SYMBOL to ENSGID and Gene name ------------------------------

chr_map <- data.frame( refseq = paste0("NC_", sprintf("%06d", 86019:86038), ".1"), chr = 1:20 )

chr_map <- rbind.data.frame(chr_map, 
                            c("NC_086039.1", "X"),
                            c("NC_086040.1", "Y"),
                            c("NC_001665.2", "MT"))


gtf_raw <- data.table::fread("data/GCF_036323735.1_GRCr8_plus_MT.gtf.gz", 
                             sep = "\t", header = FALSE, quote = "", fill = TRUE)

colnames(gtf_raw)[1:9] <- c("seqnames", "source", "type", "start", "end", "score", "strand", "phase", "attribute")

gtf_map <- gtf_raw |> 
  filter(type == "gene") |> 
  transmute(
    refseq       = seqnames,
    gene_start   = start,
    gene_end     = end,
    gene_id      = str_match( attribute, '(?:^|; )gene_id "([^"]+)"')[, 2],
    SYMBOL       = str_match( attribute, '(?:^|; )gene "([^"]+)"' )[, 2],
    ENTREZID     = str_match( attribute, 'db_xref "GeneID:([0-9]+)"')[, 2],
    GENENAME     = str_match( attribute, 'description "([^"]+)"')[, 2],
    gene_biotype = str_match( attribute, 'gene_biotype "([^"]+)"')[, 2]
  ) |>
  distinct() |> 
  data.frame() |> 
  filter(SYMBOL %in% rownames(counts)) |> 
  left_join(chr_map) |> 
  mutate(chr = ifelse(is.na(chr), "scaffold", chr)) |> 
  arrange(SYMBOL) |> 
  relocate(chr, .before = "gene_start") |> 
  relocate(SYMBOL, .before = "refseq")


gtf_map |> tabyl(chr)

gtf_map |> filter(is.na(ENTREZID))                                 # all found

gtf_map |> filter(refseq == "NC_001665.2") |> arrange(gene_start)  # found EntrezID for all expressed mitochondria encoded genes
# "mt-Rnr1" "mt-Rnr2" "ND1" "Trnm" "ND2" "COX1" "COX2" "ATP8" "ATP6" "COX3" "ND3" "ND4L" "ND4" "Trns2" "ND5" "ND6" "CYTB" "Trnp" 

table(is.na(gtf_map$ENTREZID))  # no missing value
table(table(gtf_map$SYMBOL))    # all unique
table(table(gtf_map$ENTREZID))  # all unique

rm(gtf_raw, chr_map, common); gc()



# Update mitochondrial gene symbols to match Mitocarta 3.0 -----------------

## update symbols for the counts ----
rn <- rownames(counts)
w  <- rn %in% names(mt_map)
data.frame( rn[w], mt_map[ rn[w] ] )

rn[w] <- mt_map[ rn[w] ]
rownames(counts) <- rn
rm(rn, w)

grep("^mt-", rownames(counts), v = T)


## update symbols for the gtf_map ----
gtf_map$gene_id <- NULL

rn <- gtf_map$SYMBOL
w  <- rn %in% names(mt_map)
data.frame( rn[w], mt_map[ rn[w] ] )

rn[w] <- mt_map[ rn[w] ]
gtf_map$SYMBOL <- rn
rownames(gtf_map) <- gtf_map$SYMBOL
rm(rn, w, mt_map)

grep("^mt-", rownames(gtf_map), v = T)



## Sort and check order is same
setdiff( rownames(counts),  rownames(gtf_map) )
setdiff( rownames(gtf_map), rownames(counts) )
identical( sort(rownames(counts)), sort(rownames(gtf_map) ))

counts  <- counts [ sort(rownames(counts)), ]
gtf_map <- gtf_map[ sort(rownames(gtf_map)), ]

identical( rownames(counts), rownames(gtf_map) )


# Prepare the sample annotation --------------------------------------------------
snames <- colnames(counts)

meta <- data.frame(
  sample    = snames,
  tissue    = str_sub(snames, 1, 2),
  diet      = ifelse(str_detect(snames, "HG"), "HGI", "LGI"),
  treatment = ifelse(str_detect(snames, "21D"), "CL", "S"),
  row.names = snames) |> 
  
  mutate(tissue = recode_values(tissue, "Bt" ~ "iBAT", "Mv" ~ "mVAT", "St" ~ "SAT"),
         tissue = factor(tissue, levels = c("mVAT", "SAT", "iBAT")),
         
         group = paste(diet, treatment, sep = "-"),
         group = factor(group, levels = c("LGI-S", "HGI-S", "HGI-CL"))
  )

meta |> tabyl(tissue, group)
# tissue LG-S HG-S HG-CL
#   iBAT    3    3     3
#   mVAT    3    3     3
#    SAT    3    3     3



# Create the DESeq2 object ------------------------------------------------
identical(colnames(counts), rownames(meta))

dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData   = meta,
                              design    = ~ -1 + tissue + group)

identical( rownames(counts), rownames(gtf_map) )
rowData(dds) <- gtf_map

rm(counts, meta, gtf_map, snames)

grp_colors <- c("LGI-S" = "#08519C", "LGI-CL" = "#6BAED6",  
                "HGI-S" = "#CB181D", "HGI-CL" = "#FC9272")

save(dds, grp_colors, file = "results/dds27.rda")



# Quality Control check ---------------------------------------------------


## PCA global -------------------------------------------------------------
vst_expr <- vst(dds)

ntop <- 3000
pc0 <- plotPCA(vst_expr, intgroup = "group", ntop = ntop) 
pc  <- plotPCA(vst_expr, intgroup = "group", ntop = ntop, returnData = T) 

head(pc)

ggplot(pc, aes(x = PC1, y = PC2, col = tissue, shape = group)) +
  geom_point(size = 5)  +
  theme_bw(base_size = 18) +
  labs(x=pc0@labels$x, y=pc0@labels$y) 

rm(vst_expr, ntop, pc0, pc, depot)


## PCA by tissue ----------------------------------------------------------
plots <- list()

for(depot in levels(dds$tissue)){
  
  vst_expr <- vst(dds[ , which(dds$tissue == depot)])
  
  ntop <- 3000
  pc0 <- plotPCA(vst_expr, intgroup = "group", ntop = ntop) 
  pc  <- plotPCA(vst_expr, intgroup = "group", ntop = ntop, returnData = T) 
  
  plots[[depot]] <- ggplot(pc, aes(x = PC1, y = PC2, label = name, col = group)) +
    geom_point() +
    facet_wrap( ~ tissue, scales = "free") +
    scale_color_manual(values = grp_colors) +
    ggrepel::geom_text_repel(show.legend = FALSE, size = 5) +
    theme_bw(base_size = 18) +
    theme(legend.position = "none") +
    labs(x=pc0@labels$x, y=pc0@labels$y) 
  
  rm(vst_expr, ntop, pc0, pc, depot)
}

wrap_plots(plots, nrow = 1)



# plot mitochondrial encoded protein coding genes -------------------------
keep <- c(
  "mt-Nd1", "mt-Nd2", "mt-Nd3", "mt-Nd4l", "mt-Nd4", "mt-Nd5", "mt-Nd6", # Complex I
  "mt-Cytb",                                                             # Complex III
  "mt-Co1", "mt-Co2", "mt-Co3",                                          # Complex IV
  "mt-Atp6", "mt-Atp8",                                                  # Complex V
  "mt-Rnr1", "mt-Rnr2"
)

vst_expr <- vst(dds)

mdf <- data.frame(
  tissue = vst_expr$tissue,
  group  = vst_expr$group,
  assay(vst_expr)[keep, ] |> t(),
  check.names = FALSE) |> 
  pivot_longer(cols = !c(group, tissue)) |> 
  mutate(name = factor(name, levels = keep))

ggplot(mdf, aes(x = group, y = value, fill = group)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  facet_grid(tissue ~ name) + 
  labs(x = NULL, y = "VST normalized expression") +
  theme_bw() +
  theme(axis.text.x  = element_blank(), axis.ticks.x = element_blank())

