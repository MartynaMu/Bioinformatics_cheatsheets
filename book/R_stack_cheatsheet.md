R stack for -omics data
================
Martyna Siewiera
2024-02-14

# Scope

This is a quick overview of my current R stack that I use for proteomic
data analysis. The aim of this cheatsheet is to organize packages and
their standard use as it is easy to get lost in all available
Bioconductor packages that in the end don’t work well or are missing a
lot of functions and need to be combined with another packages. This is
what has worked for me so far and I wanted a place where everything is
summed up in a somewhat sequential order of analysis.

So far I have been mainly working with two or three experimental groups
but there was an instance I had to create long functions just to
automatize DE, term enrichment and their visualization for multiple
groups. These functions can be found under “3.4 Functions for multiple
groups”.

## 1. affy

Currently used for testing normalization performance with MA plots. Blue
line represents theoretical perfect loess curve, red line - actual loess
curve. The more similar they are, the better.

<b>Note:</b> Normalization principle is that majority of features
(proteins/genes) do not change their expression under most conditions
and so their trend is used to normalize samples.

M = log2(x/y) aka log2FC

A = (log2(x) + log2(y))/2 = log2(xy)\*1/2, aka means of means=average
expression

``` r
library(affy)
# Named list with column indexes for each experimental group
groups <- list(1:3,4:6,7:9)
names(groups) <- colnames(design) #from limma pkg - search below

#take 2 experimental groups from matrix
x <- rowMeans(qnorm[,c(groups$x3D)])
y <- rowMeans(qnorm[,c(groups$x2D)])

M <- x - y
A <- (x + y)/2

#plot and save
maplot <- ma.plot(A = A, M = M, cex = 1)
title("... normalisation")
```

## 2. limma

Originally dedicated to microarray and RNA-seq data, which is also
useful for proteomics. Currently used for normalization methods and
differential expression analysis.

1.  Load required packages

    ``` r
    library(tidyverse)
    library(limma)
    ```

2.  Prepare an expression matrix - log2 transformed, normalized, tidy.

    ``` r
    mat <- read.delim("expr_matrix.tsv")
    mat <- as.matrix(mat)
    head(mat)

    #different normalization methods, assess with boxplots with jitter and MA plots
    # normalizeVSN()
    # normalizeQuantiles()
    # normalizeCyclicLoess()
    # normalizeMedianValues()
    ```

3.  Create a design table, which is a 0/1 matrix of assigning samples to
    an experimental group. There are many ways to do it, but I find the
    following one the easiest. In the <code>model.matrix()</code>
    function in <code>rep()</code> first digit represents the column
    index that will be filled with nr of replicates in a group indicated
    by second digit.

    ``` r
    design <- model.matrix(~ 0+factor(c(rep(1,3),
                                        rep(2,3),
                                        rep(3,3),
                                        rep(4,3))))
    # assign names of experimental groups
    colnames(design) <- c("ctrl", "treat1", "treat2", "treat3")
    # assign sample names from expression matrix to exp groups or manually
    row.names(design) <- colnames(mat)
    design
    ```

4.  Create an lm fit and matrix of contrasts - meaning comparisons
    between groups that will be made. To make things clear in next
    steps, assign a number in the comment as it will be needed for coef
    arguments later.

    ``` r
    fit <- lmFit(mat,design)
    cont.matrix <- makeContrasts(ctrl_treat1 = ctrl - treat1, #1
                                 ctrl_treat2 = ctrl - treat2, #2
                                 ctrl_treat3 = ctrl - treat3, #3
                                 levels = design) 
    ```

5.  Fit the contrasts and calculate Bayes statistics.

    ``` r
    fit_cont <- contrasts.fit(fit, cont.matrix)
    fit_bayes <- eBayes(fit_cont)
    ```

6.  Write down the whole differential expression matrix into a
    data.frame. ‘coef’ argument indicates which comparison is taken into
    account. See that expression matrix has to have protein/gene ids in
    rownames.

    ``` r
    DEP <- topTable(fit_bayes, 
                    coef=1, 
                    adjust.method = "BH", 
                    genelist = rownames(mat), 
                    resort.by = "logFC", 
                    number = nrow(mat))
    ```

7.  For a quick volcano plot do this:

    ``` r
    volcanoplot(fit_bayes,
                coef=1,
                highlight = 10,
                names = rownames(mat))
    abline(h = 1.2, v = c(-1.5,1.5), lty = 2)
    ```

8.  A much more preferred version:

    ``` r
    deps <- topTable(fit_bayes, coef = 1, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))

    #Retain columns for volcano plot and assign Category of expression
    volc_data <- deps
    volc_data$Category <- "Not significant"
    volc_data$Category[volc_data$logFC <= -1.3 & volc_data$adj.P.Val <= 0.05] <- "Down-regulated"
    volc_data$Category[volc_data$logFC >= 1.3 & volc_data$adj.P.Val <= 0.05] <- "Up-regulated"

    #Assign label of DEPs
    volc_data$Label <- NA
    volc_data$Label[volc_data$Category != "Not significant"] <- volc_data$ID[volc_data$Category != "Not significant"]

    #Plot
    theme_set(theme_light())
    v <- ggplot(volc_data,
                aes(x = logFC,
                    y = (-1*log10(adj.P.Val)),
                    color = Category,
                    label = Label))+
      geom_point(alpha = .5)+
      scale_color_manual(values = c("Down-regulated" = "blue", "Not significant" = "gray", "Up-regulated" = "red"))+
      labs(title = names(coefs[i]),
           y = "-log10 adj. p-value",
           x = "log2 fold change")+
      geom_hline(yintercept = 1.3, lty=2)+
      geom_vline(xintercept = c(-1.3,1.3), lty=2)+
      scale_x_continuous(limits = c(-5.5,5.5))+
      scale_y_continuous(limits = c(0,10.5))

    #Additional formatting, adding n of upregs and downregs
    v+
      theme(legend.position = "top",
            text = element_text(size=15))+
      geom_text_repel(colour = "black", max.overlaps = 15, size = 3.5)+
      annotate("text", x = -5.5, y=0, col= "blue",label=paste0("n=",length(downreg)))+
      annotate("text", x = 5.5, y=0, col= "red",label=paste0("n=",length(upreg)))
    ```

## 3. clusterProfiler

Used for Gene Ontology, KEGG, Reactome and WikiPathways GSEA or ORA
analysis. Very useful feature for simplifying term results (removing
similar/redundant terms, easing up the interpretation). Used with
‘enrichplot’ to visualize results.

**Important note:** This package likes to not cooperate, esp. if it
comes to plots and gseKEGG functions. To resolve the problem consider
these:

- make sure you have the latest version of BiocManager;

- update clusterProfiler, enrichplot, org db, dependancies etc;

- If still not working, use
  “remotes::install_github(”YuLab-SMU/clusterProfiler”))” to download
  dev version - last time worked;

- still nothing - download KEGG db and use internal db.  

1.  Load libraries

    ``` r
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    ```

2.  Prepare geneList which is a named list of genes/proteins with their
    log2 fold changes/pvalues from DE like limma’s topTable.

    ``` r
    geneList <- DEP[,2] #fc
    names(geneList) <- as.character(DEP[,1]) #gene symbols
    geneList <- sort(geneList, decreasing = TRUE) #sort ranks
    ```

3.  Some functions, like ‘gseKEGG’ don’t accept symbols as input, in
    this case it is needed to convert symbols into eg. ENTREZID.

    ``` r
    entrezid <- bitr(as.character(DEP[,1]), 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

    geneList2 <- filter(DEP, ID %in% entrezid$SYMBOL)[,2] #fc - filter out unmapped
    names(geneList2) <- as.character(entrezid$ENTREZID) #gene symbols
    geneList2 <- sort(geneList2, decreasing = TRUE) #sort ranks
    ```

4.  Sometimes you want to use external tools like WebGESTALT. Write down
    the geneList into tab-seperated file:

    ``` r
    write.table(geneList, "genelist.txt", sep = "\t", quote = FALSE)
    ```

    ### 3.1 GO GSEA

    ``` r
    gsea <- gseGO(geneList = geneList,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP", # or CC, MF
                  minGSSize = 50, # min nr of genes mapped to term
                  maxGSSize = 100, # max nr of genes mapped to term
                  pvalueCutoff = 0.05,
                  verbose = TRUE,
                  keyType = "SYMBOL")
    ```

- To view results from GSEA:

  ``` r
  gsea_results <- gsea@result
  upreg <- filter(gsea_results, NES > 0 & p.adjust < 0.05)
  downreg <- filter(gsea_results, NES < 0 & p.adjust < 0.05)
  head(gsea_results)
  ```

- Sometimes it is possible to merge similar terms through ‘simplify()’:

  ``` r
  simp <- clusterProfiler::simplify(gsea)
  simp_results <- simp@result
  ```

  ### 3.2 KEGG GSEA

  ``` r
  kegg <- clusterProfiler::gseKEGG(geneList = geneList2,
                  organism = "hsa", 
                  keyType = "ncbi-geneid", # if ENTREZID in geneList
                  minGSSize = 20,
                  maxGSSize = 500,
                  pvalueCutoff = 0.05,
                  verbose = TRUE)

  kegg_result <- kegg@result
  head(kegg_result)
  #this function converts entrezid back to readable symbols
  kegg <- setReadable(kegg)
  ```

  ### 3.3 Visualisation of GSEA results

  You can choose from barplots, dotplots, ridgeplots, umaps, cnetplots
  etc. Check clustProfiler and enrichplot vignette for possibilities.

  ``` r
  ridgeplot(kegg) +
    geom_vline(xintercept = 0, lty = 2)
  ```

  In case you’d like to make manual barplot, which is nicely ordered by
  p-value:

  ``` r
  library(ggplot2)
  theme_set(theme_minimal())

  upreg %>% 
    arrange(desc(p.adjust)) %>%
    mutate(Description = factor(Description, levels = Description)) %>%
    ggplot()+
      geom_col(aes(x = setSize, y = Description, fill = -log10(p.adjust)))+
      scale_fill_gradient(low = "blue", high = "red")
  ```

  ### 3.4 Functions for multiple groups

  ``` r
  # Create a named list with coef numbers for topTable and its name of comparison
  coefs <- seq.int(1,3,1)
  names(coefs) <- c("ctrlVtrt1", "ctrlVtrt2", "trt1Vtrt2")
  ```

  - GO GSEA

    ``` r
    compute_gsea <- function(fit, GO = c("BP|CC|MF"), coefs, simplify = TRUE) {
      ## Case 1. Simplified table results + figures---------
      if (simplify == TRUE) {
        path_figs <- paste0("figures/gsego_results/gsego_", GO, "_results/unfiltered/")
        path_data <- paste0("data/gsego_results/gsego_", GO, "_results/unfiltered/simplified/")

        if (dir.exists(path_figs) == FALSE) {
          dir.create(path_figs,showWarnings = TRUE,recursive = TRUE)
        }

        if (dir.exists(path_data) == FALSE) {
          dir.create(path_data,showWarnings = TRUE,recursive = TRUE)
        }

        for (i in coefs) {
          filename_data <- paste0("gsego_", GO, "_", names(coefs[i]), ".tsv")
          filename_fig <- paste0("gsego_", GO, "_top30_", names(coefs[i]), ".png")

          deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))

          geneList <- deps[,2] #fc
          names(geneList) <- as.character(deps[,1]) #gene symbols
          geneList <- sort(geneList, decreasing = TRUE) #sort ranks

          # define keytype for genes
          keyType = "SYMBOL"

          # GO GSEA 
          # min/max GSSize means how many genes per term should be used
          gsea <- gseGO(geneList = geneList,
                        OrgDb = org.Hs.eg.db,
                        ont = GO,
                        minGSSize = 15,
                        maxGSSize = 100,
                        pvalueCutoff = 0.05,
                        verbose = TRUE,
                        keyType = keyType)


          # remove redundant GO terms via GOSemSim methods
          simp <- clusterProfiler::simplify(gsea)
          simp_results <- simp@result
          write.table(file = paste0(path_data, filename_data), 
                      x = simp_results, 
                      quote = FALSE, 
                      sep = "\t",
                      row.names = FALSE)

          # Visualize
          ridgeplot(simp, label_format = 60)+
            geom_vline(xintercept = 0, lty = 2)+
            labs(x = "log2FC")

          ggsave(filename = filename_fig,
                 device = "png",
                 path = path_figs,
                 width = 800,
                 height = 780,
                 units = "px",
                 dpi = 100)
        }
      } else if (simplify == FALSE) {
        ## Case 2. NOT simplified table results ONLY---------
          path_data <- paste0("data/gsego_results/gsego_", GO, "_results/unfiltered/notsimplified/")

          if (dir.exists(path_data) == FALSE) {
            dir.create(path_data,showWarnings = TRUE,recursive = TRUE)
          }

          for (i in coefs) {
            filename_data <- paste0("gsego_", GO, "_", names(coefs[i]), ".tsv")

            deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))

            geneList <- deps[,2] #fc
            names(geneList) <- as.character(deps[,1]) #gene symbols
            geneList <- sort(geneList, decreasing = TRUE) #sort ranks

            # define keytype for genes
            keyType = "SYMBOL"

            # GO GSEA 
            # min/max GSSize means how many genes per term should be used
            gsea <- gseGO(geneList = geneList,
                          OrgDb = org.Hs.eg.db,
                          ont = GO,
                          minGSSize = 15,
                          maxGSSize = 100,
                          pvalueCutoff = 0.05,
                          verbose = TRUE,
                          keyType = keyType)


            gsea_results <- gsea@result
            write.table(file = paste0(path_data, filename_data), 
                        x = gsea_results, 
                        quote = FALSE, 
                        sep = "\t",
                        row.names = FALSE)
        }
      }
    }
    ```

  - KEGG GSEA

    ``` r
    compute_keggsea <- function(fit, coefs) {
      path_figs <- paste0("figures/gsekegg_results/unfiltered/")
      path_data <- paste0("data/gsekegg_results/unfiltered/")
      if (dir.exists(path_figs) == FALSE) {
        dir.create(path_figs,showWarnings = TRUE,recursive = TRUE)
      }

      if (dir.exists(path_data) == FALSE) {
        dir.create(path_data,showWarnings = TRUE,recursive = TRUE)
      }

      for (i in coefs) {
        filename_data <- paste0("gsekegg_", names(coefs[i]), ".tsv")
        filename_fig <- paste0("gsekegg_", "top30_", names(coefs[i]), ".png")

        deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))

        entrezid <- bitr(as.character(deps[,1]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        geneList <- filter(deps, ID %in% entrezid$SYMBOL)[,2] #fc
        names(geneList) <- as.character(entrezid$ENTREZID) #gene symbols
        geneList <- sort(geneList, decreasing = TRUE) #sort ranks

        kegg <- gseKEGG(geneList = geneList,                
                        organism = "hsa", 
                        keyType = "ncbi-geneid",
                        minGSSize = 20,
                        maxGSSize = 100,
                        pvalueCutoff = 0.05,
                        verbose = TRUE)

        kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

        kegg_results <- kegg@result
        write.table(file = paste0(path_data, filename_data), 
                    x = kegg_results, 
                    quote = FALSE, 
                    sep = "\t",
                    row.names = FALSE)

        # Visualize
        ridgeplot(kegg, label_format = 60)+
          geom_vline(xintercept = 0, lty = 2)+
          labs(x = "log2FC")

        ggsave(filename = filename_fig,
               device = "png",
               path = path_figs,
               width = 800,
               height = 780,
               units = "px",
               dpi = 100)
      }
    }
    ```

  - GO ORA

    ``` r
    enrich_go <- function(fit, coefs, GO="BP|CC|MF", simplify = TRUE) {
      # simplify=TRUE -------------------------
      if (simplify == TRUE) {
        path_figs <- paste0("figures/enrichgo_results/enrich_", GO, "_results/unfiltered/")
        path_data <- paste0("data/enrichgo_results/enrich_", GO, "_results/unfiltered/simplified/")

        if (dir.exists(path_figs) == FALSE) {
          dir.create(path_figs,showWarnings = TRUE,recursive = TRUE)
        }

        if (dir.exists(path_data) == FALSE) {
          dir.create(path_data,showWarnings = TRUE,recursive = TRUE)
        }

        for (i in coefs) {
          filename_data <- paste0("enrich_", GO, "_sim_", names(coefs[i]))
          filename_fig <- paste0("enrich_", GO, "_top30_sim_", names(coefs[i]))

          deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))

          geneList <- deps[,2] #fc
          names(geneList) <- as.character(deps[,1]) #gene symbols
          geneList <- sort(geneList, decreasing = TRUE) #sort ranks

          # define keytype for genes
          keyType = "SYMBOL"

          #extract ids
          upreg <- deps$ID[deps$logFC >= 1.3 & deps$adj.P.Val <= 0.05]
          downreg <- deps$ID[deps$logFC <= -1.3 & deps$adj.P.Val <= 0.05]
          message(paste0(names(coefs[i]),": ", "Nr of upregulated genes:", length(upreg)))
          message(paste0(names(coefs[i]),": ", "Nr of downregulated genes:", length(downreg)))

          #enrich upreg
          go_enr <- enrichGO(gene = upreg,
                             universe = names(geneList),
                             OrgDb = org.Hs.eg.db,
                             keyType = "SYMBOL",
                             ont = GO,
                             pvalueCutoff = 0.05,
                             pAdjustMethod = "BH")
          go_enr_results <- go_enr@result
          write.table(file = paste0(path_data, filename_data, "_upreg.tsv"), 
                      x = go_enr_results, 
                      quote = FALSE, 
                      sep = "\t",
                      dec = ",",
                      row.names = FALSE)
          #Visualize upreg
          if(min(go_enr_results$p.adjust) <= 0.05) {
            go_enr_sim <- simplify(go_enr)
            barplot(go_enr_sim, showCategory = 30, label_format = 60)
            ggsave(filename = paste0(filename_fig,"_upreg.png"),
                   device = "png",
                   path = path_figs,
                   width = 800,
                   height = 780,
                   units = "px",
                   dpi = 100)
          } else {
            message(paste0("No enrichment for upregulated genes in ", names(coefs[i])))
          }
          #enrich downreg
          go_enr <- enrichGO(gene = downreg,
                             universe = names(geneList),
                             OrgDb = org.Hs.eg.db,
                             keyType = "SYMBOL",
                             ont = GO,
                             pvalueCutoff = 0.05,
                             pAdjustMethod = "BH")
          go_enr_results <- go_enr@result
          write.table(file = paste0(path_data, filename_data, "_downreg.tsv"), 
                      x = go_enr_results, 
                      quote = FALSE, 
                      sep = "\t",
                      dec = ",",
                      row.names = FALSE)
          #Visualize downreg
          if(min(go_enr_results$p.adjust) <= 0.05) {
            go_enr_sim <- simplify(go_enr)
            barplot(go_enr_sim, showCategory = 30, label_format = 60)
            ggsave(filename = paste0(filename_fig,"_downreg.png"),
                   device = "png",
                   path = path_figs,
                   width = 800,
                   height = 780,
                   units = "px",
                   dpi = 100)
          } else {
            message(paste0("No enrichment for downregulated genes in ", names(coefs[i])))
          }
          # simplify=FALSE ----------------------------------------------
        } 
      } else if (simplify == FALSE) {
        path_data <- paste0("data/enrichgo_results/enrich_", GO, "_results/unfiltered/not_simplified/")

        if (dir.exists(path_figs) == FALSE) {
          dir.create(path_figs,showWarnings = TRUE,recursive = TRUE)
        }

        if (dir.exists(path_data) == FALSE) {
          dir.create(path_data,showWarnings = TRUE,recursive = TRUE)
        }

        for (i in coefs) {
          filename_data <- paste0("enrich_", GO, "_notsim_", names(coefs[i]))

          deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))

          geneList <- deps[,2] #fc
          names(geneList) <- as.character(deps[,1]) #gene symbols
          geneList <- sort(geneList, decreasing = TRUE) #sort ranks

          # define keytype for genes
          keyType = "SYMBOL"

          #extract ids
          upreg <- deps$ID[deps$logFC >= 1.3 & deps$adj.P.Val <= 0.05]
          downreg <- deps$ID[deps$logFC <= -1.3 & deps$adj.P.Val <= 0.05]
          message(paste0(names(coefs[i]),": ", "Nr of upregulated genes:", length(upreg)))
          message(paste0(names(coefs[i]),": ", "Nr of downregulated genes:", length(downreg)))

          #enrich upreg
          go_enr <- enrichGO(gene = upreg,
                             universe = names(geneList),
                             OrgDb = org.Hs.eg.db,
                             keyType = "SYMBOL",
                             ont = GO,
                             pvalueCutoff = 0.05,
                             pAdjustMethod = "BH")
          go_enr_results <- go_enr@result
          write.table(file = paste0(path_data, filename_data, "_notsim_upreg.tsv"), 
                      x = go_enr_results, 
                      quote = FALSE, 
                      sep = "\t",
                      dec = ",",
                      row.names = FALSE)

          if(min(go_enr_results$p.adjust) >= 0.05) {
            message(paste0("No enrichment for upregulated genes in ", names(coefs[i])))
          } 
          #enrich downreg
          go_enr <- enrichGO(gene = downreg,
                             universe = names(geneList),
                             OrgDb = org.Hs.eg.db,
                             keyType = "SYMBOL",
                             ont = GO,
                             pvalueCutoff = 0.05,
                             pAdjustMethod = "BH")
          go_enr_results <- go_enr@result
          write.table(file = paste0(path_data, filename_data, "_notsim_downreg.tsv"), 
                      x = go_enr_results, 
                      quote = FALSE, 
                      sep = "\t",
                      dec = ",",
                      row.names = FALSE)
          if(min(go_enr_results$p.adjust) >= 0.05) {
            message(paste0("No enrichment for downregulated genes in ", names(coefs[i])))
          }
        }
      }  
    }
    ```

  - Volcanoes

    ``` r
    for (i in coefs) {
      if (dir.exists("figures/volcanoes/") == FALSE) {
        dir.create("figures/volcanoes/",showWarnings = TRUE,recursive = TRUE)
      }
      if (dir.exists("data/comparison_stats/") == FALSE) {
        dir.create("data/comparison_stats/",showWarnings = TRUE,recursive = TRUE)
      }

      deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))

      volc_data <- deps
      volc_data$Category <- "Not significant"
      volc_data$Category[volc_data$logFC <= -1.3 & volc_data$adj.P.Val <= 0.05] <- "Down-regulated"
      volc_data$Category[volc_data$logFC >= 1.3 & volc_data$adj.P.Val <= 0.05] <- "Up-regulated"
      volc_data$Label <- NA
      volc_data$Label[volc_data$Category != "Not significant"] <- volc_data$ID[volc_data$Category != "Not significant"]

      upreg <- volc_data$ID[volc_data$logFC >= 1.3 & volc_data$adj.P.Val <= 0.05]
      downreg <- volc_data$ID[volc_data$logFC <= -1.3 & volc_data$adj.P.Val <= 0.05]

      write.table(file = paste0("data/comparison_stats/", names(coefs[i]), ".tsv"), 
                  x = volc_data, 
                  quote = FALSE, 
                  sep = "\t",
                  row.names = FALSE)

      theme_set(theme_light())
      v <- ggplot(volc_data,
                  aes(x = logFC,
                      y = (-1*log10(adj.P.Val)),
                      color = Category,
                      label = Label))+
        geom_point(alpha = .5)+
        scale_color_manual(values = c("Down-regulated" = "blue", "Not significant" = "gray", "Up-regulated" = "red"))+
        labs(title = names(coefs[i]),
             y = "-log10 adj. p-value",
             x = "log2 fold change")+
        geom_hline(yintercept = 1.3, lty=2)+
        geom_vline(xintercept = c(-1.3,1.3), lty=2)+
        scale_x_continuous(limits = c(-5.5,5.5))+
        scale_y_continuous(limits = c(0,10.5))
      v+
        theme(legend.position = "top",
              text = element_text(size=15))+
        geom_text_repel(colour = "black", max.overlaps = 15, size = 3.5)+
        annotate("text", x = -5.5, y=0, col= "blue",label=paste0("n=",length(downreg)))+
        annotate("text", x = 5.5, y=0, col= "red",label=paste0("n=",length(upreg)))

      ggsave(filename = paste0("volcano", names(coefs[i]), ".png"),
             path = "figures/volcanoes/",
             device = "png",
             width = 640,
             height = 660,
             units = "px",
             dpi = 100)
    }
    ```

  - ANOVA filtering

    For multiple comparisons it is best to check if removing low
    variance proteins results in better clustering/differentiation -
    e.g. with heatmap/dendrogram and/or PCA

    ``` r
    # drop 'coef' argument for anova
    anova <- topTable(fit_bayes, adjust = "BH", genelist = rownames(mat), resort.by = "adj.P.val", number = nrow(qnorm))

    anova_filtr <- anova$ProbeID[anova$adj.P.Val <= 0.05]
    temp <- filter(qnorm, !Gene %in% anova_filtr)

    pheatmap::pheatmap(temp[1:18],
                       scale = "row",
                       show_rownames = FALSE,
                       cutree_rows = 5,
                       clustering_distance_cols = "euclidean",
                       clustering_distance_rows = "correlation",
                       annotation_col = annot_col,
                       annotation_colors = annot_colors,
                       main = "ANOVA filtered")
    ```

  - Visualize terms on a heatmap

    ``` r
    term = "negative regulation of cell cycle process"
    filtr <- str_split_1(gsea_results$core_enrichment[gsea_results$Description == term], pattern = "/")

    pheatmap::pheatmap(qnorm[qnorm$Gene %in% filtr,1:18],
                       scale = "row",
                       show_rownames = TRUE,
                       clustering_distance_cols = "euclidean",
                       clustering_distance_rows = "correlation",
                       annotation_col = annot_col,
                       annotation_colors = annot_colors,
                       main = term)
    ```

## 4. pheatmap

Global heatmap with experimental groups as column annotation.

``` r
# Tidying up--------------------------------------------------------------------
#drop nas and set rownames
msqrob <- drop_na(msqrob)
msqrob <- column_to_rownames(msqrob, "Gene")
msqrob <- as.matrix(msqrob[1:18])

# Global heatmap-----------------------------------------------------------------------
#Column annotation - each column is a named grouping variable
annot_col <- data.frame(
  Condition.L1 = as.factor(
    c(rep(c("2D", "3D.young", "3D.old"), each = 4),
      rep(c("PDX.2D", "PDX.3D"), each = 3))),
  Condition.L2 = as.factor(
    c(rep("2D", 4), 
      rep("3D", 8), 
      rep("PDX", 6))))
#set sample names as rownames
row.names(annot_col) <- colnames(msqrob)[1:18]

#Colors mapping for column annotation
#A named list for each grouping with color values
library(RColorBrewer)
annot_colors <- list(Condition.L1 = rev(brewer.pal(5, "Paired")),
                     Condition.L2 = c("firebrick1", "darkgreen", "blue"))
#Map names of experimental groups to the colors
names(annot_colors$Condition.L1) <- levels(annot_col$Condition.L1)
names(annot_colors$Condition.L2) <- levels(annot_col$Condition.L2)

#Plot hm
hm <- pheatmap::pheatmap(msqrob[1:18],
           scale = "row",
           show_rownames = FALSE,
           clustering_distance_cols = "euclidean",
           clustering_distance_rows = "correlation",
           annotation_col = annot_col,
           annotation_colors = annot_colors,
           main = "Global heatmap")
```

I like to check clustering after ANOVA filtering \<0.05 adj.p-value and
then after DE (limma pkg below).

Explore row dendrogram of hm after filtering of matrix.

``` r
filtr_deps <- unlist(deps_list) |> unique()

# View row tree, decide on nr of clusters
hm$tree_row %>% 
  as.dendrogram() %>%
  plot(horiz = TRUE)

# assign genes to clusters and annotate in hm, extract and explore in string network/go
clusters <- cutree(hm$tree_row, k = 5)
clusters <- as.data.frame(clusters) |> mutate(clusters=as.factor(clusters))
clipr::write_clip(rownames(clusters)[clusters$clusters == 5])

#plot hm again with clustered rows separated
hm <- pheatmap::pheatmap(qnorm[qnorm$Gene %in% filtr_deps,1:18],
                   scale = "row",
                   show_rownames = FALSE,
                   cutree_rows = 5,
                   annotation_row = clusters,
                   clustering_distance_cols = "euclidean",
                   clustering_distance_rows = "correlation",
                   annotation_col = annot_col,
                   annotation_colors = annot_colors,
                   main = "Significant DEPs")
```

Full list of all DEPs across comparisons:

``` r
#create deps list if multiple comparisons
deps_list <- function(fit){
  comb <- list()
  for (i in coefs) {
    deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))
    
    volc_data <- deps
    volc_data$Category <- "Not significant"
    volc_data$Category[volc_data$logFC <= -1.3 & volc_data$adj.P.Val <= 0.05] <- "Down-regulated"
    volc_data$Category[volc_data$logFC >= 1.3 & volc_data$adj.P.Val <= 0.05] <- "Up-regulated"
    volc_data$Label <- NA
    volc_data$Label[volc_data$Category != "Not significant"] <- volc_data$ID[volc_data$Category != "Not significant"]
    
    upreg <- volc_data$Label[volc_data$Category == "Up-regulated"]
    downreg <- volc_data$Label[volc_data$Category == "Down-regulated"]
    
    comb <- append(comb, c(list(upreg), list(downreg)))
  }
  comb_names <- rep((names(coefs)),each=2,times=1) |> paste0(rep(c("_upreg", "_downreg"), every=2))
  names(comb) <- comb_names
  return(comb)
}
```

## 5. PCAtools

Principal component analysis. Better results after removing low variance
features (e.g. after ANOVA).

``` r
library(PCAtools)
#column to rownames for pca function
mat <- as.matrix(msqrob[1:18])
#pca
p <- pca(mat, metadata = annot_col, scale = TRUE) #annot_col same as in pheatmap

#plot PCA
biplot(p, 
       lab = colnames(mat), 
       labSize = 5, 
       showLoadings = FALSE, #show genes driving the variation
       boxedLoadingsNames = FALSE, 
       colby = "Condition.L1", #colors of points
       colkey = annot_colors$Condition.L1, #names of points, same as in pheatmap
       title = "This is PCA biplot",
       encircle = TRUE, #encircle groups (no stats applied - for that use 'eclipse' arg)
       hline = 0, #add 0,0 coordinates
       vline = 0)

#show screeplot
screeplot(p)

#loadings plot to screen genes driving the variation in components
plotloadings(p)
```

## 6. eulerr

Create area-proportional venn/euler diagrams. Most frequently used for
peptide/protein overlap between samples.

``` r
#create a named list with peptide/proteins ids for each group/sample/comparison
comb <- list(l1, l2, l3, l4)
names(comb) <- c("Sample1", "Sample2", "Sample3", "Sample4")

euler(comb) |> plot(quantities = TRUE, #area-proportional?
                    main = "title", 
                    fills = c("#AB3428", "#FFB703", "white")) #circle colors
```

## 7. corrplot

Visualize correlation matrix in a heatmap. Assess pre- and
post-normalization.

``` r
#compute correlation coefs on matrix w/o NAs
corr <- cor(log2(mat), use="complete.obs")
#visualize
corrplot(corr, method = "color", hclust.method = "ward.D2", is.corr = FALSE, col = COL1("OrRd"))
```

## 8. ggpubr & patchwork

Organize figures in ready-to-publish panels.

# Session info

``` r
sessionInfo()
```

    ## R version 4.3.2 (2023-10-31 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19045)
    ## 
    ## Matrix products: default
    ## 
    ## 
    ## locale:
    ## [1] LC_COLLATE=Polish_Poland.utf8  LC_CTYPE=Polish_Poland.utf8   
    ## [3] LC_MONETARY=Polish_Poland.utf8 LC_NUMERIC=C                  
    ## [5] LC_TIME=Polish_Poland.utf8    
    ## 
    ## time zone: Europe/Warsaw
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] corrplot_0.92              eulerr_7.0.0              
    ##  [3] PCAtools_2.14.0            ggrepel_0.9.5             
    ##  [5] pheatmap_1.0.12            enrichplot_1.22.0         
    ##  [7] clusterProfiler_4.11.0.001 limma_3.58.1              
    ##  [9] affy_1.80.0                Biobase_2.62.0            
    ## [11] BiocGenerics_0.48.1        lubridate_1.9.3           
    ## [13] forcats_1.0.0              stringr_1.5.1             
    ## [15] dplyr_1.1.4                purrr_1.0.2               
    ## [17] readr_2.1.5                tidyr_1.3.1               
    ## [19] tibble_3.2.1               ggplot2_3.4.4             
    ## [21] tidyverse_2.0.0           
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3        rstudioapi_0.15.0        
    ##   [3] jsonlite_1.8.8            magrittr_2.0.3           
    ##   [5] farver_2.1.1              rmarkdown_2.25           
    ##   [7] fs_1.6.3                  zlibbioc_1.48.0          
    ##   [9] vctrs_0.6.5               DelayedMatrixStats_1.24.0
    ##  [11] memoise_2.0.1             RCurl_1.98-1.14          
    ##  [13] ggtree_3.10.0             S4Arrays_1.2.0           
    ##  [15] htmltools_0.5.7           SparseArray_1.2.3        
    ##  [17] gridGraphics_0.5-1        plyr_1.8.9               
    ##  [19] cachem_1.0.8              igraph_2.0.1.1           
    ##  [21] lifecycle_1.0.4           pkgconfig_2.0.3          
    ##  [23] rsvd_1.0.5                Matrix_1.6-5             
    ##  [25] R6_2.5.1                  fastmap_1.1.1            
    ##  [27] gson_0.1.0                GenomeInfoDbData_1.2.11  
    ##  [29] MatrixGenerics_1.14.0     digest_0.6.34            
    ##  [31] aplot_0.2.2               colorspace_2.1-0         
    ##  [33] patchwork_1.2.0           AnnotationDbi_1.64.1     
    ##  [35] S4Vectors_0.40.2          dqrng_0.3.2              
    ##  [37] irlba_2.3.5.1             RSQLite_2.3.5            
    ##  [39] beachmat_2.18.0           fansi_1.0.6              
    ##  [41] timechange_0.3.0          abind_1.4-5              
    ##  [43] httr_1.4.7                polyclip_1.10-6          
    ##  [45] compiler_4.3.2            bit64_4.0.5              
    ##  [47] withr_3.0.0               BiocParallel_1.36.0      
    ##  [49] viridis_0.6.5             DBI_1.2.1                
    ##  [51] ggforce_0.4.1             MASS_7.3-60.0.1          
    ##  [53] DelayedArray_0.28.0       HDO.db_0.99.1            
    ##  [55] tools_4.3.2               ape_5.7-1                
    ##  [57] scatterpie_0.2.1          glue_1.7.0               
    ##  [59] nlme_3.1-164              GOSemSim_2.28.1          
    ##  [61] grid_4.3.2                shadowtext_0.1.3         
    ##  [63] reshape2_1.4.4            fgsea_1.28.0             
    ##  [65] generics_0.1.3            gtable_0.3.4             
    ##  [67] tzdb_0.4.0                preprocessCore_1.64.0    
    ##  [69] data.table_1.15.0         hms_1.1.3                
    ##  [71] ScaledMatrix_1.10.0       BiocSingular_1.18.0      
    ##  [73] tidygraph_1.3.1           utf8_1.2.4               
    ##  [75] XVector_0.42.0            pillar_1.9.0             
    ##  [77] yulab.utils_0.1.4         splines_4.3.2            
    ##  [79] tweenr_2.0.2              treeio_1.26.0            
    ##  [81] lattice_0.22-5            bit_4.0.5                
    ##  [83] tidyselect_1.2.0          GO.db_3.18.0             
    ##  [85] Biostrings_2.70.2         knitr_1.45               
    ##  [87] gridExtra_2.3             IRanges_2.36.0           
    ##  [89] stats4_4.3.2              xfun_0.41                
    ##  [91] graphlayouts_1.1.0        statmod_1.5.0            
    ##  [93] matrixStats_1.2.0         stringi_1.8.3            
    ##  [95] lazyeval_0.2.2            ggfun_0.1.4              
    ##  [97] yaml_2.3.8                evaluate_0.23            
    ##  [99] codetools_0.2-19          ggraph_2.1.0             
    ## [101] qvalue_2.34.0             BiocManager_1.30.22      
    ## [103] ggplotify_0.1.2           cli_3.6.2                
    ## [105] affyio_1.72.0             munsell_0.5.0            
    ## [107] Rcpp_1.0.12               GenomeInfoDb_1.38.5      
    ## [109] png_0.1-8                 parallel_4.3.2           
    ## [111] blob_1.2.4                DOSE_3.28.2              
    ## [113] sparseMatrixStats_1.14.0  bitops_1.0-7             
    ## [115] viridisLite_0.4.2         tidytree_0.4.6           
    ## [117] scales_1.3.0              crayon_1.5.2             
    ## [119] rlang_1.1.3               cowplot_1.1.3            
    ## [121] fastmatch_1.1-4           KEGGREST_1.42.0

# More sources

- <https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf>

- <https://yulab-smu.top/biomedical-knowledge-mining-book/index.html>

- <https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/>
