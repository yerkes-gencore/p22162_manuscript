plotCountsGKT <- function (dds, gene, intgroup = "condition", normalized = TRUE, 
                           transform = TRUE, main, xlab = "group", returnData = FALSE, 
                           replaced = FALSE, pc, color, ...) 
{
  stopifnot(length(gene) == 1 & (is.character(gene) | (is.numeric(gene) & 
                                                         (gene >= 1 & gene <= nrow(dds)))))
  if (!all(intgroup %in% names(colData(dds)))) 
    stop("all variables in 'intgroup' must be columns of colData")
  if (!returnData) {
    if (!all(sapply(intgroup, function(v) is(colData(dds)[[v]], 
                                             "factor")))) {
      stop("all variables in 'intgroup' should be factors, or choose returnData=TRUE and plot manually")
    }
  }
  if (missing(pc)) {
    pc <- if (transform) 
      0.5
    else 0
  }
  if (is.null(sizeFactors(dds)) & is.null(normalizationFactors(dds))) {
    dds <- estimateSizeFactors(dds)
  }
  cnts <- counts(dds, normalized = normalized, replaced = replaced)[gene, 
  ]
  group <- if (length(intgroup) == 1) {
    colData(dds)[[intgroup]]
  }
  else if (length(intgroup) == 2) {
    lvls <- as.vector(t(outer(levels(colData(dds)[[intgroup[1]]]), 
                              levels(colData(dds)[[intgroup[2]]]), function(x, 
                                                                            y) paste(x, y, sep = ":"))))
    droplevels(factor(apply(as.data.frame(colData(dds)[, 
                                                       intgroup, drop = FALSE]), 1, paste, collapse = ":"), 
                      levels = lvls))
  }
  else {
    factor(apply(as.data.frame(colData(dds)[, intgroup, drop = FALSE]), 
                 1, paste, collapse = ":"))
  }
  #  data <- data.frame(count = cnts + pc, group = as.integer(group)) # orig
  #  dataGKT <- data.frame(count = cnts, group = as.integer(group))
  #  dataGKT <- data.frame(count = counts(dds, normalized = normalized, replaced = replaced)[gene,], colData(dds)[intgroup])
  #  dataGKT <- data.frame(count = counts(dds, normalized = normalized, replaced = replaced)[gene,], group = as.factor(group))
  dataGKT <- data.frame(count = counts(dds, normalized = normalized, replaced = replaced)[gene,], group = group)
  #  str(dataGKT) # diagnostic
  #  str(group) # diagnostic
  #  str(xgroup) #diagnostic
  logxy <- if (transform) 
    "y"
  else ""
  if (missing(main)) {
    main <- if (is.numeric(gene)) {
      rownames(dds)[gene]
    }
    else {
      gene
    }
  }
  ylab <- ifelse(normalized, "normalized count", "count")
  if (returnData) 
    #    return(data.frame(count = data$count, colData(dds)[intgroup])) # orig
    return(data.frame(count = dataGKT$count, colData(dds)[intgroup]))
  #    plot(data$group + runif(ncol(dds), -0.05, 0.05), data$count, 
  #       xlim = c(0.5, max(data$group) + 0.5), log = logxy, xaxt = "n", 
  #       xlab = xlab, ylab = ylab, main = main, ...)
  #  axis(1, at = seq_along(levels(group)), levels(group))
  (ggplot(dataGKT, aes(group, count))
    + geom_violin()
    + geom_boxplot(width = 0.2,outlier.shape = NA)
    #    + geom_sina()
    + xlab(paste(intgroup,collapse = ":"))
    + theme(axis.text.x = element_text(angle = 90))
    + labs(title = main)
  )
}

pcaDataGKT <- function (object, intgroup = "condition", ntop = 500) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(pca$x, group = group, 
                  intgroup.df, name = colnames(object))
  attr(d, "percentVar") <- percentVar
  return(d)
}

pcaPlotGKT <- function(object, intgroup = "condition", xpc = 1, ypc = 2, ntop = 500)
{
  d <- pcaDataGKT(object, intgroup, ntop)
  percentVar <- round(100 * attr(d, "percentVar"))
  ggplot(d,aes_string(x = names(d)[xpc], y = names(d)[ypc])) + labs(x=paste0(names(d)[xpc],": ", percentVar[xpc], "% variance"), y=paste0(names(d)[ypc],": ", percentVar[ypc], "% variance"))
}

str_wrap_balance <- function(string, width = 80, indent = 0, exdent = 0,USE.NAMES = FALSE) {
  out <- str_wrap(string, width, indent, exdent)
  vapply(out, function(string,width,indent,exdent) {
    wraps <- str_count(string,"\n")
    if(wraps > 0 && width > 1) {
      bwidth <- width
      repeat {
        bwidth <- bwidth - 1
        bstring <- str_wrap(string, bwidth, indent, exdent)
        bwraps <- str_count(bstring,"\n")
        if(bwraps > wraps || bwidth <= 1) break
        string <- bstring
      }
    }
    string
  },character(1),width,indent,exdent,USE.NAMES = USE.NAMES)
}


as.mutate.data.frame <- function(object, var = "gene_id", subset = "all", relabel = "", contrastlabel = "", select = c("gene","baseMean","log2FoldChange","pvalue","padj","significant","filterThreshold","subset","contrast"))
{
  df <- rownames_to_column(as.data.frame(object), var = var) %>%
    mutate(significant = factor(if_else(is.na(padj),"removed",if_else(padj<object@metadata$alpha,"significant","not_significant")),
                                levels = c("significant","not_significant","removed")),
           subset = subset,
           contrast = contrastlabel,
           filterThreshold = as.numeric(object@metadata$filterThreshold)) %>%
    rename_at(vars(2:7), function(x) paste0(names(object)[1:6], " (", object@elementMetadata$description[1:6], ")"))
  
  if(subset != "all") {
    subset
    names(df)[2:7]
    #    df <- df %>% rename_at(vars(2:7), function(x) gsub(": ", paste0(": ",subset," "), gsub("all",subset, names(df)[2:7])))
    df <- df %>% rename_at(vars(2:7), function(x) gsub("all", subset, names(df)[2:7]))
  }
  if(relabel != "") {
    df <- df %>% rename_at(vars(2:7), function(x) gsub(relabel, "", names(df)[2:7], fixed = TRUE))
    relabel_nous <- gsub("_"," ",relabel)
    df <- df %>% rename_at(vars(2:7), function(x) gsub(relabel_nous, "", names(df)[2:7], fixed = TRUE))
    # df <- df %>% rename_at(vars(2:7), function(x) gsub("p-values", paste0("p-values: ",newlabel), names(df)[2:7], fixed = TRUE))
  }
  # if(relabel != "" & newlabel != "") {
  #   df <- df %>% rename_at(vars(2:7), function(x) gsub(relabel, newlabel, names(df)[2:7], fixed = TRUE))
  #   relabel_nous <- gsub("_"," ",relabel)
  #   df <- df %>% rename_at(vars(2:7), function(x) gsub(relabel_nous, newlabel, names(df)[2:7], fixed = TRUE))
  #   df <- df %>% rename_at(vars(2:7), function(x) gsub("p-values", paste0("p-values: ",newlabel), names(df)[2:7], fixed = TRUE))
  # }
  if(length(select) > 1) {
    df <- df %>% select(contains(select))
  }
  df
}


wrap_underscore_strings_balance <- function(string, width) {
  str_replace_all(str_wrap_balance(str_replace_all(string,"_"," "),width)," ","_")
}

greg_volplot <- function(result){
  ## Greg's volcano plot for his filter
  gVolData <- as.mutate.data.frame(result)
  gVolPlot <- gVolData %>%
    dplyr::arrange(!is.na(.data[[colnames(gVolData)[5]]]),desc(.data[[colnames(gVolData)[5]]])) %>%
    ggplot(aes(x=.data[[colnames(gVolData)[2]]],
               y=.data[[colnames(gVolData)[3]]])) +
    geom_vline(aes(xintercept = filterThreshold), color = "goldenrod") +
    geom_text(aes(x = filterThreshold, y=-5,
                  label = format(filterThreshold, digits=3)),
              color = "goldenrod", hjust = 0.5, vjust = 1) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    geom_point(aes(color = significant), size = 0.1) +
    scale_x_log10() +
    theme_bw() +
    scale_color_manual(values = c("significant" = "red",
                                  "not_significant" = "green",
                                  "removed" = "blue",
                                  "filterThreshold" = "goldenrod"), drop = FALSE) + labs(y="Log2FoldChange") + labs(y="Log2FoldChange")
  return(gVolPlot)
}

##### Derrik functions

### makes a heatmap of the given list of genes, separating samples slightly by group variable
heatmap_from_genelist <- function(geneList, data=analysis$rldDrop){
  ### makes a heatmap of the given list of genes, separating samples slightly by group variable
  ## data should be of type DESeqTransform
  hmap <- data[rownames(data) %in% geneList,
               !colnames(data) %in% analysis$config$dropSamples]
  slices <- as.numeric(hmap@colData$Group)
  labels <- levels(hmap@colData$Group)[unique(slices)]
  baseline <- rowMedians(assay(hmap[,hmap@colData$Group %in% analysis$config$designBaseline]))
  hmap <- assay(hmap) - baseline
  Heatmap(hmap, show_row_names = TRUE, heatmap_legend_param = list(title="log2 fold\ndifference\nfrom\nmedian\nweek 0\nexpression"),border="black",column_order=colnames(hmap),
          #  width = ncol(hmap)*unit(5, "mm"), 
          # height = nrow(hmap)*unit(5, "mm"),
          rect_gp = gpar(color="black"),
          column_title=" ",
          cluster_rows=FALSE,
          column_split=slices,
          top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(col=NA),
                                                              #,fill = 2:length(unique(exp.design$Group))+1),
                                                              labels = labels, 
                                                              labels_gp = gpar(col = "black", fontsize = 10), labels_rot=90, height=unit(8,"cm"))), #
          column_gap=unit(2, "mm"),
          col=colorRamp2(c(-4, 0, 4), c("blue", "white", "red")))
}

### Make dotplot from fGSEA result, showing top n pathways
gsea_dotplot <- function(res, n=20, findtoppaths = TRUE, title='Top Enriched Pathways' ){
  ### Make dotplot from fGSEA result, showing top n pathways
  res <- res %>% mutate(test=((-log10(pval))*abs(NES)))
  res <- res %>% mutate(perc=100*lengths(leadingEdge)/size) %>% mutate(name=paste0(wrap_underscore_strings_balance(pathway,36), "\nn=",size))
  if (findtoppaths){
    res <- rbind(head(res%>% arrange(desc(test)) %>% filter(ES>0), n=round(n/2)), head(res %>% arrange(desc(test)) %>% filter(ES<0), n=round(n/2)))
  }
  res$name <- factor(res$name)
  ggplot(res) + geom_point(aes(x=perc, y=name, size=-log10(pval), color=NES)) +
    scale_color_gradient2(low="blue", mid="white",high="red", midpoint=0, breaks=c(-2,-1,0,1,2), limits=c(min(res$NES,-1), max(res$NES,1))) +
    theme_classic(11) +
    theme(panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_line(colour = "grey92"),
          panel.grid.major.y = element_line(colour = "grey92"),#element_blank(),
          panel.grid.minor.y = element_line(colour = "grey92")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x="% of genes in leading edge", y="Gene set", color = "Normalized\nenrichment\nscore", size="Nom p-val", title=title, caption='GSEA p-value calculations are not continuous,\nthere may be several or many pathways with the same p-value\nn = number of genes in pathway') +  scale_radius(name="NOM p-val", range=c(1,8), breaks=-log10(c(.1,0.01, 0.001, 0.0001)), limits=c(0,4), labels=c(0.1,0.01, 0.001, 0.0001)) + scale_y_discrete(limits=res$name)
}


### Read Cls File for fGSEA
GSEA.ReadClsFile <- function(file = "NULL") {
  ### Read Cls File for fGSEA
    cls.cont <- readLines(file)
    num.lines <- length(cls.cont)
    cls.cont[[3]] <- gsub("\\t", " ", cls.cont[[3]])  #Converts any tabs to spaces
    class.list <- unlist(strsplit(cls.cont[[3]], " "))  #Splits CLS on spaces
    s <- length(class.list)
    t <- table(class.list)[c(unique(class.list))]
    l <- length(t)
    phen <- vector(length = l, mode = "character")
    phen.label <- vector(length = l, mode = "numeric")
    class.v <- vector(length = s, mode = "numeric")
    for (i in 1:l) {
      phen[i] <- noquote(names(t)[i])
      phen.label[i] <- i - 1
    }
    for (i in 1:s) {
      for (j in 1:l) {
        if (class.list[i] == phen[j]) {
          class.v[i] <- phen.label[j]
        }
      }
    }
    return(list(phen = phen, class.v = class.v))
}


### Mapping results bins plot
mapping_plot <- function(mapBins=analysis$mapBins, sort=FALSE){
  if (sort) {
    analysis$mapBins <- analysis$mapBins[,mixedorder(colnames(analysis$mapBins), decreasing = TRUE)]
  }
  rownames_to_column(as_tibble(mapBins, rownames = NA), var = "map_result") %>%
  pivot_longer(!map_result, names_to = "SampleID", values_to = "count") %>%
  #mutate(SampleID = str_sub(SampleID,end = -18)) %>%
  ggplot(aes(x = SampleID, y = count, fill = factor(map_result, levels = c("N_unmapped","N_multimapping","N_noFeature","N_ambiguous","N_identified")))) +
  geom_bar(stat = "identity") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  guides(fill=guide_legend(title="Map Result")) +
  ggtitle(paste0(analysis$config$analysis," Mapping")) + labs(x="Sample") + if (!sort) {aes(fct_inorder(SampleID)) }
}

### RLE plots
RLE_plot <- function(data, title){
  as_tibble(data, rownames = NA) %>% pivot_longer(everything(), names_to = "Sample", values_to = "RLE") %>% 
    ggplot(aes(x=Sample,y=RLE)) +
    geom_hline(yintercept = 0, color = "red") + geom_violin(draw_quantiles = c(0.25,0.75), trim = TRUE, color = "lightgreen", alpha = 0.1) +
    geom_boxplot(alpha = 0) + theme_bw() + theme(axis.text.x = element_text(vjust = 0.5,angle = 90),axis.title.x = element_blank(),aspect.ratio = 0.55) + 
    ggtitle(title) + aes(fct_inorder(Sample))#+ scale_y_continuous(limits = c(-9,3.25), expand = c(0,0))
}

### PCA plot from config, dependant on pcaPlotGKT data, just adds some overlays
PCA_plot_from_config <- function(analysis=analysis){
  pcaPlotGKT(analysis$vst,
             intgroup = names(colData(analysis$vst)),
             xpc = 1, ypc = 2) +
    #  geom_path(aes(group = SubjectID)) 
    geom_point(aes(color = (if (is.null(analysis$config$pcaMapping$color)) NULL else .data[[analysis$config$pcaMapping$color]]),
                   shape = (if (is.null(analysis$config$pcaMapping$shape)) NULL else .data[[analysis$config$pcaMapping$shape]])),
               size = 5) +
    labs(color=analysis$config$pcaMapping$color, shape=analysis$config$pcaMapping$shape) +
    (if (is.null(analysis$config$pcaMapping$label)) NULL else geom_text_repel(aes(label = .data[[analysis$config$pcaMapping$label]]),
                    size = 4, hjust = 0.5, vjust = -0.5, alpha=0.5)) +
    scale_x_continuous(expand = c(0.5,0)) +
    theme_bw() + ggtitle(paste0(analysis$config$analysis," PCA")) +
    (if (!is.null(analysis$config$pcaMapping$path)){
      geom_path(aes(linetype=.data[[analysis$config$pcaMapping$path]]))
    }) +
    # geom_path(aes(linetype=(if (is.null(analysis$config$pcaMapping$path)) NULL else .data[[analysis$config$pcaMapping$path]]))) +
    theme(legend.key.width = unit(1.2, "cm")) +
    #labs(color="Weeks") +
    (if (!is.null(analysis$config$pcaMapping$ellipse)){
      stat_ellipse(aes(color=.data[[analysis$config$pcaMapping$ellipse]]), type="norm", level=0.67)
    })+
    theme(text = element_text(size=10)) #, arrow=arrow(ends="last", type="closed", length=unit(0.1, "inches")))
  
  
}

#### DESeq2 results object
### Generate volcano plot from DESeq2 results table
### also references config GeneList to label genes of interest
generate_volplot <- function(result){
  ## Volcano plot
  volData <- result[!is.na(result$padj),]
  volData <- volData[order(volData$padj),]
  labs <- if (!is.null(analysis$config$geneList[[1]])) analysis$config$geneList else rownames(volData[1:20,])
  volplot <- (EnhancedVolcano(volData,
                              x = 'log2FoldChange',
                              y = 'padj',
                              lab = rownames(volData),
                              selectLab=labs,
                              drawConnectors = TRUE,
                              colConnectors = "lightgrey",
                              pCutoff = analysis$config$alpha,
                              FCcutoff = log2(1.3),
                              title =NULL,# paste0(numer_contrast, " vs ", denom_contrast),
                              caption=NULL,
                              subtitle=NULL,
                              labSize = 3
                              
  )) 
  return(volplot)
}

gen_text_summary <- function(result){
  ## Text summary
  volData <- result[!is.na(result$padj),]
  volData <- volData[order(volData$padj),]
  summary_out <- capture.output(summary(result))
  
}

make_worksheet <- function(result, wb){
  ## Write results to excel sheeet in an open workbook
  pair <- result@metadata$contrast
  if (nchar(pair)>31){
    pair <- substr(pair,1,31)
  }
  volData <- result[!is.na(result$padj),]
  volData <- volData[order(volData$padj),]
  addWorksheet(wb, pair)
  writeData(wb, sheet=pair, x=as.data.frame(volData), rowNames=TRUE)
}

write_output_workbook <- function(wb, results, analysis){
  change_colnames <- function(x){
    a <- rownames_to_column(as.data.frame(x), var="Gene")
    colnames(a) <- c("Gene", "baseMean (mean of normalized counts)", "log2FoldChange (Log2 fold change MLE)", "lfcSE (Log fold change standard error)", "stat","pvalue (Wald test p-value)", "padj (BH adjusted p-values)")
    a
  }
  merge_tables <- lapply(results, change_colnames)
  suffixes <- paste0(".", names(merge_tables))
  res <- merge_tables[[1]]
  for (i in head(seq_along(merge_tables),-1)){
    res <- merge(res, merge_tables[[i+1]], all=TRUE,
                 suffixes = suffixes[i:(i++1)], by="Gene")
  }
  res <- res %>% select(-contains("stat"))
  addWorksheet(wb, "all")
  writeData(wb, sheet="all", x=res, rowNames=FALSE)
  
  tmp <- analysis$rldDrop
  baseline <- rowMeans(assay(tmp[,tmp@colData$Group %in% analysis$config$designBaseline]))
  tmp <- assay(tmp) - baseline
  addWorksheet(wb, "rle")
  writeData(wb, sheet="rle", x= tmp, rowNames=TRUE)
}

my_results_transform <- function(result, wb, write_out=TRUE){
  setClass(
    ## custom class for storing results data
    "myDESRclass",
    contains="DESeqResults",
    slots=c("visualizations", "sig.genes", 'gsea')
  ) -> myDESRclass
  ## Convert results object to new class and generate all outputs
  result <- as(result, "myDESRclass")
  if (write_out){
    make_worksheet(result, wb)
  }
  result@visualizations$summary <- capture.output(summary(result))
  result@visualizations$volplot <- generate_volplot(result)
  result@visualizations$volplot_fixed <- generate_volplot(result) + ylim(c(0,25)) + xlim(c(-5,5))
  result@visualizations$greg_volplot <- greg_volplot(result)
  # result@visualizations$finalfig <- arrangeGrob(textGrob(paste(result@visualizations$summary[1:8],
  #                                                              collapse="\n")),
  #                       ggarrange(result@visualizations$volplot, result@visualizations$volplot_fixed,
  #                                 common.legend = T),
  #                       result@visualizations$greg_volplot,
  #                       ncol=2, nrow=2, layout_matrix=layout,
  #                       top=textGrob(result@metadata$contrast, gp=gpar(fontsize=25)))#, bottom="padj < 0.1, | FoldChange | > 1.3")#, fig.lab=unlist(out[[i]]$pair))
  result@gsea <- run_fgsea(result, gmt.file, nperm=10000)
  
  volData <- result[!is.na(result$padj),]
  volData <- volData[order(volData$padj),]
  result@sig.genes <-  volData[volData$padj<analysis$config$alpha,]@rownames
  return(result)
}

##### GSEA

run_fgsea <- function(result, gmt.file=gmt.file, nperm=1000){
  result <- result %>% na.omit()
  tmp1 <- result$stat
  names(tmp1) <- rownames(result)
  tmp1 <- tmp1[!is.na(tmp1)]
  res <- fgseaSimple(pathways=gmt.file,
                     stats=tmp1,
                     nperm=nperm)
  return(res)
}
