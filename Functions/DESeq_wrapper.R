DESeq_wrap = function(data = data.phenotype.clean, genes = data.gene.count.clean, effects = NULL, formula = NULL, alpha = 0.05, estimate_disp = TRUE, detail = FALSE, lfc = 1){
  
  if(is.null(effects) & is.null(formula)){stop("Must specify either a list of additive effects or a character string specifying a formula")}
  #if effects were supplied, produce a simple additive model
  if(is.null(effects)==FALSE){  
    formula = as.formula(paste0("~ ",paste(effects, collapse = " + ")))
  }
  #if formula was supplied, use that instead of simply adding the effects
  if(is.null(formula)==FALSE){
    formula = as.formula(formula)
    effects = all.vars(formula)
  }
  data.foo = data[complete.cases(dplyr::select(data,c(effects))),]
  genes.foo = genes[ , as.character(data.foo$WaspID)]
  dds.gene.model = DESeqDataSetFromMatrix(countData = genes.foo,
                                          colData = data.foo,
                                          design = formula)
  #if user requires, genewise dispersion estimates are simply taken as the mean of normalizes counts (i.e. no data sharing)
  #this exists so that user can check out how impactful information sharing actually is
  if(estimate_disp==TRUE){
    dds.gene.model.deg = DESeq(dds.gene.model, quiet = TRUE)
  }else if(estimate_disp==FALSE){
    dds.gene.model.deg = estimateSizeFactors(dds.gene.model)
    dds.gene.model.deg = estimateDispersionsGeneEst(dds.gene.model.deg)
    dispersions(dds.gene.model.deg) = mcols(dds.gene.model.deg)$dispGeneEst
    dds.gene.model.deg = nbinomWaldTest(dds.gene.model.deg, quiet = TRUE)
  }
  comparisons = resultsNames(dds.gene.model.deg)
  resultstable = data.frame()
  details = c()
  for(i in 1:length(comparisons)){
    comparison = results(dds.gene.model.deg, 
                         name = comparisons[i],     
                         alpha = alpha,
                         lfcThreshold = log2(lfc))
    DEGs = (length(na.omit(which(comparison$padj<alpha))))
    DEGs_pct = (length(na.omit(which(comparison$padj<alpha)))/length(na.omit(comparison$padj)))*100
    resultstable = rbind(resultstable, data.frame("variable" = comparisons[i],
                                                  "features" = DEGs,
                                                  "features_pct" = DEGs_pct))
    if(detail == TRUE){details = c(details,comparison)}
  }
  if(detail == TRUE){
    names(details) = comparisons
    return(details)
  } else { 
    return(resultstable)
  }
  
}