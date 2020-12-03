topGO_wrapper = function(geneScores,
                         geneScoresDE = FALSE,
                         geneScoresDirection = c("Up","Down",NA),
                         GOmapping,
                         algorithm = c("classic", "elim", "weight", "weight01", "lea", "parentchild"),
                         statistic = c("fisher", "ks", "t", "globaltest", "sum"),
                         nodeSize = 5,
                         discretisedDE = F,
                         p = 0.05) {
  
  #if user specifies that gene scores are in form of TRUE/FALSE, choose significant genes based on cutoff
  #otherwise, gene scores are assumed to be supplied in a continuous form, e.g p-values
  topDiffGenes = ifelse(discretisedDE,
                        function(x)return(x < 0.05),
                        function(x)return(x))
  
  #if gene scores were provided as a DESeq2 output, process them to get them into the right form for topGO
  if(geneScoresDE==TRUE){
    
    #if a direction is supplied, cut down to only those genes that are either up or down regulated relative to the base condition
    if(is.na(geneScoresDirection)==FALSE){
      if(geneScoresDirection=="Down"){geneScores = subset(geneScores, log2FoldChange<0)
      } else                         {geneScores = subset(geneScores, log2FoldChange>0)}
    }
    
    DEGs = na.omit(geneScores$pvalue)
    names(DEGs) = rownames(geneScores)[which(is.na(geneScores$pvalue)==FALSE)]  
    geneScores = DEGs
    
  }
  
  ontologies = c("BP", "MF", "CC")
  
  for(i in 1:length(ontologies)){
    
    # instantiate a topGO object to work with
    GOdata = new(
      "topGOdata",
      ontology = ontologies[i],
      allGenes = geneScores,
      geneSel = topDiffGenes,
      annot = annFUN.gene2GO,
      gene2GO = GOmapping,
      nodeSize = nodeSize
    )
    
    # run enrichment test
    result = runTest(GOdata, algorithm = algorithm, statistic = statistic)
    assign(paste0("result",ontologies[i]), result)
    
    table = GenTable(GOdata, result = result, topNodes = 50, numChar = 1000)
    assign(paste0("genTable",ontologies[i]), table)
    
  }
  
  #summarise results in a simple table
  summaryTable = rbind(geneData(resultBP),geneData(resultMF),geneData(resultCC))
  row.names(summaryTable) = ontologies
  
  #summarise results ina table suitable for printing
  consolid = rbind(genTableBP, genTableMF, genTableCC)
  consolid$Ontology = c(rep("BP",nrow(genTableBP)),
                        rep("MF",nrow(genTableMF)),
                        rep("CC",nrow(genTableCC)))
  consolid$result = as.numeric(consolid$result)
  consolid = consolid[order(consolid$result),] 
  consolid = consolid[order(consolid$Ontology),] %>% 
    subset(result<p) %>%
    dplyr::select(GO.ID,Ontology,Term,result) 
  consolid = dplyr::mutate(consolid, result = formatC(consolid$result, digits = 2, format = "fg"))
  consolid$result = as.numeric(consolid$result)
  consolid$result[consolid$result<0.00001] = formatC(consolid$result[consolid$result<0.00001],digits = 1,format="e")
  
  out = list(summary = summaryTable, 
             BP_result = genTableBP,
             MP_result = genTableMF,
             CC_result = genTableCC,
             consolidated_result = consolid)
  
  return(out)
  
}
