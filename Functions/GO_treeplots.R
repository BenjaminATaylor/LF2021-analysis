GO_treeplots = function(GOframe, path){
  
  require(rrvgo)
  require(org.Dm.eg.db)
  require(treemap)
  require(grid)
  require(GOSemSim)
  require(RColorBrewer)
  
  #create similarity matrices in parent environment if necessary
  if(!exists("BPsim")){BPsim <<- GOSemSim::godata(org.Dm.eg.db, ont="BP")}
  if(!exists("MFsim")){MFsim <<- GOSemSim::godata(org.Dm.eg.db, ont="MF")}
  if(!exists("CCsim")){CCsim <<- GOSemSim::godata(org.Dm.eg.db, ont="CC")}
  
  #check numbers of terms in each ontology
  BPTRUE = nrow(subset(GOframe$consolidated_result, Ontology == "BP"))>1
  MFTRUE = nrow(subset(GOframe$consolidated_result, Ontology == "MF"))>1
  CCTRUE = nrow(subset(GOframe$consolidated_result, Ontology == "CC"))>1
  TRUEsum = sum(BPTRUE,CCTRUE,MFTRUE)
  
  #prepare plotting environment
  png(filename=path, width=22, height=45*(TRUEsum/3), units = "cm", res = 800)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(TRUEsum, 1)))
  row = 1
  aspR = 1.5
  
  #BP
  simMatrix <- suppressWarnings(calculateSimMatrix(GOframe$consolidated_result$GO.ID,
                                                   orgdb="org.Dm.eg.db",
                                                   semdata = BPsim,
                                                   ont="BP",
                                                   method="Rel"))
  scores <- setNames(-log10(as.numeric(GOframe$consolidated_result$result)), GOframe$consolidated_result$GO.ID)
  if(!BPTRUE){print("Insufficient GO terms for plotting in ontology: BP")
  } else {
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.9,
                                  orgdb="org.Dm.eg.db")
  vp <- viewport(layout.pos.col=1, layout.pos.row=1)
  suppressWarnings(treemap::treemap(reducedTerms, index = c("parentTerm", "term"), 
                                    vSize = "score", type = "index", 
                                    title = "Biological Process", fontsize.title = 18,
                                    palette = colorRampPalette(brewer.pal(n = 9 ,name = "Reds")[3:7])(length(unique(reducedTerms$parent))), 
                                    fontcolor.labels = c("black", "white"), fontface.labels = c("bold","italic"), fontsize.labels = c(18,8),
                                    bg.labels = 20, overlap.labels = 0.5,
                                    border.col = "#00000080", 
                                    vp = vp, vsize = 1/3, aspRatio = aspR))
  #increment position
  row = row + 1
  }
  
  #MF
  simMatrix <- suppressWarnings(calculateSimMatrix(GOframe$consolidated_result$GO.ID,
                                                   orgdb="org.Dm.eg.db",
                                                   ont="MF",
                                                   method="Rel"))
  scores <- setNames(-log10(as.numeric(GOframe$consolidated_result$result)), GOframe$consolidated_result$GO.ID)
  if(!MFTRUE){print("Insufficient GO terms for plotting in ontology: MF")
  } else {
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.9,
                                  orgdb="org.Dm.eg.db")
  vp <- viewport(layout.pos.col=1, layout.pos.row=row)
  suppressWarnings(treemap::treemap(reducedTerms, index = c("parentTerm", "term"), 
                                    vSize = "score", type = "index", 
                                    title = "Molecular Function", fontsize.title = 18,
                                    palette = colorRampPalette(brewer.pal(n = 9 ,name = "Greens")[3:7])(length(unique(reducedTerms$parent))), 
                                    fontcolor.labels = c("black", "white"), fontface.labels = c("bold","italic"), fontsize.labels = c(18,12),
                                    bg.labels = 20, overlap.labels = 0.5,
                                    border.col = "#00000080", 
                                    vp = vp, vsize = 1/3, aspRatio = aspR))
  row = row + 1
  }
  
  #CC
  simMatrix <- suppressWarnings(calculateSimMatrix(GOframe$consolidated_result$GO.ID,
                                                   orgdb="org.Dm.eg.db",
                                                   ont="CC",
                                                   method="Rel"))
  scores <- setNames(-log10(as.numeric(GOframe$consolidated_result$result)), GOframe$consolidated_result$GO.ID)
  if(!CCTRUE){print("Insufficient GO terms for plotting in ontology: CC")
    } else {
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.9,
                                  orgdb="org.Dm.eg.db")
    
  vp <- viewport(layout.pos.col=1, layout.pos.row=row)
  suppressWarnings(treemap::treemap(reducedTerms, index = c("parentTerm", "term"),  
                                    vSize = "score", type = "index", 
                                    title = "Cellular Component", fontsize.title = 18,
                                    palette = colorRampPalette(brewer.pal(n = 9 ,name = "Blues")[3:7])(length(unique(reducedTerms$parent))), 
                                    fontcolor.labels = c("black", "white"), fontface.labels = c("bold","italic"), fontsize.labels = c(18,12),
                                    bg.labels = 20, overlap.labels = 0.5,
                                    border.col = "#00000080", 
                                    vp = vp, vsize = 1/3, aspRatio = aspR))
    }
  dev.off()
  
}

GO_treeplots(darkturquoiseGO, path = "/home/benjamin/Desktop/treeplot_temp.png")
