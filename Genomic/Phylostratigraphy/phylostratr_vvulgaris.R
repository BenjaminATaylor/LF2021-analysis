library(phylostratr)
library(tidyverse)

setwd("/home/benjamin/Documents/LF_2020_repo/Genomic/Phylostratigraphy/Vvulgaris")

phylostrata_pdom = read.csv("/home/benjamin/Documents/LF_2020_repo/Genomic/Phylostratigraphy/phylostrata_vvulgaris.csv")

#bin phylostrata
unique(phylostrata_pdom[,c("ps","mrca_name")])
phylostrata_pdom$bin = dplyr::case_when(phylostrata_pdom$ps %in% c(18:19) ~ "Aculeata",
                                   phylostrata_pdom$ps %in% c(16:17) ~ "Hymenoptera",
                                   phylostrata_pdom$ps %in% c(13:15) ~ "Insecta",
                                   #phylostrata_pdom$ps %in% c(4:9) ~ "Metazoa",
                                   phylostrata_pdom$ps %in% c(1:12) ~ "Ancient")
#phylostrata_pdom$qseqid = substr(phylostrata_pdom$qseqid,1,14)
phylostrata_pdom = phylostrata_pdom %>% dplyr::select(c("bin","V2.y")) %>%
  distinct()

# divide into groups based on DESEq2 caste-bias
phylo_up = subset(phylostrata_pdom, V2.y %in% genes_up_Pdom_caste)$bin
phylo_down = subset(phylostrata_pdom, V2.y %in% genes_down_Pdom_caste)$bin
phylo_neutral = subset(phylostrata_pdom, (V2.y %in% row.names(data.gene.count.clean)) & !(V2.y %in% c(genes_up_Pdom_caste, genes_down_Pdom_caste)))$bin

# divide into groups based on SVM caste-bias
# phylo_up = subset(phylostrata_pdom, V2.y %in% keyfeatures_up)$bin
# phylo_down = subset(phylostrata_pdom, V2.y %in% keyfeatures_down)$bin
# phylo_neutral = subset(phylostrata_pdom, (V2.y %in% row.names(data.gene.count.clean)) & !(V2.y %in% c(keyfeatures_up, keyfeatures_down)))$bin

#plot
strataframe =  rbind(data.frame(treatment = "Queen",strata = phylo_up),
                     data.frame(treatment = "Worker",strata = phylo_down),
                     data.frame(treatment = "Neutral",strata = phylo_neutral)) %>%
  #mutate(strata = factor(strata,levels = c("Lineage-specific","Hymenoptera","Arthropoda","Metazoa","Ancient"))) %>%
  mutate(strata = factor(strata,levels = c("Aculeata","Hymenoptera","Insecta","Ancient"))) %>%
  mutate(treatment = factor(treatment,levels = c("Queen","Worker","Neutral")))

strataplot = ggplot(strataframe, aes(x = treatment, fill = strata)) +
  geom_bar(position = "fill", width = 0.9, colour = "black",) +
  #scale_fill_brewer(palette = "Blues",direction = -1)+
  scale_fill_manual(values = brewer.pal(6,"Blues")[3:6], guide = guide_legend(reverse=T)) +
  scale_x_discrete(labels = c("Queen\nbiased","Worker\nbiased","Unbiased"),
                   expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  xlab("Differential expression\nwith reproductive status") +
  ylab("Proportion") +
  coord_flip() +
  theme_bw() +
  theme(aspect.ratio = 1/4, 
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        legend.position = "top",
        legend.spacing.x = unit(0.5,"cm"),
        plot.title = element_text(size=18, face="bold", hjust = 0.5),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15, face = "bold"),
        axis.title.y = element_text(size=15, face = "bold"),
        panel.grid = element_line(size = 0.2, colour = "gray95"),
        axis.ticks = element_line(size=1),
        panel.border = element_rect(size=1.5),
        strip.text.x = element_text(size=15, face = "bold"),
        strip.background = element_rect(size = 1.5))

#save plot
ggsave(strataplot, filename = "Pdom_caste_phylostrata.png", device = "png", path ="/home/benjamin/Dropbox/Ben PhD/PhD Docs/THESIS/data_chapter_3/figures",
       width = 25, height = 8, units = "cm")


#make contignency table for hypothesis testing
stratacontingency = rbind(table(phylo_up),
                          table(phylo_down),
                          table(phylo_neutral)) %>% 'row.names<-'(c("queen","worker","neutral"))
#
# first test whether genes are more likely to be ancient in postive/negative vs neutral
stratacontingency_ancient = cbind(stratacontingency[,2],rowSums(stratacontingency[,c(1,3,4)]))
#negatively-correlated genes are more likely to be ancient
fisher.test(stratacontingency_ancient[c(1,3),])
#positively-correlated genes are not more likely to be ancient
fisher.test(stratacontingency_ancient[c(2,3),])

# first test whether genes are more likely to be ancient in postive/negative vs neutral
stratacontingency_aculeata = cbind(stratacontingency[,1],rowSums(stratacontingency[,c(2:4)]))
#negatively-correlated genes are more likely to be ancient
fisher.test(stratacontingency_aculeata[c(1,3),])
#positively-correlated genes are not more likely to be ancient
fisher.test(stratacontingency_aculeata[c(2,3),])
