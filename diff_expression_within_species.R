library(edgeR)
library(dplyr)
library(magrittr)
library(stringr)
library(readr)
library(tibble)
library(VennDiagram)
library(circlize)
library(ggplot2)
library(GGally)
library(intergraph)
library(sna)
library(network)
library(qgraph)


####################################
#        DE genes analysis         #
####################################

#read in featurecounts output
con_counts<-read.table("/Users/katie/Desktop/Bg/begonia_duplicate_expression/featurecount/con/con_counts", header=TRUE)
ple_counts<-read.table("/Users/katie/Desktop/Bg/begonia_duplicate_expression/featurecount/ple/ple_counts", header=TRUE)

# gen the gene IDs
con_geneid<-con_counts$Geneid
ple_geneid<-ple_counts$Geneid

# prepend them with CON_ and PLE_
#con_names<-paste("CON", con_geneid, sep="_")
#ple_names<-paste("PLE", ple_geneid, sep="_")

# remove all the path names fro  the column names to make them shorter
con_cols<-colnames(con_counts) %>% 
  str_remove("X.home.kemelianova.to_download.trimmed_decontaminated.star.") %>% 
  str_remove("Aligned.sortedByCoord.out.bam")

ple_cols<-colnames(ple_counts) %>% 
  str_remove("X.home.kemelianova.to_download.trimmed_decontaminated.star.") %>% 
  str_remove("Aligned.sortedByCoord.out.bam")

# rename columns using the shorter version
colnames(con_counts)<-con_cols
colnames(ple_counts)<-ple_cols

con_lengths<-con_counts %>% dplyr::select(Length) %>% set_rownames(con_geneid)
ple_lengths<-ple_counts %>% dplyr::select(Length) %>% set_rownames(ple_geneid)


# select columns of the counts, excluding chr start end etc info at beginning
con_counts<-con_counts[,7:23]
row.names(con_counts)<-con_geneid
ple_counts<-ple_counts[,7:23]
row.names(ple_counts)<- ple_geneid

# assign groups to the columns in con and ple
# missing letters are for Con root 2 and Ple petiole 3
con_groups<-c('A','A','A','B','B','B','C','C','C','D','D','D','E','E','F','F','F')
ple_groups<-c('G','G','G','H','H','H','I','I','I','J','J','K','K','K','L','L','L')

# import file that has transcript IDs and the GO terms mapping to them to be used after topTags later
con_go<-read_tsv("/Users/katie/Desktop/Bg/begonia_duplicate_expression/con_go_annotations_trans.txt")
ple_go<-read_tsv("/Users/katie/Desktop/Bg/begonia_duplicate_expression/ple_go_annotations_trans.txt")
colnames(con_go)<-c("transcript_id", "GO_terms")
colnames(ple_go)<-c("transcript_id", "GO_terms")


# create DGE list
# keep only rows where counts per million is above 100 in at least 2 samples
# calculate normalisation factors
# estimate tagwise and common dispersal
con_d<-DGEList(counts=con_counts,group=factor(con_groups))
con_d.full<-con_d
con_d.keep<-rowSums(cpm(con_d)>5)>=2
con_d<-con_d[con_d.keep,]
con_d_lengths<-con_lengths[con_d.keep,] %>% data.frame()
con_d$samples$lib.size<-colSums(con_d$counts)
con_d<-calcNormFactors(con_d, method="TMM")
plotMDS(con_d, method="bcv", col=as.numeric(con_d$samples$group))
con_d<-estimateCommonDisp(con_d)
con_d<-estimateTagwiseDisp(con_d, prior.n = 4)

# CON tissues
#A = female flower
#B = leaf
#C = male flower
#D = petiole
#E = root
#F = vegetative bud

# PLE tissues
#G = female flower 
#H = leaf = coef
#I = male flower 
#J = petiole 
#K = root
#L = vegetative bud 

ple_d<-DGEList(counts= ple_counts,group=factor(ple_groups))
ple_d.full<-ple_d
ple_d.keep<-rowSums(cpm(ple_d)>5)>=2
ple_d<-ple_d[ple_d.keep,]
ple_d_lengths<-ple_lengths[ple_d.keep,] %>% data.frame()
rownames(ple_d_lengths)<-rownames(ple_d$counts)
ple_d$samples$lib.size<-colSums(ple_d$counts)
ple_d<-calcNormFactors(ple_d, method="TMM")
plotMDS(ple_d, method="bcv", col=as.numeric(ple_d$samples$group))
ple_d<-estimateCommonDisp(ple_d)
ple_d<-estimateTagwiseDisp(ple_d, prior.n = 4)
design<-model.matrix(~ple_groups)


tissue_reference<-list("A" = "female_flower",
                       "B" = "leaf",
                       "C" = "male_flower",
                       "D" = "petiole",
                       "E" = "root",
                       "F" = "vegetative_bud",
                       "G" = "female_flower",
                       "H" = "leaf",
                       "I" = "male_flower",
                       "J" = "petiole",
                       "K" = "root",
                       "L" = "vegetative_bud")
de_con<-c()
de_ple<-c()
write_de_genes<-function(species, dge_object, pair1, pair2, annotation){
  tissue1 = tissue_reference[[pair1]]
  tissue2 = tissue_reference[[pair2]]
  exact_test<-exactTest(dge_object, pair=c(pair1, pair2))
  top_tags<-exact_test %>% topTags(n=10000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(annotation, by="transcript_id")
  comparison_all<-top_tags %>% filter(FDR < 0.05 & abs(logFC) > 2)
  comparison_upregulated<- top_tags %>% filter(FDR < 0.05 & logFC > 2) %>% dplyr::select(transcript_id)
  comparison_downregulated<-top_tags %>% filter(FDR < 0.05 & logFC < -2) %>% dplyr::select(transcript_id)
  up_outfile<-paste(species, tissue1, "v", tissue2, tissue2, "up", sep="_")
  down_outfile<-paste(species, tissue1, "v", tissue2, tissue2, "down", sep="_")
  all_oufile<-paste(species, tissue1, "v", tissue2, tissue2, "all", sep="_")
  write.table(comparison_upregulated, file=up_outfile, quote=FALSE, row.names = FALSE, col.names = FALSE)
  write.table(comparison_downregulated, file=down_outfile, quote=FALSE, row.names = FALSE, col.names = FALSE)
  write.table(comparison_all, file=all_oufile, quote=FALSE, row.names = FALSE, col.names = FALSE)
  if (species == "con"){
    de_con<<-c(de_con, comparison_all$transcript_id)
  }
  if (species == "ple"){
    de_ple<<-c(de_ple, comparison_all$transcript_id)
  }
    
}
  
  
write_de_genes("con", con_d, "A", "B", con_go)
write_de_genes("con", con_d, "A", "C", con_go)
write_de_genes("con", con_d, "A", "D", con_go)
write_de_genes("con", con_d, "A", "E", con_go)
write_de_genes("con", con_d, "A", "F", con_go)
write_de_genes("con", con_d, "B", "C", con_go)
write_de_genes("con", con_d, "B", "D", con_go)
write_de_genes("con", con_d, "B", "E", con_go)
write_de_genes("con", con_d, "B", "F", con_go)
write_de_genes("con", con_d, "C", "D", con_go)
write_de_genes("con", con_d, "C", "E", con_go)
write_de_genes("con", con_d, "C", "F", con_go)
write_de_genes("con", con_d, "D", "E", con_go)
write_de_genes("con", con_d, "D", "F", con_go)
write_de_genes("con", con_d, "E", "F", con_go)
  
  
  write_de_genes("ple", ple_d, "G", "H", ple_go)
  write_de_genes("ple", ple_d, "G", "I", ple_go)
  write_de_genes("ple", ple_d, "G", "J", ple_go)
  write_de_genes("ple", ple_d, "G", "K", ple_go)
  write_de_genes("ple", ple_d, "G", "L", ple_go)
  write_de_genes("ple", ple_d, "H", "I", ple_go)
  write_de_genes("ple", ple_d, "H", "J", ple_go)
  write_de_genes("ple", ple_d, "H", "K", ple_go)
  write_de_genes("ple", ple_d, "H", "L", ple_go)
  write_de_genes("ple", ple_d, "I", "J", ple_go)
  write_de_genes("ple", ple_d, "I", "K", ple_go)
  write_de_genes("ple", ple_d, "I", "L", ple_go)
  write_de_genes("ple", ple_d, "J", "K", ple_go)
  write_de_genes("ple", ple_d, "J", "L", ple_go)
  write_de_genes("ple", ple_d, "K", "L", ple_go)
  






########################################################
#    Venn of GO terms shared and unique of DE genes    #
########################################################

con_transcript2go<-strsplit(con_go$GO_terms, ",")
names(con_transcript2go)<-con_go$transcript_id

ple_transcript2go<-strsplit(ple_go$GO_terms, ",")
names(ple_transcript2go)<-ple_go$transcript_id
  
con_de_go<-sapply(de_con, function(x) con_transcript2go[[x]])
ple_de_go<-sapply(de_ple, function(x) ple_transcript2go[[x]])

a<-unlist(con_de_go) %>% as.character()
b<-unlist(ple_de_go) %>% as.character()

intersect(a,b)
# elements in con but not in ple
setdiff(a, b)
# elements in ple but not in con
setdiff(b, a)

length(unique(c(a,b)))

myCol<-c("mistyrose1", "lightskyblue1")
v<-venn.diagram(list(cat_a=a, cat_b=b), 
                filename = "test_venn.png", 
                category.names = c("B. conchifolia" , "B. plebeja"),
                lwd = 1,
                lty = 'blank',
                fill = myCol,
                cex = .6,
                fontface = "bold",
                fontfamily = "sans",
                
                imagetype="png" ,
                height = 600 , 
                width = 600 , 
                resolution = 300,
                compression = "lzw",
                
                cat.cex = 0.5,
                cat.fontface = "bold",
                cat.default.pos = "outer",
                cat.pos = c(-30, 32),
                cat.dist = c(0.05, 0.05),
                cat.fontfamily = "sans",
                )


####################################
#    DE gene circos plots    #
####################################

#A = female flower
#B = leaf
#C = male flower
#D = petiole
#E = root
#F = vegetative bud


from = c("Female flower",
"Female flower",
"Female flower",
"Female flower",
"Female flower",
"Leaf",
"Leaf",
"Leaf",
"Leaf",
"Male flower",
"Male flower",
"Male flower",
"Petiole",
"Petiole",
"Root")

to<-c("Leaf",
"Male flower",
"Petiole",
"Root",
"Vegetative bud",
"Male flower",
"Petiole",
"Root",
"Vegetative bud",
"Petiole",
"Root",
"Vegetative bud",
"Root",
"Vegetative bud",
"Vegetative bud")


con_de_count<-c(1241, 1444, 1015, 1778, 1226, 2366, 1030, 2018, 1182, 1861, 2766, 2425, 1826, 1221, 1410)
ple_de_count<-c(1455, 1248, 1205, 2188, 1312, 2592, 1191, 2634, 1579, 1976, 2659, 2235, 1953, 948, 1667)

sum(con_de_count)
sum(ple_de_count)


con_chord_df<-data.frame(from=from, to=to, de_count=con_de_count)
ple_chord_df<-data.frame(from=from, to=to, de_count=ple_de_count)

circos_cols<-c("hotpink",
               "green4",
               "orange",
               "darkseagreen2",
               "tan4",
               "skyblue1")

pdf('con_circos_plot.pdf')
par(cex = 2, mar = c(0, 0, 0, 0))
chordDiagram(con_chord_df, grid.col = circos_cols)
dev.off()
circos.clear()
pdf('ple_circos_plot.pdf')
chordDiagram(ple_chord_df, grid.col = circos_cols)
dev.off()




####################################
#    GO enrichment circos plots    #
####################################

con_female_flower_v_leaf_go_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_female_flower_v_leaf_leaf_all.GOseq.enriched")
con_female_flower_v_male_flower_go_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_female_flower_v_male_flower_male_flower_all.GOseq.enriched")
con_female_flower_v_petiole_go_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_female_flower_v_petiole_petiole_all.GOseq.enriched")
con_female_flower_v_root_go_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_female_flower_v_root_root_all.GOseq.enriched")
con_female_flower_v_vegetative_bud_go_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_female_flower_v_vegetative_bud_vegetative_bud_all.GOseq.enriched")
con_leaf_v_male_flower_go_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_leaf_v_male_flower_male_flower_all.GOseq.enriched")
con_leaf_v_petiole_go_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_leaf_v_petiole_petiole_all.GOseq.enriched")
con_leaf_v_root_go_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_leaf_v_root_root_all.GOseq.enriched")
con_leaf_v_vegetative_bud_go_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_leaf_v_vegetative_bud_vegetative_bud_all.GOseq.enriched")
con_male_flower_v_petiole_go_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_male_flower_v_petiole_petiole_all.GOseq.enriched")
con_male_flower_v_root_go_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_male_flower_v_root_root_all.GOseq.enriched")
con_male_flower_v_vegetative_bud_go_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_male_flower_v_vegetative_bud_vegetative_bud_all.GOseq.enriched")
con_petiole_v_root_go_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_petiole_v_root_root_all.GOseq.enriched")
con_petiole_v_vegetative_bud_go_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_petiole_v_vegetative_bud_vegetative_bud_all.GOseq.enriched")
con_root_v_vegetative_bud_go_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_root_v_vegetative_bud_vegetative_bud_all.GOseq.enriched")

ple_female_flower_v_leaf_go_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_female_flower_v_leaf_leaf_all.GOseq.enriched")
ple_female_flower_v_male_flower_go_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_female_flower_v_male_flower_male_flower_all.GOseq.enriched")
ple_female_flower_v_petiole_go_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_female_flower_v_petiole_petiole_all.GOseq.enriched")
ple_female_flower_v_root_go_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_female_flower_v_root_root_all.GOseq.enriched")
ple_female_flower_v_vegetative_bud_go_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_female_flower_v_vegetative_bud_vegetative_bud_all.GOseq.enriched")
ple_leaf_v_male_flower_go_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_leaf_v_male_flower_male_flower_all.GOseq.enriched")
ple_leaf_v_petiole_go_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_leaf_v_petiole_petiole_all.GOseq.enriched")
ple_leaf_v_root_go_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_leaf_v_root_root_all.GOseq.enriched")
ple_leaf_v_vegetative_bud_go_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_leaf_v_vegetative_bud_vegetative_bud_all.GOseq.enriched")
ple_male_flower_v_petiole_go_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_male_flower_v_petiole_petiole_all.GOseq.enriched")
ple_male_flower_v_root_go_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_male_flower_v_root_root_all.GOseq.enriched")
ple_male_flower_v_vegetative_bud_go_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_male_flower_v_vegetative_bud_vegetative_bud_all.GOseq.enriched")
ple_petiole_v_root_go_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_petiole_v_root_root_all.GOseq.enriched")
ple_petiole_v_vegetative_bud_go_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_petiole_v_vegetative_bud_vegetative_bud_all.GOseq.enriched")
ple_root_v_vegetative_bud_go_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_root_v_vegetative_bud_vegetative_bud_all.GOseq.enriched")

con_enriched_list<-list(ple_female_flower_v_leaf_go_enriched,
ple_female_flower_v_male_flower_go_enriched,
ple_female_flower_v_petiole_go_enriched,
ple_female_flower_v_root_go_enriched,
ple_female_flower_v_vegetative_bud_go_enriched,
ple_leaf_v_male_flower_go_enriched,
ple_leaf_v_petiole_go_enriched,
ple_leaf_v_root_go_enriched,
ple_leaf_v_vegetative_bud_go_enriched,
ple_male_flower_v_petiole_go_enriched,
ple_male_flower_v_root_go_enriched,
ple_male_flower_v_vegetative_bud_go_enriched,
ple_petiole_v_root_go_enriched,
ple_petiole_v_vegetative_bud_go_enriched,
ple_root_v_vegetative_bud_go_enriched)

ple_enriched_list<-list(con_female_flower_v_leaf_go_enriched,
con_female_flower_v_male_flower_go_enriched,
con_female_flower_v_petiole_go_enriched,
con_female_flower_v_root_go_enriched,
con_female_flower_v_vegetative_bud_go_enriched,
con_leaf_v_male_flower_go_enriched,
con_leaf_v_petiole_go_enriched,
con_leaf_v_root_go_enriched,
con_leaf_v_vegetative_bud_go_enriched,
con_male_flower_v_petiole_go_enriched,
con_male_flower_v_root_go_enriched,
con_male_flower_v_vegetative_bud_go_enriched,
con_petiole_v_root_go_enriched,
con_petiole_v_vegetative_bud_go_enriched,
con_root_v_vegetative_bud_go_enriched)


con_enriched_count<-sapply(con_enriched_list, function(x) x %>% filter(over_represented_FDR < 0.05) %>% dplyr::select(term) %>% data.frame() %>% nrow())
ple_enriched_count<-sapply(ple_enriched_list, function(x) x %>% filter(over_represented_FDR < 0.05) %>% dplyr::select(term) %>% data.frame() %>% nrow())



con_chord_df<-data.frame(from=from, to=to, de_count=con_enriched_count)
ple_chord_df<-data.frame(from=from, to=to, de_count=ple_enriched_count)


chordDiagram(con_chord_df, grid.col = circos_cols)
chordDiagram(ple_chord_df, grid.col = circos_cols)
circos.clear()


###############################################
#   GO enrichment up and down DE barplots     #
###############################################



con_female_flower_v_leaf_leaf_up_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_female_flower_v_leaf_leaf_up.GOseq.enriched")
con_female_flower_v_male_flower_male_flower_up_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_female_flower_v_male_flower_male_flower_up.GOseq.enriched")
con_female_flower_v_petiole_petiole_up_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_female_flower_v_petiole_petiole_up.GOseq.enriched")
con_female_flower_v_root_root_up_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_female_flower_v_root_root_up.GOseq.enriched")
con_female_flower_v_vegetative_bud_vegetative_bud_up_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_female_flower_v_vegetative_bud_vegetative_bud_up.GOseq.enriched")
con_leaf_v_male_flower_male_flower_up_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_leaf_v_male_flower_male_flower_up.GOseq.enriched")
con_leaf_v_petiole_petiole_up_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_leaf_v_petiole_petiole_up.GOseq.enriched")
con_leaf_v_root_root_up_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_leaf_v_root_root_up.GOseq.enriched")
con_leaf_v_vegetative_bud_vegetative_bud_up_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_leaf_v_vegetative_bud_vegetative_bud_up.GOseq.enriched")
con_male_flower_v_petiole_petiole_up_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_male_flower_v_petiole_petiole_up.GOseq.enriched")
con_male_flower_v_root_root_up_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_male_flower_v_root_root_up.GOseq.enriched")
con_male_flower_v_vegetative_bud_vegetative_bud_up_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_male_flower_v_vegetative_bud_vegetative_bud_up.GOseq.enriched")
con_petiole_v_root_root_up_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_petiole_v_root_root_up.GOseq.enriched")
con_petiole_v_vegetative_bud_vegetative_bud_up_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_petiole_v_vegetative_bud_vegetative_bud_up.GOseq.enriched")
con_root_v_vegetative_bud_vegetative_bud_up_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_root_v_vegetative_bud_vegetative_bud_up.GOseq.enriched")

con_female_flower_v_leaf_leaf_down_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_female_flower_v_leaf_leaf_down.GOseq.enriched")
con_female_flower_v_male_flower_male_flower_down_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_female_flower_v_male_flower_male_flower_down.GOseq.enriched")
con_female_flower_v_petiole_petiole_down_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_female_flower_v_petiole_petiole_down.GOseq.enriched")
con_female_flower_v_root_root_down_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_female_flower_v_root_root_down.GOseq.enriched")
con_female_flower_v_vegetative_bud_vegetative_bud_down_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_female_flower_v_vegetative_bud_vegetative_bud_down.GOseq.enriched")
con_leaf_v_male_flower_male_flower_down_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_leaf_v_male_flower_male_flower_down.GOseq.enriched")
con_leaf_v_petiole_petiole_down_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_leaf_v_petiole_petiole_down.GOseq.enriched")
con_leaf_v_root_root_down_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_leaf_v_root_root_down.GOseq.enriched")
con_leaf_v_vegetative_bud_vegetative_bud_down_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_leaf_v_vegetative_bud_vegetative_bud_down.GOseq.enriched")
con_male_flower_v_petiole_petiole_down_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_male_flower_v_petiole_petiole_down.GOseq.enriched")
con_male_flower_v_root_root_down_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_male_flower_v_root_root_down.GOseq.enriched")
con_male_flower_v_vegetative_bud_vegetative_bud_down_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_male_flower_v_vegetative_bud_vegetative_bud_down.GOseq.enriched")
con_petiole_v_root_root_down_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_petiole_v_root_root_down.GOseq.enriched")
con_petiole_v_vegetative_bud_vegetative_bud_down_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_petiole_v_vegetative_bud_vegetative_bud_down.GOseq.enriched")
con_root_v_vegetative_bud_vegetative_bud_down_enriched <- read_tsv("go_sidb/con/enrichment_input_con/con_root_v_vegetative_bud_vegetative_bud_down.GOseq.enriched")

ple_female_flower_v_leaf_leaf_up_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_female_flower_v_leaf_leaf_up.GOseq.enriched")
ple_female_flower_v_male_flower_male_flower_up_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_female_flower_v_male_flower_male_flower_up.GOseq.enriched")
ple_female_flower_v_petiole_petiole_up_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_female_flower_v_petiole_petiole_up.GOseq.enriched")
ple_female_flower_v_root_root_up_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_female_flower_v_root_root_up.GOseq.enriched")
ple_female_flower_v_vegetative_bud_vegetative_bud_up_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_female_flower_v_vegetative_bud_vegetative_bud_up.GOseq.enriched")
ple_leaf_v_male_flower_male_flower_up_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_leaf_v_male_flower_male_flower_up.GOseq.enriched")
ple_leaf_v_petiole_petiole_up_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_leaf_v_petiole_petiole_up.GOseq.enriched")
ple_leaf_v_root_root_up_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_leaf_v_root_root_up.GOseq.enriched")
ple_leaf_v_vegetative_bud_vegetative_bud_up_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_leaf_v_vegetative_bud_vegetative_bud_up.GOseq.enriched")
ple_male_flower_v_petiole_petiole_up_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_male_flower_v_petiole_petiole_up.GOseq.enriched")
ple_male_flower_v_root_root_up_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_male_flower_v_root_root_up.GOseq.enriched")
ple_male_flower_v_vegetative_bud_vegetative_bud_up_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_male_flower_v_vegetative_bud_vegetative_bud_up.GOseq.enriched")
ple_petiole_v_root_root_up_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_petiole_v_root_root_up.GOseq.enriched")
ple_petiole_v_vegetative_bud_vegetative_bud_up_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_petiole_v_vegetative_bud_vegetative_bud_up.GOseq.enriched")
ple_root_v_vegetative_bud_vegetative_bud_up_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_root_v_vegetative_bud_vegetative_bud_up.GOseq.enriched")

ple_female_flower_v_leaf_leaf_down_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_female_flower_v_leaf_leaf_down.GOseq.enriched")
ple_female_flower_v_male_flower_male_flower_down_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_female_flower_v_male_flower_male_flower_down.GOseq.enriched")
ple_female_flower_v_petiole_petiole_down_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_female_flower_v_petiole_petiole_down.GOseq.enriched")
ple_female_flower_v_root_root_down_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_female_flower_v_root_root_down.GOseq.enriched")
ple_female_flower_v_vegetative_bud_vegetative_bud_down_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_female_flower_v_vegetative_bud_vegetative_bud_down.GOseq.enriched")
ple_leaf_v_male_flower_male_flower_down_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_leaf_v_male_flower_male_flower_down.GOseq.enriched")
ple_leaf_v_petiole_petiole_down_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_leaf_v_petiole_petiole_down.GOseq.enriched")
ple_leaf_v_root_root_down_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_leaf_v_root_root_down.GOseq.enriched")
ple_leaf_v_vegetative_bud_vegetative_bud_down_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_leaf_v_vegetative_bud_vegetative_bud_down.GOseq.enriched")
ple_male_flower_v_petiole_petiole_down_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_male_flower_v_petiole_petiole_down.GOseq.enriched")
ple_male_flower_v_root_root_down_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_male_flower_v_root_root_down.GOseq.enriched")
ple_male_flower_v_vegetative_bud_vegetative_bud_down_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_male_flower_v_vegetative_bud_vegetative_bud_down.GOseq.enriched")
ple_petiole_v_root_root_down_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_petiole_v_root_root_down.GOseq.enriched")
ple_petiole_v_vegetative_bud_vegetative_bud_down_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_petiole_v_vegetative_bud_vegetative_bud_down.GOseq.enriched")
ple_root_v_vegetative_bud_vegetative_bud_down_enriched <- read_tsv("go_sidb/ple/enrichment_input_ple/ple_root_v_vegetative_bud_vegetative_bud_down.GOseq.enriched")




con_up_enriched_list<-list(con_female_flower_v_leaf_leaf_up_enriched,
con_female_flower_v_male_flower_male_flower_up_enriched,
con_female_flower_v_petiole_petiole_up_enriched,
con_female_flower_v_root_root_up_enriched,
con_female_flower_v_vegetative_bud_vegetative_bud_up_enriched,
con_leaf_v_male_flower_male_flower_up_enriched,
con_leaf_v_petiole_petiole_up_enriched,
con_leaf_v_root_root_up_enriched,
con_leaf_v_vegetative_bud_vegetative_bud_up_enriched,
con_male_flower_v_petiole_petiole_up_enriched,
con_male_flower_v_root_root_up_enriched,
con_male_flower_v_vegetative_bud_vegetative_bud_up_enriched,
con_petiole_v_root_root_up_enriched,
con_petiole_v_vegetative_bud_vegetative_bud_up_enriched,
con_root_v_vegetative_bud_vegetative_bud_up_enriched)

con_down_enriched_list<-list(con_female_flower_v_leaf_leaf_down_enriched,
con_female_flower_v_male_flower_male_flower_down_enriched,
con_female_flower_v_petiole_petiole_down_enriched,
con_female_flower_v_root_root_down_enriched,
con_female_flower_v_vegetative_bud_vegetative_bud_down_enriched,
con_leaf_v_male_flower_male_flower_down_enriched,
con_leaf_v_petiole_petiole_down_enriched,
con_leaf_v_root_root_down_enriched,
con_leaf_v_vegetative_bud_vegetative_bud_down_enriched,
con_male_flower_v_petiole_petiole_down_enriched,
con_male_flower_v_root_root_down_enriched,
con_male_flower_v_vegetative_bud_vegetative_bud_down_enriched,
con_petiole_v_root_root_down_enriched,
con_petiole_v_vegetative_bud_vegetative_bud_down_enriched,
con_root_v_vegetative_bud_vegetative_bud_down_enriched)




ple_up_enriched_list<-list(ple_female_flower_v_leaf_leaf_up_enriched,
ple_female_flower_v_male_flower_male_flower_up_enriched,
ple_female_flower_v_petiole_petiole_up_enriched,
ple_female_flower_v_root_root_up_enriched,
ple_female_flower_v_vegetative_bud_vegetative_bud_up_enriched,
ple_leaf_v_male_flower_male_flower_up_enriched,
ple_leaf_v_petiole_petiole_up_enriched,
ple_leaf_v_root_root_up_enriched,
ple_leaf_v_vegetative_bud_vegetative_bud_up_enriched,
ple_male_flower_v_petiole_petiole_up_enriched,
ple_male_flower_v_root_root_up_enriched,
ple_male_flower_v_vegetative_bud_vegetative_bud_up_enriched,
ple_petiole_v_root_root_up_enriched,
ple_petiole_v_vegetative_bud_vegetative_bud_up_enriched,
ple_root_v_vegetative_bud_vegetative_bud_up_enriched)

ple_down_enriched_list<-list(ple_female_flower_v_leaf_leaf_down_enriched,
ple_female_flower_v_male_flower_male_flower_down_enriched,
ple_female_flower_v_petiole_petiole_down_enriched,
ple_female_flower_v_root_root_down_enriched,
ple_female_flower_v_vegetative_bud_vegetative_bud_down_enriched,
ple_leaf_v_male_flower_male_flower_down_enriched,
ple_leaf_v_petiole_petiole_down_enriched,
ple_leaf_v_root_root_down_enriched,
ple_leaf_v_vegetative_bud_vegetative_bud_down_enriched,
ple_male_flower_v_petiole_petiole_down_enriched,
ple_male_flower_v_root_root_down_enriched,
ple_male_flower_v_vegetative_bud_vegetative_bud_down_enriched,
ple_petiole_v_root_root_down_enriched,
ple_petiole_v_vegetative_bud_vegetative_bud_down_enriched,
ple_root_v_vegetative_bud_vegetative_bud_down_enriched)



con_enriched_up_count<-sapply(con_up_enriched_list, function(x) x %>% filter(over_represented_FDR < 0.05) %>% dplyr::select(term) %>% data.frame() %>% nrow())
con_enriched_down_count<-sapply(con_down_enriched_list, function(x) x %>% filter(over_represented_FDR < 0.05) %>% dplyr::select(term) %>% data.frame() %>% nrow())

ple_enriched_up_count<-sapply(ple_up_enriched_list, function(x) x %>% filter(over_represented_FDR < 0.05) %>% dplyr::select(term) %>% data.frame() %>% nrow())
ple_enriched_down_count<-sapply(ple_down_enriched_list, function(x) x %>% filter(over_represented_FDR < 0.05) %>% dplyr::select(term) %>% data.frame() %>% nrow())



comparisons<-c("fem_flower vs leaf",
"female flower vs male flower",
"female flower vs petiole",
"female flower vs root",
"female flower vs veg bud",
"leaf vs male flower",
"leaf vs petiole",
"leaf vs root",
"leaf vs veg bud",
"male flower vs petiole",
"male flower vs root",
"male flower vs veg bud",
"petiole vs root",
"petiole vs veg bud",
"root vs veg bud")



con_up_df<-data.frame(comparison=comparisons, 
           count=con_enriched_up_count, 
           species="B. conchifolia",
           direction = "up")

con_down_df<-data.frame(comparison=comparisons, 
           count=con_enriched_down_count, 
           species="B. conchifolia",
           direction = "down")


ple_up_df<-data.frame(comparison=comparisons, 
           count=ple_enriched_up_count, 
           species="B. plebeja",
           direction = "up")

ple_down_df<-data.frame(comparison=comparisons, 
           count=ple_enriched_down_count, 
           species="B. plebeja",
           direction = "down")


de_df<-rbind(con_up_df, con_down_df, ple_up_df, ple_down_df)

ggplot(de_df, aes(x=comparison, y=count, colour=species, fill=species)) + 
  geom_bar(stat="identity", position = "dodge", width=0.5) + facet_grid(rows = vars(direction)) + 
  coord_flip() +
  labs(x = "GO term", y = "Number DE in Category") +
  theme(legend.text=element_text(size=10),
                                 legend.title=element_text(size=15),
                                 axis.title.x=element_text(size=15), 
                                 axis.title.y=element_text(size=15),
                                 axis.text.x= element_text(size=10),
                                 axis.text.y= element_text(size=10))





brks <- seq(-120, 120, 20)
lbls = paste0(as.character(c(seq(120, 0, -20), seq(20, 120, 20))), "")
de_df_mirror<-de_df %>% 
  mutate(count_mirror=ifelse(de_df$direction == "down", -(de_df$count), de_df$count))
  
de_df_mirror<-de_df_mirror %>%
  mutate(colour_label = case_when(
    species == "B. conchifolia" & direction == 'up' ~ "B. conchifolia up",
    species == "B. conchifolia" & direction == 'down' ~ "B. conchifolia down",
    species == "B. plebeja" & direction == 'up' ~ "B. plebeja up",
    species == "B. plebeja" & direction == 'down' ~ "B. plebeja down"
  ))


ggplot(de_df_mirror, aes(x = comparison, y = count_mirror, fill = colour_label)) +   # Fill column
  geom_bar(stat = "identity", width = .6) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls) + 
  coord_flip() +
  scale_fill_manual(values=c("firebrick2", "lightpink2", "skyblue2", "dodgerblue2")) + 
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title.x=element_text(size=15), 
        axis.title.y=element_text(size=15),
        axis.text.x= element_text(size=10),
        axis.text.y= element_text(size=15)) +
  labs(x = "Comparison", y = "Number of enriched GO terms", fill='')

de_df_mirror %>% filter(comparison == "fem_flower vs leaf")



#################################################################
#   enriched GO that are unique to a species, per comparison   #
#################################################################


# For each tissue, get GO terms that are unique to each species (up AND down), get the GO term mappings, plot the top


con_up_enriched_list_go<-list(con_female_flower_v_leaf_leaf_up_enriched %>% dplyr::select(category) %>% pull(),
                           con_female_flower_v_male_flower_male_flower_up_enriched %>% dplyr::select(category) %>% pull(),
                           con_female_flower_v_petiole_petiole_up_enriched %>% dplyr::select(category) %>% pull(),
                           con_female_flower_v_root_root_up_enriched %>% dplyr::select(category) %>% pull(),
                           con_female_flower_v_vegetative_bud_vegetative_bud_up_enriched %>% dplyr::select(category) %>% pull(),
                           con_leaf_v_male_flower_male_flower_up_enriched %>% dplyr::select(category) %>% pull(),
                           con_leaf_v_petiole_petiole_up_enriched %>% dplyr::select(category) %>% pull(),
                           con_leaf_v_root_root_up_enriched %>% dplyr::select(category) %>% pull(),
                           con_leaf_v_vegetative_bud_vegetative_bud_up_enriched %>% dplyr::select(category) %>% pull(),
                           con_male_flower_v_petiole_petiole_up_enriched %>% dplyr::select(category) %>% pull(),
                           con_male_flower_v_root_root_up_enriched %>% dplyr::select(category) %>% pull(),
                           con_male_flower_v_vegetative_bud_vegetative_bud_up_enriched %>% dplyr::select(category) %>% pull(),
                           con_petiole_v_root_root_up_enriched %>% dplyr::select(category) %>% pull(),
                           con_petiole_v_vegetative_bud_vegetative_bud_up_enriched %>% dplyr::select(category) %>% pull(),
                           con_root_v_vegetative_bud_vegetative_bud_up_enriched %>% dplyr::select(category) %>% pull())

con_down_enriched_list_go<-list(con_female_flower_v_leaf_leaf_down_enriched %>% dplyr::select(category) %>% pull(),
                             con_female_flower_v_male_flower_male_flower_down_enriched %>% dplyr::select(category) %>% pull(),
                             con_female_flower_v_petiole_petiole_down_enriched %>% dplyr::select(category) %>% pull(),
                             con_female_flower_v_root_root_down_enriched %>% dplyr::select(category) %>% pull(),
                             con_female_flower_v_vegetative_bud_vegetative_bud_down_enriched %>% dplyr::select(category) %>% pull(),
                             con_leaf_v_male_flower_male_flower_down_enriched %>% dplyr::select(category) %>% pull(),
                             con_leaf_v_petiole_petiole_down_enriched %>% dplyr::select(category) %>% pull(),
                             con_leaf_v_root_root_down_enriched %>% dplyr::select(category) %>% pull(),
                             con_leaf_v_vegetative_bud_vegetative_bud_down_enriched %>% dplyr::select(category) %>% pull(),
                             con_male_flower_v_petiole_petiole_down_enriched %>% dplyr::select(category) %>% pull(),
                             con_male_flower_v_root_root_down_enriched %>% dplyr::select(category) %>% pull(),
                             con_male_flower_v_vegetative_bud_vegetative_bud_down_enriched %>% dplyr::select(category) %>% pull(),
                             con_petiole_v_root_root_down_enriched %>% dplyr::select(category) %>% pull(),
                             con_petiole_v_vegetative_bud_vegetative_bud_down_enriched %>% dplyr::select(category) %>% pull(),
                             con_root_v_vegetative_bud_vegetative_bud_down_enriched %>% dplyr::select(category) %>% pull())




ple_up_enriched_list_go<-list(ple_female_flower_v_leaf_leaf_up_enriched %>% dplyr::select(category) %>% pull(),
                           ple_female_flower_v_male_flower_male_flower_up_enriched %>% dplyr::select(category) %>% pull(),
                           ple_female_flower_v_petiole_petiole_up_enriched %>% dplyr::select(category) %>% pull(),
                           ple_female_flower_v_root_root_up_enriched %>% dplyr::select(category) %>% pull(),
                           ple_female_flower_v_vegetative_bud_vegetative_bud_up_enriched %>% dplyr::select(category) %>% pull(),
                           ple_leaf_v_male_flower_male_flower_up_enriched %>% dplyr::select(category) %>% pull(),
                           ple_leaf_v_petiole_petiole_up_enriched %>% dplyr::select(category) %>% pull(),
                           ple_leaf_v_root_root_up_enriched %>% dplyr::select(category) %>% pull(),
                           ple_leaf_v_vegetative_bud_vegetative_bud_up_enriched %>% dplyr::select(category) %>% pull(),
                           ple_male_flower_v_petiole_petiole_up_enriched %>% dplyr::select(category) %>% pull(),
                           ple_male_flower_v_root_root_up_enriched %>% dplyr::select(category) %>% pull(),
                           ple_male_flower_v_vegetative_bud_vegetative_bud_up_enriched %>% dplyr::select(category) %>% pull(),
                           ple_petiole_v_root_root_up_enriched %>% dplyr::select(category) %>% pull(),
                           ple_petiole_v_vegetative_bud_vegetative_bud_up_enriched %>% dplyr::select(category) %>% pull(),
                           ple_root_v_vegetative_bud_vegetative_bud_up_enriched %>% dplyr::select(category) %>% pull())

ple_down_enriched_list_go<-list(ple_female_flower_v_leaf_leaf_down_enriched %>% dplyr::select(category) %>% pull(),
                             ple_female_flower_v_male_flower_male_flower_down_enriched %>% dplyr::select(category) %>% pull(),
                             ple_female_flower_v_petiole_petiole_down_enriched %>% dplyr::select(category) %>% pull(),
                             ple_female_flower_v_root_root_down_enriched %>% dplyr::select(category) %>% pull(),
                             ple_female_flower_v_vegetative_bud_vegetative_bud_down_enriched %>% dplyr::select(category) %>% pull(),
                             ple_leaf_v_male_flower_male_flower_down_enriched %>% dplyr::select(category) %>% pull(),
                             ple_leaf_v_petiole_petiole_down_enriched %>% dplyr::select(category) %>% pull(),
                             ple_leaf_v_root_root_down_enriched %>% dplyr::select(category) %>% pull(),
                             ple_leaf_v_vegetative_bud_vegetative_bud_down_enriched %>% dplyr::select(category) %>% pull(),
                             ple_male_flower_v_petiole_petiole_down_enriched %>% dplyr::select(category) %>% pull(),
                             ple_male_flower_v_root_root_down_enriched %>% dplyr::select(category) %>% pull(),
                             ple_male_flower_v_vegetative_bud_vegetative_bud_down_enriched %>% dplyr::select(category) %>% pull(),
                             ple_petiole_v_root_root_down_enriched %>% dplyr::select(category) %>% pull(),
                             ple_petiole_v_vegetative_bud_vegetative_bud_down_enriched %>% dplyr::select(category) %>% pull(),
                             ple_root_v_vegetative_bud_vegetative_bud_down_enriched %>% dplyr::select(category) %>% pull())



# for each tissue comparison, get the unique GO enriched terms between con and ple
# also shared ones
in_con_notin_ple_down<-mapply(function(one, two) unique(setdiff(one, two)), con_down_enriched_list_go, ple_down_enriched_list_go)
in_ple_notin_con_down<-mapply(function(one, two) unique(setdiff(two, one)), con_down_enriched_list_go, ple_down_enriched_list_go)

in_con_notin_ple_up<-mapply(function(one, two) unique(setdiff(one, two)), con_up_enriched_list_go, ple_up_enriched_list_go)
in_ple_notin_con_up<-mapply(function(one, two) unique(setdiff(two, one)), con_up_enriched_list_go, ple_up_enriched_list_go)

in_ple_and_con_down<-mapply(function(one, two) intersect(one, two), con_down_enriched_list_go, ple_down_enriched_list_go)
in_ple_and_con_up<-mapply(function(one, two) intersect(one, two), con_up_enriched_list_go, ple_up_enriched_list_go)


comparison_names<-list("Female flower vs Leaf",
                    "Female flower vs Male flower",
                    "Female flower vs Petiole",
                    "Female flower vs Root",
                    "Female flower vs Vegetative bud",
                    "Leaf vs Male flower",
                    "Leaf vs Petiole",
                    "Leaf vs Root",
                    "Leaf vs Vegetative bud",
                    "Male flower vs Petiole",
                    "Male flower vs Root",
                    "Male flower vs Vegetative bud",
                    "Petiole vs Root",
                    "Petiole vs Vegetative bud",
                    "Root vs Vegetative bud")


# conch specific Go terms enriched in con but not in ple per tissue comparison (upregulaed reference tissue)
con_up<-mapply(function(one, two) two %>% dplyr::filter(category %in% one & over_represented_FDR < 0.005) %>% dplyr::select(category, numDEInCat, ontology, go_term), in_con_notin_ple_up, con_up_enriched_list, SIMPLIFY=FALSE)
con_up_labelled<-mapply(function(names, tibbles) tibbles %>% mutate(comparison=names, species="B. conchifolia", direction="up"), comparison_names, con_up, SIMPLIFY = FALSE)

# conch specific Go terms enriched in con but not in ple per tissue coparison (downregulated reference tissue)
con_down<-mapply(function(one, two) two %>% dplyr::filter(category %in% one & over_represented_FDR < 0.005) %>% dplyr::select(category, numDEInCat, ontology, go_term), in_con_notin_ple_down, con_down_enriched_list, SIMPLIFY=FALSE)
con_down_labelled<-mapply(function(names, tibbles) tibbles %>% mutate(comparison=names, species="B. conchifolia", direction="down"), comparison_names, con_down, SIMPLIFY = FALSE)

# as above but for ple
ple_up<-mapply(function(one, two) two %>% dplyr::filter(category %in% one & over_represented_FDR < 0.005) %>% dplyr::select(category, numDEInCat, ontology, go_term), in_ple_notin_con_up, ple_up_enriched_list, SIMPLIFY=FALSE)
ple_up_labelled<-mapply(function(names, tibbles) tibbles %>% mutate(comparison=names, species="B. plebeja", direction="up"), comparison_names, ple_up, SIMPLIFY = FALSE)

ple_down<-mapply(function(one, two) two %>% dplyr::filter(category %in% one & over_represented_FDR < 0.005) %>% dplyr::select(category, numDEInCat, ontology, go_term), in_ple_notin_con_down, ple_down_enriched_list, SIMPLIFY=FALSE)
ple_down_labelled<-mapply(function(names, tibbles) tibbles %>% mutate(comparison=names, species="B. plebeja", direction="down"), comparison_names, ple_down, SIMPLIFY = FALSE)


# bind all the data frames together
# edit the go_term column to remove the ontology in front of go term definitions
all_up<-bind_rows(con_up_labelled, con_down_labelled, ple_up_labelled, ple_down_labelled) %>% filter(direction == "up") %>% dplyr::select(-direction) %>% data.frame()
all_down<-bind_rows(con_up_labelled, con_down_labelled, ple_up_labelled, ple_down_labelled) %>% filter(direction == "down") %>% dplyr::select(-direction) %>% data.frame()
all_up$go_term<-substring(all_up$go_term, 4)
all_down$go_term<-substring(all_down$go_term, 4)


all_con<-bind_rows(con_up_labelled, con_down_labelled, ple_up_labelled, ple_down_labelled) %>% filter(species == "B. conchifolia") %>% data.frame()
all_ple<-bind_rows(con_up_labelled, con_down_labelled, ple_up_labelled, ple_down_labelled) %>% filter(species == "B. plebeja") %>% data.frame()
all_con$go_term<-substring(all_con$go_term, 4)
all_ple$go_term<-substring(all_ple$go_term, 4)

# order the go terms by number DE genes in category for plotting
all_up$go_term<-factor(all_up$go_term, levels = unique(all_up$go_term[order(all_up$numDEInCat)]))
all_down$go_term<-factor(all_down$go_term, levels = unique(all_down$go_term[order(all_down$numDEInCat)]))

all_con$go_term<-factor(all_con$go_term, levels = unique(all_con$go_term[order(all_con$numDEInCat)]))
all_ple$go_term<-factor(all_ple$go_term, levels = unique(all_ple$go_term[order(all_ple$numDEInCat)]))


# do it by species, remove direction altogether as it doesnt add much and mnakes the graphs look weird
ggplot(all_con, aes(x=go_term, y=numDEInCat, fill=comparison)) + 
  labs(x = "GO term", y = "Number DE in Category") +
  geom_bar(stat="identity", width = 0.7) + coord_flip() + theme(legend.text=element_text(size=10),
                                                                      legend.title=element_text(size=15),
                                                                      axis.title.x=element_text(size=15), 
                                                                      axis.title.y=element_text(size=15),
                                                                      axis.text.x= element_text(size=10),
                                                                      axis.text.y= element_text(size=12)) +
  scale_fill_brewer(palette = "Paired")


ggplot(all_ple, aes(x=go_term, y=numDEInCat, fill=comparison)) + 
  geom_bar(stat="identity", position = "dodge") + coord_flip() +
  labs(x = "GO term", y = "Number DE in Category") +
  geom_bar(stat="identity") + coord_flip() + theme(legend.text=element_text(size=10),
                                                                       legend.title=element_text(size=15),
                                                                       axis.title.x=element_text(size=15), 
                                                                       axis.title.y=element_text(size=15),
                                                                       axis.text.x= element_text(size=10),
                                                                       axis.text.y= element_text(size=12)) +
  scale_fill_manual(values=c("chocolate1", 
                              "aquamarine3", 
                              "cornflowerblue", 
                              "firebrick2", 
                              "gold2",
                              "darkorchid1",
                              "darkolivegreen2",
                              "deeppink1",
                              "lightgoldenrod1",
                              "lightsteelblue3",
                              "maroon",
                              "pink1",
                              "springgreen2",
                              "tan1",
                              "wheat3"))



colours_ple<-c("chocolate1", 
               "aquamarine3", 
               "cornflowerblue", 
               "firebrick2", 
               "gold2",
               "darkorchid1",
               "darkolivegreen2",
               "deeppink1",
               "lightgoldenrod1",
               "lightsteelblue3",
               "maroon",
               "pink1",
               "springgreen2",
               "tan1",
               "wheat3")










pdf("con_comparison_network.pdf")
con_comparison_numbers<-all_con[,c("numDEInCat", "comparison")]
con_comparison_numbers<-aggregate(numDEInCat ~ comparison, con_comparison_numbers, sum)
con_comparison_numbers<-tidyr::separate(con_comparison_numbers, comparison, into = c("test1", "test2"), sep = " vs ",)
colnames(con_comparison_numbers)<-c("from", "to", "d")
m <- dodgr_dists(con_comparison_numbers)
m <- as.matrix(m)
qgraph(m,edge.labels=TRUE)
dev.off()

pdf("ple_comparison_network.pdf")
ple_comparison_numbers<-all_ple[,c("numDEInCat", "comparison")]
ple_comparison_numbers<-aggregate(numDEInCat ~ comparison, ple_comparison_numbers, sum)
ple_comparison_numbers<-tidyr::separate(ple_comparison_numbers, comparison, into = c("test1", "test2"), sep = " vs ",)
colnames(ple_comparison_numbers)<-c("from", "to", "d")
m <- dodgr_dists(ple_comparison_numbers)
m <- as.matrix(m)
qgraph(m,edge.labels=TRUE)
dev.off()




# check how many GO terms are unique to one tissue comparison, and before that print how mnay unqiue GO terms there are
length(unique(all_con$go_term))
length(which(sort(summary(all_con$go_term)) == 1))

length(unique(all_con$go_term))
length(which(sort(summary(all_ple$go_term)) == 1))


tissue_graph<-unique(c(test3$test1, test3$test2))


g4 <- graph(x)
ggnet2(g4)            
             
install.packages("dodgr")
library(dodgr)

aggregate(numDEInCat ~ row, test2, sum)


test$comparison

test<-all$go_term[1]

test[1]

substring(all$go_term, 4)


fflower_vs<-c(0, 611, 671, 415, 832, 689)
leaf_vs<-c(716, 0, 1071, 447, 871, 625)
mflower_vs<-c(658, 1245, 0, 873, 1391, 1307)
petiole_vs<-c(675, 573, 1005, 0, 950, 788)
root_vs<-c(1085, 1216, 1331, 902, 0, 761)
vegbud_vs<-c(799, 761, 1215, 556, 913, 0)


de_upreg<-data.frame(fflower_vs, leaf_vs, mflower_vs, petiole_vs, root_vs, vegbud_vs)


PLE_TRINITY_DN10424_c0_g2_i3<-c(146.510082987706,	2508.90811274799,	26.7097241343816,	591.039566176487,	8.86216943978959,	323.823009957749)
PLE_TRINITY_DN10702_c0_g1_i6<-c(410.783033582042,	2476.58077169722,	59.7230672964667,	530.094683714974,	7.41027437540702,	275.11157132182)
PLE_TRINITY_DN1234_c0_g1_i1<-	c(29.5028385950133,	150.642755197154,	10.675726223133,	40.7030102291845,	0.617789897782811,	23.9577421204015)
PLE_TRINITY_DN1394_c0_g1_i1<-	c(575.282545625639,	5267.25356802233,	118.018142652291,	1074.89751304855,	6.48322070175823,	571.211324511977)
PLE_TRINITY_DN1453_c0_g1_i1<-	c(21.6959964621876,	56.3566708106038,	3.20822054580883,	16.0957822120655,	3.24223505545978,	28.5339195986808)
PLE_TRINITY_DN2018_c0_g1_i2<-	c(38.7911994845226,	97.3584127178022,	4.89241384036555,	31.7182753919664,	4.53468162451518,	28.192437904522)
PLE_TRINITY_DN4425_c0_g1_i1<-	c(34.0463717444646,	68.6382835205116,	17.8969548663024,	18.1565841307301,	14.7051993076339,	18.6446103929537)
PLE_TRINITY_DN6516_c2_g2_i7<-	c(2013.00277213284,	9186.21879254671,	289.135916297716,	3302.34114531942,	12.1226282031063,	1870.76113260932)
PLE_TRINITY_DN9177_c0_g1_i4<-	c(19.7930155850893,	68.4810972773486,	14.2084412716347,	101.149191426205,	37.2026445832983,	89.8678653700418)
PLE_TRINITY_DN9399_c0_g1_i1<-	c(66.5089090693356,	188.760072835038,	34.065036427242,	65.7403447856097,	13.6504968936632,	43.3695622195009)
PLE_TRINITY_DN9399_c0_g3_i1<-	c(16.7637807607527,	97.6832412484189,	3.67351638476552,	14.6028223131999,	0.943263641264809,	13.0337048038031)
PLE_TRINITY_DN9990_c1_g1_i1<-	c(95.0061154105097,	255.6537302859,	12.8952088777428,	58.1396374637146,	0.850974253167925,	37.1973712371114)





PLE_TRINITY_DN10702_c0_g1_i6 <-	c(410.783033582042,	2476.58077169722,	59.7230672964667,	530.094683714974,	7.41027437540702,	275.11157132182)
PLE_TRINITY_DN1234_c0_g1_i1	<- c(29.5028385950133,	150.642755197154,	10.675726223133,	40.7030102291845,	0.617789897782811,	23.9577421204015)
PLE_TRINITY_DN1394_c0_g1_i1	<- c(575.282545625639,	5267.25356802233,	118.018142652291,	1074.89751304855,	6.48322070175823,	571.211324511977)
PLE_TRINITY_DN2018_c0_g1_i2	<- c(38.7911994845226,	97.3584127178022,	4.89241384036555,	31.7182753919664,	4.53468162451518,	28.192437904522)
PLE_TRINITY_DN4425_c0_g1_i1	<- c(34.0463717444646,	68.6382835205116,	17.8969548663024,	18.1565841307301,	14.7051993076339,	18.6446103929537)
PLE_TRINITY_DN6516_c2_g2_i7	<- c(2013.00277213284,	9186.21879254671,	289.135916297716,	3302.34114531942,	12.1226282031063,	1870.76113260932)
PLE_TRINITY_DN9177_c0_g1_i4	<- c(19.7930155850893,	68.4810972773486,	14.2084412716347,	101.149191426205,	37.2026445832983,	89.8678653700418)
PLE_TRINITY_DN9399_c0_g1_i1	<- c(66.5089090693356,	188.760072835038,	34.065036427242,	65.7403447856097,	13.6504968936632,	43.3695622195009)

CON_TRINITY_DN6879_c1_g3_i3	<- c(277.112050969756,	1126.41063031514,	41.7909253733399,	116.906805925568,	95.881496783914,	314.310403657318)
CON_TRINITY_DN1904_c0_g1_i1	<- c(49.0200288060174,	294.296919234734,	16.30566133939,	66.104384345018,	3.60434303861895,	12.8854498439226)
CON_TRINITY_DN6223_c0_g3_i3	<- c(508.403565573151,	3496.3256154533,	79.9175466309462,	491.33310325017,	100.870622457619,	481.592539485971)
CON_TRINITY_DN1734_c0_g1_i3	<- c(21.2640097686692,	49.9893625047817,	6.07949906495243,	21.6591378449562,	9.9789960937048,	19.5930395322612)
CON_TRINITY_DN862_c0_g1_i2	<- c(24.0562779199061,	37.0783275479418,	18.640343689413,	29.9518036363914,	19.0132345892191,	23.0084679056202)
CON_TRINITY_DN6158_c0_g2_i4	<- c(1050.1812024465,	5228.28860093757,	237.317257124915,	621.429047451046,	420.28481213993,	1445.78170043933)
CON_TRINITY_DN7064_c3_g3_i1	<- c(12.7578153940939,	8.32960667291691,	0.241362928891985,	7.58597185372785,	14.1326710658355,	12.8228679780768)
CON_TRINITY_DN7269_c1_g2_i2	<- c(42.7687704145285,	93.04670463628,	32.9993641146538,	83.0000073393078,	19.8664892138165,	40.8692874935894)



CON_TRINITY_DN6879_c1_g3_i3,
CON_TRINITY_DN1904_c0_g1_i1,
CON_TRINITY_DN6223_c0_g3_i3,
CON_TRINITY_DN1734_c0_g1_i3,
CON_TRINITY_DN862_c0_g1_i2,
CON_TRINITY_DN6158_c0_g2_i4,
CON_TRINITY_DN7064_c3_g3_i1,
CON_TRINITY_DN7269_c1_g2_i2



test_light_con<-data.frame(CON_TRINITY_DN6879_c1_g3_i3,
                  CON_TRINITY_DN1904_c0_g1_i1,
                  CON_TRINITY_DN6223_c0_g3_i3,
                  CON_TRINITY_DN1734_c0_g1_i3,
                  CON_TRINITY_DN862_c0_g1_i2,
                  CON_TRINITY_DN6158_c0_g2_i4,
                  CON_TRINITY_DN7064_c3_g3_i1,
                  CON_TRINITY_DN7269_c1_g2_i2)


test_light_ple<-data.frame(PLE_TRINITY_DN10702_c0_g1_i6,
                  PLE_TRINITY_DN1234_c0_g1_i1,
                  PLE_TRINITY_DN1394_c0_g1_i1,
                  PLE_TRINITY_DN2018_c0_g1_i2,
                  PLE_TRINITY_DN4425_c0_g1_i1,
                  PLE_TRINITY_DN6516_c2_g2_i7,
                  PLE_TRINITY_DN9177_c0_g1_i4,
                  PLE_TRINITY_DN9399_c0_g1_i1)






test_light<-data.frame(PLE_TRINITY_DN10702_c0_g1_i6,
                       CON_TRINITY_DN6879_c1_g3_i3,
                       PLE_TRINITY_DN1234_c0_g1_i1,
                       CON_TRINITY_DN1904_c0_g1_i1,
                       PLE_TRINITY_DN1394_c0_g1_i1,
                       CON_TRINITY_DN6223_c0_g3_i3,
                       PLE_TRINITY_DN2018_c0_g1_i2,
                       CON_TRINITY_DN1734_c0_g1_i3,
                       PLE_TRINITY_DN4425_c0_g1_i1,
                       CON_TRINITY_DN862_c0_g1_i2,
                       PLE_TRINITY_DN6516_c2_g2_i7,
                       CON_TRINITY_DN6158_c0_g2_i4,
                       PLE_TRINITY_DN9177_c0_g1_i4,
                       CON_TRINITY_DN7064_c3_g3_i1,
                       PLE_TRINITY_DN9399_c0_g1_i1,
                       CON_TRINITY_DN7269_c1_g2_i2)

test_light<-as.matrix(test_light)

test_light_con<-as.matrix(test_light_con)
test_light_ple<-as.matrix(test_light_ple)

pheatmap(log2(test_light), cluster_rows = TRUE, cluster_cols = TRUE)

pheatmap(log2(test_light_con), cluster_rows = TRUE, cluster_cols = TRUE)
pheatmap(log2(test_light_ple), cluster_rows = TRUE, cluster_cols = TRUE)




#########

CON_TRINITY_DN6879_c1_g3_i3<-c(277.112050969756, 1126.41063031514, 41.7909253733399, 116.906805925568, 95.881496783914, 314.310403657318)
CON_TRINITY_DN1904_c0_g1_i1<-c(49.0200288060174, 294.296919234734, 16.30566133939, 66.104384345018, 3.60434303861895, 12.8854498439226)
CON_TRINITY_DN6223_c0_g3_i3<-c(508.403565573151, 3496.3256154533, 79.9175466309462, 491.33310325017, 100.870622457619, 481.592539485971)
CON_TRINITY_DN862_c0_g1_i2 <-c(24.0562779199061, 37.0783275479418, 18.640343689413, 29.9518036363914, 19.0132345892191, 23.0084679056202)

PLE_TRINITY_DN10702_c0_g1_i6<-c(410.783033582042, 2476.58077169722, 59.7230672964667, 530.094683714974, 7.41027437540702, 275.11157132182)
PLE_TRINITY_DN1234_c0_g1_i1 <-c(29.5028385950133, 150.642755197154, 10.675726223133, 40.7030102291845, 0.617789897782811, 23.9577421204015)
PLE_TRINITY_DN1394_c0_g1_i1 <-c(575.282545625639, 5267.25356802233, 118.018142652291, 1074.89751304855, 6.48322070175823, 571.211324511977)
PLE_TRINITY_DN4425_c0_g1_i1 <-c(34.0463717444646, 68.6382835205116, 17.8969548663024, 18.1565841307301, 14.7051993076339, 18.6446103929537)





test_light2<-data.frame(CON_TRINITY_DN6879_c1_g3_i3,
                       CON_TRINITY_DN1904_c0_g1_i1,
                       CON_TRINITY_DN6223_c0_g3_i3,
                       CON_TRINITY_DN862_c0_g1_i2,
                       PLE_TRINITY_DN10702_c0_g1_i6,
                       PLE_TRINITY_DN1234_c0_g1_i1,
                       PLE_TRINITY_DN1394_c0_g1_i1,
                       PLE_TRINITY_DN4425_c0_g1_i1)
test_light2<-as.matrix(test_light2)

pheatmap(log2(test_light2), cluster_rows = TRUE, cluster_cols = TRUE)

########### root vs veg but response to chitin



CON_TRINITY_DN8704_c0_g1_i1<- c(0.375058801011086, 6.56683275048367, 1.5898260214002, 2.87972321927564, 0.663110220871657, 0.84456267290407)
CON_TRINITY_DN6566_c2_g2_i2<- c(5.70603277799682, 2.82956624570527, 9.29086592280381, 2.55425505328203, 8.39125201897979, 6.60118208667222)
CON_TRINITY_DN6428_c3_g6_i1<- c(0, 0, 0.0630941689296661, 0, 5.7783576991365, 4.19080630211211)
CON_TRINITY_DN6428_c3_g8_i1<- c(0, 1.50483775555447, 0.0521695365869513, 0.319569796927699, 2.95927410832759, 0.466872369331984)
CON_TRINITY_DN7052_c1_g1_i4<- c(53.4979541159031, 6.16272134893453, 326.504344758021, 10.8425086240611, 12.8226485591497, 6.05636450549996)
CON_TRINITY_DN8003_c0_g1_i1<- c(0, 0, 0, 0, 10.7328453352164, 0.136421110424267)
CON_TRINITY_DN4507_c0_g2_i1<- c(1.71788310405484, 0, 0.0929429215691094, 0.136320048302985, 6.51700808284003, 3.30771015869177)
CON_TRINITY_DN4507_c0_g1_i1<- c(4.76579625790262, 2.52254369989016, 3.65999587394966, 2.99750482475297, 19.2552739264758, 15.4999592917205)
CON_TRINITY_DN7520_c0_g6_i1<- c(24.4922996704888, 0.475525086302726, 15.8069397150784, 0.89402359732081, 0.328310524223138, 0.599666411817255)
CON_TRINITY_DN5026_c0_g2_i1<- c(1.85911571145037, 0.1910226741438, 0.827872527348919, 0.638738064490399, 4.80720122355383, 4.68212679613554)
CON_TRINITY_DN3461_c0_g1_i1<- c(0, 0, 0, 0.119078895482203, 12.1466848202506, 0)
CON_TRINITY_DN3461_c0_g2_i1<- c(0, 0, 0, 0, 30.6105400410935, 0.0826977162628845)
CON_TRINITY_DN6729_c0_g9_i1<- c(33.1579283135185, 107.204603595753, 21.610272911086, 62.5779131973003, 39.6314447570796, 113.328978336146)
CON_TRINITY_DN6993_c4_g1_i1<- c(3.25250919441379, 0.77959347478842, 11.0159476123005, 1.05641788968733, 14.7882419009452, 3.22892174213397)
CON_TRINITY_DN6882_c4_g4_i5<- c(5.35560812504056, 2.27491677339933, 3.29328870567352, 3.92635868536667, 11.0232769075161, 12.7642218157982)
CON_TRINITY_DN5957_c1_g2_i1<- c(3.67405917859885, 4.10408133204951, 9.30093618045327, 0.829197515061433, 10.0054073080056, 3.79750451190861)
CON_TRINITY_DN5307_c1_g2_i1<- c(0, 0.389688368792383, 1.03952615506951, 0.484401178443507, 3.89122119475069, 8.19586812227167)
PLE_TRINITY_DN1325_c0_g1_i3 <- c(0.659323384622814, 0.975400044696627, 2.88298269327096, 2.21088792589068, 2.08527663374336, 0.3035289682697)
PLE_TRINITY_DN1545_c0_g1_i1 <- c(1.49236746342907, 2.83376413696846, 4.90836152914577, 3.21371794386438, 0.568093675426948, 5.95890766378839)
PLE_TRINITY_DN1776_c0_g2_i1 <- c(0, 0, 0, 0, 13.7630436526205, 0)
PLE_TRINITY_DN1776_c0_g3_i2 <- c(0.299786783054826, 0.259186674888013, 0.616049168801522, 0.406397421908351, 7.09116234114136, 0)
PLE_TRINITY_DN4742_c0_g4_i1 <- c(6.87420737763873, 1.19862120433423, 18.5641080977738, 4.87168383371093, 29.9779705988592, 3.43091966477076)
PLE_TRINITY_DN5247_c0_g1_i2 <- c(0.085223083299511, 0.436437409833707, 0, 0.562132793875544, 22.430430041409, 0.360371695969352)
PLE_TRINITY_DN5341_c0_g3_i2 <- c(0.712456821630132, 0, 0.118077678276962, 0.381829525578476, 9.60224878919893, 0.462122895091024)
PLE_TRINITY_DN5341_c0_g4_i1 <- c(2.35804907281021, 0.369516382874314, 3.90270133662589, 4.70608771775706, 14.9247294934707, 3.0716684997354)
PLE_TRINITY_DN5411_c3_g2_i2 <- c(6.6152814132776, 0.970730644891229, 59.3716546329681, 0.876233429871122, 0, 1.81823589864504)
PLE_TRINITY_DN6010_c0_g1_i4 <- c(4.44486710756923, 0.866118357215631, 2.34420927487125, 0.820796137324788, 11.4850381231546, 0.646378297289276)
PLE_TRINITY_DN6573_c0_g2_i1 <- c(0.904464031287028, 0.968749760436948, 0.637527914057625, 1.19880490747864, 15.0950351844776, 1.55040491206715)
PLE_TRINITY_DN6573_c0_g3_i2 <- c(0, 0, 0, 0, 49.5314506255628, 0.378813608752038)
PLE_TRINITY_DN7191_c2_g1_i1 <- c(47.2312084958709, 73.1355968752219, 29.6828267493179, 75.2331288842908, 17.0049897184004, 86.5433552625001)
PLE_TRINITY_DN7629_c1_g1_i1 <- c(3.34276522512501, 0.785915767258281, 10.0099609353064, 2.69320707680029, 33.9828741992774, 7.14002766676061)
PLE_TRINITY_DN8063_c1_g1_i1 <- c(1.68292583732554, 2.92635711661422, 2.17352329159728, 0.774740001986499, 12.4658706785786, 1.87939975163101)
PLE_TRINITY_DN8779_c1_g4_i2 <- c(6.02361490031432, 2.13692100644751, 10.3729769131239, 0.497831632604912, 22.8104338179194, 2.77413419278609)
PLE_TRINITY_DN9848_c1_g2_i3 <- c(0.157057187688541, 0, 1.62160467011718, 3.61990164480147, 0.872834021308513, 5.66567136672498)



test_chitin<-data.frame(CON_TRINITY_DN8704_c0_g1_i1,
                        CON_TRINITY_DN6566_c2_g2_i2,
                        CON_TRINITY_DN6428_c3_g6_i1,
                        CON_TRINITY_DN6428_c3_g8_i1,
                        CON_TRINITY_DN7052_c1_g1_i4,
                        CON_TRINITY_DN8003_c0_g1_i1,
                        CON_TRINITY_DN4507_c0_g2_i1,
                        CON_TRINITY_DN4507_c0_g1_i1,
                        CON_TRINITY_DN7520_c0_g6_i1,
                        CON_TRINITY_DN5026_c0_g2_i1,
                        CON_TRINITY_DN3461_c0_g1_i1,
                        CON_TRINITY_DN3461_c0_g2_i1,
                        CON_TRINITY_DN6729_c0_g9_i1,
                        CON_TRINITY_DN6993_c4_g1_i1,
                        CON_TRINITY_DN6882_c4_g4_i5,
                        CON_TRINITY_DN5957_c1_g2_i1,
                        CON_TRINITY_DN5307_c1_g2_i1,
                        PLE_TRINITY_DN1325_c0_g1_i3,
                        PLE_TRINITY_DN1545_c0_g1_i1,
                        PLE_TRINITY_DN1776_c0_g2_i1,
                        PLE_TRINITY_DN1776_c0_g3_i2,
                        PLE_TRINITY_DN4742_c0_g4_i1,
                        PLE_TRINITY_DN5247_c0_g1_i2,
                        PLE_TRINITY_DN5341_c0_g3_i2,
                        PLE_TRINITY_DN5341_c0_g4_i1,
                        PLE_TRINITY_DN5411_c3_g2_i2,
                        PLE_TRINITY_DN6010_c0_g1_i4,
                        PLE_TRINITY_DN6573_c0_g2_i1,
                        PLE_TRINITY_DN6573_c0_g3_i2,
                        PLE_TRINITY_DN7191_c2_g1_i1,
                        PLE_TRINITY_DN7629_c1_g1_i1,
                        PLE_TRINITY_DN8063_c1_g1_i1,
                        PLE_TRINITY_DN8779_c1_g4_i2,
                        PLE_TRINITY_DN9848_c1_g2_i3
)
rownames(test_chitin) <- c("female flower", "leaf", "male flower", "petiole", "root", "veg bud")


test_chitin<-as.matrix(test_chitin)

pheatmap(log2(test_chitin+1), cluster_rows = TRUE, cluster_cols = FALSE)
