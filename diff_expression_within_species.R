library(edgeR)
library(dplyr)
library(magrittr)
library(stringr)
library(readr)
library(tibble)
library(VennDiagram)
library(circlize)
library(ggplot2)

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

con_lengths<-con_counts %>% select(Length) %>% set_rownames(con_geneid)
ple_lengths<-ple_counts %>% select(Length) %>% set_rownames(ple_geneid)


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

con_chord_df<-data.frame(from=from, to=to, de_count=con_de_count)
ple_chord_df<-data.frame(from=from, to=to, de_count=ple_de_count)

pdf('con_circos_plot.pdf')
chordDiagram(con_chord_df, grid.col = circos_cols)
dev.off()
circos.clear()
pdf('ple_circos_plot.pdf')
chordDiagram(ple_chord_df, grid.col = circos_cols)
dev.off()

circos_cols<-c("hotpink",
"green4",
"orange",
"darkseagreen2",
"tan4",
"skyblue1")


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
  geom_bar(stat="identity", position = "dodge") + facet_grid(rows = vars(direction)) + coord_flip() 
  
#  theme(axis.text.y=element_text(size=12), 
#        axis.title.x=element_blank(), 
#        axis.title.y=element_blank())







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

# order the go terms by number DE genes in category for plotting
all_up$go_term<-factor(all_up$go_term, levels = unique(all_up$go_term[order(all_up$numDEInCat)]))
all_down$go_term<-factor(all_down$go_term, levels = unique(all_down$go_term[order(all_down$numDEInCat)]))


ggplot(all_up, aes(x=go_term, y=numDEInCat, colour=comparison, fill=comparison)) + 
  geom_bar(stat="identity", position = "dodge") + coord_flip() + facet_grid(cols = vars(species))

ggplot(all_down, aes(x=go_term, y=numDEInCat, colour=comparison, fill=comparison)) + 
  geom_bar(stat="identity", position = "dodge") + coord_flip() + facet_grid(cols = vars(species))

+ facet_grid(rows = vars(direction)) + coord_flip() 

test<-all$go_term[1]

test[1]

substring(all$go_term, 4)
