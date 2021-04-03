library(edgeR)
library(dplyr)
library(magrittr)
library(stringr)
library(readr)

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


#A = female flower
#B = leaf
#C = male flower
#D = petiole
#E = root
#F = vegetative bud

con_fflower_v_leaf_et <- exactTest(con_d, pair=c("A","B"))
con_fflower_v_mflower_et <- exactTest(con_d, pair=c("A","C"))
con_fflower_v_petiole_et <- exactTest(con_d, pair=c("A","D"))
con_fflower_v_root_et <- exactTest(con_d, pair=c("A","E"))
con_fflower_v_vbud_et <- exactTest(con_d, pair=c("A","F"))
con_leaf_v_mflower_et <- exactTest(con_d, pair=c("B","C"))
con_leaf_v_petiole_et <- exactTest(con_d, pair=c("B","D"))
con_leaf_v_root_et <- exactTest(con_d, pair=c("B","E"))
con_leaf_v_vbud_et <- exactTest(con_d, pair=c("B","F"))
con_mflower_v_petiole_et <- exactTest(con_d, pair=c("C","D"))
con_mflower_v_root_et <- exactTest(con_d, pair=c("C","E"))
con_mflower_v_vbud_et <- exactTest(con_d, pair=c("C","F"))
con_petiole_v_root_et <- exactTest(con_d, pair=c("D","E"))
con_petiole_v_vbud_et <- exactTest(con_d, pair=c("D","F"))
con_root_v_vbud_et <- exactTest(con_d, pair=c("E","F"))

con_fflower_v_leaf_toptags<-con_fflower_v_leaf_et %>% topTags(n=1000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(con_go, by="transcript_id") 
con_fflower_v_mflower_toptags<-con_fflower_v_mflower_et %>% topTags(n=1000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(con_go, by="transcript_id") 
con_fflower_v_petiole_toptags<-con_fflower_v_petiole_et %>% topTags(n=1000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(con_go, by="transcript_id") 
con_fflower_v_root_toptags<-con_fflower_v_root_et %>% topTags(n=1000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(con_go, by="transcript_id") 
con_fflower_v_vbud_toptags<-con_fflower_v_vbud_et %>% topTags(n=1000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(con_go, by="transcript_id") 
con_leaf_v_mflower_toptags<-con_leaf_v_mflower_et %>% topTags(n=1000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(con_go, by="transcript_id") 
con_leaf_v_petiole_toptags<-con_leaf_v_petiole_et %>% topTags(n=1000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(con_go, by="transcript_id") 
con_leaf_v_root_toptags<-con_leaf_v_root_et %>% topTags(n=1000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(con_go, by="transcript_id") 
con_leaf_v_vbud_toptags<-con_leaf_v_vbud_et %>% topTags(n=1000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(con_go, by="transcript_id") 
con_mflower_v_petiole_toptags<-con_mflower_v_petiole_et %>% topTags(n=1000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(con_go, by="transcript_id") 
con_mflower_v_root_toptags<-con_mflower_v_root_et %>% topTags(n=1000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(con_go, by="transcript_id") 
con_mflower_v_vbud_toptags<-con_mflower_v_vbud_et %>% topTags(n=1000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(con_go, by="transcript_id") 
con_petiole_v_root_toptags<-con_petiole_v_root_et %>% topTags(n=1000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(con_go, by="transcript_id") 
con_petiole_v_vbud_toptags<-con_petiole_v_vbud_et %>% topTags(n=1000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(con_go, by="transcript_id") 
con_root_v_vbud_toptags<-con_root_v_vbud_et %>% topTags(n=1000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(con_go, by="transcript_id") 




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



#G = female flower = intercept
#H = leaf = coef 2
#I = male flower = coef 3
#J = petiole = coef 4
#K = root = coef 5
#L = vegetative bud = coef 6

ple_fflower_v_leaf_et <- exactTest(ple_d, pair=c("G","H"))
ple_fflower_v_mflower_et <- exactTest(ple_d, pair=c("G","I"))
ple_fflower_v_petiole_et <- exactTest(ple_d, pair=c("G","J"))
ple_fflower_v_root_et <- exactTest(ple_d, pair=c("G","K"))
ple_fflower_v_vbud_et <- exactTest(ple_d, pair=c("G","L"))
ple_leaf_v_mflower_et <- exactTest(ple_d, pair=c("H","I"))
ple_leaf_v_petiole_et <- exactTest(ple_d, pair=c("H","J"))
ple_leaf_v_root_et <- exactTest(ple_d, pair=c("H","K"))
ple_leaf_v_vbud_et <- exactTest(ple_d, pair=c("H","L"))
ple_mflower_v_petiole_et <- exactTest(ple_d, pair=c("I","J"))
ple_mflower_v_root_et <- exactTest(ple_d, pair=c("I","K"))
ple_mflower_v_vbud_et <- exactTest(ple_d, pair=c("I","L"))
ple_petiole_v_root_et <- exactTest(ple_d, pair=c("J","K"))
ple_petiole_v_vbud_et <- exactTest(ple_d, pair=c("J","L"))
ple_root_v_vbud_et <- exactTest(ple_d, pair=c("K","L"))

ple_fflower_v_leaf_toptags<-ple_fflower_v_leaf_et %>% topTags(n=10000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(ple_go, by="transcript_id") 
ple_fflower_v_mflower_toptags<-ple_fflower_v_mflower_et %>% topTags(n=10000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(ple_go, by="transcript_id") 
ple_fflower_v_petiole_toptags<-ple_fflower_v_petiole_et %>% topTags(n=10000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(ple_go, by="transcript_id") 
ple_fflower_v_root_toptags<-ple_fflower_v_root_et %>% topTags(n=10000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(ple_go, by="transcript_id") 
ple_fflower_v_vbud_toptags<-ple_fflower_v_vbud_et %>% topTags(n=10000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(ple_go, by="transcript_id") 
ple_leaf_v_mflower_toptags<-ple_leaf_v_mflower_et %>% topTags(n=10000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(ple_go, by="transcript_id") 
ple_leaf_v_petiole_toptags<-ple_leaf_v_petiole_et %>% topTags(n=10000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(ple_go, by="transcript_id") 
ple_leaf_v_root_toptags<-ple_leaf_v_root_et %>% topTags(n=10000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(ple_go, by="transcript_id") 
ple_leaf_v_vbud_toptags<-ple_leaf_v_vbud_et %>% topTags(n=10000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(ple_go, by="transcript_id") 
ple_mflower_v_petiole_toptags<-ple_mflower_v_petiole_et %>% topTags(n=10000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(ple_go, by="transcript_id") 
ple_mflower_v_root_toptags<-ple_mflower_v_root_et %>% topTags(n=10000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(ple_go, by="transcript_id") 
ple_mflower_v_vbud_toptags<-ple_mflower_v_vbud_et %>% topTags(n=10000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(ple_go, by="transcript_id") 
ple_petiole_v_root_toptags<-ple_petiole_v_root_et %>% topTags(n=10000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(ple_go, by="transcript_id") 
ple_petiole_v_vbud_toptags<-ple_petiole_v_vbud_et %>% topTags(n=10000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(ple_go, by="transcript_id") 
ple_root_v_vbud_toptags<-ple_root_v_vbud_et %>% topTags(n=10000) %>% data.frame() %>% rownames_to_column("transcript_id") %>% inner_join(ple_go, by="transcript_id") 



nrow(ple_fflower_v_leaf_toptags[ple_fflower_v_leaf_toptags$FDR < 0.05,])
nrow(ple_fflower_v_mflower_toptags[ple_fflower_v_mflower_toptags$FDR < 0.05,])
nrow(ple_fflower_v_petiole_toptags[ple_fflower_v_petiole_toptags$FDR < 0.05,])
nrow(ple_fflower_v_root_toptags[ple_fflower_v_root_toptags$FDR < 0.05,])
nrow(ple_fflower_v_vbud_toptags[ple_fflower_v_vbud_toptags$FDR < 0.05,])
nrow(ple_leaf_v_mflower_toptags[ple_leaf_v_mflower_toptags$FDR < 0.05,])
nrow(ple_leaf_v_petiole_toptags[ple_leaf_v_petiole_toptags$FDR < 0.05,])
nrow(ple_leaf_v_root_toptags[ple_leaf_v_root_toptags$FDR < 0.05,])
nrow(ple_leaf_v_vbud_toptags[ple_leaf_v_vbud_toptags$FDR < 0.05,])
nrow(ple_mflower_v_petiole_toptags[ple_mflower_v_petiole_toptags$FDR < 0.05,])
nrow(ple_mflower_v_root_toptags[ple_mflower_v_root_toptags$FDR < 0.05,])
nrow(ple_mflower_v_vbud_toptags[ple_mflower_v_vbud_toptags$FDR < 0.05,])
nrow(ple_petiole_v_root_toptags[ple_petiole_v_root_toptags$FDR < 0.05,])
nrow(ple_petiole_v_vbud_toptags[ple_petiole_v_vbud_toptags$FDR < 0.05,])
nrow(ple_root_v_vbud_toptags[ple_root_v_vbud_toptags$FDR < 0.05,])




test<-ple_fflower_v_leaf_toptags %>% filter(FDR < 0.05 & logFC > 2) %>% dplyr::select(transcript_id)
write.table(test, file="ple_fflower_v_leaf_DE", quote=FALSE, row.names = FALSE, col.names = FALSE)



tail(ple_petiole_v_root_toptags$FDR, n=200)

colnames(ple_fflower_v_leaf_toptags)




