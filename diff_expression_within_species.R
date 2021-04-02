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

# create DGE list
# keep only rows where counts per million is above 100 in at least 2 samples
# calculate normalisation factors
# estimate tagwise and common dispersal
con_d<-DGEList(counts=con_counts,group=factor(con_groups))
con_d.full<-con_d
con_d.keep<-rowSums(cpm(con_d)>1)>=2
con_d<-con_d[con_d.keep,]
con_d_lengths<-con_lengths[con_d.keep,] %>% data.frame()
con_d$samples$lib.size<-colSums(con_d$counts)
con_d<-calcNormFactors(con_d, method="TMM")
plotMDS(con_d, method="bcv", col=as.numeric(con_d$samples$group))
con_d<-estimateCommonDisp(con_d)
con_d<-estimateTagwiseDisp(con_d, prior.n = 4)


ple_d<-DGEList(counts= ple_counts,group=factor(ple_groups))
ple_d.full<-ple_d
ple_d.keep<-rowSums(cpm(ple_d)>1)>=2
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

fflower_v_leaf_et <- exactTest(ple_d, pair=c("G","H"))
fflower_v_mflower_et <- exactTest(ple_d, pair=c("G","I"))
fflower_v_petiole_et <- exactTest(ple_d, pair=c("G","J"))
fflower_v_root_et <- exactTest(ple_d, pair=c("G","K"))
fflower_v_vbud_et <- exactTest(ple_d, pair=c("G","L"))

leaf_v_mflower_et <- exactTest(ple_d, pair=c("H","I"))
leaf_v_petiole_et <- exactTest(ple_d, pair=c("H","J"))
leaf_v_root_et <- exactTest(ple_d, pair=c("H","K"))
leaf_v_vbud_et <- exactTest(ple_d, pair=c("H","L"))

mflower_v_petiole_et <- exactTest(ple_d, pair=c("I","J"))
mflower_v_root_et <- exactTest(ple_d, pair=c("I","K"))
mflower_v_vbud_et <- exactTest(ple_d, pair=c("I","L"))

petiole_v_root_et <- exactTest(ple_d, pair=c("J","K"))
petiole_v_vbud_et <- exactTest(ple_d, pair=c("J","L"))

root_v_vbud_et <- exactTest(ple_d, pair=c("K","L"))

fflower_v_leaf_toptags<-topTags(fflower_v_leaf_et)
fflower_v_mflower_toptags<-topTags(fflower_v_mflower_et)
fflower_v_petiole_toptags<-topTags(fflower_v_petiole_et)
fflower_v_root_toptags<-topTags(fflower_v_root_et)
fflower_v_vbud_toptags<-topTags(fflower_v_vbud_et)
leaf_v_mflower_toptags<-topTags(leaf_v_mflower_et)
leaf_v_petiole_toptags<-topTags(leaf_v_petiole_et)
leaf_v_root_toptags<-topTags(leaf_v_root_et)
leaf_v_vbud_toptags<-topTags(leaf_v_vbud_et)
mflower_v_petiole_toptags<-topTags(mflower_v_petiole_et)
mflower_v_root_toptags<-topTags(mflower_v_root_et)
mflower_v_vbud_toptags<-topTags(mflower_v_vbud_et)
petiole_v_root_toptags<-topTags(petiole_v_root_et)
petiole_v_vbud_toptags<-topTags(petiole_v_vbud_et)
root_v_vbud_toptags<-topTags(root_v_vbud_et)







con_annot<-read_tsv("/Users/katie/Desktop/Bg/begonia_duplicate_expression/con_trinotate_annotation_report.xls")
ple_annot<-read_tsv("/Users/katie/Desktop/Bg/begonia_duplicate_expression/ple_trinotate_annotation_report.xls")



ple_annot$transcript_id


fflower_v_leaf_toptags<-
fflower_v_mflower_toptags<-
fflower_v_petiole_toptags<-
fflower_v_root_toptags<-
fflower_v_vbud_toptags<-
leaf_v_mflower_toptags<-
leaf_v_petiole_toptags<-
leaf_v_root_toptags<-
leaf_v_vbud_toptags<-
mflower_v_petiole_toptags<-
mflower_v_root_toptags<-
mflower_v_vbud_toptags<-
petiole_v_root_toptags<-
petiole_v_vbud_toptags<-
root_v_vbud_toptags<-


topTags(fflower_v_leaf_et)





con 
(77736/200849)*100

(32848/42614)*100

ple
(86926/251473)*100

(43119/59106)*100


con_d$counts %>% data.frame() %>% 
  select(CONfemaleFlower1, CONfemaleFlower2, CONfemaleFlower3) %>% 
  rpkmByGroup(gene.length=con_d_lengths)

CONfemaleFlower<-con_d$counts %>% data.frame() %>% select(CONfemaleFlower1, CONfemaleFlower2, CONfemaleFlower3) %>% rpkmByGroup(gene.length=con_d_lengths)  %>% data.frame()
CONleaf<-con_d$counts %>% data.frame() %>% select(CONleaf1, CONleaf2, CONleaf3) %>% rpkmByGroup(gene.length=con_d_lengths) %>% data.frame()
CONmaleFlower<-con_d$counts %>% data.frame() %>% select(CONmaleFlower1, CONmaleFlower2, CONmaleFlower3) %>% rpkmByGroup(gene.length=con_d_lengths) %>% data.frame()
CONpetiole<-con_d$counts %>% data.frame() %>% select(CONpetiole1, CONpetiole2, CONpetiole3) %>% rpkmByGroup(gene.length=con_d_lengths) %>% data.frame()
CONroot<-con_d$counts %>% data.frame() %>% select(CONroot1, CONroot3, CONvegBud1) %>% rpkmByGroup(gene.length=con_d_lengths) %>% data.frame()
CONvegBud<-con_d$counts %>% data.frame() %>% select(CONvegBud1, CONvegBud2, CONvegBud3) %>% rpkmByGroup(gene.length=con_d_lengths) %>% data.frame()

PLEfemaleFlower<-ple_d$counts %>% data.frame() %>% select(PLEfemaleFlower1, PLEfemaleFlower2, PLEfemaleFlower3) %>% rpkmByGroup(gene.length=ple_d_lengths)  %>% data.frame()
PLEleaf<-ple_d$counts %>% data.frame() %>% select(PLEleaf1, PLEleaf2, PLEleaf3) %>% rpkmByGroup(gene.length=ple_d_lengths)  %>% data.frame()
PLEmaleFlower<-ple_d$counts %>% data.frame() %>% select(PLEmaleFlower1, PLEmaleFlower2, PLEmaleFlower3) %>% rpkmByGroup(gene.length=ple_d_lengths)  %>% data.frame()
PLEpetiole<-ple_d$counts %>% data.frame() %>% select(PLEpetiole1, PLEpetiole2) %>% rpkmByGroup(gene.length=ple_d_lengths)  %>% data.frame()
PLEroot<-ple_d$counts %>% data.frame() %>% select(PLEroot1, PLEroot2, PLEroot3) %>% rpkmByGroup(gene.length=ple_d_lengths) %>% data.frame()
PLEvegBud<-ple_d$counts %>% data.frame() %>% select(PLEvegBud1, PLEvegBud2, PLEvegBud3) %>% rpkmByGroup(gene.length=ple_d_lengths)  %>% data.frame()


# bind con and ple together per species
con_avg_cpm<-cbind(CONfemaleFlower, CONleaf, CONmaleFlower, CONpetiole, CONroot, CONvegBud)
ple_avg_cpm<-cbind(PLEfemaleFlower, PLEleaf, PLEmaleFlower, PLEpetiole, PLEroot, PLEvegBud)

# transfer row names
rownames(con_avg_cpm)<-rownames(con_d$counts)
rownames(ple_avg_cpm)<-rownames(ple_d$counts)

# assign colnames
colnames(con_avg_cpm)<-c("CONfemaleFlower", "CONleaf", "CONmaleFlower", "CONpetiole", "CONroot", "CONvegBud")
colnames(ple_avg_cpm)<-c("PLEfemaleFlower", "PLEleaf", "PLEmaleFlower", "PLEpetiole", "PLEroot", "PLEvegBud")

# write to table
write.table(con_avg_cpm, file="con_avg_cpm", sep="\t", quote = FALSE, row.names = TRUE)
write.table(ple_avg_cpm, file="ple_avg_cpm", sep="\t", quote = FALSE, row.names = TRUE)



