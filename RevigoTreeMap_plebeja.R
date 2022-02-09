# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006857","oligopeptide transport",0.069,8.000,0.986,0.000,"oligopeptide transport"),
c("GO:0010150","leaf senescence",0.006,15.000,0.785,0.000,"leaf senescence"),
c("GO:0010233","phloem transport",0.000,5.000,0.888,0.321,"leaf senescence"),
c("GO:0010014","meristem initiation",0.001,9.000,0.865,0.363,"leaf senescence"),
c("GO:0010050","vegetative phase change",0.000,5.000,0.823,0.364,"leaf senescence"),
c("GO:0035019","somatic stem cell population maintenance",0.007,4.000,0.848,0.382,"leaf senescence"),
c("GO:0006949","syncytium formation",0.010,7.000,0.849,0.460,"leaf senescence"),
c("GO:0051301","cell division",1.474,31.000,0.992,0.000,"cell division"),
c("GO:0055114","(obsolete) oxidation-reduction process",0.130,67.000,1.000,0.000,"(obsolete) oxidation-reduction process"),
c("GO:0009834","plant-type secondary cell wall biogenesis",0.010,16.000,0.965,0.007,"plant-type secondary cell wall biogenesis"),
c("GO:0009873","ethylene-activated signaling pathway",0.013,23.000,0.852,0.007,"ethylene-activated signaling pathway"),
c("GO:2000024","regulation of leaf development",0.001,8.000,0.956,0.121,"ethylene-activated signaling pathway"),
c("GO:2000012","regulation of auxin polar transport",0.002,7.000,0.915,0.124,"ethylene-activated signaling pathway"),
c("GO:0009637","response to blue light",0.020,17.000,0.883,0.204,"ethylene-activated signaling pathway"),
c("GO:2001141","regulation of RNA biosynthetic process",9.906,6.000,0.930,0.231,"ethylene-activated signaling pathway"),
c("GO:0048016","inositol phosphate-mediated signaling",0.010,5.000,0.887,0.343,"ethylene-activated signaling pathway"),
c("GO:0010362","negative regulation of anion channel activity by blue light",0.000,5.000,0.940,0.438,"ethylene-activated signaling pathway"),
c("GO:0010200","response to chitin",0.004,16.000,0.910,0.477,"ethylene-activated signaling pathway"),
c("GO:0060919","auxin influx",0.000,6.000,0.928,0.515,"ethylene-activated signaling pathway"),
c("GO:0009864","induced systemic resistance, jasmonic acid mediated signaling pathway",0.001,6.000,0.867,0.607,"ethylene-activated signaling pathway"),
c("GO:0009644","response to high light intensity",0.005,12.000,0.889,0.679,"ethylene-activated signaling pathway"),
c("GO:0015979","photosynthesis",0.221,23.000,0.969,0.009,"photosynthesis"),
c("GO:0018298","protein-chromophore linkage",0.140,14.000,0.939,0.051,"protein-chromophore linkage"),
c("GO:0015995","chlorophyll biosynthetic process",0.050,12.000,0.909,0.108,"protein-chromophore linkage"),
c("GO:0009813","flavonoid biosynthetic process",0.008,10.000,0.945,0.129,"protein-chromophore linkage"),
c("GO:0005992","trehalose biosynthetic process",0.069,11.000,0.932,0.152,"protein-chromophore linkage"),
c("GO:0006270","DNA replication initiation",0.124,8.000,0.851,0.166,"protein-chromophore linkage"),
c("GO:0010206","photosystem II repair",0.002,6.000,0.950,0.186,"protein-chromophore linkage"),
c("GO:0015936","coenzyme A metabolic process",0.187,5.000,0.926,0.208,"protein-chromophore linkage"),
c("GO:0006799","polyphosphate biosynthetic process",0.029,4.000,0.926,0.231,"protein-chromophore linkage"),
c("GO:0009809","lignin biosynthetic process",0.005,12.000,0.884,0.253,"protein-chromophore linkage"),
c("GO:0006267","pre-replicative complex assembly involved in nuclear cell cycle DNA replication",0.007,4.000,0.866,0.637,"protein-chromophore linkage"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
jpeg( file="figure9_revigo_treemap_plebeja.jpeg", width = 1280, height = 880  ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none",
  align.labels = list(c("left", "top"), c("center", "center")),
  fontsize.labels = c(19, 17)
)

dev.off()

