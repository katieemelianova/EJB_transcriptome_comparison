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
revigo.data <- rbind(c("GO:0009737","response to abscisic acid",0.056,54.000,0.902,0.000,"response to abscisic acid"),
c("GO:0009414","response to water deprivation",0.020,28.000,0.888,0.644,"response to abscisic acid"),
c("GO:0051301","cell division",1.474,37.000,0.989,0.000,"cell division"),
c("GO:0055114","(obsolete) oxidation-reduction process",0.130,48.000,1.000,0.000,"(obsolete) oxidation-reduction process"),
c("GO:2000652","regulation of secondary cell wall biogenesis",0.001,11.000,1.000,0.000,"regulation of secondary cell wall biogenesis"),
c("GO:0009834","plant-type secondary cell wall biogenesis",0.010,21.000,0.816,0.007,"plant-type secondary cell wall biogenesis"),
c("GO:0030261","chromosome condensation",0.088,5.000,0.783,0.267,"plant-type secondary cell wall biogenesis"),
c("GO:0010207","photosystem II assembly",0.008,7.000,0.799,0.305,"plant-type secondary cell wall biogenesis"),
c("GO:0010417","glucuronoxylan biosynthetic process",0.002,6.000,0.762,0.472,"plant-type secondary cell wall biogenesis"),
c("GO:0008283","cell population proliferation",0.101,6.000,0.991,0.009,"cell population proliferation"),
c("GO:0006270","DNA replication initiation",0.124,10.000,0.811,0.009,"DNA replication initiation"),
c("GO:0006782","protoporphyrinogen IX biosynthetic process",0.218,8.000,0.882,0.185,"DNA replication initiation"),
c("GO:0010345","suberin biosynthetic process",0.002,7.000,0.903,0.263,"DNA replication initiation"),
c("GO:0000727","double-strand break repair via break-induced replication",0.011,6.000,0.830,0.350,"DNA replication initiation"),
c("GO:0006267","pre-replicative complex assembly involved in nuclear cell cycle DNA replication",0.007,4.000,0.684,0.637,"DNA replication initiation"),
c("GO:0007018","microtubule-based movement",0.347,25.000,0.990,0.010,"microtubule-based movement"),
c("GO:0005975","carbohydrate metabolic process",5.986,35.000,0.959,0.073,"carbohydrate metabolic process"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap_conchifolia.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  #index = c("description"),
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
  fontsize.labels = c(18, 15),
)

dev.off()

