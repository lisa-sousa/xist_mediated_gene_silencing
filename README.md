# Kinetics of Xist-induced gene silencing depend on combinations of epigenetic and genomic features

### Abstract
To initiate X-chromosome inactivation (XCI), the long non-coding RNA Xist mediates chromosome-wide gene silencing of one X chromosome in female mammals to equalize gene dosage between the sexes. The efficiency of gene silencing, however is highly variable across genes, with some genes even escaping XCI in somatic cells. A gene’s susceptibility to Xist-mediated silencing appears to be determined by a complex interplay of epigenetic and genomic features; however, the underlying rules remain poorly understood. We have quantified chromosome-wide gene silencing kinetics at the level of the nascent transcriptome using allele-specific Precision nuclear Run-On sequencing (PRO-seq). We have developed a Random Forest machine learning model that can predict the measured silencing dynamics based on a large set of epigenetic and genomic features and tested its predictive power experimentally. The genomic distance to the Xist locus, followed by gene density and distance to LINE elements are the prime determinants of the speed of gene silencing, Moreover, we find two distinct gene clusters associated with different silencing pathways: a cluster that require Xist repeat A for silencing, which is known to recruit the Spen-pathway and a second cluster where genes are pre-marked by Polycomb complexes and tend to rely on the B- and C repeats in Xist for silencing, known to ecruit polycomb complexes during XCI. Moreover, a series of features associated with active transcriptional elongation and chromatin 3D structure are enriched at rapidly silenced genes. Our machine learning approach can thus uncover the complex combinatorial rules underlying gene silencing during X inactivation.

### Structure of the repository
The repository is structured into four parts according to four main analysis steps, each analysis is step structured into a seperate folder.

#### fit_silencing_halftime
This folder contains two scripts that are used to compute half-times for each gene. 

a) *fit_halftimes.R*: computes half-times for the PRO-seq and mRNA-seq (undifferentiated and differentiated) data set; takes as input the allelic ratio from the PRO-seq/mRNA-seq experiment and the Gencode M9 gene annotation

b) *fit_halftimes_pyrosequencing.R*: computes half-times for the 11 pyrosequencing candidate genes; takes as input the list of candidate genes and the allelic ratios from the pyrosequencing experiment

#### chip_seq_analysis
This folder contains all scripts that were used to pre-process all chip-seq data sets.

a) *1_mapping.py*: this script maps chip-seq data to the given genome; takes as input a sra or fastq file; the path to the used tools as well as the genome assembly that should be used must be specified in the script; an example call is given in the header of the script

b) *2_plot_coverage.R*: this script claculates the library coverage for a list of given chip-seq experiments; it takes as input the input directories for chip and control experiments as well as a metadata file where the chip-seq name (e.g. CTCF), the GEO number, file name of experiment and control are defined 

c) *3_plot_fingerprints.R*: this script produced the fingerprint plots (e.g. used to create Figure SF5) for each chip-seq data set; it takes as input the input directories for chip and control experiments as well as a metadata file where the chip-seq name (e.g. CTCF), the GEO number, file name of experiment and control are defined

d) *4_get_active_gene_promoter.R*: this script assigns each gene its active TSS based on the overlap of a gene's TSS with a regulatory region (caculated as described in the manuscript); it takes as input a Gencode gene annotation and the file with the processed regulatory regions

e) *5_plot_heatmap_deepTools_control_vs_experiment.R*: this script creates the deepTools heatmap plots (e.g. used to create Figure SF4/SF5); its takes as input a bed file with gene annotation, the input directories for chip and control experiments as well as a metadata file where the chip-seq name (e.g. CTCF), the GEO number, file name of experiment and control are defined; the plotting window (+/- x bp around the TSS) is defined in the script (a/b)

f) *6_normalize_ChIP_data_with_normr.R*: this script is used to calculate the normalized enrichemnt scores in a window specified for each chip-seq feature (e.g. 500 bp around TSS or gene body); its takes as input a bed file with gene annotation, the input directories for chip and control experiments as well as a metadata file where the chip-seq name (e.g. CTCF), the GEO number, file name of experiment and control, the # of bp upstream of TSS (a), the # of bp downstream of TSS (b) (for gene body specify either a=NA and b=x which mean from TSS to x bp downstream OR gene end if x is longer then the gene; or specify a=NA and b=NA for whole gene body); binsize for normR normalization (needs to be adapted to the width of the chip-seq peak, e.g. 1000-1500 for broad peaks and 250-750 for narrow peaks) are defined

#### modelling
This folder contain three scripts which are used to create the feature matrix and to build the Random Forest, do predictions and plot the forest-guided clustering.
a) *create_feature_matrix.R*: this script is used to compute the feature matrix, it loads the pre-processed chip-seq data (normalized enrichment for pre-defined regions by *6_normalize_ChIP_data_with_normr.R*) and computes each genomic feature as described in the manuscript
b) *model_functions.R*: this script contains all functions used by *model.R* to load the feature matrix and half-times, create boxplots for each feature seperated by class, optimize Random Forest parameters, run stability test of Random Forest model, get the number of top features that optimize the error rate, make predictions for genes used to build the Random Forst and new genes, plot error rate and feature importance of the RF model, perform forest-guided (proximity) clustering and optimize number of clusters k, plot heatmap and boxplots for clustering
c) *model.R*: this script is used to build the RF model, do prediction and do the forest-guided clustering for a set of threshold combination specified in the script and outputs all the plots for each threshold into a subfolder and the statistics for all threshold combinations into a text document in the main folder; this script was used to create all figures related to the RF model (e.g. model error SF7, feature importance SF8/9, clustering stability SF 10, cluster heatmaps Fig. 5/6)


