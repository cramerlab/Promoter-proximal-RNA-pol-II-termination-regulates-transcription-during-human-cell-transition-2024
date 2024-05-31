# Promoter-proximal RNA polymerase II termination regulates transcription during human cell type transition
Codes for the analysis presented in the article : Promoter-proximal RNA polymerase II termination regulates transcription during human cell type transition

### Data availability
Fastq files and final data CSV file used on this analysis are available from GEO series GSE235181.

### Data processing

TT-seq and mNET-seq reads were aligned using STAR and filtered with Samtools. 
ChIP-seq and ChIP-nexus reads were aligned using Bowtie2 and filtered with Samtools. 

### Annotations and other utilities

Annotation:Major isoform annotation used in the study are available in data/

Other scripts and routines to run the analysis are available in data/GatheringRoutines and data/Utility which can be sourced by sourcing the data/setup.R script
Source the data/processing_fucntions.R script if it is not sourced in the previous step.


### Analysis

The analysis.R script contains step by step codes to perform all the analysis described in the article. The script generates all the necessary datasets required to derive conclusions of the study independently.

#### Main functionalities of the analysis.R script

1. Can define necessary information about the location of input BAM files, location of output directory etc.

2. Create RLE tracks for each datasets used in the study.

3. Run the routines for counting on the different datasets.

4. Run the differential expression routine to identify the genes with significant changes in TT-seq expression

5. Calculate the promoter-pproximal occupancy from mNET-seq data (routines available for calling pause positions and defining promoter-proximal regions)

6. Cluster genes based on kinetic parameter trends

7. Exponential fitting of ChIP-nexus data

8. Estimation of productive initiation frequency, apparent pause duration, half-life, total turnover rate, termination fraction.

