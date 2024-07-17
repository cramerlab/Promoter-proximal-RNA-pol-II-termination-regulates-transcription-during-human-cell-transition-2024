# Promoter-proximal RNA polymerase II termination regulates transcription during human cell type transition
Scripts for the analysis presented in the article: Promoter-proximal RNA polymerase II termination regulates transcription during human cell type transition

### Data availability
Fastq files used for this analysis are available from GEO series GSE235181.

### Data processing

TT-seq and mNET-seq reads were aligned using STAR.
ChIP-seq and ChIP-nexus reads were aligned using Bowtie2.

### Annotations and other utilities

Annotation: The major isoform annotation used in the study is available in data/

Other scripts and routines for performing the analysis are available in data/GatheringRoutines and data/Utility, which can be obtained by sourcing the script data/setup.R.  
If not sourced in the previous step, source the script data/processing_functions.R.

### Analysis

The analysis.R script provides step-by-step code to perform all of the analyses described in the article. The script generates all the necessary datasets to derive the conclusions of the study.

#### Main functionalities of the analysis.R script

1. Define the necessary information about the location of the input BAM files, the location of the output directory, etc.
2. Create RLE tracks for each datasets used in the study.
3. Run the routines for read counting on the different datasets.
4. Perform differential expression analysis to identify genes with significant changes in RNA synthesis (based on TT-seq data).
5. Calculate the promoter-proximal occupancy of RNA polymerase (Pol) II from mNET-seq data (incl. routines for determining pause positions and defining promoter-proximal regions).
6. Classification of genes based on expression pattern and gene clustering based on transcripiton kinetic parameters.
7. Exponential fitting of ChIP-nexus data (following a time course of transcription initiation inhibition).
8. Estimation of productive initiation frequency, apparent pause duration, promoter-proximal Pol II half-life, total turnover rate and termination fraction.
  
