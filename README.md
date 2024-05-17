# Gene Regulatory Networks Reveal Sex Differences in Lung Adenocarcinoma

This repository contains code for replicating the analysis presented in the paper "Gene Regulatory Networks Reveal Sex Differences in Lung Adenocarcinoma": doi: https://doi.org/10.1101/2023.09.22.559001

## Discovery Data
Uniformly processed RNA-Seq data were downloaded from the Recount3 database for two discovery datasets on May 26, 2022: (i) healthy lung tissue samples from the Genotype Tissue Expression (GTEx) Project (version 8) and (ii) lung adenocarcinoma samples from The Cancer Genome Atlas (TCGA). Clinical data for GTEx were obtained from the dbGap website (https://dbgap.ncbi.nlm.nih.gov/) under study accession phs000424.v8.p2. Clinical data for TCGA were accessed from Recount3.

## Validation Data
Two independent studies from the Gene Expression Omnibus (GEO) were used as validation datasets: GSE47460 (or, LGRC) and GSE68465. From the LGRC (downloaded on Feb 12, 2023) data, we used 108 samples (59 female and 49 male) annotated as “control” samples for validation.

## Constructing Individual-specific Gene Regulatory Networks
PANDA and LIONESS algorithms were used to construct individual sample-specific gene regulatory networks from the discovery and validation datasets, using Python package netzooPy version 0.9.10. In addition to the gene expression data, two other data sources were utilized to construct the regulatory networks: Transcription factor/target gene regulatory prior (derived by mapping Transcription factor motifs from the Catalog of Inferred Sequence Binding Preferences (CIS-BP) to the promoter of target genes) and protein-protein interaction (using the interaction scores from StringDb v11.5 between all Transcription factor in the regulatory prior).

The networks are publicly available on the GRAND database: https://grand.networkmedicine.org/downloads/

## Code
R code for replicating the analysis are documented in README.txt

An executable notebook with the same title as the paper is available on netbooks (https://netbooks.networkmedicine.org) demonstrating a small representative differential targeting analysis on TCGA data.
