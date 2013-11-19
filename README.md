flaveria_pipeline
=================

Pipeline for processing RNAseq data from Flaveria.
Four species, two C3 and two C4.

Flaveria species:
 - Flaveria pringlei
 - Flaveria trinervia
 - ?
 - ?

Data contains samples from six sections of leaves from each species with three replicates of each section.

## Steps in pipeline

Fastqc
 - Quality checking raw reads

Trimmomatic
 - Removing adapters?
 - Trimming low quality sequences at the ends of reads. Cutoff set to 10? 15? [1]

Reptile
 - Read error correction using hamming trees
 - Chosen because it came out best in this paper [2]

Khmer
 - preparing reads for assembly
 - reduce coverage

Assemblotron
 - optimal transcriptome assembly

Annotation
 - Reciprocal best usearch transcripts against Arabidopsis
 - Round-robin rbusearch between flaveria species to update annotation

Expression quantification
 - mapping reads output from trimmomatic to transcripts
 - using Express[3] or Sailfish[4]

Differential Expression Analysis
 - using EBSeq or possibly BaySeq

## How to Run

Type `rake`




[1] Paper about optimising trimmomatic cutoff
[2] Paper comparing read error correction methods
[3] Paper about Express
[4] Paper about Sailfish