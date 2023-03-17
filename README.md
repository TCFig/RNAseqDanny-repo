# RNAseqDanny-repo
This repository has code from an online video class focused on creating a pipeline for RNA sequencing - available from https://www.youtube.com/playlist?list=PLhR2Go-lh6X63hnyBzwWNvsaw1R79ESPI.

The data we worked on is from an experiment in S. cerevisiae where they were used to study the mechanism of endogenous H2S that promotes the growth rate of yeast - available from https://www.ebi.ac.uk/ena/browser/view/PRJNA714796.

The present code is only from the second class, since the first class is only for setting up the experiment, by installing the prerequired software.

The pipeline is composed of:
- Trimmomatric, for read trimming and adapter removal
- STAR, to perform the splice aware alignment
- picard, to remove duplicate reads (due to pcr)
- gatk, to do INDEL realignment and base recalibration

The script entitled buildGenome.R was used to create a primary assembly from the chromossomes available from ensembl (https://ftp.ensembl.org/pub/release-109/fasta/saccharomyces_cerevisiae/dna/). This was done because the database only had the genome TopLevel fasta file.
