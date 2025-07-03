#!/bin/bash

# download miniconda
# manually there will be some options here. 
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# install certain conda environments

# BWA
conda create --name "for_bwa"
conda activate for_bwa
conda install bioconda::bwa
conda deactivate

# SAMtools
conda create --name "for_samtools"
conda activate for_samtools
conda install bioconda::samtools
conda deactivate

# BCFtools
# this one required a specific C lib be installed too on ubuntu
conda create --name "for_bcftools"
conda activate for_bcftools
conda install -c conda-forge libopenblas
conda install bioconda::bcftools
conda deactivate

# STAR
conda create --name "for_star"
conda activate for_star
conda install bioconda::star
conda deactivate

# BEDtools
conda create --name "for_bedtools"
conda activate for_bedtools
conda install bioconda::bedtools
conda deactivate

## SANTA CRUZ
## https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/
## ref genome
wget https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
# chrom sizes
wget https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
## get the md5sum file
wget https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/md5sum.txt


# GENCODE
# primary assembly fasta
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
# genes
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
## chrom sizes
samtools faidx GRCh38.primary_assembly.genome.fa
cut -f1,2 GRCh38.primary_assembly.genome.fa.fai > chrom.sizes.txt








