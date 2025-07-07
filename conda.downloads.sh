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


### UCSC uses chr
### gencode, refseq, ensmeble do not. 

## SANTA CRUZ
## santa cruz doesnt exactly have a GTF file. you would have to extract
## it from their genes file separately from table browser.
## https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/
## ref genome
wget https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
# chrom sizes
wget https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
## get the md5sum file
wget https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/md5sum.txt
## dbsnp
wget http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155Common.vcf.gz
## GTF equivolent
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ucscGenes.gtf.gz


# GENCODE
# primary assembly fasta
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
# genes GTF
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
## chrom sizes file generate
samtools faidx GRCh38.primary_assembly.genome.fa
cut -f1,2 GRCh38.primary_assembly.genome.fa.fai > chrom.sizes.txt
## dbSNP 
wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi

#Download the GRCh38 Assembly Report:
## fix the naming conventions
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40/GCF_000001405.40_assembly_report.txt
awk -F '\t' 'BEGIN{OFS="\t"} !/^#/ && $5 != "na" {print $7, $5}' GCF_000001405.40_assembly_report.txt > chr_map.txt
bcftools annotate --rename-chrs chr_map.txt -O z -o dbsnp_grch38_chr.vcf.gz GCF_000001405.40.gz
tabix -p vcf dbsnp_grch38_chr.vcf.gz


## RefSeq
# Genomic sequence (FASTA)
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
# Annotation (GFF)
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz







