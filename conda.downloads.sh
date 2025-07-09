#!/bin/bash

# download miniconda
# manually there will be some options here. 
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# install certain conda environments

# minimap2 for long reads
conda create --name "for_minimap2"
conda activate for_minimap2
conda install bioconda::minimap2
conda deactivate
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
##  For mapping Ensembl IDs to gene symbols:
grep -v "^#" gencode.v44.annotation.gtf | awk '$3 == "gene"' | \
awk -F'\t|; ' '{print $1, $9, $10, $11}' > ensembl_gene_symbols.txt
## make the star file
STAR --runMode genomeGenerate \
     --runThreadN 8 \
     --genomeDir STAR_genome \
     --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile gencode.v44.annotation.gtf \
     --sjdbOverhang 100


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

## Download docker (for nextflow)
# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "${UBUNTU_CODENAME:-$VERSION_CODENAME}") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update

sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
sudo docker run hello-world
## install nextflow
java -version                           # Check that Java v11+ is installed. i think java 17+ needs to be here
curl -s https://get.nextflow.io | bash  # Download Nextflow
chmod +x nextflow                       # Make executable
mv nextflow ~/bin/                      # Add to user's $PATH. had to make ~/bin. do not use /bin/
