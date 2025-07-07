## #!/bin/bash
set -euo pipefail

# Create directories
mkdir -p ucsc gencode refseq dbsnp

# -------- UCSC --------
echo "Downloading UCSC reference files..."
cd ucsc

# Genome FASTA and chrom sizes
wget -N https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
wget -N https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
wget -N https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/md5sum.txt

# Optional UCSC Genes GTF
wget -N http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ucscGenes.gtf.gz

cd ..

# -------- GENCODE --------
echo "Downloading GENCODE reference files..."
cd gencode

# Genome FASTA and GTF
wget -N https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
wget -N https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz

# Generate chrom sizes
echo "Indexing genome to create chrom sizes..."
samtools faidx GRCh38.primary_assembly.genome.fa.gz
cut -f1,2 GRCh38.primary_assembly.genome.fa.fai > chrom.sizes.txt

cd ..

# -------- RefSeq --------
echo "Downloading RefSeq reference files..."
cd refseq

# FASTA and GFF
wget -N https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
wget -N https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz

cd ..

# -------- dbSNP (shared) --------
echo "Downloading dbSNP for GRCh38 and renaming chromosomes..."
cd dbsnp

# Download raw dbSNP VCF
wget -N https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
wget -N https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi

# Download and build chr map
wget -N https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40/GCF_000001405.40_assembly_report.txt
awk -F '\t' 'BEGIN{OFS="\t"} !/^#/ && $5 != "na" {print $7, $5}' GCF_000001405.40_assembly_report.txt > chr_map.txt

# Convert to UCSC-style
bcftools annotate --rename-chrs chr_map.txt -O z -o dbsnp_grch38_chr.vcf.gz GCF_000001405.40.gz
tabix -p vcf dbsnp_grch38_chr.vcf.gz

cd ..

echo "✅ All reference files downloaded and prepared."
