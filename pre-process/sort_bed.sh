echo start at `date`

module load miniconda
conda activate snapatac

# Make a bed file from output file of EB1_gene_peak.R and use bedtools intersect to get overlapped motifs in these peaks
cat ../output/EB1/test_gene_peak.txt | awk -v FS='\t' -v OFS='\t' '{print $1,$2,$3,$7,0,".",$4,$5,$6,$8}' | sort -k1,1 -k2,2n > ../output/EB1/test_gene_peak.bed
bedtools intersect -wa -wb -sorted -a ../data/database/homer.KnownMotifs.hg38.gene.sorted.bed -b ../output/EB1/test_gene_peak.bed > ../output/EB1/test_gene_TF.bed

echo end at `date`
