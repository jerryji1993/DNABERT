# sed 's/N//g' $1 | sed 's/^(>.*)\n/^(>.*)\t/g' | tr -d '\n' | sed 's/\t/\n/g' > GRCh38.genome.1.fa

# use a loop to split whole genome to individual chromosomes, then remove N and newline char
for i in {1..22}
do
    samtools faidx GRCh38.primary_assembly.genome.fa chr$i > chr$i.fa
    sed '/^>/d' chr$i.fa | sed 's/N//g' | tr -d '\n\t\r' > GRCh38.chr$i.fa
    rm chr$i.fa
done