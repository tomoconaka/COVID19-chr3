hgdp_tgp_data=/home/richards/tomoko.nakanishi/scratch/data/hgdp_tgp
DATADIR=/home/richards/tomoko.nakanishi/scratch/09.COVID19/03.ICDA/03.GWAS_summary/summary_rev/
OUTDIR=/home/richards/tomoko.nakanishi/scratch/09.COVID19/03.ICDA/04.PRS

zcat $DATADIR/B2_meta_filtered.txt.gz | tail -n+2 | awk '{OFS="\t"; print $1, $2, $2+length($4)-1, toupper($4)"/"toupper($5), "+"}' > $DATADIR/meta.input.txt 
vep -i $DATADIR/meta.input.txt \
--plugin LoF,loftee_path:/scratch/richards/sirui.zhou/ensembl-vep/lof.v38/loftee/,human_ancestor_fa:/home/richards/sirui.zhou/scratch/WES_200K/scratch/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
--plugin CADD,/scratch/richards/sirui.zhou/ensembl-vep/whole_genome_SNVs.tsv.gz,/scratch/richards/sirui.zhou/ensembl-vep/gnomad.genomes.r3.0.indel.tsv.gz \
--offline --dir_plugins ~/scratch/data/VEP/loftee --dir_cache ~/scratch/data/VEP/ --force_overwrite \
--cache -o $DATADIR/meta.lof.cadd.gnomad.txt \
--fork 20 --everything --assembly GRCh38

less -S $DATADIR/meta.lof.cadd.gnomad.txt | grep -v "#" | cut -f1,13 | grep "rs" | sort | uniq > $DATADIR/marker.rsid.map

awk '{print $0, "SNP"}' <(zcat $DATADIR/B2_meta_filtered.txt.gz | head -1 ) > $DATADIR/B2_meta_filtered.rsid.txt
awk '(FNR == NR){m[$1]=$2; next}($1"_"$2"_"toupper($4)"/"toupper($5) in m){print $0,m[$1"_"$2"_"toupper($4)"/"toupper($5)]}($1"_"$2"_"toupper($5)"/"toupper($4) in m){print $0,m[$1"_"$2"_"toupper($5)"/"toupper($4)]}' \
$DATADIR/marker.rsid.map \
<(zcat $DATADIR/B2_meta_filtered.txt.gz | tail -n+2) >> $DATADIR/B2_meta_filtered.rsid.txt

for p in 5e-8 5e-7 5e-6 5e-5 5e-4 5e-3 5e-2
do
for r in 0.1 0.3 0.5 0.7
do
plink \
    --bfile $DATADIR/hgdp_1kg \
    --clump-p1 ${p} \
    --clump-r2 ${r} \
    --clump-kb 250 \
    --clump $DATADIR/B2_meta_filtered.rsid.txt \
    --clump-snp-field SNP \
    --clump-field P-value \
    --out ${OUTDIR}/${p}_${r}
done
done

