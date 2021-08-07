conda activate hail
export PYSPARK_SUBMIT_ARGS="--driver-memory 400G pyspark-shell"
export JAVA_OPTS="-Xss128k -Xms256m -Xmx512m -XX:PermSize=64m -XX:MaxPermSize=1024m"

DATADIR=/home/richards/tomoko.nakanishi/scratch/09.COVID19/03.ICDA/03.GWAS_summary
PNDDIR=/home/richards/tomoko.nakanishi/09.COVID19/scratch/03.ICDA/03.GWAS_summary/qqplot
for i in $(ls $DATADIR/*)
do
j=$(basename $i)
python GWASplot.py --input $i --outdir $PNDDIR --pheno --qqplot
done

