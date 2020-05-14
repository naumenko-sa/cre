#PBS -l walltime=20:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=40g,mem=40g

if [ -z $family ]
then
	family=$1
fi

/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin/bcbio-variation-recall ensemble --cores=7 --numpass 2 --names gatk-haplotype,platypus,freebayes,samtools --nofiltered ${family}.ensemble.vcf.gz /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa ${family}-gatk-haplotype-annotated-decomposed.vcf.gz ${family}-platypus-annotated-decomposed.vcf.gz ${family}-freebayes-annotated-decomposed.vcf.gz ${family}-samtools-annotated-decomposed.vcf.gz

