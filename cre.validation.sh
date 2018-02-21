#!/bin/bash

export PATH=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin:$PATH &&  export RTG_JAVA_OPTS='-Xms750m' && export RTG_MEM=10000m && rtg vcfeval --threads 1  \
    -b /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/validation/giab-NA12878/truth_small_variants.vcf.gz \
    --bed-regions /hpf/largeprojects/ccmbio/naumenko/validation/lynette/work/validate/NA12878-1/ensemble/NA12878-1-sort-callable_sample-NA12878-1-ensemble-wrm.bed \
    -c $1 \
    -t /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/rtg/GRCh37.sdf -o rtg --vcf-score-field='GQ'

for f in {tp,fp,fn};
do
    echo snp $f `bcftools view --types snps rtg/$f.vcf.gz | grep -vc "^#"` >> $1.stat
    echo indels $f `bcftools view --exclude-types snps rtg/$f.vcf.gz | grep -vc "^#"` >> $1.stat
done
