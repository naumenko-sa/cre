for i in `ls *.bam`; do qsub ~/bioscripts/scripts/bam.coverage.sh -v bam=$(pwd)/$i,bed=/hpf/largeprojects/ccmbio/dennis.kao/gene_data/protein_coding_genes.exons.fixed.sorted.bed; done
