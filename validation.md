# Latest NA12878 WES and WGS
<table>
    <tr><td rowspan="2">Data</td><td rowspan="2">Tool</td><td colspan="2">SNP</td><td colspan="2">Indel</td></tr>
    <tr><td>FDR</td><td>FNR</td><td>FDR</td><td>FNR</td></tr>
    <tr><td>WES</td><td>Ensembl 2of4 (bcbio 1.1.2)</td><td>0.22%</td><td>0.23%</td><td>7.41%</td><td>5.98%</td></tr>
    <tr><td>WGS</td><td>Gatk 4.0.12.0 (bcbio)</td><td>0.10%</td><td>0.12%</td><td>0.95%</td><td>1.19%</td></tr>
</table>

# WES - 2019-01-23

**Data:**

|Sample name|NA12878|
|------------------|-------------|
|Cycles|2x100|
|Avg coverage for Ensembl protein coding exons|156|
|Median bp coverage for Ensembl protein coding exons|162|
|Reads (single)|142,641,214|
|Insert size|190|

**Methods:**
I intersect callable regions from bcbio x giab bed file x capture kit bed file = 50,436,223 bp in my case. It includes some non-coding sites, coding only should be 30M. Then I use RTG to report baseline-TP, FP, FN. 
My bcbio config:
```
details:
- algorithm:
    aligner: bwa
    effects: false
    ensemble:
      numpass: 2
      use_filtered: false
    mark_duplicates: true
    realign: false
    recalibrate: false
    remove_lcr: false
    save_diskspace: true
    tools_off:
    - vqsr
    tools_on:
    - noalt_calling
    variantcaller:
    - gatk-haplotype
    - freebayes
    - platypus
    - samtools
    validate: giab-NA12878/truth_small_variants.vcf.gz
    validate_regions: giab-NA12878/truth_regions.bed
```

**Results:**

**SNPs:**

variant caller | gatk  4.0.12.0 | gatk  4.0.1.2 | freebayes  1.1.0.46 | ensemble
-- | -- | -- | -- | --
tp | 38,614 | 38,627 | 38,649 | 38,625
fp | 203 | 211 | 285 | 266
fn | 129 | 116 | 94 | 91
FDR | 0.52% | 0.54% | 0.73% | 0.68%
FNR | 0.33% | 0.30% | 0.24% | 0.23%
Target | 38743 | 38743 | 38743 | 38743

**Indels:**

variant caller | gatk  4.0.12.0 | gatk  4.0.1.2 | freebayes  1.1.0.46 | ensemble
-- | -- | -- | -- | --
tp | 2874 | 2883 | 2776 | 2855
fp | 364 | 384 | 238 | 228
fn | 170 | 161 | 268 | 189
FDR | 11.24% | 11.75% | 7.90% | 7.40%
FNR | 5.58% | 5.29% | 8.80% | 6.21%
