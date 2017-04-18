# cre
clinical research exome - excel report generation using data from [bcbio variant2](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#germline-variant-calling) 
germline variant calling pipeline.


# 1. Create a project (projects) to run with bcbio.

## 1a. If you start from bam files.
Suppose you have a WES trio, or a cohort of trios, each sample is a bam file. 
Then create a file table.txt where each row is 
sample_name	family_name	file.bam, i.e.
```
CH0517	1004	/somewhere/1004_CH0517.bam
CH0518	1004	/somewhere/1004_CH0518.bam
CH0519	1004	/somewhere/1004_CH0519.bam
CH0133	1013	/somewhere/1013_CH0133.bam
CH0135	1013	/somewhere/1013_CH0135.bam
CH0136	1013	/somewhere/1013_CH0136.bam
```
