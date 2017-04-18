# cre
clinical research exome - excel report generation using results from [bcbio variant2](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#germline-variant-calling) 
germline variant calling pipeline.


# 1. Create a project (projects) to run with bcbio.

## 1a. If you start with bam files.
Suppose you have a WES trio, or a cohort of trios, each sample is a bam file. 
Then create a file table.txt where each line is sample_name<tab>family_name<tab>file.bam, i.e.
```
CH0517	1004	/somewhere/1004_CH0517.bam
CH0518	1004	/somewhere/1004_CH0518.bam
CH0519	1004	/somewhere/1004_CH0519.bam
CH0133	1013	/somewhere/1013_CH0133.bam
CH0135	1013	/somewhere/1013_CH0135.bam
CH0136	1013	/somewhere/1013_CH0136.bam
```

Then run 
```
[bcbio.prepare_families.sh](../master/bcbio.prepare_families.sh) [table.txt] 
```
or use qsub, if you have a large cohort:
```
bcbio.prepare_families.sh	-v project_list=table.txt
```

The script will prepare a folder for each project(family) using (bcbio.templates.exome.yaml)[../master/bcbio.templates.exome.yaml] template.
Using this template is important, because later the final report will be produced from 4 callers.

The details of the template:
* 4 callers, the order is important, because variant metrics in ensemble calling (like AD) are picked up from the first caller in the list first
* ensemble calling
* realignment and recalibration. I know that there is no much sense in it for precision/sensitivity, but people are still asking for realignment and recalibration.
* no bed file. Let callers call every variant which has coverage, we will filter poorly covered variants later. Modern exome capture kits are so perfect, that
we can discover a useful non-coding variant. No sense to filter them out during that stage.
* effects: VEP. There is a holywar VEP/snpEff/Annovar. My choice is VEP. You will see later both Ensembl and Refseq in the report, so no reason for using Annovar.
* effects_transcripts: all. We want all effects of a variant on all transcripts to be reported.
* aligner: bwa. Even staring with bam files, bwa is used. Sometimes input bam files aligned against older reference, or different (chr) naming scheme. It is better to have a bam file consistent with calls made.

# 2. Run bcbio
