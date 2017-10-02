# cre
Simple clinically oriented research exome - excel report generation using results from [bcbio variant2](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#germline-variant-calling) 
germline variant calling pipeline.

#-1. TODO
* do all additional annotations with vcfanno
* no need in gemini db - now bcbio makes all gemini annotations with vcfanno

# 0. Prerequisites and credits

## 0.1 Prerequisites
* Install **Bcbio** (HPC or server) and add it to the PATH. bcbio installs many other useful tools through bionconda.
* Clone **cre** to ~/cre and add it to the PATH (HPC).
* Install R and packages: stringr,data.table,plyr (HPC or laptop, if you'd like to use it for the report generation).
* Install OMIM (HPC or laptop).
  * Goto https://omim.org/downloads/ and request the latest database. It makes sense to renew it once a year.
  * In a couple of days you will get genemap2.txt,genemap.txt,mim2gene.txt,mimTitles.percent.txt,mimTitles.txt,morbidmap.txt. Put them into OMIM_DIR where you want.
  * Preprocess OMIM with [omim.sh](../master/omim.sh). It creates two tables: omim.txt with omim description of diseases related to 3700 genes, and omim_inheritance.txt - the inheritance mode
  for genes in OMIM.
```
    cd OMIM_DIR
    ~/cre/omim.sh
```
* Install Orphanet (HPC or laptop) with [orphanet.sh](../master/orphanet.sh) ~ 3600 genes:
```
    cd ~/cre
    wget http://www.orphadata.org/data/xml/en_product6.xml
    orphanet.sh
```
* [Optional, for now uses old scores from ~/cre/exac_scores.txt] Install EXaC scores.
* [Optional, for now uses old table from ~/cre/imprinting.txt] Install imprinting annotation.

If you already have bcbio project results, you may start from step 3. However, note that resulting file names
may have changed in bcbio since you had run the project, and cre follows the latest naming scheme, like project-ensemble-annotated-decomposed.vcf.gz.

## 0.2 Credits

This work was inspired by 
* [bcbio](https://github.com/chapmanb/bcbio-nextgen/) and [gemini](https://github.com/arq5x/gemini) communities. Thank you all!
* Kristin Kernohan from Children Hospital of Eastern Ontario (CHEO), who generated most ideas about the report contents. Thank you, Kristin, for all of the discussions!

Thank you colleagues at [CCM](https://ccm.sickkids.ca/), for seminars and personal discussions.

# 1. Create a project (projects) to run with bcbio.

## 1a. If you start with a cohort of samples/families.
Suppose you have a cohort of trios and individual samples, each sample represented by a bam file.
Then create a file table.txt where each line is sample_name[tab]family_name[tab]file.bam, i.e.
```
Ashkenazim_HG002	Ashkenazim	/data/Ashkenazim_HG002.bam
Ashkenazim_HG003	Ashkenazim	/data/Ashkenazim_HG003.bam
Ashkenazim_HG004	Ashkenazim	/data/Ashkenazim_HG004.bam
CH0517	1004	/somewhere/1004_CH0517.bam
CH0518	1004	/somewhere/1004_CH0518.bam
CH0519	1004	/somewhere/1004_CH0519.bam
CH0133	1013	/somewhere/1013_CH0133.bam
CH0135	1013	/somewhere/1013_CH0135.bam
CH0136	1013	/somewhere/1013_CH0136.bam
```

Then run [bcbio.prepare_families.sh](../master/bcbio.prepare_families.sh) [table.txt] 
or use qsub, if you have a large cohort:
```
qsub ~/cre/bcbio.prepare_families.sh	-v project_list=table.txt
```

The script will prepare a folder for each project(family) using [bcbio.templates.exome.yaml](../master/bcbio.templates.exome.yaml) template.
Using this template is important, because later the report will be produced from 4 callers.

The details of the template:
* 4 callers, the order is important, because variant metrics in ensemble calling (like AD) are picked up from the first caller in the list first
* ensemble calling
* realignment and recalibration. I know that there is no much sense in it for precision/sensitivity, but people are still asking for realignment and recalibration.
* no bed file. Let callers call every variant which has coverage, we will filter poorly covered variants later. Modern exome capture kits are so perfect, that
we can discover a useful non-coding variant. No sense to filter them out during that stage.
* effects: VEP. There is a holywar VEP/snpEff/Annovar. My choice is VEP. You will see later both Ensembl and Refseq in the report, so no reason for using Annovar.
* effects_transcripts: all. We want all effects of a variant on all transcripts to be reported.
* aligner: bwa. Even staring with bam files, bwa is used. Sometimes input bam files aligned against older reference, or different (chr) naming scheme. It is better to have a bam file consistent with calls made.

## 1b. If you start with a single sample/family.
Supporse you have fastq files or bam files for a single sample/single family. You may try NIST Ashkenazim trio: ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio (download OsloUniversityHospital exomes).
* Rename fastq files in the format family_sample_1.fq.gz, family_sample_2.fq.gz, or family_sample.bam and place them into family/input folder.
* run [cre.prepare_bcbio_run.sh](../master/cre.prepare_bcbio_run.sh) [family] [no_align]. By default, it uses the template with alignment, if you set $2, it will use another template - without alignment,
calling only.

## 1c. If you start with Illumina basespace.
* use basespace-cli to dump bcl files to HPF.

## 1d. If you start with bcl files.
* create a samples sheet and run [bcl2fq.sh](../master/bcl2fq.sh).

# 2. Run bcbio

It really depends on your HPC job management system and policies. Our HPC uses torque. It is possible to run jobs from bcbio (it automatically submits jobs to the queue,
and then those jobs communicate with each other via network). I've tried this (parallel execution) for some time. I was stuck with some problems, maybe HPC-related.
So I decided to keep it simple, and to run one multicore job for one node per family(project) using torque.

To run one project: [bcbio.pbs](../master/bcbio.pbs):
```
cd [project_name]/work
qsub ~/cre/bcbio.pbs -v project=[project_name],[threads=[number_of_threads]]
```
Project should have a folder project_name in the current directory.

To run many project (N) as job array: [bcbio.array.pbs](../master.bcbio.array.pbs) - requires projects.txt (list of project=list of folders)
in the current directory and folders for projects.
```
qsub -t 1-N ~/cre/bcbio.array.pbs
```
Use a number instead of N, i.e. 100. I'm using 5 cores x 50G of RAM per project. It is HPC-specific. Our policies encourage submitting small jobs.
I can wait for 2-5 days for a project when working with cohorts. Faster processing is possible using more memory and cores, or with bcbio parallel execution.

# 3.clean project directory  and create family.csv report for import to excel
[cre.sh](../master/cre.sh) [family]
or 
```
qsub cre.sh -v family=[family]
```
by default it does both cleanup and make_report steps.
If you want to do one of those steps:
```
qsub cre.sh -v family=[family],cleanup=0,make_report=1
```

If you run a large cohort, some families maybe finished while some are still running.
To delete unnecessary filed, cleanup finished projects:
```
mkdir ready4report
for f in `cat all_dirs.txt`;do if [ -d ${f}/final/ ];then mv $f ready4report/;fi;done;

```
And then run cre for them:
```
cd ready4report
ls | grep -v projects.txt > projects.txt
for f in `cat projects.txt`;do qsub ~/cre/cre.sh -v family=$f;done;
```

What it does.
During the cleanup step:
* moves project results and sample bam files to family dir
* removes work and final dirs from bcbio project
* removes gemini databases for individual callers (we need only ensemble gemini database)

During the report generation step:
* dumps variants from gemini database to tab text file
* dumps variant impacts from gemini database to tab text file
* annotates variants with refseq in addition to ensembl
* gets coverage from GATK Haplotype calls, freebayes, and platypus
* build excel report based on gemini variants table, variant impacts, coverage information and some other fields.

# 4. Step 3 in detail

## 4.0 [Report description](https://docs.google.com/document/d/1zL4QoINtkUd15a0AK4WzxXoTWp2MRcuQ9l_P9-xSlS4/edit?usp=sharing).
## 4.1 [Report example for Ashkenazim trio from NIST](https://drive.google.com/open?id=0B_bLL10GwDnsN29vY3RRdGlXMWM).
## 4.2 [gemini.gemini2txt.sh](../master/gemini.gemini2txt.sh) [project-ensembl.db] - dumps a gemini database into text file.
I prefer to do report generation in R, and first I tried to access gemini database from R using sqlite interface. It turned out impossible, because
of packaging of genotype BLOB fields. I ended up with gemini utility query to dump fields I need from variants database. Filters are described in the doc.
## 4.3 [gemini.variant_impacts.sh](../master/gemini.variant_impacts.sh) [project-ensembl.db] dumps variant impacts from gemini.
## 4.4 creates a vcf file with rare and potentially harmful variants, the same set of variants will be shown in the excel report 
```
cat ${family}-ensemble.db.txt | cut -f 23,24  | sed 1d | sed s/chr// > ${family}-ensemble.db.txt.positions
    bcftools view -R ${family}-ensemble.db.txt.positions -o ${family}.vcf.gz -O z ${family}-ensemble-annotated-decomposed.vcf.gz

```
## 4.5 gets coverage from VCFs produced by GATK, platypus, and freebayes - requires gatk wrapper from bcbio.
```
vcf.freebayes.getAO.sh ${family}-freebayes-annotated-decomposed.vcf.gz
vcf.gatk.get_depth.sh ${family}-gatk-haplotype-annotated-decomposed.vcf.gz
vcf.platypus.getNV.sh ${family}-platypus-annotated-decomposed.vcf.gz

```
## 4.6 Rscript ~/cre/[cre.R](../master/cre.R) [family] - creates report family.csv.

# 5. How to create a database of variants

- cre.database.sh [input_dir] [output_dir] - creates sample-wise and variant-wise reports, which are necessary for annotation with cre.R.
- cre.database.pull_gene.sh [database_prefix] [gene_name] - pull a gene report from the database.

# 6. Coverage plots

- cre.coverage.bamstats05.sh - calculate coverage
- cheo.R - plot coverage pictures

# 7. List of all scripts

* bcbio.array.pbs
* bcbio.pbs
* bcbio.prepare_families.sh
* cre.prepare_bcbio_run.sh
* bcbio.rename_old_names.sh
* bcl2fastq.sh
* cre.bam.validate.sh
* cre.coverage.bamstats05.sh - calculate coverage
* cre.fixit.sh - fixes sample names
* cre.gemini_load.sh loads vep-annotated vcf to gemini db.
* cre.package.sh
* cre.sh
* gemini.gemini2txt.sh
* gemini.refseq.sh
* gemini.variant_impacts.sh
* gemini.vep.refseq.sh
* gemini.vep.sh
* omim.sh
* orphanet.sh
* vcf.freebayes.getAO.sh
* vcf.gatk.get_depth.sh
* vcf.platypus.getNV.sh
* vcf.samtools.get_depth.sh
* vcf.split_multi.sh
* vep4seqr_hg38.sh
* vep4seqr.sh
