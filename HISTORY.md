## 0.0.4 (23 May 2019)
- fixed a bug in clinvar column

## 0.0.3 (21 May 2019)
- bug fixes, WES report generated for NA12878 looks ok

## 0.0.1 (15 May 2019)
- initial release, compatible with bcbio_1.1.5 and gnomad2.1

## previous notes
- 2018-11-07: new gnomad obseverved/expected scores instead of pLi and exac_missense. cre.gnomad_scores.R
- 2018-11-02: added back Gerp_score, updated OMIM
- 2017-09-22: added Info_refseq and Maf_exac to the database report
- 2017-09-14: improved cre.database.sh: it creates databases for cre.R
- 2017-07-13: added cre.package.sh. It packages reports to send.
- 2017-07-07: added bam.validate.sh, it produces bam.check file with picard's validation status.
- 2017-07-04: bcbio uses vep-merged, a merged cache of refseq and vep, no need to annotate with refseq anymore, removed refseq from cre.sh.
