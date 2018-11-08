-- mysql -u root -p hgmd_pro < ~/cre/cre.hgmd2csv.sql
select 
	v.chrom, v.pos, v.id, v.ref, v.alt,
	a.gene, a.tag, a.author, a.allname, a.vol, a.page, a.year, a.pmid,a.dbsnp
into outfile 'hgmd.csv' 
	fields terminated by ',' 
	optionally enclosed by '"' lines terminated by '\n' 
from hgmd_hg19_vcf v, allmut a
where v.id = a.acc_num;