#!/bin/bash

# HLH gene panel for MAS diagnostics

gemini query --header -q "select s.name as sample,v.chrom,v.start+1 as pos,v.ref,v.alt,v.gene,v.rs_ids,v.impact,v.max_aaf_all 
			  from variants v, samples s
			  where v.max_aaf_all <= 0.05 and 
			  v.gene in ('AP3B1','BLOC1S6','CD27','GATA2','ITK','LYST','NLRC4','PRF1','RAB27A','SH2D1A','SLC7A7','STX11','STXBP2','UNC13D','XIAP')" $1
