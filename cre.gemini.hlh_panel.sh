#!/bin/bash

# HLH gene panel for MAS diagnostics

gemini query --header -q "select chrom,end,ref,alt,gene,rs_ids,impact,max_aaf_all 
			  from variants 
			  where max_aaf_all <= 0.05 and 
			  gene in ('AP3B1','BLOC1S6','CD27','GATA2','ITK','LYST','NLRC4','PRF1','RAB27A','SH2D1A','SLC7A7','STX11','STXBP2','UNC13D','XIAP')" $1
