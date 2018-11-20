# http://bioconductor.org/packages/release/data/annotation/html/MafDb.TOPMed.freeze5.hg19.html 
# https://bioconductor.org/packages/release/data/annotation/html/SNPlocs.Hsapiens.dbSNP144.GRCh37.html
# http://gnomad.broadinstitute.org/
# https://bravo.sph.umich.edu/freeze3a/hg19/

install = function()
{
    source("https://bioconductor.org/biocLite.R")
    biocLite("MafDb.TOPMed.freeze5.hg19")
    biocLite("SNPlocs.Hsapiens.dbSNP144.GRCh37")
}

init = function()
{
    library("SNPlocs.Hsapiens.dbSNP144.GRCh37")
    library("MafDb.TOPMed.freeze5.hg19")
}

af_by_rsid = function(rs_id)
{
    # DEBUG:
    # from Johar 2016
    rs_id = "rs202193903"
    # rs_id = "rs903331232"
    #- not found
    #ls("package:MafDb.TOPMed.freeze5.hg19")
    mafdb = MafDb.TOPMed.freeze5.hg19
    #mafdb
    #citation(mafdb)
    #populations(mafdb)
    ## lookup allele frequencies for rs1129038, a SNP associated to blue and brown eye colors
    ## as reported by Eiberg et al. Blue eye color in humans may be caused by a perfectly associated
    ## founder mutation in a regulatory element located within the HERC2 gene inhibiting OCA2 expression.
    ## Human Genetics, 123(2):177-87, 2008 [http://www.ncbi.nlm.nih.gov/pubmed/18172690]
    snpdb = SNPlocs.Hsapiens.dbSNP144.GRCh37
    #print(rs_id)
    rng = snpsById(snpdb, ids=rs_id,ifnotfound ="drop")
    #rng
    seqlevelsStyle(rng) = seqlevelsStyle(mafdb)
    if (length(rng) > 0){
        scores = gscores(mafdb, rng)
        result = scores$AF
    }else{
        result = NA
    }
    return(result)
    
    #gscores(mafdb, GRanges("chr15:28356859"))
}

analysis = function()
{
    setwd("~/Desktop/work")
    hgmd = read.csv("hgmd.csv",header = T,stringsAsFactors = F)
    
    print("Total variants")
    print(nrow(hgmd))
    
    hgmd$AF = NULL
    
    hgmd = hgmd[hgmd$tag == "DM" | hgmd$tag == "DM?",]
    for (i in 1:nrow(hgmd))
    {
        print(i)
        hgmd[i,"AF"] = af_by_rsid(hgmd[i,"dbsnp"])
    }
    
    write.csv(hgmd,"hgmd_with_topmed.csv",row.names = F)
    
    hgmd = hgmd[!is.na(hgmd$AF),]
    
    print("With AF:")
    print(nrow(hgmd))
    
    hgmd = hgmd[hgmd$AF>0.01,]
    
    print("AF > 1%")
    print(nrow(hgmd))
}
