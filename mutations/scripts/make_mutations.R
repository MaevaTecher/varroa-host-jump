library(tidyverse)
library(vcfR)

acgt = c("A","C","G","T","A","C","G","T")
sons = c("Son1","Son15","Son17","Son19","Son2","Son7")

b = read_tsv("random_sites.bed",col_names=c("CHROM","POS0","POS","ID","W","STRAND"))
b = b %>% select("CHROM","POS","ID") %>% arrange(ID) %>% mutate(DNL = rep_len(sons, length.out=nrow(b)), ID=NULL)

a0 = read.vcfR("recall_0_genotype.vcf.gz")
CHROM = getCHROM(a0)
POS = getPOS(a0)
N = tibble(CHROM=CHROM,POS=POS,N=1:length(CHROM))

gt_ad = NULL

for(dnl in sons) {
    o = b %>% filter(DNL == dnl) %>% inner_join(.,N) %>% pull(N) %>% sort()
    v = a0[o,] %>% vcfR2tidy()
    gt = v$gt %>% mutate(Indiv = str_replace(Indiv,"^LB/","")) %>% filter(Indiv == dnl)

    # n is original; m is mutant
    n = gt$gt_GT_alleles
    u = sample(c(1,2,3),nrow(gt),replace=TRUE)
    m = acgt[match(n, acgt) + u]

    refalt = ifelse(is.na(v$fix$ALT), v$fix$REF, str_c(v$fix$REF, v$fix$ALT,sep=",")) %>% str_split(",")
    ad = gt$gt_AD %>% str_split(",")
    AD = v$gt %>% mutate(Indiv = str_replace(Indiv,"^LB/","")) %>% select(ChromKey, POS, Indiv,gt_AD) %>%
        spread(Indiv, gt_AD) %>% select(-ChromKey, -POS)

    for(j in 1:length(m)) {
        x = ad[[j]]
        pm = match(m[j],refalt[[j]],nomatch=0)
        pn = match(n[j],refalt[[j]],nomatch=0)
        if(pm == 0) {
            refalt[[j]] = c(refalt[[j]],m[j])
            AD[j,] = AD[j,] %>% str_c(",0")
            x = c(x,"0")
            pm = length(x)
        }
        ad[[j]][pn] = x[pm]
        ad[[j]][pm] = x[pn]
    }

    ad = sapply(ad,str_c,collapse=",")
    v$fix$ALT = refalt %>% lapply("[", -1) %>% sapply(str_c,collapse=",")
    AD[[dnl]] = ad

    records = v$fix %>% select(CHROM,POS,REF,ALT) %>% mutate(DNL=dnl)
    records = bind_cols(records, AD)
    gt_ad = bind_rows(gt_ad,records)
}
gt_ad = gt_ad %>% arrange(CHROM,POS)

a = a0

a@fix[,"ALT"] = gt_ad$ALT
a@fix[,"INFO"] = str_c("DNL=",gt_ad$DNL)
a@gt[,1] = "AD"
a@gt[,-1] = as.matrix(gt_ad[,-(1:5)])
colnames(a@gt) = colnames(a@gt) %>% str_replace("LB/","")

write.vcf(a,file='recall_1_mutated.vcf.gz')
