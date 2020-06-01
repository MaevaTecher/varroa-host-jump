library(vcfR)
library(tidyverse)
library(DescTools)
library(ADGofTest)

genome = read_tsv("contigs.txt",col_names=FALSE) %>% 
    rename(CHROM = X1, LEN = X2) %>% mutate(PROP = LEN/sum(LEN))
w = setNames(genome$PROP,genome$CHROM)
glen = sum(genome$LEN)

sons = c("Son1","Son15","Son17","Son19","Son2","Son7")

gtab = genome %>% select(CHROM,POS=LEN)
gtab = gtab %>% mutate(POS=0) %>% bind_rows(.,gtab) %>% arrange(CHROM,POS)
gtab = bind_cols(gtab[rep(1:nrow(gtab),times=length(sons)),],Indiv=rep(sons,each=nrow(gtab)))

a0 = read.vcfR("0_variants.vcf.gz")

gt = extract_gt_tidy(a0, format_fields="DP",alleles=FALSE)

tab = gt %>% mutate(gt_DP = replace_na(gt_DP,0)) %>%
    group_by(Indiv) %>% summarise(AvgDP=mean(gt_DP))

a2 = read.vcfR("2_denovos.libs.vcf.gz")

y = vcfR2tidy(a2,info_fields="DNL",info_only=TRUE)$fix %>%
    mutate(Indiv = str_replace(DNL,"^LB/",""))

y = y %>% select(CHROM,POS,Indiv)

gg = function(x,p) {
    if(length(x) == 1) {
        return(1)
    }
    GTest(x,p=p)$p.value
}

x = y %>% group_by(Indiv,CHROM) %>% summarise(n=n()) %>% group_by(Indiv) %>% summarise(GTest=gg(n,w[CHROM]))
tab = left_join(tab,x)

yy = bind_rows(y,gtab) %>% arrange(CHROM,POS)

space = yy %>% group_by(Indiv,CHROM) %>% mutate(Dist = c(NA,POS[-1] - POS[-length(POS)])) %>% filter(POS > 0)

tab = space %>% group_by(Indiv) %>% summarise(MeanSpacing = mean(Dist)) %>% left_join(tab,.)
tab = space %>% group_by(Indiv) %>% summarise(MedianSpacing = median(Dist)) %>% left_join(tab,.)
tab = space %>% group_by(Indiv) %>% summarise(AdTest = ad.test(Dist,pexp,rate=length(Dist)/glen)$p.value) %>% left_join(tab,.)

dnl = y$Indiv
x = table(dnl)
tab$Count2 = replace_na(x[tab$Indiv],0)

a3 = read.vcfR("3_deduped.vcf.gz")
dnl = extract_info_tidy(a3, info_fields="DNL") %>% mutate(DNL = str_replace(DNL,"^LB/",""))
x = table(dnl$DNL)
tab$Count3 = replace_na(x[tab$Indiv],0)

a4 = read.vcfR("4_filtered.vcf.gz")
dnl = extract_info_tidy(a4, info_fields="DNL") %>% mutate(DNL = str_replace(DNL,"^LB/",""))
x = table(dnl$DNL)
tab$Count4 = replace_na(x[tab$Indiv],0)


write_csv(tab,"stats.csv")
