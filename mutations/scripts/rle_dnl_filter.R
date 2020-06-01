library(vcfR)
library(tidyverse)

a = read.vcfR("2_denovos.libs.vcf.gz")

db = vcfR2tidy(a,info_fields="DNL",info_only=TRUE)

r = rle(db$fix$DNL)

r = rep(r$lengths,r$lengths)

d = db$fix[r <= 6,]

g = d %>% mutate(POS0 = POS-1L) %>% select(CHROM, POS0, POS)

write_tsv(g,"2_deruns.bed",col_names = FALSE)
