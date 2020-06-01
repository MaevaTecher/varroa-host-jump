library(vcfR)
library(tidyverse)

a = read.vcfR("1_denovos.vcf.gz")

dnl = extract_info_tidy(a, info_fields="DNL")

r = rle(dnl$DNL)

out = tibble(DNL=r$values,Count=r$lengths)

write_csv(out,"stats_rle_dnl.csv")