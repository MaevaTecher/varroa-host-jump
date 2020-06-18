library(vcfR)
library(tidyverse)

a = read.vcfR("0_variants.vcf.gz")

gt = extract_gt_tidy(a, format_fields="DP",alleles=FALSE)

out = gt %>% mutate(gt_DP = replace_na(gt_DP,0)) %>% group_by(Indiv) %>% summarise(DP=mean(gt_DP))

write_csv(out,"stats_avg_dp.csv")
