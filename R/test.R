rm(list = ls())


library(gdsfmt)
library(SNPRelate)
library(stringi)
library(tidyverse)


setwd('/user/doau0129/work/chapter1_2/outputs/09_1_snpfiltration/')

vcf_file <- 'filterd_bi-allelic_changed_casz1.vcf.gz'

ld_threshold <- 1

out_prefix <- 'test'
gds_file <- str_c(out_prefix,"_casz1.gds")

snpgdsVCF2GDS(vcf.fn=vcf_file, out.fn=gds_file, method="biallelic.only")
genofile <- snpgdsOpen(gds_file)
snpset <- snpgdsLDpruning(genofile, ld.threshold = ld_threshold, autosome.only = FALSE)
snpset.id <- unlist(snpset)

pca <- snpgdsPCA(genofile, snp.id = snpset.id, num.thread = 2, autosome.only = FALSE)

pca_tib <- pca$eigenvect %>%
  as_tibble() %>%
  set_names(nm = str_pad(1:length(names(.)),2,pad = '0') %>%
              str_c("PC",.)) %>%
  mutate(sample = pca$sample.id,
         spec = str_sub(sample,-6,-4),
         loc = str_sub(sample,-3,-1),
         pop = str_c(spec,loc))

explained_var <- tibble(variation = pca$varprop*100) %>%
  mutate(PC = str_c("PC", str_pad(row_number(),
                                  width = 2,
                                  pad = "0"))) %>%
  select(PC, variation) %>%
  filter(!is.na(variation))

save(pca,file = str_c(out_prefix,"_casz1_pca.RData"))
pca_tib %>% write_tsv(path = str_c(out_prefix,"_casz1_eigenvectors.tsv"))
explained_var %>% write_tsv(path = str_c(out_prefix,"_casz1_eixplained_var.tsv"))


png("casz1_pca.png")
plot(pca$eigenvect[,1],pca$eigenvect[,2] ,col=as.numeric(substr(pca$sample, 1,3) == 'CCM')+3, pch=2)
dev.off()


snpgdsClose(genofile)
system(str_c("rm ", gds_file))
showfile.gds(closeall=TRUE)
