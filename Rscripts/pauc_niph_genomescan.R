## Fst scans all chromosomes

library(plyr)
library(ggplot2)

# First run:
scans = ldply(list.files(pattern = "_biasnps_mis80_maf01_mDP5-25_fst_niph_pauc.weir.fst"), read.delim, header=TRUE)
# Second run:
scans = ldply(list.files(pattern = "biasnps_DP5-25_mac4_mis70_fst_niph_pauc_5kb.windowed.weir.fst"), read.delim, header=TRUE) %>%
  filter(N_VARIANTS > 1)
# Third run:
scans = ldply(list.files(pattern = "biasnps_DP5-25_mac4_mis70_fst_niph33_pauc233_5kb.windowed.weir.fst"), read.delim, header=TRUE) %>%
  filter(N_VARIANTS > 5)
# fourth run:
scans = ldply(list.files(pattern = "_biasnps_DP5-25_mac4_mis70.windowed.weir.fst"), read.delim, header=TRUE) %>%
  filter(N_VARIANTS > 5)

hist(scans$N_VARIANTS,
     breaks = 100)

table(scans$CHROM)

p = ggplot(scans, aes(x = BIN_START, y = WEIGHTED_FST)) +
  geom_point(size = 0.1) + 
  ylim(0,max(scans$MEAN_FST)) +
  facet_wrap(~CHROM, nrow = 1, scales = 'free_x') + 
  labs(x = "Chr", y = "FST")
p

### scans without windowing:

# Third run:
scans = ldply(list.files(pattern = "_biasnps_DP5-25_mac4_mis70_fst_niph33_pauc233.weir.fst"), read.delim, header=TRUE)
scans$pos_mb = (as.numeric(scans$POS))/1000000
scans$CHROM <- as.character(scans$CHROM)

hist(scans$WEIR_AND_COCKERHAM_FST,
     xlim = c(0,1),
     breaks = 10)

p2 = ggplot(scans, aes(x = pos_mb, y = WEIR_AND_COCKERHAM_FST)) +
  geom_point(size = 0.1) + 
  ylim(0,1) +
  facet_wrap(~CHROM, nrow = 1, scales = 'free_x') + 
  labs(x = "Genomic position (Mbp)", y = "FST") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(path = figstub, 
       filename = "biasnps_DP5-25_mac4_mis70_fst_niph33_pauc233.weir.fst.pdf",
       p2, 
       width = 8, height = 3)

detach("package:plyr", unload = TRUE)

scans %>% 
  group_by(CHROM) %>%
  dplyr::summarise(mean_fst = mean(WEIR_AND_COCKERHAM_FST, na.rm = T),
                   median_fst = median(WEIR_AND_COCKERHAM_FST, na.rm = T))






