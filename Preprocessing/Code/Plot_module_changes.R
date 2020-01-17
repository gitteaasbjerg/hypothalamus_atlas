library("ggalluvial")
library('ggplot2')
library('tidyverse')

matrix = read.csv('/raid5/home/cbmr/kzd307/gitte/hypothalamus_atlas/new_data/matrix.csv', header = TRUE)
colnames(matrix)[1]  <- "genes"
arrange(matrix,big1)

keep_modules <- unique(matrix$big9)[2:11]
keep_modules <- matrix$big9 %in% keep_modules
keep_modules <- which(keep_modules)

new_matrix <- matrix[keep_modules,]

ggplot(data = new_matrix,
       aes(axis1 = big1, axis2 = big3, axis3 = big5, axis4 = big7, axis5=big9,
           y = genes)) +
  scale_x_discrete(limits = c("big1", "big2", "big3","big5","big8"), expand = c(.1, .05)) +
  xlab("") +
  geom_alluvium(aes(fill = big8),aes.bind=TRUE) +
  geom_stratum() +
  #geom_stratum() + geom_text(stat = "stratum", label.strata = TRUE) +
  theme_minimal() +
  ggtitle("")
ggsave('/raid5/home/cbmr/kzd307/gitte/hypothalamus_atlas/new_data/diagram1.png', height = 21, width = 14)



