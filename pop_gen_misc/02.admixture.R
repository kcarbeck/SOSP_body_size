# Plot structure from ADMIXTURE 
# April 5, 2021
# Author: Katherine Carbeck

library(tidyverse)

samplelist <- read_delim("plink.list", delim=" ",
                       col_names = c("sample", "population"))

all_data <- tibble(sample=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

for (k in 1:4){
  data <- read_delim(paste0("plink.",k,".Q"),
                     col_names = paste0("Q",seq(1:k)),
                     delim=" ")
  data$sample <- samplelist$sample
  data$k <- k
  
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-k) -> data
  all_data <- rbind(all_data,data)
}
all_data



# Plot
all_data %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Greys",name="K",
                    labels=seq(1:8)) +
  facet_wrap(~k,ncol=1)
