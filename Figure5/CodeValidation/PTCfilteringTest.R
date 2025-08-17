library(dplyr)

# Create a mock dataframe
set.seed(123)  # For reproducibility
df <- data.frame(
  ensembl_gene_id = rep(c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"), each = 3),
  transcript_id = paste0("TX", 1:15),
  PTC.Status = c(TRUE, FALSE, TRUE,  # GENE1: Has both TRUE & FALSE
                 TRUE, TRUE, TRUE,   # GENE2: Only TRUE
                 FALSE, FALSE, FALSE, # GENE3: Only FALSE
                 TRUE, FALSE, TRUE,  # GENE4: Has both TRUE & FALSE
                 TRUE, TRUE, FALSE)  # GENE5: Has both TRUE & FALSE
)

# View the mock dataframe
print(df)



filtered_df <- df %>%
  group_by(ensembl_gene_id) %>%
  filter(any(PTC.Status == TRUE) & any(PTC.Status == FALSE)) %>%
  ungroup()

# View the filtered result
print(filtered_df)
