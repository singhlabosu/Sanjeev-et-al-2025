# Prepare variables
data <- SMGnorm
data <- UPF3norm
data <- E117RvsWT
group_var <- "quartile"
value_var <- "baseMean"
value_var <- "log2foldchange"
value_var <- "UTR3lengthfromcds"
value_var <- "baseMeanHCT"
value_var <- "log2FoldChange"

# Get all pairwise combinations
group_levels <- unique(data[[group_var]])
pairs <- combn(group_levels, 2, simplify = FALSE)

# Perform tests and store results
results <- lapply(pairs, function(pair) {
  group1 <- pair[1]
  group2 <- pair[2]
  
  x <- data[[value_var]][data[[group_var]] == group1]
  y <- data[[value_var]][data[[group_var]] == group2]
  
  test <- wilcox.test(x, y, exact = FALSE)
  
  data.frame(
    group1 = group1,
    group2 = group2,
    p_value = test$p.value,
    W_statistic = test$statistic,
    alternative = test$alternative,
    method = test$method,
    stringsAsFactors = FALSE
  )
})

# Combine into a dataframe
results_df <- bind_rows(results)

# Adjust p-values (Bonferroni, same as in your test)
results_df <- results_df %>%
  mutate(p_adjusted = p.adjust(p_value, method = "bonferroni"),
         hypothesis = paste(group1, "vs", group2))
results_df
write_clip(results_df)
