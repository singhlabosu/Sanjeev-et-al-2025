data <- PYMkdDESeq
group_var <- "Class1"
value_var <- "log2FoldChange"

data <- StabilityTable
group_var <- "Compartment"
value_var <- "RelativeStability"

data <- E117RvsWT
group_var <- "Class1"
value_var <- "log2FoldChange"

data <- SMG67DESeq
group_var <- "Class1"
value_var <- "log2FoldChange"

data <- DESeqFritz
group_var <- "Class1"
value_var <- "log2FoldChange"


data <- StabilityTable
group_var <- "Compartment"
value_var <- "RelativeStability"


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






#ERTG slcn######

# Initialize results dataframe
comparison_results <- data.frame(
  comparison = character(),
  test = character(),
  p_value = numeric(),
  statistic = numeric(),
  n_PTCpos = integer(),
  n_PTCneg = integer(),
  stringsAsFactors = FALSE
)

# Loop over each unique comparison group
for (comp in unique(E117RvsWTslcn$Region)) {
  print(comp)
  subdf <- E117RvsWTslcn %>% filter(Region == comp)
  # Grouped values
  ER_plus  <- subdf %>% filter(Class1 == "ER+") %>% pull(log2FoldChange)
  TG_plus <- subdf %>% filter(Class1 == "TG+") %>% pull(log2FoldChange)
  # Perform Wilcoxon test (non-parametric)
  test <- wilcox.test(ER_plus, TG_plus, alternative = "two.sided")
  # Append results
  comparison_results <- rbind(comparison_results, data.frame(
    comparison = comp,
    test = "Wilcoxon",
    p_value = test$p.value,
    statistic = unname(test$statistic),
    n_PTCpos = length(ER_plus),
    n_PTCneg = length(TG_plus)))
}

# Apply Bonferroni correction
comparison_results$p_adj <- p.adjust(comparison_results$p_value, method = "bonferroni")

# View results
print(comparison_results)

#Compartments PYM1kdFoldchange######

# Initialize results dataframe
comparison_results <- data.frame(
  Group = character(),
  Comparison = character(),
  n_Ann_yes = integer(),
  n_Ann_no = integer(),
  test = character(),
  p_value = numeric(),
  W_statistic = numeric(),
  Alternate_Hypothesis = character(),
  stringsAsFactors = FALSE
)

# Loop over each unique comparison group
for (comp in unique(PYMkdDESeqlong$Classification)) {
  print(comp)
  subdf <- PYMkdDESeqlong %>% filter(Classification == comp)
  # Grouped values
  Ann_yes  <- subdf %>% filter(Annotation == "Yes") %>% pull(log2FoldChange)
  Ann_no <- subdf %>% filter(Annotation == "No") %>% pull(log2FoldChange)
  # Perform Wilcoxon test (non-parametric)
  test <- wilcox.test(Ann_yes, Ann_no)
  # Append results
  comparison_results <- rbind(comparison_results, data.frame(
    Group = comp,
    Comparison = "Annotation:Yes vs No",
    n_Ann_yes = length(Ann_yes),
    n_Ann_no = length(Ann_no),
    test = "Wilcoxon rank sum test with continuity correction",
    p_value = test$p.value,
    W_statistic = unname(test$statistic),
    Alternate_Hypothesis = "true location shift is not equal to 0"
    ))
}

# Apply Bonferroni correction
comparison_results$p_adj_bonferroni <- p.adjust(comparison_results$p_value, method = "bonferroni")

# View results
print(comparison_results)
write_clip(comparison_results, object_type = "auto" )
write.table(comparison_results, "clipboard", sep="\t", row.names=T)

