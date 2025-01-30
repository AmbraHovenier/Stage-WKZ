y_columns <- colnames(CVD_Inflammation_and_IO_data)[grepl("\\.y$", colnames(CVD_Inflammation_and_IO_data))]
x_columns <- colnames(CVD_Inflammation_and_IO_data)[grepl("\\.x$", colnames(CVD_Inflammation_and_IO_data))]
x_and_y_proteins <- CVD_Inflammation_and_IO_data %>%
  select(all_of(x_columns), all_of(y_columns), Sample) %>%
  mutate(Sample = gsub("A", "", Sample)) %>%
  mutate(Sample = gsub("^(\\D*)0+", "\\1", Sample))  %>%
  mutate(across(everything(), as.numeric))

protein_names <- gsub(".y", "", y_columns)
print(protein_names)

# Loop for each column and calculate R2 squared
for(protein_name in protein_names) {
  x_protein <- paste0(protein_name, ".x")
  y_protein <- paste0(protein_name, ".y")
  r_squared <- r_squared_function(y_protein, x_protein, x_and_y_proteins)
  print(paste(protein_name, r_squared))
}

# Function to calculate the R2 value
r_squared_function <- function(y_protein, x_protein, x_and_y_proteins) {
  formula <- as.formula(paste0("`", y_protein, "` ~ `", x_protein, "`"))
  model <- lm(formula, data = x_and_y_proteins)
  r_squared <- summary(model)$r.squared
  if(r_squared < 1) {
    plot_function(y_protein, x_protein, x_and_y_proteins, r_squared)
  }
  return(r_squared)
}

# Function to create the plot
plot_function <- function(y_protein, x_protein, x_and_y_proteins, r_squared) {
  xprotein <- gsub(".x", " (CVD)", x_protein)
  yprotein <- gsub(".y", " (Inflammation)", y_protein)
  protein <- gsub(".x", "", x_protein)
  plot <- ggplot(
    x_and_y_proteins,
    aes(x = !!sym(x_protein), y = !!sym(y_protein))) + 
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.5) +  # Add regression line
    labs(title = paste(protein, ": CVD panel vs Inflammation panel with R²"),
         x = xprotein, y = yprotein) +
    annotate("text", x = Inf, y = -Inf, label = paste("R² =", round(r_squared, 3)),
             hjust = 1.1, vjust = -1.1, size = 5, color = "red")
  print(plot)
}
