# Extract sample nr and patient IDs from the IO data
Samples <- IO_data$Sample
Patient_IDs <- Clinical_data_IO_panel$PatientID

# Filter the CVD and Inflammation data to include only the samples present in the IO data
CVD_and_Inflammation_IO_panel <- CVD_and_Inflammation_data %>%
  filter(Sample %in% Samples)
# Filter the clinical data from CVD and Inflammation to include only the patient IDs present in the IO data
Clinical_data_CVD_Inflammation_IO_panel <- Clinical_data_CVD_and_Inflammation %>%
  filter(PatientID %in% Patient_IDs)

# Get the age groups, samples with age groups and database with age groups from for the CVD and Inflammation data
age_group_CVD_results <- get_age_groups(Clinical_data_CVD_Inflammation_IO_panel)
sample_with_age_group_CVD <- age_group_CVD_results$sample_with_age_group
database_with_age_groups_CVD <- age_group_CVD_results$database_with_age_groups

# Get the age groups, samples with age groups and database with age groups from for the IO data
age_group_IO_results <- get_age_groups(Clinical_data_IO_panel)
sample_with_age_group_IO <- age_group_IO_results$sample_with_age_group
database_with_age_groups_IO <- age_group_IO_results$database_with_age_groups

# Filter the data using the function from Analysis.R
filtered_data_CVD_result <- get_filtered_data(CVD_and_Inflammation_IO_panel, "CVD_and_Inflammation_data", sample_with_age_group_CVD)
filtered_data_IO_result <- get_filtered_data(IO_data, "IO_data", sample_with_age_group_IO)

# Extract filtered data
filtered_data_CVD <- filtered_data_CVD_result$data_filtered
filtered_data_IO <- filtered_data_IO_result$data_filtered

# Extract protein data
protein_data_CVD <- filtered_data_CVD_result$protein_data_filt
protein_data_IO <- filtered_data_IO_result$protein_data_filt

# Perform PCA on the protein data
pca_CVD <- get_PCA_data(protein_data_CVD, sample_with_age_group_CVD)
pca_IO <- get_PCA_data(protein_data_IO, sample_with_age_group_IO)

# Create PCA plots and show these plots
pca_plot_CVD <- get_PCA_plot(pca_CVD, filtered_data_CVD)
pca_plot_CVD
pca_plot_IO <- get_PCA_plot(pca_IO, filtered_data_IO)
pca_plot_IO

# Perform cluster PCA on the protein data
cluster_pca_CVD <- get_PCA_kmeans(pca_CVD, centers = 3, sample_with_age_group = sample_with_age_group_CVD)
cluster_pca_IO <- get_PCA_kmeans(pca_IO, centers = 3, sample_with_age_group = sample_with_age_group_IO)

# Create cluster PCA plots and show these plots
cluster_pca_plot_CVD <- get_cluster_pca_plot(cluster_pca_CVD)
cluster_pca_plot_CVD
cluster_pca_plot_IO <- get_cluster_pca_plot(cluster_pca_IO)
cluster_pca_plot_IO

# Extract the sample nr with the samples for both datasets
sample_with_clusters_CVD <- data.frame(Sample = cluster_pca_CVD$data$sample, Cluster_CVD = cluster_pca_CVD$data$cluster)
sample_with_clusters_IO <- data.frame(Sample = cluster_pca_IO$data$sample, Cluster_IO = cluster_pca_IO$data$cluster)

# Merge the sample with cluster datasets
merged_sample_with_clusters <- sample_with_clusters_CVD %>%
  left_join(sample_with_clusters_IO, by = "Sample")

# Create a confusion matrix to compare clusters
confusion_matrix <- table(merged_sample_with_clusters$Cluster_IO, merged_sample_with_clusters$Cluster_CVD)
confusion_matrix_df <- as.data.frame(confusion_matrix)
confusion_matrix_df <- spread(confusion_matrix_df, Var2, Freq)

# Set row names and remove the first column
rownames(confusion_matrix_df) <- confusion_matrix_df$Var1
confusion_matrix_df <- confusion_matrix_df[, -1]
confusion_matrix_df <- cbind("IO panel" = rownames(confusion_matrix_df), confusion_matrix_df)

# Create a gt table for the confusion matrix
gt_table <- gt(confusion_matrix_df) %>%
  tab_header(
    title = "Confusion Matrix",
    subtitle = "Clusters Comparison"
  ) %>%
  tab_spanner(
    label = "CVD & Inflammation panel",
    columns = c("1", "2", "3")
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 0
  )

# Show the gt table
gt_table
