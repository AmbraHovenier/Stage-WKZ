####################################################################
##################### FUNCTIONS ####################################
####################################################################

# Load all needed packages
load_all_packages <- function(){
  packages <- c("ggplot2", "ggfortify", "dplyr", "readxl", "haven", "tidyverse", 
                "factoextra", "gt", "NbClust", "fmsb", "RColorBrewer", "GISTools", 
                "plotly", "crayon", "ggrepel")
  
  # Check for missing packages and install them
  missing_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(missing_packages)) install.packages(missing_packages)
  
  # Load all packages and suppress startup messages
  invisible(lapply(packages, function(pkg) suppressPackageStartupMessages(library(pkg, character.only = TRUE))))
}

# Function to perform k-means clustering on PCA scores and add clusters to the PCA object
get_PCA_kmeans <- function(PCAobj, varExplained = 90, centers, sample_with_age_group){
  
  # Set seed for reproducibility
  set.seed(123)
  
  # Determine the number of Principal Components (PCs) that explain the specified percentage of variance
  # and extract the PCA scores for the selected PCs
  PCs <- get_number_of_PCs(PCAobj$pca, varExplained = varExplained)
  scores <- PCAobj$pca$x[,PCs]
  
  # Perform k-means clustering on the selected PCA scores
  km <- kmeans(scores, centers = centers, nstart = 10)
  
  # Add the resulting cluster assignments and sample IDs to the PCA object under the 'data' field
  PCAobj$data$cluster <- km$cluster
  PCAobj$data$sample <- sample_with_age_group$Sample
  PCAobj$data$age_group <- sample_with_age_group$Sampling_age_group
  
  # Return the updated PCA object with cluster and sample ID
  return(PCAobj)
}

# Function to determine the optimal number of clusters for PCA data using various methods
get_optimal_number_of_clusters <- function(PCAobj, varExplained = 90){
  
  # Determine the number of Principal Components (PCs) that explain the specified percentage of variance
  # and extract the PCA scores for the selected PCs
  PCs <- get_number_of_PCs(PCAobj$pca, varExplained)
  scores <- PCAobj$pca$x[,PCs]
  
  # Standardize the PCA scores, by scaling the data
  scaled_data <- scale(scores)
  
  # Generate plots for the gap statistic, silhouette and WSS methods to determine the optimal number of clusters
  number_of_clusters_gap_stat_plot <- fviz_nbclust(scaled_data, FUNcluster = kmeans, method = "gap_stat", nstart = 25, iter.max = 100)
  number_of_clusters_silhouette_plot <- fviz_nbclust(scaled_data, FUNcluster = kmeans, method = "silhouette", nstart = 25, iter.max = 100)
  number_of_clusters_wss_plot <- fviz_nbclust(scaled_data, FUNcluster = kmeans, method = "wss", nstart = 25, iter.max = 100)
  
  # Use the NbClust package to determine the optimal number of clusters using the silhouette and WSS methods
  number_of_clusters_silhouette <- NbClust(scaled_data, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 20, method = "kmeans", index = "silhouette")
  number_of_clusters_wss <- NbClust(scaled_data, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 20, method = "kmeans", index = "ch")
  
  # Return a list containing the results from both NbClust methods and the three cluster evaluation plots
  return(list(number_of_clusters_silhouette = number_of_clusters_silhouette,
              number_of_clusters_wss = number_of_clusters_wss,
              number_of_clusters_gap_stat_plot = number_of_clusters_gap_stat_plot,
              number_of_clusters_silhouette_plot = number_of_clusters_silhouette_plot,
              number_of_clusters_wss_plot = number_of_clusters_wss_plot))
}

# Function to determine the number of Principal Components (PCs) that explain 
# a specified percentage of variance in a PCA object
get_number_of_PCs <- function(PCAobj, varExplained = 90){
  
  # Extract standard deviations of the PCs from the PCA object
  sd <- PCAobj$sdev
  
  # Calculate the variance for each PC
  var <- sd^2
  
  # Calculate the percentage of total variance explained by each PC
  # and calculate the cumulative percentage of variance explained by the PCs
  var.percent <- var/sum(var) * 100
  var.percent.cum <- cumsum(var.percent)
  
  # Return the indices of PCs that explain less than the specified variance
  # and the first PC that exceeds the specified variance
  return(c(which(var.percent.cum < varExplained), 
           which(var.percent.cum > varExplained)[1]))  
}

# Function to perform PCA on protein data
get_PCA_data <- function(protein_data, sample_with_age_group){
  
  # Store the input protein data in a variable for later use
  data <- protein_data
  data <- data %>%
    mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
  
  # Perform PCA on the protein data using the prcomp function
  # This function standardizes the data by default before performing PCA
  pca_data <- prcomp(data)
  
  # Add the sample nr and age group to the data
  data$sample <- sample_with_age_group$Sample
  data$age_group <- sample_with_age_group$Sampling_age_group
  
  # Return a list containing the PCA results and the original data
  # 'pca' holds the PCA output, and 'data' retains the original data
  return(list(pca = pca_data, data = data))
}

# Function to generate a PCA plot using the autoplot function
get_PCA_plot <- function(pca_data, data_filtered){
  
  # 'pca_data$pca' contains the PCA results
  # 'data_filtered' is the dataset used for coloring and labeling points in the plot
  # 'label = FALSE' indicates that labels for points should not be displayed
  # 'scale = TRUE' standardizes the data for plotting
  # 'colour' specifies the variable used for coloring points
  pca_plot <- autoplot(pca_data$pca, 
                       data = pca_data$data, 
                       label = FALSE, 
                       scale = TRUE,
                       colour = "age_group") +
    scale_color_brewer(palette = "Set2")
  
  # Return the PCA plot object
  return(pca_plot)
}

# Function to generate a cluster PCA plot using the autoplot function
get_cluster_pca_plot <- function(pca_data_with_kmeans, sample_with_age_group){
  
  # Transform the cluster column to factors
  pca_data_with_kmeans$data$cluster <- as.factor(pca_data_with_kmeans$data$cluster)
  pca_data_with_kmeans$data$age_group <- as.factor(pca_data_with_kmeans$data$age_group)
  
  # Create the base PCA plot
  cluster_pca_plot <- autoplot(pca_data_with_kmeans$pca,
                               data = pca_data_with_kmeans$data,
                               colour = "cluster",
                               shape = "age_group",
                               label = FALSE,
                               scale = TRUE,
                               frame = TRUE,
                               x = 1,
                               y = 2,
                               # loadings = TRUE,
                               # loadings.label = TRUE,
                               # loadings.label.size= 3,
                               frame.type = "norm") +
    scale_color_brewer(palette = "Set2") +
    scale_fill_brewer(palette = "Set2") +
    coord_fixed(ratio = 1)  # Ensure equal scaling of x and y axes
  
  # Return the cluster PCA plot object
  return(cluster_pca_plot)
}

# Function to generate a 3D PCA plot
get_cluster_pca_plot_3d <- function(pca_data_with_kmeans) {
  
  # Transform the cluster column to factors
  pca_data_with_kmeans$data$cluster <- as.factor(pca_data_with_kmeans$data$cluster)
  
  # Extract PCA scores
  pca_scores <- as.data.frame(pca_data_with_kmeans$pca$x)
  pca_scores$cluster <- pca_data_with_kmeans$data$cluster
  
  # Create the 3D PCA plot
  plot_3d <- plot_ly(pca_scores, 
                     x = ~PC1, 
                     y = ~PC2, 
                     z = ~PC3, 
                     color = ~cluster, 
                     colors = "Set2",
                     type = "scatter3d", 
                     mode = "markers") %>%
    layout(title = "3D PCA Plot with Clusters",
           scene = list(xaxis = list(title = 'PC1'),
                        yaxis = list(title = 'PC2'),
                        zaxis = list(title = 'PC3')))
  
  # Return the 3D PCA plot object
  return(plot_3d)
}

# Function to extract PCA loadings and identify the top 5 loadings for the first two Principal Components (PC1 & PC2)
get_loadings <- function(pca_data){
  
  # Extract the PCA loadings from the PCA results
  # Convert the loadings for PC1 and PC2 into a data frame
  loadings <- as.data.frame(pca_data$rotation[, 1:2])
  
  loadings <- cbind(Protein = rownames(loadings), loadings)
  
  # Identify the top 10 loadings for PC1
  # Arrange the loadings in descending order based on the absolute values of PC1
  # Select only the PC1 column and take the top 5 loadings
  top_10_loadings_pc1 <- loadings %>%
    arrange(desc(abs(PC1))) %>%
    select(PC1) %>%
    head(10)
  
  # Identify the top 10 loadings for PC2 following the same steps for PC1
  top_10_loadings_pc2 <- loadings %>%
    arrange(desc(abs(PC2))) %>%
    select(PC2) %>%
    head(10)
  
  # Return a list containing the full loadings and the top 5 loadings for PC1 and PC2
  return(list(loadings = loadings, 
              top_10_loadings_pc1 = top_10_loadings_pc1,
              top_10_loadings_pc2 = top_10_loadings_pc2))
}

# Function to retrieve the top 5 proteins based on the PCA loadings
get_top_5_proteins <- function(pca_data){
  
  # Get PCA variable information (cos2 and contributions)
  pca_var <- get_pca_var(pca_data)
  
  # Extract cos2 values for the first two Principal Components (PCs)
  # and calculate the total cos2 values for each variable across the first two PCs
  cos2_values <- pca_var$cos2[, 1:2]
  total_cos2 <- rowSums(cos2_values)
  
  # Identify the top 5 variables based on total cos2 values
  top_5_variables <- names(sort(total_cos2, decreasing = TRUE)[1:5])
  
  # Extract contributions for the first PC
  # and identify the top 5 contributions for the first PC
  contributions_pc1 <- pca_var$contrib[, 1]
  top_5_contributions_pc1 <- sort(contributions_pc1, decreasing = TRUE)[1:5]
  
  # Extract contributions for the second PC
  # and identify the top 5 contributions for the second PC
  contributions_pc2 <- pca_var$contrib[, 2]
  top_5_contributions_pc2 <- sort(contributions_pc2, decreasing = TRUE)[1:5]
  
  # Return a list containing the top variables and contributions for both PCs
  return(list(top_5_total = top_5_variables, 
              top_5_pc1 = top_5_contributions_pc1,
              top_5_pc2 = top_5_contributions_pc2))
}

# Function to transform and filter data
get_filtered_data <- function(data, dataset_name, sample_with_age_group){
  data_name <- dataset_name
  if(data_name == "CVD_and_Inflammation_data" || data_name == "CVD_Inflammation_and_IO_data"){
    # List with all cytokines (column names)
    cytokines <<- colnames(data)[-1]
    
    # Filter out rows where the 'Sample' column contains 'HC'
    # Remove 'A' from the 'Sample' column values
    # Remove leading zeros from the 'Sample' column values, retaining the numbers
    # Convert the 'Sample' column values to numeric values
    # Join with 'sample_with_age_group' dataframe on the 'Sample' column
    # Convert the 'Sampling_age_group' column values to factor values
    data_filtered <- data %>%
      filter(!grepl("HC", Sample)) %>%
      mutate(Sample = gsub("A", "", Sample)) %>%
      mutate(Sample = gsub("^(\\D*)0+", "\\1", Sample))  %>%
      mutate(Sample = as.numeric(Sample))  %>% 
      mutate(across(everything(), as.numeric)) %>%
      left_join(sample_with_age_group, by =  "Sample") %>%
      mutate(as.factor(Sampling_age_group))
    
    # Identify columns ending with '.y'
    # These columns are the overlapping proteins originated from the Inflammation panel
    y_columns <- colnames(data_filtered)[grepl("\\.y$", colnames(data_filtered))]
    
    # Select only the cytokine columns, excluding '.y' columns, for the protein_data_filt dataset
    protein_data_filt <- data_filtered[, cytokines] %>%
      select(-all_of(y_columns))
    
    # Remove the 'Sample' and '.y' columns from the filtered dataset
    data_filtered <- data_filtered %>%
      select(-Sample, -all_of(y_columns))
    
    # Return a list containing the filtered data and the relevant protein data
    return(list(data_filtered = data_filtered, 
                protein_data_filt = protein_data_filt))
  } else if(data_name == "IO_data"){
    # Remove columns Plate ID and QC Warning
    # Remove 'A' from the Sample column values
    # Convert all values into numeric values
    data_filtered <- data %>%
      select(-`Plate ID`, -`QC Warning`) %>%
      mutate(Sample = gsub("A", "", Sample)) %>%
      mutate(across(everything(), as.numeric))
    
    # Impute NA values with the LOD values of that protein
    # Convert all values into numeric values
    data_filtered <- data_filtered %>%
      mutate(across(everything(), ~ ifelse(is.na(.), lod_row[[cur_column()]], .))) %>%
      mutate(across(everything(), as.numeric))

    # Select only the cytokine columns for the protein_data_filt dataset
    protein_data_filt <- data_filtered %>%
      select(-Sample)
    
    # Add sample nr to the filtered dataset
    data_filtered <- data_filtered %>%
      left_join(sample_with_age_group, by =  "Sample") %>%
      mutate(as.factor(Sampling_age_group))
    
    # Return a list containing the filtered data and the relevant protein data
    return(list(data_filtered = data_filtered,
                protein_data_filt = protein_data_filt))
  } else{
    print("Invalid dataset")
  }
  
}

# Function to summarize cluster data from PCA results
get_table_values <- function(table, group_var){
  variables <- colnames(table)
  
  # Create an empty list for the table values
  table_values <- list()
  
  if(group_var == "Sampling_age_group") {
    cohort_table_values <- list()
    # Count the number of patients in the data and add this number to cohort_count
    cohort_count <- table %>%
      summarise(count = n()) %>%
      select(count_cohort = count)
    cohort_table_values$cohort_count <- cohort_count
    
    if("Gender" %in% variables) {
      # Count the number of males and females and calculate the percentages
      gender_count <- table %>%
        count(Gender) %>%
        mutate(Total_count = sum(n),
               Percentage = sprintf("%.1f", (n / Total_count) * 100)) %>%
        select(Gender_cohort = Gender, n_gender_cohort = n, Percentage_gender_cohort = Percentage)
      cohort_table_values$gender_count <- gender_count
    }
    
    if("Age_onset_group" %in% variables) {
      # Count the number of patients in each age of onset group and calculate the percentages
      onset_count <- table %>%
        count(Age_onset_group) %>%
        mutate(Total_count = sum(n),
               Percentage = sprintf("%.1f", (n / Total_count) * 100)) %>%
        select(Age_onset_group_cohort = Age_onset_group, n_onset_cohort = n, Percentage_onset_cohort = Percentage)
      cohort_table_values$onset_count <- onset_count
    }
    
    if("Age_sampling" %in% variables & "EASI" %in% variables) {
      # Calculate the averages and standard deviations for the age and EASI score
      averages_cohort <- table %>%
        summarise(avg_age_cohort = mean(Age_sampling),
                  sd_age_cohort = sd(Age_sampling),
                  avg_EASI_cohort = mean(EASI, na.rm = TRUE),
                  sd_EASI_cohort = sd(EASI, na.rm = TRUE)) %>%
        mutate(across(everything(), ~ sprintf("%.1f", .)))
      cohort_table_values$averages <- averages_cohort
    }
    
    if("Atopy_yes_or_no" %in% variables) {
      # Count the number of patients with comorbidities and calculate the percentages
      comorbidities <- table %>%
        count(Atopy_yes_or_no) %>%
        mutate(Total_count = sum(n),
               Percentage = sprintf("%.1f", (n / Total_count) * 100)) %>%
        select(comorbidities_cohort = Atopy_yes_or_no, n_comorbidities_cohort = n, Percentage_comorbidities_cohort = Percentage)
      cohort_table_values$comorbidities <- comorbidities
    }
    
    if("AA" %in% variables){
      # Count the number of cases of Allergic Asthma (AA) and calculate their percentages
      aa_cohort <- table %>%
        count(AA) %>%
        mutate(Total_count = sum(n),
               Percentage = sprintf("%.1f", (n / Total_count) * 100)) %>%
        select(AA_cohort = AA, n_AA_cohort = n, Percentage_AA_cohort = Percentage)
      cohort_table_values$aa <- aa_cohort
    }
    
    if("ARC" %in% variables){
      # Count the number of cases of Allergic Rhinoconjunctivitis (ARC) and calculate their percentages
      arc_cohort <- table %>%
        count(ARC) %>%
        mutate(Total_count = sum(n),
               Percentage = sprintf("%.1f", (n / Total_count) * 100)) %>%
        select(ARC_cohort = ARC, n_ARC_cohort = n, Percentage_ARC_cohort = Percentage)
      cohort_table_values$arc <- arc_cohort
    }
    
    if("FA" %in% variables){
      # Count the number of cases of Food Allergy (FA) and calculate their percentages
      fa_cohort <- table %>%
        count(FA) %>%
        mutate(Total_count = sum(n),
               Percentage = sprintf("%.1f", (n / Total_count) * 100)) %>%
        select(FA_cohort = FA, n_FA_cohort = n, Percentage_FA_cohort = Percentage)
      cohort_table_values$fa <- fa_cohort
    }
  } 
  # Count the number of entries for each age group
  count_per_cluster <- table %>%
    group_by(!!sym(group_var)) %>%
    summarise(count = n())
  table_values$count_per_cluster <- count_per_cluster
  
  if("Gender" %in% variables) {
    # Count the number of males and females per cluster/age group and calculate their percentages
    gender_per_cluster <- table %>%
      group_by(!!sym(group_var), Gender) %>%
      count(Gender) %>%
      group_by(!!sym(group_var)) %>%
      mutate(Total_count = sum(n),
             Percentage = sprintf("%.1f", (n / Total_count) * 100)) %>%
      select(!!sym(group_var), Gender, n_gender = n, Percentage_gender = Percentage)
    table_values$gender_per_cluster <- gender_per_cluster
  }
  
  if("Age_onset_group" %in% variables) {
    # Count the number of cases per age-onset group for each cluster/age group and calculate their percentages
    onset_per_cluster <- table %>%
      group_by(!!sym(group_var), Age_onset_group) %>%
      count(Age_onset_group) %>%
      group_by(!!sym(group_var)) %>%
      mutate(Total_count = sum(n),
             Percentage = sprintf("%.1f", (n / Total_count) * 100)) %>%
      select(!!sym(group_var), Age_onset_group, n_onset = n, Percentage_onset = Percentage)
    table_values$onset_per_cluster <- onset_per_cluster
  }
  
  if("Age_sampling" %in% variables & "EASI" %in% variables){
    # Calculate the mean and standard deviation of age and EASI score for each cluster/age group
    averages_per_cluster <- table %>%
      group_by(!!sym(group_var)) %>%
      summarise(avg_age = mean(Age_sampling),
                sd_age = sd(Age_sampling),
                avg_EASI = mean(EASI, na.rm = TRUE), 
                sd_EASI = sd(EASI, na.rm = TRUE)) %>%
      mutate(across(everything(), ~ as.numeric(sprintf("%.1f", .)))) 
    if(group_var == "Sampling_age_group") {
      averages_per_cluster <- averages_per_cluster %>%
        mutate(Sampling_age_group = recode(Sampling_age_group,
                                           `1` = "0-4 years",
                                           `2` = "5-11 years",
                                           `3` = "12-17 years")) %>%
        mutate(!!sym(group_var) := as.factor(!!sym(group_var)))
    }
    table_values$averages_per_cluster <- averages_per_cluster
  }
  
  if("Atopy_yes_or_no" %in% variables){
    # Count comorbidities for each cluster/age group and calculate their percentages
    comorbidities_per_cluster <- table %>%
      group_by(!!sym(group_var), Atopy_yes_or_no) %>%
      count(Atopy_yes_or_no) %>%
      group_by(!!sym(group_var)) %>%
      mutate(Total_count = sum(n),
             Percentage = sprintf("%.1f", (n / Total_count) * 100)) %>%
      select(!!sym(group_var), comorbidities = Atopy_yes_or_no, n_comorbidities = n, Percentage_comorbidities = Percentage)
    table_values$comorbidities_per_cluster <- comorbidities_per_cluster
  }
  
  if("AA" %in% variables){
    # Count the number of cases of Allergic Asthma (AA) for each cluster/age group and calculate their percentages
    aa_per_cluster <- table %>%
      group_by(!!sym(group_var), AA) %>%
      count(AA) %>%
      group_by(!!sym(group_var)) %>%
      mutate(Total_count = sum(n),
             Percentage = sprintf("%.1f", (n / Total_count) * 100)) %>%
      select(!!sym(group_var), AA, n_AA = n, Percentage_AA = Percentage)
    table_values$aa_per_cluster <- aa_per_cluster
  }
  
  if("ARC" %in% variables){
    # Count the number of cases of Allergic Rhinoconjunctivitis (ARC) for each cluster/age group and calculate their percentages
    arc_per_cluster <- table %>%
      group_by(!!sym(group_var), ARC) %>%
      count(ARC) %>%
      group_by(!!sym(group_var)) %>%
      mutate(Total_count = sum(n),
             Percentage = sprintf("%.1f", (n / Total_count) * 100)) %>%
      select(!!sym(group_var), ARC, n_ARC = n, Percentage_ARC = Percentage)
    table_values$arc_per_cluster <- arc_per_cluster
  }
  
  if("FA" %in% variables){
    # Count the number of cases of Food Allergy (FA) for each cluster/age group and calculate their percentages
    fa_per_cluster <- table %>%
      group_by(!!sym(group_var), FA) %>%
      count(FA) %>%
      group_by(!!sym(group_var)) %>%
      mutate(Total_count = sum(n),
             Percentage = sprintf("%.1f", (n / Total_count) * 100)) %>%
      select(!!sym(group_var), FA, n_FA = n, Percentage_FA = Percentage)
    table_values$fa_per_cluster <- fa_per_cluster
  }
  
  if(group_var == "Sampling_age_group") {
    # Return a list of all summarized tables and the values for the total cohort
    return(list(table_values = table_values,
                cohort_table_values = cohort_table_values))
  } else {
    # Return a list of all summarized tables
    return(table_values)
  }
}

# Function to create a table of cluster data from PCA analysis
get_cluster_table <- function(pca_data_with_kmeans, database_with_age_groups){
  
  samples <- pca_data_with_kmeans$data$sample
  clusters <- pca_data_with_kmeans$data$cluster
  pca_cluster_data <<- data.frame(Sample = samples, Cluster = clusters) %>%
    left_join(database_with_age_groups, join_by("Sample"))
  
  # Retrieve table values from PCA cluster data
  table_values <- get_table_values(pca_cluster_data, "Cluster")
  # Prepare the the cluster table based on the table values
  table <- prepare_cluster_data(table_values)
  # Mutate the cluster data to include necessary transformations
  table <- mutate_cluster_data(table, "Cluster")
  # Reshape the cluster data into a wide format suitable for summarization
  wide_table <- reshape_cluster_data(table, "Cluster")
  # Generate a final GT table from the reshaped cluster data and original PCA data
  gt_table <- generate_cluster_table(wide_table, pca_cluster_data, "Cluster_name")
  
  # Return the final GT table for display
  return(gt_table)
}

# Function to bind the cluster data
prepare_cluster_data <- function(table_values){
  
  # Bind all rows from 'table_values' into a single dataframe
  table <- bind_rows(table_values)
  
  # Return the combined table with additional data
  return(table)
}

# Function to mutate and format data for the cluster table
mutate_cluster_data <- function(table, group_var){
  variables <- colnames(table)
  
  if(group_var == "Cluster") {
    table <- table %>%
      # Create a new column for cluster names with sample counts
      # Remove the count column
      mutate("Cluster_name" = paste0("Cluster ", Cluster, " (n=", count, ")")) %>%
      select(-count)
  } else if(group_var == "Sampling_age_group") {
    table <- table %>%
      mutate("Sampling_age" = paste0(Sampling_age_group, " (n=", count, ")")) %>%
      select(-count)
  } else if(group_var == "Cohort") {
    table <- table %>%
      mutate("Cohort" = paste0("Total cohort (n=", count_cohort, ")")) %>%
      select(-count_cohort)
    colnames(table) <- gsub("_cohort", "", colnames(table))
    variables <- colnames(table)
  }

  
  if("avg_age" %in% variables){
    table <- table %>%
      # Format age data as mean and standard deviation
      # Remove the original age data
      mutate("Age (y), mean (SD)" = paste0(avg_age, " (", sd_age, ")")) %>%
      select(-avg_age, -sd_age)
  }
  
  if("Gender" %in% variables) {
    table <- table %>%
      # Format gender counts as number and percentage for males/females
      # Remove the original gender data
      mutate("Male, no. (%)" = ifelse(Gender == 0, paste0(n_gender, " (", Percentage_gender, ")"), NA)) %>%
      select(-Gender, -n_gender, -Percentage_gender)
  }
  
  if("avg_EASI" %in% variables){
    table <- table %>%
      # Format EASI score as mean and standard deviation
      # Remove the original EASI score data
      mutate("EASI score, mean (SD)" = paste0(avg_EASI, " (", sd_EASI, ")")) %>%
      select(-avg_EASI, -sd_EASI)
  }
  
  if ("Age_onset_group" %in% variables) {
    table <- table %>%
      # Format onset age group data for different age ranges
      # Remove the original onset age group data
      mutate("0-1 years" = ifelse(Age_onset_group == "0-1 years", paste0(n_onset, " (", Percentage_onset, ")"), NA)) %>%
      mutate("2-4 years" = ifelse(Age_onset_group == "2-4 years", paste0(n_onset, " (", Percentage_onset, ")"), NA)) %>%
      mutate("5-11 years" = ifelse(Age_onset_group == "5-11 years", paste0(n_onset, " (", Percentage_onset, ")"), NA)) %>%
      mutate("12-17 years" = ifelse(Age_onset_group == "12-17 years", paste0(n_onset, " (", Percentage_onset, ")"), NA)) %>%
      mutate("Missing" = ifelse(is.na(Age_onset_group), paste0(n_onset, " (", Percentage_onset, ")"), NA)) %>%
      select(-Age_onset_group, -n_onset, -Percentage_onset)
  }
  
  if("comorbidities" %in% variables) {
    table <- table %>%
      # Format comorbidity data as number and percentages
      # Remove the original comorbidity data
      mutate("Comorbidities, no. (%)" = ifelse(comorbidities == 1, paste0(n_comorbidities, " (", Percentage_comorbidities, ")"), NA)) %>%
      select(-comorbidities, -n_comorbidities, -Percentage_comorbidities)
  }
  
  if("AA" %in% variables & "ARC" %in% variables & "FA" %in% variables) {
    table <- table %>%
      # Format allergy data as number and percentages asthma/rhinoconjunctivitis/food
      # Remove the orinial allergy data
      mutate("Allergic asthma, no. (%)" = ifelse(AA == 1, paste0(n_AA, " (", Percentage_AA, ")"), NA)) %>%
      mutate("Allergic rhinoconjunctivitis, no. (%)" = ifelse(ARC == 1, paste0(n_ARC, " (", Percentage_ARC, ")"), NA)) %>%
      mutate("Food allergy, no. (%)" = ifelse(FA == 1, paste0(n_FA, " (", Percentage_FA, ")"), NA)) %>%
      select(-AA, -n_AA, -Percentage_AA, -ARC, -n_ARC, -Percentage_ARC, -FA, -n_FA, -Percentage_FA)
  }
  
  # Convert all columns to character type for uniformity
  table <- table %>%
    mutate(across(everything(), as.character))

  # Return the mutated and formatted table
  return(table)
}

# Function to pivot data from wide to long format and back to wide
reshape_cluster_data <- function(table, group_var){
  if(group_var == "Cohort") {
    # Detect the value of n that is not NA
    valid_n <- table %>%
      filter(!is.na(Cohort) & grepl("Total cohort \\(n=", Cohort)) %>%
      pull(Cohort) %>%
      unique() %>%
      .[1]
    
    # Sub "Total cohort (n=NA)" with the valid value of n
    table <- table %>%
      mutate(Cohort = ifelse(Cohort == "Total cohort (n=NA)", valid_n, Cohort))
  }
  
  # Convert the table from wide to long format
  long_table <- table %>%
    pivot_longer(cols = -!!sym(group_var), names_to = "Variable", values_to = "Value")
  # Filter out rows with NA values in the Value column
  long_table <- long_table[long_table$Value != "NA (NA)", ]
  if(group_var == "Cluster") {
    # Filter out rows where Cluster name contains NA
    long_table <- long_table[!grepl("^Cluster \\d+ \\(n=NA\\)$", long_table$Value), ]
  } else if(group_var == "Sampling_age_group") {
    # Filter out rows where Age of sampling contains NA
    long_table <- long_table[!grepl("(0-4 years|5-11 years|12-17 years) \\(n=NA\\)", long_table$Value), ]
  } else if(group_var == "Cohort") {
    long_table <- long_table[!grepl("Total cohort \\(n=NA\\)", long_table$Value), ]
  }

  # Convert the long table back to wide format 
  # The values from the Cluster column will become the new column names in the wide table
  wide_table <- long_table %>%
    pivot_wider(names_from = !!sym(group_var), values_from = Value) %>%
    # Remove any NA columns that may have been created
    select(-`NA`) %>%
    # Remove rows with NA in the Variable column
    drop_na(Variable)
  
  # Return the reshaped wide table
  return(wide_table)
}

# Function to finalize and style the cluster table
generate_cluster_table <- function(table, data, group_var){
  if(group_var == "Sampling_age") {
    wide_table <- table$wide_table
    wide_table_cohort <- table$wide_table_cohort
  } else {
    wide_table <- table
  }
  
  variables <- colnames(data)
  if("Age_onset" %in% variables) {
    # Create a new empty row for the Age of onset variable
    new_row <- data.frame(matrix("", nrow = 1, ncol = ncol(wide_table))) 
    colnames(new_row) <- colnames(wide_table)  
    new_row[1, "Variable"] <- "Age of onset, no. (%)"  
    
    # Append the new row to the wide table
    wide_table <- rbind(wide_table, new_row)
    
    if(group_var == "Sampling_age") {
      # Create a new empty row for the Age of onset variable
      new_row <- data.frame(matrix("", nrow = 1, ncol = ncol(wide_table_cohort))) 
      colnames(new_row) <- colnames(wide_table_cohort)  
      new_row[1, "Variable"] <- "Age of onset, no. (%)"  
      
      # Append the new row to the wide table
      wide_table_cohort <- rbind(wide_table_cohort, new_row)
    }
  }

  # Update column names to use the first row as headers (Cluster name)
  colnames(wide_table) <- unlist(wide_table[1, ])
  # Remove the first row , now used as header
  wide_table <- wide_table[-1, ]
  
  if(group_var == "Cluster_name") {
    # Count the number of clusters, 
    # using the number columns in the table - 1 (Clinical characteristics column)
    clusters <- ncol(wide_table) - 1
    # Calculate the p-values
    p_values <- calculate_p_values(data, clusters, variables, "Cluster_name")
  } else if (group_var == "Sampling_age") {
    groups <- 3
    p_values <- calculate_p_values(data, groups, variables, "Sampling_age")
  }

  # Add the p-values to the table
  wide_table <- wide_table %>%
    left_join(p_values, by = group_var) %>%
    mutate(
      p_value = ifelse(p_value < 0.001, "<.001", sprintf("%.3f", p_value)),  
      p_value = as.character(p_value),
      p_value = replace_na(p_value, "") 
    ) 
  
  if(group_var == "Sampling_age") {
    wide_table <- wide_table %>%
      rename(Variable = Sampling_age)
    wide_table <- wide_table_cohort %>%
      left_join(wide_table, by = "Variable") %>%
      rename(Sampling_age = Variable)
  }
    
  # Define the desired order of variables for the final table
  variable_order <- c("Age (y), mean (SD)",
                      "Male, no. (%)",
                      "EASI score, mean (SD)",
                      "Age of onset, no. (%)",
                      "0-1 years",
                      "2-4 years",
                      "5-11 years",
                      "12-17 years",
                      "Missing",
                      "Comorbidities, no. (%)",
                      "Allergic asthma, no. (%)",
                      "Allergic rhinoconjunctivitis, no. (%)",
                      "Food allergy, no. (%)")
  
  # Ensure "12-17 years" row is present
  if(!"12-17 years" %in% wide_table[[group_var]]) {
    new_row <- data.frame(matrix("", nrow = 1, ncol = ncol(wide_table)))
    colnames(new_row) <- colnames(wide_table)
    new_row[[group_var]] <- "12-17 years"
    wide_table <- rbind(wide_table, new_row)
  }
  
  # Convert Cluster_name to a factor with specified levels and arrange by this order
  wide_table <- wide_table %>%
    mutate(!!sym(group_var) := factor(!!sym(group_var), levels = variable_order)) %>%
    arrange(!!sym(group_var))
  
  # Create a formatted table using the gt package
  gt_table <- wide_table %>%
    gt() %>%
    cols_label(!!sym(group_var) := "Clinical characteristic") %>%
    tab_footnote(
      footnote = "Age at the moment of sample collection",
      locations = cells_body(
        columns = c(!!sym(group_var)), 
        rows = !!sym(group_var) == "Age (y), mean (SD)"
      )) %>%
    tab_footnote(
      footnote = "Categoric variables are presented as counts and percentages; continuous variables are presented as means with SDs."
    ) %>%
    tab_footnote(
      footnote = "Clinical characteristics between age groups were compared by using a one-way ANOVA test for continuous variables and chi-sqaure test for categoric variabels."
    ) %>%
    tab_style(
      style = list(
        cell_text(weight = "bold")
      ),
      locations = cells_column_labels(
        everything()
      )
    ) %>%
    opt_stylize(style = 1, color = "gray", add_row_striping = TRUE)
  
  # Return the finilized table
  return(gt_table)
}

# Function to calculate the p-values for various variables across the clusters
calculate_p_values <- function(table, clusters, variables, group_var){
  if(group_var == "Cluster_name") {
    relevant_columns <- c("Age_sampling", "Cluster", "EASI", "Gender", "Age_onset", "Atopy_yes_or_no", "AA", "ARC", "FA")
  } else if(group_var == "Sampling_age") {
    relevant_columns <- c("Age_sampling", "Sampling_age_group", "EASI", "Gender", "Age_onset", "Atopy_yes_or_no", "AA", "ARC", "FA")
  }
  
  available_columns <- intersect(variables, relevant_columns)
  
  # Select relevant columns 
  table <- table %>%
    select(all_of(available_columns))
  # Replace 0 with 2 in the Atopy_yes_or_no column for future use  
  if("Atopy_yes_or_no" %in% available_columns) {
    table <- table %>%
      mutate(Atopy_yes_or_no = replace(Atopy_yes_or_no, Atopy_yes_or_no == 0, 2))
  }
  
  # Initialize p-values with NA
  age_p_value <- NA
  easi_p_value <- NA
  comorbidities_p_value <- NA
  AA_p_value <- NA
  ARC_p_value <- NA
  FA_p_value <- NA
  gender_p_value <- NA
  age_onset_p_value <- NA
  
  if(group_var == "Cluster_name") {
    # Calculate p-values for continuous variables
    if("Age_sampling" %in% available_columns) {
      age_p_value <- calculate_ttest_anova_p_values(table, clusters, "Age_sampling", "Cluster")
    }
    if("EASI" %in% available_columns) {
      easi_p_value <- calculate_ttest_anova_p_values(table, clusters, "EASI", "Cluster")
    }
    
    # Calculate p-values for categorial values
    if("Atopy_yes_or_no" %in% available_columns) {
      comorbidities_p_value <- calculate_chisq_p_value(table, "Atopy_yes_or_no", clusters, "Cluster")
    }
    if("AA" %in% available_columns) {
      AA_p_value <- calculate_chisq_p_value(table, "AA", clusters, "Cluster")
    }
    if("ARC" %in% available_columns) {
      ARC_p_value <- calculate_chisq_p_value(table, "ARC", clusters, "Cluster")
    }
    if("FA" %in% available_columns) {
      FA_p_value <- calculate_chisq_p_value(table, "FA", clusters, "Cluster")
    }
    if("Gender" %in% available_columns) {
      gender_p_value <- calculate_chisq_p_value(table, "Gender", clusters, "Cluster")
    }
    if("Age_onset" %in% available_columns) {
      age_onset_p_value <- calculate_chisq_p_value(table, "Age_onset", clusters, "Cluster")
    }
  } else if (group_var == "Sampling_age") {
    # Calculate p-values for continuous variables
    if("Age_sampling" %in% available_columns) {
      age_p_value <- calculate_ttest_anova_p_values(table, clusters, "Age_sampling", "Sampling_age_group")
    }
    if("EASI" %in% available_columns) {
      easi_p_value <- calculate_ttest_anova_p_values(table, clusters, "EASI", "Sampling_age_group")
    }
    
    # Calculate p-values for categorial values
    if("Atopy_yes_or_no" %in% available_columns) {
      comorbidities_p_value <- calculate_chisq_p_value(table, "Atopy_yes_or_no", clusters, "Sampling_age_group")
    }
    if("AA" %in% available_columns) {
      AA_p_value <- calculate_chisq_p_value(table, "AA", clusters, "Sampling_age_group")
    }
    if("ARC" %in% available_columns) {
      ARC_p_value <- calculate_chisq_p_value(table, "ARC", clusters, "Sampling_age_group")
    }
    if("FA" %in% available_columns) {
      FA_p_value <- calculate_chisq_p_value(table, "FA", clusters, "Sampling_age_group")
    }
    if("Gender" %in% available_columns) {
      gender_p_value <- calculate_chisq_p_value(table, "Gender", clusters, "Sampling_age_group")
    }
    if("Age_onset" %in% available_columns) {
      age_onset_p_value <- calculate_chisq_p_value(table, "Age_onset", clusters, "Sampling_age_group")
    }
  }
  
  # Combine p-values into a dataframe with corresponding characteristic names
  p_values <- data.frame(
    variable = c("Age (y), mean (SD)", 
                     "EASI score, mean (SD)", 
                     "Male, no. (%)", 
                     "Age of onset, no. (%)", 
                     "Comorbidities, no. (%)",
                     "Allergic asthma, no. (%)",
                     "Allergic rhinoconjunctivitis, no. (%)",
                     "Food allergy, no. (%)"),
    p_value = c(age_p_value, easi_p_value, gender_p_value, age_onset_p_value, comorbidities_p_value, AA_p_value, ARC_p_value, FA_p_value)
  )
  names(p_values)[1] <- group_var
  
  # Return the dataframe with p-values
  return(p_values)
}

# Function to calculate p-values for continuous variables
calculate_ttest_anova_p_values <- function(table, clusters, var_name, group_var) {
  # If the number of clusters is below 3, perform a t.test
  if (clusters < 3) {
    p_value <- t.test(as.formula(paste(var_name, "~", group_var)), data = table)$p.value
  } else {
    # If the number of clusters is 3 or more, perform a one-way ANOVA
    aov <- aov(as.formula(paste(var_name, "~ as.factor(", group_var, ")")), data = table)
    p_value <- summary(aov)[[1]][["Pr(>F)"]][1]
  }
  
  # Return the calculated p-value
  return(p_value)
}

# Function to calculate p-values for categorical variables
calculate_chisq_p_value <- function(data, var_name, clusters, group_var) {
  # Special case for variable 'Gender'
  if(var_name == "Gender"){
    # Count males and females in each cluster and reshape table for Chi-squared test
    table <- data %>%
      filter(!is.na(!!sym(var_name))) %>%
      mutate(!!sym(var_name) := as.character(!!sym(var_name))) %>%
      group_by(!!sym(group_var)) %>%
      summarise(Male = sum(!!sym(var_name) == 0, na.rm = TRUE),
                Female = sum(!!sym(var_name) == 1, na.rm = TRUE)) %>%
      pivot_longer(cols = c(Male, Female), names_to = var_name, values_to = "Count")
    
    # Perform Chi-squared test
    chisq_test <- chisq.test(matrix(table$Count, nrow = 2))
    
    # Return the calculated p-value
    return(chisq_test$p.value)
    
  } # Special case for variable 'Age_onset' 
  else if(var_name == "Age_onset"){
    age_onset_table <- data %>%
      # Create age group categories
      mutate(Age_group = case_when(
        is.na(!!sym(var_name)) ~ "Missing",
        !!sym(var_name) >= 0 & !!sym(var_name) <= 1 ~ "0-1",
        !!sym(var_name) >= 2 & !!sym(var_name) <= 4 ~ "2-4",
        !!sym(var_name) >= 5 & !!sym(var_name) <= 11 ~ "5-11",
        !!sym(var_name) >= 12 & !!sym(var_name) <= 17 ~ "12-17"
      )) %>%
      group_by(!!sym(group_var)) %>%
      # Count occurrences in each cluster
      summarise(
        `0-1` = sum(Age_group == "0-1", na.rm = TRUE),
        `2-4` = sum(Age_group == "2-4", na.rm = TRUE),
        `5-11` = sum(Age_group == "5-11", na.rm = TRUE),
        `12-17` = sum(Age_group == "12-17", na.rm = TRUE),
        Missing = sum(Age_group == "Missing", na.rm = TRUE)
      ) %>%
      # Reshape table for Chi-squared test
      pivot_longer(cols = c(`0-1`, `2-4`, `5-11`, `12-17`, `Missing`), 
                   names_to = "Age_group", values_to = "Count")
    
    # Create matrix for Chi-squared test
    age_onset_matrix <- matrix(age_onset_table$Count, nrow = 5, byrow = TRUE)
    # Perform Chi-squared test
    age_onset_chisq <- chisq.test(age_onset_matrix)
    
    # Return the calculated p-value
    return(age_onset_chisq$p.value)
    
  } # For all other variables
  else {
    table <- data %>%
      # Count the occurrences of 'Yes' for the variable and reshape table for Chi-squared test
      filter(!is.na(!!sym(var_name))) %>%
      mutate(!!sym(var_name) := as.character(!!sym(var_name))) %>%
      group_by(!!sym(group_var)) %>%
      summarise(Yes = sum(!!sym(var_name) == 1, na.rm = TRUE),
                No = sum(!!sym(var_name) == 2, na.rm = TRUE)) %>%
      pivot_longer(cols = c(Yes, No), names_to = var_name, values_to = "Count")

    # Perform Chi-squared test
    chisq_test <- chisq.test(matrix(table$Count, nrow = 2))
    
    # Return the calculated p-value
    return(chisq_test$p.value)
  }
}

# Function to perform ANOVA and post-hoc Tukey test on protein data
anova_test <- function(protein_data, clusters) {
  
  protein_data_with_clusters <- cbind(protein_data, Cluster = as.factor(clusters))
  
  # Perform ANOVA on each protein column in the dataset
  anova_models <- apply(protein_data, 2, function(protein_column) {
    aov(protein_column ~ protein_data_with_clusters$Cluster)
  })
  
  # Extract p-values from the ANOVA models
  anova_results <- sapply(anova_models, function(model) {
    summary(model)[[1]][["Pr(>F)"]][1]
  })
  
  # Adjust the p-values for multiple comparisons using False Discovery Rate method
  p_adjusted <- p.adjust(anova_results, method = "fdr")
  
  # Identify significant proteins based on the adjusted p-values
  significant_proteins <- which(p_adjusted < 0.05)

  # Perform Tukey's HSD post-hoc test for significant proteins
  tukey_results_list <- lapply(anova_models[significant_proteins], TukeyHSD)

  # Combine the results of the Tukey tests into a dataframe
  tukey_results <- do.call(rbind, lapply(names(tukey_results_list), function(protein) {
    tukey_result <- tukey_results_list[[protein]]$`protein_data_with_clusters$Cluster`
    data.frame(Protein = protein,
               Comparison = rownames(tukey_result),
               diff = tukey_result[,"diff"],
               lwr = tukey_result[,"lwr"],
               upr = tukey_result[,"upr"],
               p_adj = tukey_result[,"p adj"],
               stringsAsFactors = FALSE)
  }))
  
  # Filter the Tukey results to only get the significant comparisons
  significant_tukey_results <- tukey_results %>%
    filter(p_adj < 0.05)

  # Return a list containing ANOVA results, adjusted p-values, significant protein names and Tukey results
  return(list(anova_results = anova_results, 
       p_adjusted = p_adjusted, 
       significant_protein_names = names(significant_proteins),
       tukey_results = tukey_results,
       protein_data_with_clusters = protein_data_with_clusters,
       significant_tukey_results = significant_tukey_results))
}

# Function to calculate the mean protein expressions per cluster and
# identify the cluster with the highest expression per protein
get_mean_expressions <- function(proteins, data){
  
  # Select only the relevant protein columns and cluster
  data <- data %>%
    select(all_of(proteins), Cluster)
  
  # Calculate the mean value for each protein within a cluster
  average_per_cluster <- data %>%
    group_by(Cluster) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
  
  # Reshape the data to a long format and calculate the highest value per protein
  highest_values_with_clusters <- average_per_cluster %>%
  pivot_longer(cols = -Cluster, names_to = "Protein", values_to = "Value") %>%
  group_by(Protein) %>% 
  filter(Value == max(Value)) %>%
  slice(1) %>%
  ungroup()
  
  # Summarize the highest protein expressions for each cluster and concatenate the protein names
  print_highest_values_per_cluster <- highest_values_with_clusters %>%
    group_by(Cluster) %>%
    summarise(Proteins = paste(Protein, collapse = ", ")) 
  
  # Print the proteins with the highest values per cluster 
  print_highest_values_per_cluster %>%
    mutate(Cluster_info = paste("Cluster", Cluster, ":", Proteins)) %>%
    pull(Cluster_info) %>%
    cat(sep = "\n\n")
  
  # Return a list containing the dataframe with the average values per cluster and
  # the dataframe with the highest values with the corresponding cluster
  return(list(average_per_cluster = average_per_cluster,
              highest_values_with_clusters = highest_values_with_clusters))
}

# Function to calculate the highest percentage differences per protein and
# select top proteins
get_highest_diffs <- function(data, highest_expressions, loadings){
  # Remove unnecessary columns and convert diff to absolute values
  # Join the highest expressions table by Protein, then convert Value to absolute
  # values and calculate percentage difference
  data <- data %>%
    select(-Comparison, -lwr, -upr, -p_adj) %>%
    mutate(diff = abs(diff)) %>%
    left_join(highest_expressions, by = "Protein") %>%
    left_join(loadings, by = "Protein") %>%
    mutate(Value = abs(Value)) %>%
    mutate(PC1 = abs(PC1)) %>%
    mutate(PC2 = abs(PC2))
  
  # Get top 10 proteins with highest diff percentages per cluster, 
  # ensuring unique Protein names
  top_proteins_per_cluster_PC1 <- data %>%
    group_by(Cluster) %>%
    arrange(desc(PC1)) %>%
    filter(!duplicated(Protein)) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    pull(Protein)
  
  top_proteins_per_cluster_PC2 <- data %>%
    group_by(Cluster) %>%
    arrange(desc(PC2)) %>%
    filter(!duplicated(Protein)) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    pull(Protein)
  
  # Return list containing the full data with calculated columns and  
  # the top proteins per cluster
  return(list(data = data,
              top_proteins_per_cluster_PC1 = top_proteins_per_cluster_PC1,
              top_proteins_per_cluster_PC2 = top_proteins_per_cluster_PC2))
}

# Function to generate a radar chart for protein expressions
get_radar_chart <- function(significant_protein_names, protein_data_with_clusters, tukey_results_df, loadings){
  
  # Calculate mean expressions and highest values per cluster for significant proteins
  mean_data <- get_mean_expressions(significant_protein_names, protein_data_with_clusters)
  mean <- mean_data$average_per_cluster
  highest <- mean_data$highest_values_with_clusters
  
  # Get top proteins by highest loadings per principal component from Tukey test results
  sorted_diffs <- get_highest_diffs(tukey_results_df, highest, loadings)
  sorted_diffs_data <- sorted_diffs$data
  top_proteins_PC1 <- sorted_diffs$top_proteins_per_cluster_PC1
  top_proteins_PC2 <- sorted_diffs$top_proteins_per_cluster_PC2
  
  top_proteins <- unique(c(top_proteins_PC1, top_proteins_PC2))
  
  # Make one picture with both radar charts
  par(mfrow = c(1, 2), mar = c(1, 1, 1, 1)) 
  # Generate the radar charts
  generate_radar_chart(mean, top_proteins_PC1, "Radar plot PC1")
  generate_radar_chart(mean, top_proteins_PC2, "Radar plot PC2")
  
  # Return a list containing the mean expressions, the highest values with clusters
  # and both radar charts
  return(list(mean = mean,
              highest = highest,
              sorted_diffs_data = sorted_diffs_data,
              top_proteins = top_proteins))
}

# Function to generate a radar chart
generate_radar_chart <- function(data, top_proteins, title){
  # Select expression data for the top proteins, excluding the Cluster column
  data <- data %>%
    select(-Cluster) %>%
    select(all_of(top_proteins))
  names(data) <- gsub(".x", "", names(data))
  
  # Determine overall min and max values for scaling
  overall_max <- max(data, na.rm = TRUE)
  overall_min <- min(data, na.rm = TRUE)
  
  # Add the overall min and max values for each variable (ensures consistent scaling across axes)
  max_values <- rep(overall_max, ncol(data))
  min_values <- rep(overall_min, ncol(data))
  data <- rbind(max_values, min_values, data)
  
  # Define axis labels (adjust these to your desired labels)
  axis_labels <- seq(overall_min, overall_max, length.out = 7)
  
  # Create the radar chart
  radarchart(data,
             axistype = 1,  # Uses the same scale for all variables
             caxislabels = round(axis_labels, 1),
             pcol = brewer.pal(n = 8, name = "Set2"),
             pfcol = add.alpha(brewer.pal(n = 8, name = "Set2"), 0.3),
             plwd = 2,
             cglcol = "grey",
             cglty = 1,
             axislabcol = "grey",
             vlcex = 0.7,
             title = title
  )
  
  # Add a legend indicating cluster colors
  legend(x = "top",
         legend = paste("Cluster", 1:(nrow(data)-2)),
         horiz = TRUE,
         col = brewer.pal(n = 8, name = "Set2"),
         pch = 20,
         pt.cex = 1.5,
         cex = 0.8,
         bty = "n"
  )
}

# Function to generate a list of boxplots for the specified protein expressions
get_boxplots <- function(data, proteins){
  
  # Use lapply to iterate over the list of proteins
  boxplots <- lapply(proteins, function(protein) {
    # Create a boxplot for each protein
    ggplot(data, aes(x = as.factor(Cluster), y = !!sym(protein), fill = Cluster)) +
      geom_boxplot() +
      labs(x = "Cluster", y = "ADM Expression") +
      theme_minimal() +
      scale_fill_brewer(palette = "Set2") +
      ggtitle(paste("Boxplot", gsub(".x", "", protein)))
  })
  
  # Return the list of boxplots
  return(boxplots)
}

# Function to retrieve the age groups from clinical data
get_age_groups <- function(Clinical_data){

  # Select relevant columns and rename PatientID to Sample
  database_with_age_groups <- Clinical_data %>%
    select(PatientID, Age_onset, Age_sampling, Gender, EASI, 
           Atopy_yes_or_no, AA, ARC, FA) %>%
    rename(Sample = PatientID) %>%
    # Create age groups for the sampling age
    mutate(Sampling_age_group = cut(Age_sampling,
                                    breaks = c(0, 5, 12, 18),
                                    labels = c("0-4 years", "5-11 years", "12-17 years"),
                                    right = FALSE)) %>%
    # Create age groups for the age of onset
    mutate(Age_onset_group = cut(Age_onset,
                                 breaks = c(0, 2, 5, 12, 18),
                                 labels = c("0-1 years", "2-4 years", "5-11 years", "12-17 years"),
                                 right = FALSE))

  # Select the sample and its corresponding sampling age group
  sample_with_age_group <- database_with_age_groups %>%
    select(Sampling_age_group, Sample)
  
  # Return a list containing the complete database with age groups and the sample with age group
  return(list(database_with_age_groups = database_with_age_groups,
              sample_with_age_group = sample_with_age_group))
}

# Function to choose a dataset and clinical dataset
choose_dataset <- function(){
  # List of available datasets
  datasets <- list(
    "IO_data" = IO_data,
    "CVD_and_Inflammation_data" = CVD_and_Inflammation_data,
    "CVD_Inflammation_and_IO_data" = CVD_Inflammation_and_IO_data
  )
  
  # List of available clinical datasets
  clinical_datasets <- list(
    "Clinical_data_CVD_and_Inflammation" = Clinical_data_CVD_and_Inflammation,
    "Clinical_data_IO_panel" = Clinical_data_IO_panel,
    "Clinical_data_CVD_Inflammation_IO" = Clinical_data_CVD_Inflammation_IO
  )
  
  # Loop to choose a dataset
  repeat{
    cat(white("Available datasets:\n"))
    # Display available datasets
    for (dataset_name in names(datasets)) {
      cat(white("-", dataset_name, "\n"))
    }
    
    # Prompt user to enter the name of the dataset
    choice_dataset <- readline(prompt = "Enter the name of the dataset you want to use: ")
    
    # Check if the entered dataset name is valid
    if (choice_dataset %in% names(datasets)) {
      chosen_dataset <- datasets[[choice_dataset]]
      cat(green("Chosen dataset:", choice_dataset, "\n"))
      break
    } else {
      cat(red("Invalid choice. Try again.\n"))
    }
  }
  
  # Loop to choose clinical dataset
  repeat{
    cat(white("Available clinical datasets:\n"))
    # Display available clinical datasets
    for (dataset_name in names(clinical_datasets)) {
      cat(white("-", dataset_name, "\n"))
    }
    
    # Prompt user to enter the name of the clinical dataset
    choice_clinical <- readline(prompt = "Enter the name of the clinical dataset you want to use: ")
    
    # Check if the entered clinical dataset name is valid
    if (choice_clinical %in% names(clinical_datasets)) {
      chosen_clinical_dataset <- clinical_datasets[[choice_clinical]]
      cat(green("Chosen dataset:", choice_clinical, "\n"))
      break
    } else {
      cat(red("Invalid choice. Try again.\n"))
    }
  }
  
  # Return a list containing the chosen datasets and their names
  return(list(chosen_dataset = chosen_dataset,
              choice_dataset = choice_dataset,
              chosen_clinical_dataset = chosen_clinical_dataset,
              choice_clinical = choice_clinical))
}

# Function to generate a volcano plot using the Tukey's test results
get_volcano_plot <- function(tukey_results){
  
  min_p_adj <- tukey_results %>%
    filter(p_adj != 0) %>%
    summarise(min_pvalue = min(p_adj)) %>%
    pull(min_pvalue)
  
  tukey_results$p_adj[tukey_results$p_adj == 0] <- min_p_adj / 10
  # Calculate the negative log10 of the adjusted p-values
  tukey_results <- tukey_results %>%
    mutate(log10p = -log10(p_adj))
  
  # Determine significance based on the adjusted p-values and a threshold for the difference
  tukey_results <- tukey_results %>%
    mutate(Significant = ifelse(p_adj < 0.05 & abs(diff) > 1, "Significant", "Not Significant"))
  
  max_log_10_pvalue <- max(tukey_results$log10p)
  
  # Create a volcano plot using ggplot2
  volcano_plot <- ggplot(tukey_results, aes(x = diff, y = log10p, color = Significant)) +
    geom_point(alpha = 0.8, size = 2) +
    theme_minimal() +
    labs(title = "Volcano Plot",
         x = "Difference (diff)",
         y = "-Log10 P-value") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    scale_color_brewer(palette = "Set2") +
    # Add labels to the significant points
    geom_text_repel(data = subset(tukey_results, Significant == "Significant"), 
                    aes(label = gsub(".x", "", Protein)), 
                    size = 3, color = "black", 
                    box.padding = 0.3, point.padding = 0.3, 
                    segment.color = 'grey50',
                    max.overlaps = 15) + 
    ylim(0, max_log_10_pvalue + 1) +
    facet_wrap(~Comparison)
  
  # Return the volcano plot
  return(volcano_plot)
  
}

# Function to generate a characteristics table for the total cohort
get_clinical_characteristics_table <- function(database_with_age_groups){
  table_values <- get_table_values(database_with_age_groups, "Sampling_age_group")
  
  table <- prepare_cluster_data(table_values$table_values)
  table_cohort <- prepare_cluster_data(table_values$cohort_table_values)
  
  table <- mutate_cluster_data(table, "Sampling_age_group")
  table_cohort <- mutate_cluster_data(table_cohort, "Cohort")
  
  # Reshape the cluster data into a wide format suitable for summarization
  wide_table <- reshape_cluster_data(table, "Sampling_age_group")
  wide_table_cohort <- reshape_cluster_data(table_cohort, "Cohort")
  
  wide_table <- list(wide_table = wide_table,
                     wide_table_cohort = wide_table_cohort)
  
  # Generate a final GT table from the reshaped cluster data and original PCA data
  gt_table <- generate_cluster_table(wide_table, database_with_age_groups, "Sampling_age")
  
  return(gt_table)
}

####################################################################
######################### DATA #####################################
####################################################################

import_datasets <- function(){
  # Import the CVD, Inflammation and the clinical data
  CVD_data <<- read_excel("Data/Replaced_LOD_olink_Dermat_CVD_II_2019_10_24_v1.xlsx")
  Inflammation_data <<- read_excel("Data/Replaced_LOD_olink_Dermat_Inflammation_2019-10-24_v1.xlsx")
  Clinical_data_CVD_and_Inflammation <<- read_sav("data/Klinische data kind v6.sav")
  
  # Import the IO data
  IO_data <<- read_excel("data/AD_children_BAMBOO_2.xlsx")
  
  # Import the clinical data from the IO panel
  Clinical_data_IO_panel <<- read_excel("data/AD_Biomarker kind_selectie ziekteoverstijgend cohort_n80.xlsx")
}

merge_datasets <- function(){
  # Join the CVD and Inflammation datasets into one dataset
  CVD_and_Inflammation_data <<- CVD_data %>% 
    full_join(Inflammation_data, by = "Sample") 
  
  # Transform IO data to the same format as CVD and Inflammation data
  lod_row <- IO_data %>% filter(`...1` == "LOD")
  
  # Impute NA values with the LOD value of that protein
  dataset <- IO_data %>%
    mutate(across(everything(), ~ ifelse(is.na(.), lod_row[[cur_column()]], .)))
  
  dataset <- dataset %>%
    filter(grepl("^A\\d{3}$", `...1`) | `...1` == "Assay")
  
  new_colnames <- as.character(na.omit(unlist(dataset[dataset$...1 == "Assay", ])))
  colnames(dataset) <- new_colnames
  dataset <- dataset[-1, ] 
  IO_data <<- dataset %>%
    select(-`QC Deviation from median`) %>%
    rename(Sample = Assay)
  
  # Identify proteins with ".y" and remove these from the dataset
  # Remove ".x" from the protein names
  y_columns <- colnames(CVD_and_Inflammation_data)[grepl("\\.y$", colnames(CVD_and_Inflammation_data))]
  data_without_y_columns <- CVD_and_Inflammation_data %>%
    select(-all_of(y_columns)) %>%
    rename_with(~ gsub("\\.x$", "", .), ends_with(".x"))
  
  # Join the CVD_and_Inflammation dataset with the IO dataset
  CVD_Inflammation_and_IO_data <<- left_join(IO_data, data_without_y_columns, by = "Sample") %>%
    select(-`Plate ID`, -`QC Warning`)
  
  # Join the Clinical dataset of CVD_and_Inflammation with the Clinical dataset of IO
  # to ensure all clinical data is included
  clinical_dataset <- Clinical_data_IO_panel %>%
    left_join(Clinical_data_CVD_and_Inflammation, by = "PatientID") %>%
    rename(EASI = EASI.x)
  Clinical_data_IO_panel <<- clinical_dataset
  
  Clinical_data_CVD_Inflammation_IO <<- Clinical_data_IO_panel 
}

####################################################################
############### Analysis ###########################################
####################################################################

execute_workflow <- function() {
  suppressWarnings({
    repeat{
      cat("Do you want to load all needed packages?")
      load_packages <- readline(prompt = "yes / no: ")
      if(tolower(load_packages) == "yes") {
        cat("Loading packages. \n")
        # Load all packages
        load_all_packages()
        cat(green("All packages are loaded. \n"))
        break
      } else if(tolower(load_packages) == "no") {
        cat("Loading packages was skipped. \n")
        break
      } else {
        cat(yellow("Invalid input, try again. \n"))
      }
    }
    
    cat(green("Importing all datasets. \n"))
    import_datasets()
    merge_datasets()
    
    # Choose datasets
    chosen_datasets <- choose_dataset()
    dataset <- chosen_datasets$chosen_dataset
    dataset_name <- chosen_datasets$choice_dataset
    clinical_dataset <- chosen_datasets$chosen_clinical_dataset
    clinical_name <- chosen_datasets$choice_clinical
    
    # Make a new dataset containing the sample number and age groups
    age_group_results <- get_age_groups(clinical_dataset)
    sample_with_age_group <- age_group_results$sample_with_age_group
    database_with_age_groups <- age_group_results$database_with_age_groups
    
    repeat{
      cat(white("Do you want to view the characteristics table?"))
      characteristics_table <- readline(prompt = "yes / no: ")
      if(tolower(characteristics_table) == "yes") {
        cat(green("Generating characteristics table. \n"))
        # Generate clinical characteristics table
        clinical_characteristics_table <- get_clinical_characteristics_table(database_with_age_groups)
        print(clinical_characteristics_table)
        break
      } else if(tolower(characteristics_table) == "no"){
        cat(red("Characteristics table will not be shown. \n"))
        break
      } else {
        cat(yellow("Invalid input, try again. \n"))
      }
    }
    
    
    cat(green("Data is being prepared. \n"))
    # Filter data
    data_filtered_result <- get_filtered_data(dataset, dataset_name, sample_with_age_group)
    data_filtered <- data_filtered_result$data_filtered
    protein_data <- data_filtered_result$protein_data_filt
    
    repeat{
      cat(white("Do you want to perform a PCA?"))
      pca <- readline(prompt = "yes / no: ")
      if(tolower(pca) == "yes") {
        cat(green("PCA analysis is in progress. \n"))
        # Perform PCA
        pca_data <- get_PCA_data(protein_data, sample_with_age_group)
        
        repeat{
          cat(white("Do you want to view the PCA plot?"))
          pca_plot <- readline(prompt = "yes / no: ")
          if(tolower(pca_plot) == "yes") {
            cat(green("Generating PCA plot. \n"))
            # Generate PCA plot
            pca_plot <- get_PCA_plot(pca_data, data_filtered)
            print(pca_plot)
            break
          } else if(tolower(pca_plot) == "no"){
            cat(red("PCA plot will not be shown. \n"))
            break
          } else {
            cat(yellow("Invalid input, try again. \n"))
          }
        }
        break
        
      } else if(tolower(pca) == "no"){
        cat(red("PCA analysis was skipped. \n"))
        break
      } else {
        cat(yellow("Invalid input, try again. \n"))
      }
    }
    
    
    # Get top 5 proteins
    top_5_result <- get_top_5_proteins(pca_data$pca)
    
    repeat{
      cat(white("Do you want to calculate the optimal number of clusters?"))
      nr_of_clusters <- readline(prompt = "yes / no: ")
      if(tolower(nr_of_clusters) == "yes") {
        cat(green("Calculating the optimal number of clusters. \n"))
        # Determine optimal number of clusters
        number_of_clusters <- get_optimal_number_of_clusters(pca_data)
        print("Number of clusters (Silhouette method): \n")
        print(number_of_clusters$number_of_clusters_silhouette$Best.nc)
        print("Number of clusters (WSS method): \n")
        print(number_of_clusters$number_of_clusters_wss$Best.nc)
        
        repeat{
          cat(white("Do you want to see the plots for the number of clusters?"))
          plots_nr_of_clusters <- readline(prompt = "yes / no: ")
          
          if (tolower(plots_nr_of_clusters) == "yes") {
            print(number_of_clusters$number_of_clusters_wss_plot)
            print(number_of_clusters$number_of_clusters_gap_stat_plot)
            print(number_of_clusters$number_of_clusters_silhouette_plot)
            break
          } else if(tolower(plots_nr_of_clusters) == "no"){
            cat(red("The plots for the optimal number of clusters will not be shown. \n"))
            break
          } else {
            cat(yellow("Invalid input, try again. \n"))
          }
        }
        break
        
      } else if(tolower(nr_of_clusters) == "no"){
        cat(red("Calculating the optimal number of clusters was skipped. \n"))
        break
      } else {
        cat(yellow("Invalid input, try again. \n"))
      }
    }
    
    repeat{
      cat(white("Do you want to perform a cluster PCA analysis?"))
      cluster_pca <- readline(prompt = "yes / no: ")
      if(tolower(cluster_pca) == "yes") {
        repeat{
          cat(white("How many clusters do you want to use?"))
          centers_input <- readline(prompt = "Number of clusters: ")
          centers_input <- as.integer(centers_input)
          if(!is.na(centers_input)) {
            break
          } else {
            cat(yellow("Invalid input, try again. \n"))
          }
        }
        
        cat(green("Cluster PCA in progress. \n"))
        # Perform k-means clustering on PCA data
        pca_data_with_kmeans <- get_PCA_kmeans(pca_data, centers = centers_input, sample_with_age_group = sample_with_age_group)
        
        # Get loadings
        loadings_result <- get_loadings(pca_data_with_kmeans$pca)
        loadings <- loadings_result$loadings
        repeat{
          cat(white("Do you want to view the cluster PCA plot?"))
          cluster_plot <- readline(prompt = "yes / no: ")
          if(tolower(cluster_plot) == "yes") {
            cat(green("Generating cluster PCA plot. \n"))
            # Generate cluster PCA plot
            cluster_pca_plot <- get_cluster_pca_plot(pca_data_with_kmeans, sample_with_age_group)
            print(cluster_pca_plot)
            break
          } else if(tolower(cluster_plot) == "no"){
            cat(red("Cluster PCA plot will not be shown. \n"))
            break
          } else {
            cat(yellow("Invalid input, try again. \n"))
          }
        }
        
        repeat{
          cat(white("Do you want to view the 3D cluster PCA plot?"))
          cluster_plot_3d <- readline(prompt = "yes / no: ")
          if(tolower(cluster_plot_3d) == "yes") {
            cat(green("Generating 3D cluster PCA plot. \n"))
            # Generate 3D cluster PCA plot
            cluster_pca_plot_3d <- get_cluster_pca_plot_3d(pca_data_with_kmeans)
            print(cluster_pca_plot_3d)
            break
          } else if(tolower(cluster_plot_3d) == "no") {
            cat(red("3D cluster PCA plot will not be shown. \n"))
            break
          } else {
            cat(yellow("Invalid input, try again. \n"))
          }
        }
        
        break
      } else if(tolower(cluster_pca) == "no") {
        cat(red("Cluster PCA analysis was skipped. \n"))
        break
      } else {
        cat(yellow("Invalid input, try again. \n"))
      }
    }
    
    repeat{
      cat(white("Do you want to view the cluster table?"))
      table <- readline(prompt = "yes / no: ")
      if(tolower(table) == "yes") {
        cat(green("Generating cluster table. \n"))
        # Generate cluster table
        cluster_table <- get_cluster_table(pca_data_with_kmeans, database_with_age_groups)
        print(cluster_table)
        break
      } else if(tolower(table) == "no") {
        cat(red("Cluster table will not be shown. \n"))
        break
      } else {
        cat(yellow("Invalid input, try again. \n"))
      }
    }
    
    cat(green("Performing ANOVA test. \n"))
    # Perform ANOVA test
    clusters <- pca_data_with_kmeans$data$cluster
    anova_results <- anova_test(protein_data, clusters)
    significant_protein_names <- anova_results$significant_protein_names
    protein_data_with_clusters <- anova_results$protein_data_with_clusters
    
    repeat{
      cat(white("Do you want to print all significant protein names?"))
      protein_names <- readline(prompt = "yes / no: ")
      if(tolower(protein_names) == "yes"){
        cat(white("Significant proteins: \n"))
        print(sort(anova_results$significant_protein_names))
        break
      } else if(tolower(protein_names) == "no" ) {
        cat(red("Significant protein names will not be printed. \n"))
        break
      } else {
        cat(yellow("Invalid input, try again. \n"))
      }
    }
    
    cat(green("Performing Tukey's test. \n"))
    # Get Tukey's test results
    tukey_results_df <- anova_results$tukey_results
    sign_tukey_results <- anova_results$significant_tukey_results
    
    repeat{
      cat(white("Do you want to view the Tukey's test results?"))
      tukey_test <- readline(prompt = "yes / no: ")
      if(tolower(tukey_test) == "yes"){
        View(tukey_results_df)
        break
      } else if (tolower(tukey_test) == "no") {
        cat(red("Tukey's test results will not be shown. \n"))
        break
      } else {
        cat(yellow("Invalid input, try again. \n"))
      }
    }
    
    repeat{
      cat(white("Do you want to see the radar chart with all significant proteins?"))
      radar_plot <- readline(prompt = "yes / no: ")
      if(tolower(radar_plot) == "yes"){
        cat(green("Generating radar chart. \n"))
        # Generate radar chart
        radar_chart <- get_radar_chart(significant_protein_names, protein_data_with_clusters, tukey_results_df, loadings)
        sorted_data <- radar_chart$sorted_diffs_data
        break
      } else if(tolower(radar_plot) == "no") {
        cat(red("Radar chart will not be shown. \n"))
        break
      } else {
        cat(yellow("Invalid input, try again. \n"))
      }
    }
    
    repeat{
      cat(white("Do you want to view the boxplots of the proteins in the radar chart?"))
      box_plots <- readline(prompt = "yes / no: ")
      if(tolower(box_plots) == "yes") {
        cat(green("Generating boxplots. \n"))
        # Generate box plots
        box_plot <- get_boxplots(protein_data_with_clusters, radar_chart$top_proteins)
        print(box_plot)
        break
      } else if(tolower(box_plots) == "no") {
        cat(red("Boxplots will not be shown. \n"))
        break
      } else {
        cat(yellow("Invalid input, try again. \n"))
      }
    }
    
    repeat{
      cat(white("Do you want to view the volcano plots?"))
      volcano_plots <- readline(prompt = "yes / no: ")
      if(tolower(volcano_plots) == "yes") {
        cat(green("Generating volcano plots. \n"))
        # Generate volcano plot
        volcano_plot <- get_volcano_plot(tukey_results_df)
        print(volcano_plot)
        break
      } else if(tolower(volcano_plots) == "no") {
        cat(red("Volcano plots will not be shown. \n"))
        break
      } else {
        cat(yellow("Invalid input, try again. \n"))
      }
    }
    
    cat(green("End analysis."))
  })

}


# Execute the workflow
execute_workflow()
# 
