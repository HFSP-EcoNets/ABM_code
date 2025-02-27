### This is an R-script that contains utility funcitons for experimental design ###

### Load the packages
# library(tidyverse) # for managing data
library(igraph) # for visualizing the matrix 
library(Matrix) # for generating matrix
library(plot.matrix) # for plotting matrix
library(gplots) # for plotting matrix
library(ggplot2) # for plotting matrix
library(reshape2) # for reshaping matrix into an object that ggplot() can recognize


### Functions

# Function to generate a random square matrix (non-zero value = 1) with desired matrix size and number of modules (symmetric)
modular_sq_matrix <- function(size_matrix, n_module, p_within, p_between, matrix_type) {
  # Sample the size of each module (>= 2)
  size_modules <- integer()
  size_left <- size_matrix

  repeat {
    for (i in 1:n_module) {
      size_module <- ifelse(i == n_module, size_left, sample(2:(size_left-1), 1))
      size_modules <- c(size_modules, size_module)
      size_left <- size_left - size_module
    }
  
    if (size_modules[length(size_modules)] == 1) {
      size_modules <- integer()  # Reset size_modules
      size_left <- size_matrix  # Reset size_left
    } else {
      break  # Exit the loop if the last element is not 1
    }
}
  
 # Arrange the size of the modules from high to low
 # size_modules = size_modules[order(-size_modules)]
  
  # Generate probability matrix for each of the module, and generate elements using SBM
  m <- size_modules
  k <- length(size_modules)
  c <- matrix(0, nrow = k, ncol = k)
  
  # Set intra-block (i.e. intra-module) connection probabilities
  diag(c) <- p_within

  # Set inter-block connection probabilities
  c[lower.tri(c)] <- p_between
  c[upper.tri(c)] <- p_between

  g <- sample_sbm(sum(m), pref.matrix=c, block.sizes=m)
  matrix <- as_adjacency_matrix(g)

  # Set diagonal elements to 1 for matrix H
  if (matrix_type == "H") {diag(matrix) <- 1}
  
  
  # Convert sparse matrix to a dense matrix for visualization
  dense_matrix <- as.matrix(matrix)



  # Covert the dense matrix into dataframe for plotting with ggplot()
  matrix_df = melt(dense_matrix)
  matrix_df$Var1 <- as.factor(matrix_df$Var1)
  matrix_df$Var2 <- as.factor(matrix_df$Var2)
  matrix_df$value <- factor(matrix_df$value, levels = sort(unique(matrix_df$value)))


  if (matrix_type == "H"){
    #p <- Matrix::image(dense_matrix[,ncol(dense_matrix):1], useRaster = TRUE, axes = F, col = c("white", "skyblue"), asp = 1)
    p <- ggplot(matrix_df, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("white", "skyblue")) +  # Customize colors if needed
    labs(x = "P", y = "B") +      # Add axis labels
    coord_fixed(ratio = 1) + # set aspect ratio
    theme_minimal() + 
    guides(fill = FALSE) +
    theme(axis.text.x = element_text(size = 24),  # Set x-axis label font size
        axis.text.y = element_text(size = 24),  # Set y-axis label font size
        axis.title = element_text(size = 24))
  }
  else {
    #p <- Matrix::image(dense_matrix[,ncol(dense_matrix):1], useRaster = TRUE, axes = F, col = c("white", "darkgrey"), asp = 1)
    p <- ggplot(matrix_df, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("white", "darkgrey")) +  # Customize colors if needed
    labs(x = "P", y = "P") +      # Add axis labels
    coord_fixed(ratio = 1) + # set aspect ratio
    theme_minimal() + 
    guides(fill = FALSE) +
    theme(axis.text.x = element_text(size = 24),  # Set x-axis label font size
        axis.text.y = element_text(size = 24),  # Set y-axis label font size
        axis.title = element_text(size = 24))
    
    }
  

  # Return the matrix and plot
  return(list(matrix = dense_matrix, plot = p))
}

# Function that generates a square matrix (non-zero value = 1) based on the specified size and modularity (full 2 modules, with modularity as the proportion of size for the top module)
full_bimodular_sq_matrix <- function(size_matrix, modularity, row_lab, col_lab, nz_color) {
  # Create an empty matrix
  matrix <- matrix(0, nrow = size_matrix, ncol = size_matrix)
  
  # Calculate the number of elements in each module
  module_size <- floor(size_matrix * modularity)
  
  # Populate the matrix
  for (i in 1:size_matrix) {
    for (j in 1:size_matrix) {
      if (floor((i - 1) / module_size) == floor((j - 1) / module_size)) {
        matrix[i, j] <- 1
      }
    }
  }
    # Covert the dense matrix into dataframe for plotting with ggplot()
    matrix_df = melt(matrix)
    matrix_df$Var1 <- as.factor(matrix_df$Var1)
    matrix_df$Var2 <- as.factor(matrix_df$Var2)
    matrix_df$value <- factor(matrix_df$value, levels = sort(unique(matrix_df$value)))
    #p <- Matrix::image(dense_matrix[,ncol(dense_matrix):1], useRaster = TRUE, axes = F, col = c("white", "darkgrey"), asp = 1)
    p <- ggplot(matrix_df, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("white", nz_color)) +  # Customize colors if needed
    labs(x = col_lab, y = row_lab) +      # Add axis labels
    coord_fixed(ratio = 1) + # set aspect ratio
    theme_minimal() +
    guides(fill = FALSE) +
    theme(axis.text.x = element_text(size = 24),  # Set x-axis label font size
        axis.text.y = element_text(size = 24),  # Set y-axis label font size
        axis.title = element_text(size = 24))

  # Return the matrix and plot
  return(list(matrix = matrix, plot = p))
}



# Function to generate a non-structured, sparse, square matrix (non-zero value = 1) with desired matrix size and sparsity (i.e. density of non-zero elements; between 0 & 1)
sparse_sq_matrix <- function(size_matrix, sparsity, row_lab, col_lab, nz_color, symmetric = TRUE) {

  # Generate a sparse matrix with random non-zero entries
  if (symmetric){
  matrix <- rsparsematrix(size_matrix, size_matrix, density = sparsity, symmetric = TRUE)
  } else {
  matrix <- rsparsematrix(size_matrix, size_matrix, density = sparsity, symmetric = FALSE)  
  }
  
  
  # Replace non-zero entries with 1s
  matrix@x <- rep(1, nnzero(matrix))

  # Convert sparse matrix to a dense matrix for visualization
  dense_matrix <- as.matrix(matrix)

  # Set diagonal element to 0 for matrix P
  if (col_lab == "P"){
    diag(dense_matrix) <- 0
  }


  # Covert the dense matrix into dataframe for plotting with ggplot()
    matrix_df = melt(dense_matrix)
    matrix_df$Var1 <- as.factor(matrix_df$Var1)
    matrix_df$Var2 <- as.factor(matrix_df$Var2)
    matrix_df$value <- factor(matrix_df$value, levels = sort(unique(matrix_df$value)))
    
    

  # Plot the dense matrix
  # p <- Matrix::image(dense_matrix[,ncol(dense_matrix):1], useRaster = TRUE, axes = F, col = c("white", "darkgrey"), asp = 1)
  p <- ggplot(matrix_df, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("white", nz_color)) +  # Customize colors if needed
    labs(x = col_lab, y = row_lab) +      # Add axis labels
    coord_fixed(ratio = 1) + # set aspect ratio
    theme_minimal() + 
    guides(fill = FALSE) +
    theme(axis.text.x = element_text(size = 24),  # Set x-axis label font size
        axis.text.y = element_text(size = 24),  # Set y-axis label font size
        axis.title = element_text(size = 24))

  # Return the matrix and plot
  return(list(matrix = dense_matrix, plot = p))
}

# Function to generate a non-structured, sparse matrix (non-zero value = 1) with desired matrix size and sparsity (i.e. density of non-zero elements)
sparse_matrix_custom <- function(row_matrix, col_matrix, sparsity, row_lab, col_lab, nz_color) {
  # Generate a sparse matrix with random non-zero entries
  matrix <- rsparsematrix(row_matrix, col_matrix, density = sparsity, symmetric = FALSE)
  
  # Replace non-zero entries with 1s
  matrix@x <- rep(1, nnzero(matrix))

  # Convert sparse matrix to a dense matrix for visualization
  dense_matrix <- as.matrix(matrix)

  # Covert the dense matrix into dataframe for plotting with ggplot()
    matrix_df = melt(dense_matrix)
    matrix_df$Var1 <- as.factor(matrix_df$Var1)
    matrix_df$Var2 <- as.factor(matrix_df$Var2)
    matrix_df$value <- factor(matrix_df$value, levels = sort(unique(matrix_df$value)))

  p <- ggplot(matrix_df, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("white", nz_color)) +  # Customize colors if needed
    labs(x = col_lab, y = row_lab) +      # Add axis labels
    coord_fixed(ratio = 1) + # set aspect ratio
    theme_minimal() + 
    guides(fill = FALSE) +
    theme(axis.text.x = element_text(size = 24),  # Set x-axis label font size
        axis.text.y = element_text(size = 24),  # Set y-axis label font size
        axis.title = element_text(size = 24))

  # Return both the matrix and the plot
  return(list(matrix = dense_matrix, plot = p))
}

# Function to check modularity of a given matrix (using package igraph)
ajm_modularity <- function(adjacency_matrix) {
# Convert adjacency matrix to directed graph
edge_list <- which(adjacency_matrix == 1, arr.ind = TRUE)
edges <- edge_list[, 2:1]  # Swap columns to have from-to format
graph <- graph_from_edgelist(as.matrix(edges), directed = TRUE)  

# Compute modularity
wtc <- cluster_walktrap(graph) # find densely connected subgraphs
modularity <- modularity(wtc)
return(modularity)
}


# Function to detect the largest value of a column named "key" among all .csv files in a given folder
last_key <- function(folder_path) {
  # List all CSV files in the folder
  csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  # Check if there are any CSV files in the folder
  if (length(csv_files) == 0) {
    max_key <- 0
  }
  
  else{
  # Initialize a variable to store the maximum key value
  max_key <- 0
  
  # Read each CSV file and extract the maximum value from the "key" column
  for (file in csv_files) {
    # Read the CSV file
    data <- read.csv(file)
    
    # Extract the maximum value from the "key" column
    max_key_in_file <- max(data$key, na.rm = TRUE)
    
    # Update the overall maximum if needed
    if (max_key_in_file > max_key) {
      max_key <- max_key_in_file
    }
  }}
  
  # Return the overall maximum key value
  return(max_key)
}


# Function to plot a adjacency matrix using ggplot2
# matrix_df = melt(adjacency_matrix)
# matrix_df$Var1 <- as.factor(matrix_df$Var1)
# matrix_df$Var2 <- as.factor(matrix_df$Var2)
# matrix_df$value <- factor(matrix_df$value, levels = sort(unique(matrix_df$value)))

# ggplot(matrix_df, aes(x = Var2, y = Var1, fill = value)) +
#   geom_tile(color = "white") +
#   scale_fill_manual(values = c("white", "skyblue")) +  # Customize colors if needed
#   labs(x = "P", y = "B") +      # Add axis labels
#   coord_fixed(ratio = 1) + # set aspect ratio
#   theme_minimal() + 
#   guides(fill = FALSE) +
#   theme(axis.text.x = element_text(size = 24),  # Set x-axis label font size
#         axis.text.y = element_text(size = 24),  # Set y-axis label font size
#         axis.title = element_text(size = 24))   # Set axis title font size
#         #legend.text = element_text(size = 24),  # Set legend text font size
#         #legend.title = element_text(size = 24)) # Set legend title font size