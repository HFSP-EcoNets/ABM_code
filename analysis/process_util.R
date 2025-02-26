### This is the R-script that holds util funcions for data processing ###
## Load packages and customized scripts 
library(RSQLite) # connect to SQLite DB
library(tidyverse) # for data processing dplyr, readr, ggplot2, etc.
library(scales)
library(dbplyr) # dplyr backend for DBs
#library(sqldf) # provides an easy way to perform SQL selects on R data frames
library(lemon) # facet_grid and facet_wrap, but with axis lines preserved on all panels
library(knitr) # integrate computitng and reporting (see https://sachsmc.github.io/knit-git-markr-guide/knitr/knit.html)
library(RColorBrewer) # create desired color palette
library(reshape) # reshape data into long format
library(reshape2)
library(plotly) # package for heatmap
library(reticulate) # package for exporting plots generated using package plotly
library(ggplot2)
library(grid)
library(cowplot) # this automatically incluse package grid
library(microViz) # for discrete color palette of n > 20
library(gganimate)# for animated plots
library(transformr) # for animated plots
library(gifski) # for gif output
library(magick)  # Load magick package for working with GIF files
library(ggplotify)
#library(viridis) # for heatmap colors
#library(hrbrthemes) # for heatmap theme
library(av) # for generating auto video
library(base64enc)
library(covr)
library(htmltools)
library(knitr)
library(ragg)
library(rmarkdown)
library(sf)
library(testthat)
library(tibble)
library(infomapecology)
library(bipartite)
library(DT)
library(vegan)
library(Ternary)

# run the following command to install Kaleido library on your OS (also required for exporting plots generated using package plotly) 
# reticulate::install_miniconda()
# reticulate::conda_install('r-reticulate', 'python-kaleido')
# reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')


## Function that convert p_profile from string to vector ##
string_to_factor <- function(stringcolumn){
  stringcolumn <- gsub("\\[|\\]", "", stringcolumn)  # Remove square brackets
  factor_vector <- as.factor(gsub(",", "", stringcolumn))# Remove the commas and keep the element as factors
  #numeric_vector <- as.integer(factor_vector) # Convert factors to integers; return this for infection state
  return(factor_vector)
}


## Function that create interaction table ##
interaction_table <- function(result, n_p){
  # construct vectors for the table of interaction
  vector_t <- c()
  vector_p_profile <- c()
  vector_strain <- c()
  vector_plasmid_id <- c()
  vector_copy <- c()
  
  # fill in the table of interaction
  for (i in seq_along(result$t)){
    plasmid_profile <- as.numeric(strsplit(as.character(result$p_profile[i]), " ")[[1]]) # convert plasmid profile from factor (e.g. 01) to numeric vector c(0,1))
    for (j in 1:n_p){
      # create new row  
      vector_t <- c(vector_t, result$t[i])
      vector_p_profile <- c(vector_p_profile, result$p_profile[i])
      vector_strain <- c(vector_strain, result$strain[i])
      vector_plasmid_id <- c(vector_plasmid_id, j)
      vector_copy <- c(vector_copy, as.integer(plasmid_profile[j]))
    }
  }
  
  # construct the table of interaction
  interaction <- data.frame(t = vector_t, strain = vector_strain, p_profile = vector_p_profile, plasmid_id = vector_plasmid_id, copy = vector_copy)
  return(interaction)
}

## Function creating data subsets of data based on sampled time ##
data_sub <- function(time, data){
  subset(data, t == time)  
}

## Functions for to calculate means of desired variables from replicates for each of the treatments

# Function for time and numbers of community collapse
community_collapse <- function(list_of_results, t_final){
  I <- character(0)
  P <- character(0)
  collapse <- character(0)
  time_collapse <- numeric(0)
  replicate <- numeric(0)
  for (i in 1:length(list_of_results)){
    # set file name for saving
    I <- c(I, list_of_results[[i]]$group_I[1])
    P <- c(P, list_of_results[[i]]$group_P[1])
    replicate <- c(replicate, list_of_results[[i]]$replicate[1])
    time_collapse <- c(time_collapse, max(list_of_results[[i]]$t))
    
    if (max(list_of_results[[i]]$t) >= t_final){
      collapse <- c(collapse, "false")} else {
        collapse <- c(collapse, "true")}
  }
  df <- data.frame(group_I = I, group_P = P, rep = replicate, full_extiction = collapse, extinction_time = time_collapse)
  df_extinct <- subset(df, full_extiction == "true")
  
  if (nrow(df_extinct) == 0){
    return("No system collapse.")
  } else{summary_stats <- df_extinct %>%
    group_by(group_I, group_P) %>%
    summarise(mean = mean(extinction_time), se = sd(extinction_time)/sqrt(n()), count = n())
  
  # Plot means with standard deviations and notations of counts for each combination of factors
  p <- ggplot(summary_stats, aes(group_I, mean, color = group_P)) +
    geom_point(position=position_dodge(.9)) +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) +
    geom_text(aes(label = count), position=position_dodge(.9), vjust = -1, hjust = 0, size = 5, show.legend = FALSE) +
    ylim(0, t_final) +
    labs(x="I", y = "Extinction time (mean ± 1 se)", color = "P") + 
    paper_figs_theme
  return(p)}
    
}

# Function for overall community size at the end of simulation
end_community_size <- function(list_of_results, t_output){
  # combine the list of results
  df <- bind_rows(list_of_results)
  
  # collect number of replicates
  # n_rep = max(df$replicate)
  
  # unify time unit
  # df$t <- round(df$t)
  df$t <- floor(df$t/t_output) * t_output
  
  
  # extract the final t
  t_max = max(df$t)
  df_sub <- subset(df, t == t_max)
  
  # group dependent variable by treatment groups and strain_id, then then calculate mean strain abundance
  community_abundance <- df_sub %>%
    group_by(group_I, group_P, replicate) %>%
    summarize(abundance_sum = sum(abundance)) 
  
  # iterate over combinations of group_I and group_P to add missing replicate (with sum abundance = 0) due to full extinction
  for (i in unique(community_abundance$group_I)) {
    for (j in unique(community_abundance$group_P)) {
      for (k in 1:n_rep){
        # Check if the missing rep value exists in the combination
        if (!(k %in% community_abundance$replicate[community_abundance$group_I == i & community_abundance$group_P == j])) {
          # Add a row with the missing rep value
          new_row <- data.frame(group_I = i,
                                group_P = j,
                                replicate = k,
                                abundance_sum = 0)
          # Append the new row to the new_rows dataframe
          community_abundance <- bind_rows(community_abundance, new_row)
        }
      }
    }
  }
  
  # do the stats
  summary_stats <- community_abundance %>%
    group_by(group_I, group_P) %>%
    summarize(sum_abundance_mean = mean(abundance_sum), se = sd(abundance_sum)/sqrt(n()))
  
  # plot means with standard deviations and notations of counts for each combination of factors
  p <- ggplot(summary_stats, aes(group_I, sum_abundance_mean, color = group_P)) +
    geom_point(position=position_dodge(.9)) +
    geom_errorbar(aes(ymin=sum_abundance_mean-se, ymax=sum_abundance_mean+se), width=.2, position=position_dodge(.9)) +
    #ylim(17000, 30000) +
    labs(x="I", y = "End community abundance (mean ± 1 se)", color = "P") + 
    paper_figs_theme
  return(p)  
}

# Function for community count richness at a given time (count of subpopulation: S)
community_richness_count<- function(list_of_results, time, t_output){
  # combine the list of results
  df <- bind_rows(list_of_results)
  
  # collect number of replicates
  n_rep = max(df$replicate)
  
  # unify time unit
  # df$t <- round(df$t)
  df$t <- floor(df$t/t_output) * t_output # note: 5 is the t_output
  
  # extract data subset of at given time
  df_sub <- subset(df, t == time & abundance > 0)
  
  # group dependent variable by treatment groups, then then calculate mean subpopulation richness
  # note: we take the collapse replicate into consideration, with a richness of 0
  subpop_richness <- df_sub %>%
    group_by(group_I, group_P, replicate) %>%
    mutate(richness = length(abundance)) %>%
    slice(1) %>%
    select(-experiment, -t, -subpop_id, -strain, -p_profile, -abundance) %>%
    ungroup()
  
  for (i in unique(subpop_richness$group_I)) {
    for (j in unique(subpop_richness$group_P)) {
      for (k in 1:n_rep){
        # Check if the missing rep value exists in the combination
        if (!(k %in% subpop_richness$replicate[subpop_richness$group_I == i & subpop_richness$group_P == j])) {
          # Add a row with the missing rep value
          new_row <- data.frame(group_I = i,
                                group_P = j,
                                replicate = k,
                                richness = 0)
          # Append the new row to the new_rows dataframe
          subpop_richness <- bind_rows(subpop_richness, new_row)
        }
      }
    }
  }
  
  # do the stats
  summary_stats <- subpop_richness %>%
    group_by(group_I, group_P) %>%
    summarize(mean_richness = mean(richness), se = sd(richness)/sqrt(n()))

  
  # plot means with standard deviations and notations of counts for each combination of factors
  title <- sprintf("t = %d", time)  # generate title
  p <- ggplot(summary_stats, aes(group_I, mean_richness, color = group_P)) +
    geom_point(position=position_dodge(.9)) +
    geom_errorbar(aes(ymin=mean_richness-se, ymax=mean_richness+se), width=.2, position=position_dodge(.9)) +
    ylim(0, max(summary_stats$mean_richness)+1) +
    labs(title = title, x="I", y = "Richness (mean ± 1 se)", color = "P") + 
    paper_figs_theme
  return(p)  

}

# Function for community margalef richness at a given time (count & abundance of subpopulation: S/logN)
community_richness_margalef<- function(list_of_results, time, t_output){
  # combine the list of results
  df <- bind_rows(list_of_results)
  
  # collect number of replicates
  n_rep = max(df$replicate)
  
  # unify time unit
  # df$t <- round(df$t)
  df$t <- floor(df$t/t_output) * t_output # note: 5 is the t_output
  
  # extract data subset of at given time
  df_sub <- subset(df, t == time & abundance > 0)
  
  # group dependent variable by treatment groups, then then calculate mean subpopulation richness
  # note: we take the collapse replicate into consideration, with a richness of 0
  subpop_richness_margalef <- df_sub %>%
    group_by(group_I, group_P, replicate) %>%
    mutate(richness = length(abundance)/log(sum(abundance))) %>%
    slice(1) %>%
    select(-experiment, -t, -subpop_id, -strain, -p_profile, -abundance) %>%
    ungroup()
  
  for (i in unique(subpop_richness_margalef$group_I)) {
    for (j in unique(subpop_richness_margalef$group_P)) {
      for (k in 1:n_rep){
        # Check if the missing rep value exists in the combination
        if (!(k %in% subpop_richness_margalef$replicate[subpop_richness_margalef$group_I == i & subpop_richness_margalef$group_P == j])) {
          # Add a row with the missing rep value
          new_row <- data.frame(group_I = i,
                                group_P = j,
                                replicate = k,
                                richness = 0)
          # Append the new row to the new_rows dataframe
          subpop_richness_margalef <- bind_rows(subpop_richness_margalef, new_row)
        }
      }
    }
  }
  
  # do the stats
  summary_stats <- subpop_richness_margalef %>%
    group_by(group_I, group_P) %>%
    summarize(mean_richness = mean(richness), se = sd(richness)/sqrt(n()))
  
  
  # plot means with standard deviations and notations of counts for each combination of factors
  title <- sprintf("Community subpopulation Margalef richness at t = %d", time)  # generate title
  p <- ggplot(summary_stats, aes(group_I, mean_richness, color = group_P)) +
    geom_point(position=position_dodge(.9)) +
    geom_errorbar(aes(ymin=mean_richness-se, ymax=mean_richness+se), width=.2, position=position_dodge(.9)) +
    #ylim(17000, 30000) +
    labs(title = title, x="I", y = "Richness (mean ± 1 se)", color = "P") + 
    paper_figs_theme
  return(p)  
  
}


# Function for community composition at end of simulation
community_composition <- function(list_of_results, t_output){
  # combine the list of results
  df <- bind_rows(list_of_results)
  
  # collect number of replicates
  n_rep = max(df$replicate)
  
  # unify time unit
  # df$t <- round(df$t)
  df$t <- floor(df$t/t_output) * t_output
  
  # in case the floor does not belong to the fold of t_output 
  
  # extract the data subset from final t
  t_max = max(df$t)
  df_sub <- subset(df, t == t_max)
  
  # group dependent variable by treatment groups and strain_id, then then calculate mean strain abundance
  strain_abundance <- df_sub %>%
    group_by(group_I, group_P, strain) %>%
    summarize(abundance_sum = sum(abundance)/n_rep)
  
  # add a combined group (group I; group P) 
  strain_abundance$group = paste(strain_abundance$group_I, "I", "*", strain_abundance$group_P, "P")
  
  # rename group name for experiment part II
  if("strong" %in% strain_abundance$group_I){
    for(i in 1:nrow(strain_abundance)){
      strain_abundance$group[i] = ifelse(strain_abundance$group_P[i] == "strong", "negative", 
                                         ifelse(strain_abundance$group_P[i] == "weak", "positive")) 
    }
  }
  
  # check the data frame (note that strain data of a treatment group could be missing due to community extinction)
  #strain_abundance
  
  
  # calculate mean strain percentage
  strain_percentage <- strain_abundance %>%
    group_by(group_I, group_P) %>%
    mutate(strain_percentage = round(abundance_sum / sum(abundance_sum) * 100))
  
  # set column type before plotting
  strain_abundance$strain <- as.factor(strain_abundance$strain) 
  strain_abundance$group <- as.factor(strain_abundance$group)
  colnames(strain_abundance)[colnames(strain_abundance) == "strain"] <- "Strain" # change column name
  
  # create pie charts
  p <- strain_abundance %>%
    ggplot(aes(x=" ", y=abundance_sum, group=group, fill=Strain)) + 
    geom_bar(width = 1, stat = "identity", position = position_fill()) + 
    # geom_text(aes(label = abundance_sum,), position = position_fill(vjust = 0.5)) +
    coord_polar("y", start=0) +  
    facet_wrap(~ group, ncol = 2) + 
    labs(title = "Community Composition at End of Simulation") + 
    paper_figs_pie_theme +
    theme(strip.text = element_text(size = 12))
  
  # or
  #strain_percentage %>%
  #    ggplot(aes(x = "", y=strain_percentage, group=group, colour=Strain, fill=Strain)) + 
  #    geom_bar(stat="identity", width=1) +
  #    coord_polar("y", start=0) +
  #    facet_wrap(~ group, ncol = 2) + 
  #    geom_text(aes(label = paste0(round(strain_percentage), "%")), 
  #              position = position_stack(vjust = 0.5), 
  #              color = "white", size = 5) +
  #    labs(title = "Community Composition at End of Simulation") + 
  #    paper_figs_pie_theme
  
  return(p)
  
}

# Function for population composition at end of simulation
population_composition <- function(list_of_results, t_output){
  # combine the list of results
  df <- bind_rows(list_of_results)
  
  # collect number of replicates
  n_rep = max(df$replicate)
  
  # unify time unit
  # df$t <- round(df$t)
  df$t <- floor(df$t/t_output) * t_output
  
  # extract the data subset from final t
  t_max = max(df$t)
  df_sub <- subset(df, t == t_max)
  
  df_sub$p_profile <- string_to_factor(df_sub$p_profile)
  
  # remove rows with value 0
  df_sub = subset(df_sub, abundance > 0)
  
  # group dependent variable by treatment groups, strain and p_profile, then calculate mean subpopulation abundance
  subpopulation_abundance <- df_sub %>%
    group_by(group_I, group_P, strain, p_profile) %>%
    summarize(abundance_sub = sum(abundance)/n_rep)
  
  # add a combined group (group I; group P) 
  subpopulation_abundance$group = paste(subpopulation_abundance$group_I, "I", "*", subpopulation_abundance$group_P, "P")
  
  # rename group name for experiment part II
  if("strong" %in% subpopulation_abundance$group_I){
    for(i in 1:nrow(subpopulation_abundance)){
      subpopulation_abundance$group[i] = ifelse(subpopulation_abundance$group_P[i] == "strong", "negative", 
                                                ifelse(subpopulation_abundance$group_P[i] == "weak", "positive")) 
    }
  }
  
  # set column type before plotting
  subpopulation_abundance$group = as.factor(subpopulation_abundance$group)
  subpopulation_abundance$strain = as.factor(subpopulation_abundance$strain) # convert strain ID from type integer to type string
  subpopulation_abundance$p_profile = as.factor(subpopulation_abundance$p_profile) # convert infection state from type integer to type string
  colnames(subpopulation_abundance)[colnames(subpopulation_abundance) == "strain"] <- "Strain" # change column name
  colnames(subpopulation_abundance)[colnames(subpopulation_abundance) == "abundance_sub"] <- "Abundance" # change column name
  colnames(subpopulation_abundance)[colnames(subpopulation_abundance) == "p_profile"] <- "Profile" # change column name
  
  # create grouped stacked bar plots
  p <- subpopulation_abundance %>%
    ggplot(aes(x = Strain, y = Abundance, fill = Profile)) + 
    geom_bar(stat = 'identity', position = 'fill') + 
    facet_wrap(~ group, ncol = 2) +
    scale_fill_brewer(palette = "Dark2") + 
    labs(title = "Population Composition at End of Simulation", y = "Dominance") +
    paper_figs_theme +
    theme(strip.text = element_text(size = 14))
  
 return(p)
  
}

# Function for population dynamics
population_dynamics <- function(list_of_results, t_output){
  # combine the list of results
  df <- bind_rows(list_of_results)
  
  # collect number of replicates
  n_rep = max(df$replicate)
  
  # unify time unit
  # df$t <- round(df$t)
  df$t <- floor(df$t/t_output) * t_output
  
  df$p_profile <- string_to_factor(df$p_profile)
  
  
  # group dependent variable by group_I, group_P, t and strain_id, then calculate mean abundance
  strain_t <- df %>%
    group_by(group_I, group_P, t, strain) %>%
    summarize(abundance_sum = sum(abundance)/n_rep)
  
  # add a combined group (group I; group P) 
  strain_t$group = paste(strain_t$group_I, "I", "*", strain_t$group_P, "P")
  
  # rename group name for experiment part II
  if("strong" %in% strain_t$group_I){
    for(i in 1:nrow(strain_t)){
      strain_t$group[i] = ifelse(strain_t$group_P[i] == "strong", "negative", 
                                 ifelse(strain_t$group_P[i] == "weak", "positive")) 
    }
  }
  
  # set column type before plotting
  strain_t$strain = as.character(strain_t$strain) # convert strain ID from type integer to type string
  colnames(strain_t)[colnames(strain_t) == "strain"] <- "Strain" # change column name
  
  # create line plots
  p <- strain_t %>%
    ggplot(aes(x = t, y = abundance_sum, color = Strain)) +
    geom_line() +
    facet_wrap(~ group, scales = "fixed") +
    labs(title="Bacterial Population Dynamics", x="t", y = "Abundance") +
    paper_figs_theme +
    theme(strip.text = element_text(size = 14))
  
  return(p)
}

# Function for population dynamics (stacked area graph, treatment group-specific)
population_dynamics_stacked <- function(list_of_results, trt_group, t_output){
  # combine the list of results
  df <- bind_rows(list_of_results)
  
  # collect number of replicates
  n_rep = max(df$replicate)
  
  # unify time unit
  # df$t <- round(df$t)
  df$t <- floor(df$t/t_output) * t_output
  
  df$p_profile <- string_to_factor(df$p_profile)
  
  
  # group dependent variable by group_I, group_P, t and strain_id, then calculate mean abundance
  strain_t <- df %>%
    group_by(group_I, group_P, t, strain) %>%
    summarize(abundance_sum = sum(abundance)/n_rep)
  
  # add a combined group (group I; group P) 
  strain_t$group = paste(strain_t$group_I, "I", "*", strain_t$group_P, "P")
  
  # rename group name for experiment part II
  if("strong" %in% strain_t$group_I){
    for(i in 1:nrow(strain_t)){
      strain_t$group[i] = ifelse(strain_t$group_P[i] == "strong", "negative", 
                                 ifelse(strain_t$group_P[i] == "weak", "positive")) 
    }
  }
  # extract from the treatment group of choice
  strain_t_sub = subset(strain_t, group == trt_group)
  
  # set column type before plotting
  strain_t_sub$strain = as.character(strain_t_sub$strain) # convert strain ID from type integer to type string
  colnames(strain_t_sub)[colnames(strain_t_sub) == "strain"] <- "Strain" # change column name
  
  # create stacked area plots
  title <- sprintf("Bacterial Population Dynamics: %s", trt_group)  # generate title
  
  p <- strain_t_sub %>%
    ggplot(aes(x = t, y = abundance_sum, fill = Strain)) +
    geom_area() +
    labs(title=title, x="t", y = "Abundance") +
    paper_figs_theme
  
  return(p)
}

# Function for population dynamics with animated vertical line moving with time (stacked area graph, treatment group-specific)
population_dynamics_stacked_ani <- function(list_of_results, trt_group, t_final, t_output){
  # combine the list of results
  df <- bind_rows(list_of_results)
  
  # collect number of replicates
  n_rep = max(df$replicate)
  
  # unify time unit
  # df$t <- round(df$t)
  df$t <- floor(df$t/t_output) * t_output
  
  df$p_profile <- string_to_factor(df$p_profile)
  
  
  # group dependent variable by group_I, group_P, t and strain_id, then calculate mean abundance
  strain_t <- df %>%
    group_by(group_I, group_P, t, strain) %>%
    summarize(abundance_sum = sum(abundance)/n_rep)
  
  # add a combined group (group I; group P) 
  strain_t$group = paste(strain_t$group_I, "I", "*", strain_t$group_P, "P")
  
  # rename group name for experiment part II
  if("strong" %in% strain_t$group_I){
    for(i in 1:nrow(strain_t)){
      strain_t$group[i] = ifelse(strain_t$group_P[i] == "strong", "negative", 
                                 ifelse(strain_t$group_P[i] == "weak", "positive")) 
    }
  }
  # extract from the treatment group of choice
  strain_t_sub = subset(strain_t, group == trt_group)
  
  # set column type before plotting
  strain_t_sub$strain = as.character(strain_t_sub$strain) # convert strain ID from type integer to type string
  colnames(strain_t_sub)[colnames(strain_t_sub) == "strain"] <- "Strain" # change column name
  
  # create stacked area plots
  title <- sprintf("Bacterial Population Dynamics: %s", trt_group)  # generate title
  
  # create a vector for the vertical line positions
  line_positions <- data.frame(time = unique(strain_t_sub$t))
  
  # set file name for saving
  name <- sprintf("population_dynamics_%s.gif", trt_group)
  
  # set variables that will later be used in transition_states()
  s_length = t_final/t_output + 1
  
  p <- strain_t_sub %>%
    ggplot(aes(x = t, y = abundance_sum, fill = Strain), stat_indentity()) +
    geom_area() +
    labs(title=title, x="t", y = "Abundance") +
    paper_figs_theme +
    theme(text = element_text(size = 18),
          plot.title = element_text(size = 16),
          plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_vline(data = line_positions, aes(xintercept = time), col = "black", linetype = "dotted", size = 0.5) + 
    # below the gganimate specific bits
    transition_states(line_positions$time, state_length = s_length, transition_length = c(rep(5, s_length))) +
    shadow_mark(past = TRUE, future = TRUE, exclude_layer = 2)
  
  #animated_p <- animate(p, fps=5)
  # anim_save(name, animation = animated_p) # save plot into .gif file; optional
  return(animate(p, fps=5, width = 760, height = 300))
}

# Function for population dominance dynamics (stacked area graph, treatment group-specific)
population_dominance_dynamics_stacked <- function(list_of_results, trt_group, t_output){
  # combine the list of results
  df <- bind_rows(list_of_results)
  
  # collect number of replicates
  n_rep = max(df$replicate)
  
  # unify time unit
  # df$t <- round(df$t)
  df$t <- floor(df$t/t_output) * t_output
  
  df$p_profile <- string_to_factor(df$p_profile)
  
  
  # group dependent variable by group_I, group_P, t and strain_id, then calculate mean abundance
  strain_t <- df %>%
    group_by(group_I, group_P, t, strain) %>%
    summarize(abundance_sum = sum(abundance)/n_rep) %>%
    mutate(percentage = abundance_sum / sum(abundance_sum))
  
  # add a combined group (group I; group P) 
  strain_t$group = paste(strain_t$group_I, "I", "*", strain_t$group_P, "P")
  
  # rename group name for experiment part II
  if("strong" %in% strain_t$group_I){
    for(i in 1:nrow(strain_t)){
      strain_t$group[i] = ifelse(strain_t$group_P[i] == "strong", "negative", 
                                 ifelse(strain_t$group_P[i] == "weak", "positive")) 
    }
  }
  # extract from the treatment group of choice
  strain_t_sub = subset(strain_t, group == trt_group)
  
  # set column type before plotting
  strain_t_sub$strain = as.character(strain_t_sub$strain) # convert strain ID from type integer to type string
  colnames(strain_t_sub)[colnames(strain_t_sub) == "strain"] <- "Strain" # change column name
  
  # create stacked area plots
  title <- sprintf("Bacterial Population Dynamics: %s", trt_group)  # generate title
  
  p <- strain_t_sub %>%
    ggplot(aes(x = t, y = percentage, fill = Strain)) +
    geom_area() +
    labs(title=title, x="t", y = "Dominance") +
    paper_figs_theme
  
  return(p)
}

# Function for population dominance dynamics (stacked area graph, all groups)
population_dominance_dynamics_stacked_all <- function(list_of_results, t_output){
   
  # combine the list of results
  df <- bind_rows(list_of_results)
  
  # collect number of replicates
  n_rep = max(df$replicate)
  
  # unify time unit
  # df$t <- round(df$t)
  df$t <- floor(df$t/t_output) * t_output
  
  df$p_profile <- string_to_factor(df$p_profile)
  
  
  # group dependent variable by group_I, group_P, t and strain_id, then calculate mean abundance
  strain_t <- df %>%
    group_by(group_I, group_P, t, strain) %>%
    summarize(abundance_sum = sum(abundance)/n_rep) %>%
    mutate(percentage = abundance_sum / sum(abundance_sum))
  
  # add a combined group (group I; group P) 
  strain_t$group = paste(strain_t$group_I, "I", "*", strain_t$group_P, "P")
  
  # rename group name for experiment part II
  if("strong" %in% strain_t$group_I){
    for(i in 1:nrow(strain_t)){
      strain_t$group[i] = ifelse(strain_t$group_P[i] == "strong", "negative", 
                                 ifelse(strain_t$group_P[i] == "weak", "positive")) 
    }
  }
  
  list_groups = c("modular I * modular P", "modular I * uniform P", "uniform I * modular P", "uniform I * uniform P")
  
  # generate plots for each treatment group
  list_plots = list()
  for(i in seq_along(list_groups)){
    # extract from the treatment group of choice
    strain_t_sub = subset(strain_t, group == list_groups[i])
    
    # set column type before plotting
    strain_t_sub$strain = as.character(strain_t_sub$strain) # convert strain ID from type integer to type string
    colnames(strain_t_sub)[colnames(strain_t_sub) == "strain"] <- "Strain" # change column name
    
    # create stacked area plots
    title <- sprintf("%s", list_groups[i])  # generate title
    
    assign(paste0("p", i), strain_t_sub %>%
             ggplot(aes(x = t, y = percentage, fill = Strain)) +
             geom_area() +
             labs(title=title, x="t", y = "Dominance") +
             paper_figs_theme)
  }
  
  return(cowplot::plot_grid(p1, p2, p3, p4, labels = c('A', 'B', 'C', 'D'), label_size = 12))
  
}

population_dominance_dynamics_stacked_all2 <- function(list_of_results, t_output){
  # combine the list of results
  df <- bind_rows(list_of_results)
  
  # collect number of replicates
  n_rep = max(df$replicate)
  
  # unify time unit
  # df$t <- round(df$t)
  df$t <- floor(df$t/t_output) * t_output
  
  df$p_profile <- string_to_factor(df$p_profile)
  
  
  # group dependent variable by group_I, group_P, t and strain_id, then calculate mean abundance
  strain_t <- df %>%
    group_by(group_I, group_P, t, strain) %>%
    summarize(abundance_sum = sum(abundance)/n_rep) %>%
    mutate(percentage = abundance_sum / sum(abundance_sum))
  
  # add a combined group (group I; group P) 
  strain_t$group = paste(strain_t$group_I, "I", "*", strain_t$group_P, "P")
  
  # rename group name for experiment part II
  if("strong" %in% strain_t$group_I){
    for(i in 1:nrow(strain_t)){
      strain_t$group[i] = ifelse(strain_t$group_P[i] == "strong", "negative", 
                                 ifelse(strain_t$group_P[i] == "weak", "positive")) 
    }
  }
  
  # set column type before plotting
  strain_t$strain = as.character(strain_t$strain) # convert strain ID from type integer to type string
  colnames(strain_t)[colnames(strain_t) == "strain"] <- "Strain" # change column name
  
  
  # create line plots
  p <- strain_t %>%
    ggplot(aes(x = t, y = percentage, fill = Strain)) +
    geom_area() +
    facet_wrap(~ group, scales = "fixed") +
    labs(title="Bacterial Population Dynamics (Dominance)", x="t", y = "Dominance") +
    paper_figs_theme +
    theme(strip.text = element_text(size = 14))
  
  return(p)
  
}

# Function for subpopulation dynamics (strain-specific)
strain_subpopulation_dynamics <- function(list_of_results, strain_id, design, t_output, time){
  # combine the list of results
  df <- bind_rows(list_of_results)
  
  # collect number of replicates
  n_rep = max(df$replicate)
  
  # unify time unit
  # df$t <- round(df$t)
  df$t <- floor(df$t/t_output) * t_output
  
  df$p_profile <- string_to_factor(df$p_profile)
  
  # extract data subset of the specific strain 
  df_sub = subset(df, strain == strain_id & t <= time)
  
  # group dependent variable by group, t, strain id and p_profile, then calculate mean subpopulation abundance
  subpopulation_t <- df_sub %>%
    group_by(group_I, group_P, t, p_profile) %>%
    summarize(abundance_sub = sum(abundance)/n_rep)
  
  # add a combined group (group I; group P) 
  subpopulation_t$group = paste(subpopulation_t$group_I, "I", "*", subpopulation_t$group_P, "P")
  
  # rename group name for experiment part II
  if("strong" %in% subpopulation_t$group_I){
    for(i in 1:nrow(subpopulation_t)){
      subpopulation_t$group[i] = ifelse(subpopulation_t$group_P[i] == "strong", "negative", 
                                        ifelse(subpopulation_t$group_P[i] == "weak", "positive")) 
    }
  }
  
  # set column type before plotting
  subpopulation_t$p_profile = as.factor(subpopulation_t$p_profile) # convert infection state from type integer to type string
  colnames(subpopulation_t)[colnames(subpopulation_t) == "p_profile"] <- "Profile" # change column name
  
  # create line plots
  
  title <- sprintf("Bacterial Population Dynamics: Strain %s", strain_id)  # generate title
  
  
  # create elements for a fixed color bar
  binary_combinations <- expand.grid(replicate(design$n_p[1], 0:1, simplify = FALSE))
  all_categories <- apply(binary_combinations, 1, function(x) paste(strsplit(paste(x, collapse = ""), "")[[1]], collapse = " "))
  all_categories_factor <- factor(all_categories, levels = all_categories)
  
  brewerPlus <- distinct_palette() # generate color palette using function in package microViz
  fixed_colors <- brewerPlus[1:length(all_categories)]
  
  
  
  
  p <- subpopulation_t %>%
    ggplot(aes(x = t, y = abundance_sub)) +
    facet_rep_wrap(vars(group), ncol = 3, scales = "fixed", repeat.tick.labels = FALSE) +
    geom_line(aes(color = Profile)) +
    #scale_color_manual(values = brewerPlus) +
    #ylim(0,2000) +
    scale_color_manual(values = setNames(fixed_colors, all_categories)) +
    labs(title=title, x="t", y = "Abundance") +
    paper_figs_theme +
    theme(strip.text = element_text(size = 14))
  
  return(p)
  
}

# Function for subpopulation abundances across time (strain-separated, treatment-specific, at a given time period)
subpopulation_dynamics <- function(list_of_results, trt_group, time1, time2, design, t_output){
  # combine the list of results
  df <- bind_rows(list_of_results)
  
  # collect number of replicates
  n_rep = max(df$replicate)
  
  # unify time unit
  # df$t <- round(df$t)
  df$t <- floor(df$t/t_output) * t_output
  
  df$p_profile <- string_to_factor(df$p_profile)
  
  # add a combined group (group I; group P) 
  df$group = paste(df$group_I, "I", "*", df$group_P, "P")
  
  # rename group name for experiment part II
  if("strong" %in% df$group_I){
    for(i in 1:nrow(df)){
      df$group[i] = ifelse(df$group_P[i] == "strong", "negative", 
                                        ifelse(df$group_P[i] == "weak", "positive")) 
    }
  } 
  
  # extract data subset of the specific treatment 
  df_sub = subset(df, group == trt_group & time1 <= t & t < time2+1)
  
  # group dependent variable by t, strain id, and p_profile, then calculate mean subpopulation abundance
  subpopulation_t <- df_sub %>%
    group_by(t, strain, p_profile) %>%
    summarize(abundance_sub = sum(abundance)/n_rep)
  
  # set column type before plotting
  subpopulation_t$strain = as.character(subpopulation_t$strain) # convert strain ID from type integer to type string
  colnames(subpopulation_t)[colnames(subpopulation_t) == "strain"] <- "Strain" # change column name
  subpopulation_t$p_profile = as.factor(subpopulation_t$p_profile) # convert infection state from type integer to type string
  colnames(subpopulation_t)[colnames(subpopulation_t) == "p_profile"] <- "Profile" # change column name
  
  # create stacked area plots
  title <- sprintf("%s", trt_group)  # generate title
  
  # create elements for a fixed color bar
  binary_combinations <- expand.grid(replicate(design$n_p[1], 0:1, simplify = FALSE))
  all_categories <- apply(binary_combinations, 1, function(x) paste(strsplit(paste(x, collapse = ""), "")[[1]], collapse = " "))
  all_categories_factor <- factor(all_categories, levels = all_categories)
  
  brewerPlus <- distinct_palette() # generate color palette using function in package microViz
  fixed_colors <- brewerPlus[1:length(all_categories)]
  
  p <- subpopulation_t %>%
    ggplot(aes(x = t, y = abundance_sub, fill = Profile)) +
    facet_rep_wrap(vars(Strain), ncol = 1, scales = "fixed", repeat.tick.labels = FALSE) +
    geom_area() +
    #scale_fill_brewer(palette = "Dark2") +
    #scale_fill_manual(values = brewerPlus) +
    scale_fill_manual(values = setNames(fixed_colors, all_categories)) +
    labs(title=title, x="t", y = "Abundance") +
    #ylim(0, 1000) +
    paper_figs_theme +
    theme(strip.text = element_text(size = 14))
  
  dev.new()
  print(p)

  }

# Function for subpopulation abundance at a given t (strain- and treatment(experiment)-specific, with mean and se)
strain_subpopulation_t <- function(list_of_results, strain_id, design, exp, t_output, time){
  # combine the list of results
  df <- bind_rows(list_of_results)
  
  # collect number of replicates
  n_rep = max(df$replicate)
  
  # unify time unit
  # df$t <- round(df$t)
  df$t <- floor(df$t/t_output) * t_output
  
  df$p_profile <- string_to_factor(df$p_profile)
  
  # extract data subset of the specific strain & time
  df_sub = subset(df, experiment == exp & strain == strain_id & t == time)
  
  # remove unused column for simplicity
  df_sub = subset(df_sub, select = -subpop_id)
  
  # note: we take the collapse replicate & never-emerged subpoppulation into consideration, with a subpopulation abundance of 0 (to make all sample size equal)
  I <- df_sub$group_I[1]
  P <- df_sub$group_P[1]  
  for (i in unique(df_sub$p_profile)) {
      for (k in 1:n_rep){
        # Check if the missing rep value exists in the combination
        if (!(k %in% df_sub$replicate[df_sub$p_profile == i])) {
          # Add a row with the missing rep value
          new_row <- data.frame(experiment = exp,
                                group_I = I,
                                group_P = P,
                                replicate = k,
                                t = time,
                                strain = as.numeric(strain_id),
                                p_profile = i,
                                abundance = 0)
          # Append the new row to the new_rows dataframe
          df_sub <- bind_rows(df_sub, new_row)
        }
      }
    }
  
  # create data frame of statistics (mean & se of subpopulations) to plot
  summary_stats <- df_sub %>%
    group_by(p_profile) %>%
    summarise(mean = mean(abundance), se = sd(abundance/sqrt(n())), count = n())
  
  # set column type before plotting
  summary_stats$p_profile = as.factor(summary_stats$p_profile) # convert infection state from type integer to type string
  colnames(summary_stats)[colnames(summary_stats) == "p_profile"] <- "Profile" # change column name
  
  # add a column for category label
  summary_stats$Category <- c(1:length(summary_stats$Profile))
  
  # create line plots
  
  title <- sprintf("Bacterial Population Dynamics at t = %d: Strain %s", time, strain_id)  # generate title
  
  
  # create elements for a fixed color bar
  binary_combinations <- expand.grid(replicate(design$n_p[1], 0:1, simplify = FALSE))
  all_categories <- apply(binary_combinations, 1, function(x) paste(strsplit(paste(x, collapse = ""), "")[[1]], collapse = " "))
  all_categories_factor <- factor(all_categories, levels = all_categories)
  
  brewerPlus <- distinct_palette() # generate color palette using function in package microViz
  fixed_colors <- brewerPlus[1:length(all_categories)]
  
  # Plot means with standard deviations and notations of counts for each combination of factors
  p <- ggplot(summary_stats, aes(Category, mean, color = Profile)) +
    scale_color_manual(values = setNames(fixed_colors, all_categories)) +
    geom_point(position=position_dodge(.9)) +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) +
    geom_text(aes(label = count), position=position_dodge(.9), vjust = -1, hjust = 0, size = 5, show.legend = FALSE) +
    labs(title=title, x="Profile Index", y = "Abundance (mean ± 1 se)") +
    paper_figs_theme
  
  
  # check histogram (zero-cases included)
  # hist(df_sub$abundance[df_sub$p_profile == "1 0 0 0"])
  
  
  # check histogram (zero-cases excluded)
  #hist(df_sub$abundance[df_sub$p_profile == "1 0 0 0" & df_sub$abundance < 5000 & df_sub$abundance > 0])
  # hist(df_sub$abundance[df_sub$p_profile == "1 0 0 0" & df_sub$abundance > 5000])
  
  return(p)  
  
  
}



# Function for plasmid dynamics (strain-specific)
plasmid_dynamics <- function(list_of_results, strain_id, t_output){
  # combine the list of results
  df <- bind_rows(list_of_results)
  
  # collect number of replicates
  n_rep = max(df$replicate)
  
  # unify time unit
  # df$t <- round(df$t)
  df$t <- floor(df$t/t_output) * t_output
  
  df$p_profile <- string_to_factor(df$p_profile)
  
  # extract data subset of the specific strain 
  df_sub = subset(df, strain == strain_id)
  
  # group dependent variable by treatment groups, t, strain and p_profile, then calculate mean subpopulation abundance
  subpopulation_abundance <- df_sub %>%
    group_by(group_I, group_P, t, p_profile) %>%
    summarize(abundance_sub = sum(abundance)/n_rep)
  
  
  
  # add columns to report plasmid abundance of each plasmid strain
  column_names <- paste("P", 1:n_p, sep = "") # create a vector of names of new columns
  
  for (column_name in column_names){ 
    subpopulation_abundance[, column_name] <- rep(0, nrow(subpopulation_abundance))
  }# add the new columns into the dataframe
  
  
  for (i in 1:nrow(subpopulation_abundance)) {  # Go through rows 
    plasmid_profile <- as.numeric(strsplit(as.character(subpopulation_abundance$p_profile[i]), " ")[[1]]) # convert plasmid profile from factor (e.g. 01) to numeric vector c(0,1))
    p_abundance <- subpopulation_abundance$abundance_sub[i] * plasmid_profile
    subpopulation_abundance[i, 6:ncol(subpopulation_abundance)] <- as.list(p_abundance)
  } # fill in plasmid abundance
  
  
  # change column name of p_profile (for columns with name stared with P will be used later)
  colnames(subpopulation_abundance)[colnames(subpopulation_abundance) == "p_profile"] <- "fac_p_profile" 
  
  # group data frame by t and strain_id, then calculate plasmid prevalence
  strain_t_prevalence <- subpopulation_abundance %>%
    group_by(group_I, group_P, t) %>%
    summarize(
      across(starts_with("P"), ~ sum(.) / sum(abundance_sub), .names = "{.col}")
    )
  
  # add a combined group (group I; group P) 
  strain_t_prevalence$group = paste(strain_t_prevalence$group_I, "I", "*", strain_t_prevalence$group_P, "P")
  
  # rename group name for experiment part II
  if("strong" %in% strain_t_prevalence$group_I){
    for(i in 1:nrow(strain_t_prevalence)){
      strain_t_prevalence$group[i] = ifelse(strain_t_prevalence$group_P[i] == "strong", "negative", 
                                            ifelse(strain_t_prevalence$group_P[i] == "weak", "positive")) 
    }
  }
  
  
  # create line plot
  title <- sprintf("Dynamics of Plasmid Prevalence: Strain %s", strain_id)  # generate title
  
  p <- strain_t_prevalence %>%
    pivot_longer(colnames(strain_t_prevalence)[4]:colnames(strain_t_prevalence)[4+n_p-1], names_to = "Plasmid", values_to = "prevalence") %>%
    ggplot(aes(x = t, y = prevalence, color = Plasmid)) +
    geom_line() +
    facet_wrap(~ group, scales = "fixed") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_brewer(palette = "Dark2") +
    labs(title = title, x = "t", y = "Prevalence") + 
    paper_figs_theme +
    theme(strip.text = element_text(size = 14))
  
  return(p)
  
}

# Function for plasmid dynamics (of the community)
community_plasmid_dynamics <- function(list_of_results, t_output){
  # combine the list of results
  df <- bind_rows(list_of_results)
  
  # collect number of replicates
  n_rep = max(df$replicate)
  
  # unify time unit
  # df$t <- round(df$t)
  df$t <- floor(df$t/t_output) * t_output
  
  df$p_profile <- string_to_factor(df$p_profile)
  
  
  # group dependent variable by treatment groups, t, strain and p_profile, then calculate mean subpopulation abundance
  subpopulation_abundance <- df %>%
    group_by(group_I, group_P, t, strain, p_profile) %>%
    summarize(abundance_sub = sum(abundance)/n_rep)
  
  
  # add columns to report plasmid abundance of each plasmid strain
  column_names <- paste("P", 1:n_p, sep = "") # create a vector of names of new columns
  
  for (column_name in column_names){ 
    subpopulation_abundance[, column_name] <- rep(0, nrow(subpopulation_abundance))
  }# add the new columns into the dataframe
  
  
  for (i in 1:nrow(subpopulation_abundance)) {  # Go through rows 
    plasmid_profile <- as.numeric(strsplit(as.character(subpopulation_abundance$p_profile[i]), " ")[[1]]) # convert plasmid profile from factor (e.g. 01) to numeric vector c(0,1))
    p_abundance <- subpopulation_abundance$abundance_sub[i] * plasmid_profile
    subpopulation_abundance[i, 7:ncol(subpopulation_abundance)] <- as.list(p_abundance)
  } # fill in plasmid abundance
  
  
  # change column name of p_profile (for columns with name stared with P will be used later)
  colnames(subpopulation_abundance)[colnames(subpopulation_abundance) == "p_profile"] <- "fac_p_profile" 
  
  # group data frame by t and strain_id, then calculate plasmid prevalence
  t_prevalence <- subpopulation_abundance %>%
    group_by(group_I, group_P, t) %>%
    summarize(
      across(starts_with("P"), ~ sum(.) / sum(abundance_sub), .names = "{.col}")
    )
  
  # add a combined group (group I; group P) 
  t_prevalence$group = paste(t_prevalence$group_I, "I", "*", t_prevalence$group_P, "P")
  
  # rename group name for experiment part II
  if("strong" %in% t_prevalence$group_I){
    for(i in 1:nrow(t_prevalence)){
      t_prevalence$group[i] = ifelse(t_prevalence$group_P[i] == "strong", "negative", 
                                     ifelse(t_prevalence$group_P[i] == "weak", "positive")) 
    }
  }
  
  
  # create line plot
  p <- t_prevalence %>%
    pivot_longer(colnames(t_prevalence)[4]:colnames(t_prevalence)[4+n_p-1], names_to = "Plasmid", values_to = "prevalence") %>%
    ggplot(aes(x = t, y = prevalence, color = Plasmid)) +
    geom_line() +
    facet_wrap(~ group, scales = "fixed") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_brewer(palette = "Dark2") +
    labs(title = "Dynamics of Plasmid Prevalence (Community Level)", x = "t", y = "Prevalence") + 
    paper_figs_theme +
    theme(strip.text = element_text(size = 14)) 
  
  return(p)
  
}


# Function for plasmid dynamics (of the community; stacked)
community_plasmid_dynamics_stacked <- function(list_of_results, t_output){
  # combine the list of results
  df <- bind_rows(list_of_results)
  
  # collect number of replicates
  n_rep = max(df$replicate)
  
  # unify time unit
  # df$t <- round(df$t)
  df$t <- floor(df$t/t_output) * t_output
  
  df$p_profile <- string_to_factor(df$p_profile)
  
  
  # group dependent variable by treatment groups, t, strain and p_profile, then calculate mean subpopulation abundance
  subpopulation_abundance <- df %>%
    group_by(group_I, group_P, t, strain, p_profile) %>%
    summarize(abundance_sub = sum(abundance)/n_rep) 
  
  
  # add columns to report plasmid abundance of each plasmid strain
  column_names <- paste("P", 1:n_p, sep = "") # create a vector of names of new columns
  
  for (column_name in column_names){ 
    subpopulation_abundance[, column_name] <- rep(0, nrow(subpopulation_abundance))
  }# add the new columns into the dataframe
  
  
  for (i in 1:nrow(subpopulation_abundance)) {  # Go through rows 
    plasmid_profile <- as.numeric(strsplit(as.character(subpopulation_abundance$p_profile[i]), " ")[[1]]) # convert plasmid profile from factor (e.g. 01) to numeric vector c(0,1))
    p_abundance <- subpopulation_abundance$abundance_sub[i] * plasmid_profile
    subpopulation_abundance[i, 7:ncol(subpopulation_abundance)] <- as.list(p_abundance)
  } # fill in plasmid abundance
  
  
  # change column name of p_profile (for columns with name stared with P will be used later)
  colnames(subpopulation_abundance)[colnames(subpopulation_abundance) == "p_profile"] <- "fac_p_profile" 
  
  # group data frame by t and strain_id, then calculate plasmid prevalence
  t_prevalence <- subpopulation_abundance %>%
    group_by(group_I, group_P, t) %>%
    summarize(
      across(starts_with("P"), ~ sum(.) / sum(abundance_sub), .names = "{.col}")
    )
  
  # add a combined group (group I; group P) 
  t_prevalence$group = paste(t_prevalence$group_I, "I", "*", t_prevalence$group_P, "P")
  
  # rename group name for experiment part II
  if("strong" %in% t_prevalence$group_I){
    for(i in 1:nrow(t_prevalence)){
      t_prevalence$group[i] = ifelse(t_prevalence$group_P[i] == "strong", "negative", 
                                     ifelse(t_prevalence$group_P[i] == "weak", "positive")) 
    }
  }
  
  
  # create line plot
  p <- t_prevalence %>%
    pivot_longer(colnames(t_prevalence)[4]:colnames(t_prevalence)[4+n_p-1], names_to = "Plasmid", values_to = "prevalence") %>%
    ggplot(aes(x = t, y = prevalence, fill = Plasmid)) +
    geom_area() +
    facet_wrap(~ group, scales = "fixed") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_brewer(palette = "Dark2") +
    labs(title = "Dynamics of Plasmid Prevalence (Community Level)", x = "t", y = "Prevalence") + 
    paper_figs_theme +
    theme(strip.text = element_text(size = 14)) 
  
  return(p)
  
}




# # Function for demonstrating population-level bipartite network dynamics of a treatment group (stacked heat map): un-weighted
# uw_bp_network_dynamic <- function(list_of_results, trt_group, t_output){
#   # combine the list of results
#   df <- bind_rows(list_of_results)
#   
#   # collect number of replicates
#   n_rep = max(df$replicate)
#   
#   # unify time unit
#   # df$t <- round(df$t)
#   df$t <- floor(df$t/t_output) * t_output
#   
#   df$p_profile <- string_to_factor(df$p_profile)
#   
#   
#   # group dependent variable by treatment groups, t, strain and p_profile, then calculate mean subpopulation abundance
#   subpopulation_abundance <- df %>%
#     group_by(group_I, group_P, t, strain, p_profile) %>%
#     summarize(abundance = sum(abundance)/n_rep)
#   
#   # add a combined group (group I; group P) 
#   subpopulation_abundance$group = paste(subpopulation_abundance$group_I, "I", "*", subpopulation_abundance$group_P, "P")
#   
#   # rename group name for experiment part II
#   if("strong" %in% subpopulation_abundance$group_I){
#     for(i in 1:nrow(subpopulation_abundance)){
#       subpopulation_abundance$group[i] = ifelse(subpopulation_abundance$group_P[i] == "strong", "negative", 
#                                                 ifelse(subpopulation_abundance$group_P[i] == "weak", "positive")) 
#     }
#   }
#   
#   # extract from the treatment group of choice
#   subpopulation_abundance_sub = subset(subpopulation_abundance, group == trt_group)
#   
#   # create the interaction table
#   interaction <- interaction_table(subpopulation_abundance_sub, n_p)
#   
#   # table of p_profile levels
#   levels_p_profile = levels(subpopulation_abundance_sub$p_profile)
#   
#   # convert content of column p_profile from the level of p_profile to character p_profile
#   interaction$p_profile = levels_p_profile[interaction$p_profile]
#   
#   # create a vector of subpopulatoin abundance, replicate n_p times for each element
#   subp_abundance <- rep(subpopulation_abundance_sub$abundance, each = n_p)
#   
#   # fill in the vector to interaction table as a new column
#   interaction <- cbind(interaction, subp_abundance)
#   
#   # replace the copy with 0 for rows with subp_abundance = 0
#   for (i in 1:nrow(interaction)){
#     if(interaction$subp_abundance[i] == 0){
#       interaction$copy[i] = 0
#     }
#   }
#   
#   # define a vector of 6 layers (times), each with t_max/5 steps
#   times = seq(0, max(interaction$t), by = max(interaction$t)/5)
#   
#   # create layers of bipartite interaction matrix
#   matrix_list <- list()
#   for (i in 1:length(times)){
#     interaction_sub <- data_sub(times[i], interaction)
#     # check if all strain_id exist, if no (due to extinction), add the blank rows for the missing strain
#     for(j in 1:n_b){
#       if(!(j %in% interaction_sub$strain)){
#         for(k in 1:n_p){
#           new_row <- data.frame(t = times[i], subpop_id = 0, strain = j, plasmid_id = k, copy = 0)
#           interaction_sub <- rbind(interaction_sub, new_row)
#         }
#       }
#       else{interaction_sub <- interaction_sub}  
#     }
#     sum_interaction <- interaction_sub %>%
#       group_by(strain,plasmid_id) %>%
#       summarize(sum_link = sum(copy))
#     matrix_list[[i]] <- matrix(sum_interaction$sum_link, nrow = n_b, ncol = n_p, byrow = TRUE)
#     colnames(matrix_list[[i]]) <- paste0("P", 1:n_p)                             # Column names
#     rownames(matrix_list[[i]]) <- paste0("B", 1:n_b)
#   }
#   
#   # create serial heatmaps across time
#   p <- list() # create a list of plots
#   for (i in 1:length(matrix_list)){
#     data <- matrix_list[[i]]
#     p[[i]] <- plot_ly(x = colnames(data), y = rownames(data), z = data, type = "heatmap", coloraxis = "coloraxis") %>%
#       layout(font = list(size = 14))
#   }
#   
#   
#   subplot_titles = c() 
#   for (i in 1:length(times)){
#     subplot_titles[i] = paste0("t = ", times[i])
#   } # create titles of subplots
#   
#   title <- sprintf("<b>Bipartite Interaction Links: %s<b>", trt_group)  # generate title
#   
#   fig <- subplot(p, nrows = 2, margin = c(0.03,0.03,0.06,0.08)) %>% 
#     layout(title = list(text = title, font = list(size = 14, color = "#666666")),
#            coloraxis = list(colorscale = "Jet"),
#            margin = list(t = 50),
#            annotations = list(
#              list(x = 0.15 , y = 1.0, text = subplot_titles[1], xanchor = "center", yanchor = "bottom", showarrow = F, xref="paper", yref="paper"), 
#              list(x = 0.50 , y = 1.0, text = subplot_titles[2], xanchor = "center", yanchor = "bottom", showarrow = F, xref="paper", yref="paper"),
#              list(x = 0.85 , y = 1.0, text = subplot_titles[3], xanchor = "center", yanchor = "bottom", showarrow = F, xref="paper", yref="paper"),
#              list(x = 0.15 , y = 0.45, text = subplot_titles[4], xanchor = "center", yanchor = "bottom", showarrow = F, xref="paper", yref="paper"),
#              list(x = 0.50 , y = 0.45, text = subplot_titles[5], xanchor = "center", yanchor = "bottom", showarrow = F, xref="paper", yref="paper"),
#              list(x = 0.85 , y = 0.45, text = subplot_titles[6], xanchor = "center", yanchor = "bottom", showarrow = F, xref="paper", yref="paper")
#            )) # create the plot
#   
#   # note: subplots with different width...
#   
#   return(fig)
# }
# 
# 
# # Function for demonstrating population-level bipartite network dynamics of a treatment group (stacked heat map): weighted
# bp_network_dynamic <- function(list_of_results, trt_group, t_output){
#   # combine the list of results
#   df <- bind_rows(list_of_results)
#   
#   # collect number of replicates
#   n_rep = max(df$replicate)
#   
#   # unify time unit
#   # df$t <- round(df$t)
#   df$t <- floor(df$t/t_output) * t_output
#   
#   df$p_profile <- string_to_factor(df$p_profile)
#   
#   
#   # group dependent variable by treatment groups, t, strain and p_profile, then calculate mean subpopulation abundance
#   subpopulation_abundance <- df %>%
#     group_by(group_I, group_P, t, strain, p_profile) %>%
#     summarize(abundance = sum(abundance)/n_rep)
#   
#   # add a combined group (group I; group P) 
#   subpopulation_abundance$group = paste(subpopulation_abundance$group_I, "I", "*", subpopulation_abundance$group_P, "P")
#   
#   # rename group name for experiment part II
#   if("strong" %in% subpopulation_abundance$group_I){
#     for(i in 1:nrow(subpopulation_abundance)){
#       subpopulation_abundance$group[i] = ifelse(subpopulation_abundance$group_P[i] == "strong", "negative", 
#                                                 ifelse(subpopulation_abundance$group_P[i] == "weak", "positive")) 
#     }
#   }
#   
#   # extract from the treatment group of choice
#   subpopulation_abundance_sub = subset(subpopulation_abundance, group == trt_group)
#   
#   # create the interaction table
#   interaction <- interaction_table(subpopulation_abundance_sub, n_p)
#   
#   # table of p_profile levels
#   levels_p_profile = levels(subpopulation_abundance_sub$p_profile)
#   
#   # convert content of column p_profile from the level of p_profile to character p_profile
#   interaction$p_profile = levels_p_profile[interaction$p_profile]
#   
#   # create a vector of subpopulatoin abundance, replicate n_p times for each element
#   subp_abundance <- rep(subpopulation_abundance_sub$abundance, each = n_p)
#   
#   # fill in the vector to interaction table as a new column
#   interaction <- cbind(interaction, subp_abundance)
#   
#   # replace the copy with 0 for rows with subp_abundance = 0
#   for (i in 1:nrow(interaction)){
#     if(interaction$subp_abundance[i] == 0){
#       interaction$copy[i] = 0
#     }
#   }
#   
#   # define a vector of 6 layers (times), each with t_max/5 steps
#   times = seq(0, max(interaction$t), by = max(interaction$t)/5)
#   
#   # create layers of bipartite interaction matrix
#   matrix_list <- list()
#   for (i in 1:length(times)){
#     interaction_sub <- data_sub(times[i], interaction)
#     # check if all strain_id exist, if no (due to extinction), add the blank rows for the missing strain
#     for(j in 1:n_b){
#       if(!(j %in% interaction_sub$strain)){
#         for(k in 1:n_p){
#           new_row <- data.frame(t = times[i], subpop_id = 0, strain = j, plasmid_id = k, copy = 0, subp_abundance = 0)
#           interaction <- rbind(interaction_sub, new_row)
#         }
#       }
#       else{interaction_sub <- interaction_sub}  
#     }
#     sum_interaction <- interaction_sub %>%
#       group_by(strain,plasmid_id) %>%
#       summarize(sum_link = sum(copy * subp_abundance))
#     matrix_list[[i]] <- matrix(sum_interaction$sum_link, nrow = n_b, ncol = n_p, byrow = TRUE)
#     colnames(matrix_list[[i]]) <- paste0("P", 1:n_p)                             # Column names
#     rownames(matrix_list[[i]]) <- paste0("B", 1:n_b)
#   }
#   
#   # create serial heatmaps across time
#   p <- list() # create a list of plots
#   for (i in 1:length(matrix_list)){
#     data <- matrix_list[[i]]
#     p[[i]] <- plot_ly(x = colnames(data), y = rownames(data), z = data, type = "heatmap", coloraxis = "coloraxis") %>%
#       layout(font = list(size = 14))
#   }
#   
#   
#   subplot_titles = c() 
#   for (i in 1:length(times)){
#     subplot_titles[i] = paste0("t = ", times[i])
#   } # create titles of subplots
#   
#   title <- sprintf("<b>Bipartite Interaction Links (Weighted): %s<b>", trt_group)  # generate title
#   
#   fig <- subplot(p, nrows = 2, margin = c(0.03,0.03,0.06,0.08)) %>% 
#     layout(title = list(text = title, font = list(size = 14, color = "#666666")),
#            coloraxis = list(colorscale = "Jet"),
#            margin = list(t = 50),
#            annotations = list(
#              list(x = 0.15 , y = 1.0, text = subplot_titles[1], xanchor = "center", yanchor = "bottom", showarrow = F, xref="paper", yref="paper"), 
#              list(x = 0.50 , y = 1.0, text = subplot_titles[2], xanchor = "center", yanchor = "bottom", showarrow = F, xref="paper", yref="paper"),
#              list(x = 0.85 , y = 1.0, text = subplot_titles[3], xanchor = "center", yanchor = "bottom", showarrow = F, xref="paper", yref="paper"),
#              list(x = 0.15 , y = 0.45, text = subplot_titles[4], xanchor = "center", yanchor = "bottom", showarrow = F, xref="paper", yref="paper"),
#              list(x = 0.50 , y = 0.45, text = subplot_titles[5], xanchor = "center", yanchor = "bottom", showarrow = F, xref="paper", yref="paper"),
#              list(x = 0.85 , y = 0.45, text = subplot_titles[6], xanchor = "center", yanchor = "bottom", showarrow = F, xref="paper", yref="paper")
#            )) # create the plot
#   
#   # note: subplots with different width...
#   
#   return(fig)
# }
# 
# # Function for population-level network structure across time of a treatment group (animated heat map): un-weighted
# uw_bp_network_dynamic_ani = function(list_of_results, trt_group, t_final, t_output){
#   
#   df <- bind_rows(list_of_results)
#   
#   # collect number of replicates
#   n_rep = max(df$replicate)
#   
#   # unify time unit
#   # df$t <- round(df$t)
#   df$t <- floor(df$t/t_output) * t_output
#   
#   df$p_profile <- string_to_factor(df$p_profile)
#   
#   
#   # group dependent variable by treatment groups, t, strain and p_profile, then calculate mean subpopulation abundance
#   subpopulation_abundance <- df %>%
#     group_by(group_I, group_P, t, strain, p_profile) %>%
#     summarize(abundance = sum(abundance)/n_rep)
#   
#   # add a combined group (group I; group P) 
#   subpopulation_abundance$group = paste(subpopulation_abundance$group_I, "I", "*", subpopulation_abundance$group_P, "P")
#   
#   # rename group name for experiment part II
#   if("strong" %in% subpopulation_abundance$group_I){
#     for(i in 1:nrow(subpopulation_abundance)){
#       subpopulation_abundance$group[i] = ifelse(subpopulation_abundance$group_P[i] == "strong", "negative", 
#                                                 ifelse(subpopulation_abundance$group_P[i] == "weak", "positive")) 
#     }
#   }
#   
#   # extract from the treatment group of choice
#   subpopulation_abundance_sub = subset(subpopulation_abundance, group == trt_group)
#   
#   # create the interaction table, and fix value range for the color bar
#   interaction <- interaction_table(subpopulation_abundance_sub, n_p)
#   
#   # table of p_profile levels
#   levels_p_profile = levels(subpopulation_abundance_sub$p_profile)
#   
#   # convert content of column p_profile from the level of p_profile to character p_profile
#   interaction$p_profile = levels_p_profile[interaction$p_profile]
#   
#   # create a vector of subpopulatoin abundance, replicate n_p times for each element
#   subp_abundance <- rep(subpopulation_abundance_sub$abundance, each = n_p)
#   
#   # fill in the vector to interaction table as a new column
#   interaction <- cbind(interaction, subp_abundance)
#   
#   # replace the copy with 0 for rows with subp_abundance = 0
#   for (i in 1:nrow(interaction)){
#     if(interaction$subp_abundance[i] == 0){
#       interaction$copy[i] = 0
#     }
#   }
#   
#   # group data by t, strain, plsmid_id, and sum up plasmid copies
#   interaction_sum <- interaction %>%
#     group_by(t, strain, plasmid_id) %>%
#     summarize(copy_sum = sum(copy))
#   
#   # for each t, check if all strain_id exist, if no (due to extinction), add the blank rows for the missing strain
#   times = unique(interaction_sum$t)
#   for (i in times){
#     interaction_sum_sub = subset(interaction_sum, t == i)
#     for(j in 1:n_b){
#       if(!(j %in% interaction_sum_sub$strain)){
#         for(k in 1:n_p){
#           new_row <- data.frame(t = i, strain = j, plasmid_id = k, copy_sum = 0)
#           interaction_sum <- rbind(interaction_sum, new_row)
#         }
#       }
#       else{interaction_sum <- interaction_sum}  
#     }
#   }
#   
#   
#   # add characters B & P before strain & plasmid ids
#   interaction_sum$strain <- paste0("B", interaction_sum$strain)
#   interaction_sum$plasmid_id <- paste0("P", interaction_sum$plasmid_id)
#   
#   # convert dada frame into form tibble
#   interaction_sum = as_tibble(interaction_sum)
#   
#   # set title for the animated heatmap
#   title <- sprintf("Bipartite Interaction: %s", trt_group)  # generate title
#   
#   # set file name for saving
#   name <- sprintf("uw_bp_network_%s.gif", trt_group)
#   
#   # set value range (note: adjust the maximum value based on the largest link of all experimental treatment to unify the color bar across treatments)
#   value_range <- 0:8
#   # color_palette <- colorRampPalette(c("navyblue", "white"))(length(value_range))
# 
#   # breaks = as.character(value_range), limits = as.character(value_range)
#   
#   # set variables that will later be used in transition_states()
#   s_length = t_final/t_output + 1
#   
#   # create animated heatmap across time (using package ggplot2)
#   p <- ggplot(interaction_sum, aes(plasmid_id, strain, fill= copy_sum)) + 
#     geom_tile() +
#     #scale_fill_manual(drop=FALSE, values=color_palette, limits = rev(factor(value_range)), breaks =  rev(factor(value_range)), labels = rev(value_range), name="Links") +
#     scale_fill_gradient2(low="white", high="navyblue", #colors in the scale
#                          breaks=seq(min(value_range), max(value_range), 1), #breaks in the scale bar
#                          limits=c(min(value_range), max(value_range)), guide = "legend", name = "Links") +
#     coord_fixed() +
#     paper_figs_theme +
#     theme(text = element_text(size = 20),
#           plot.title = element_text(size = 16),
#           plot.subtitle=element_text(size=14),
#           panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(), 
#           panel.background = element_blank(),
#           plot.margin = unit(c(0.35,0.35,0.35,0.35), "cm")) +
#           #guides(size = guide_legend(reverse = FALSE)) +
#     # Here comes the gganimate specific bits
#     transition_states(as.integer(t), state_length = s_length, transition_length = c(rep(5, s_length)), wrap = TRUE) +
#     labs(title = title, subtitle ='t: {next_state}', x = 'Plasmids', y = 'Bacteria') 
# 
#   
#   # for the R pipeline
#   animate(p, fps=5, height = 380, width = 380)
#   # animate(p, renderer = ffmpeg_renderer()) # show in mpeg format
#   # anim_save(name, animation = animate(p, fps=5, height = 380, width = 380)) # save plot into .gif file
#   
#   }
# 
# 
# # Function for population-level network structure across time of a treatment group (animated heat map): weighted
# w_bp_network_dynamic_ani = function(list_of_results, trt_group, t_final, t_output){
#   df <- bind_rows(list_of_results)
#   
#   # collect number of replicates
#   n_rep = max(df$replicate)
#   
#   # unify time unit
#   # df$t <- round(df$t)
#   df$t <- floor(df$t/t_output) * t_output
#   
#   df$p_profile <- string_to_factor(df$p_profile)
#   
#   
#   # group dependent variable by treatment groups, t, strain and p_profile, then calculate mean subpopulation abundance
#   subpopulation_abundance <- df %>%
#     group_by(group_I, group_P, t, strain, p_profile) %>%
#     summarize(abundance = sum(abundance)/n_rep)
#   
#   # add a combined group (group I; group P) 
#   subpopulation_abundance$group = paste(subpopulation_abundance$group_I, "I", "*", subpopulation_abundance$group_P, "P")
#   
#   # rename group name for experiment part II
#   if("strong" %in% subpopulation_abundance$group_I){
#     for(i in 1:nrow(subpopulation_abundance)){
#       subpopulation_abundance$group[i] = ifelse(subpopulation_abundance$group_P[i] == "strong", "negative", 
#                                                 ifelse(subpopulation_abundance$group_P[i] == "weak", "positive")) 
#     }
#   }
#   
#   # extract from the treatment group of choice
#   subpopulation_abundance_sub = subset(subpopulation_abundance, group == trt_group)
#   
#   # create the interaction table
#   interaction <- interaction_table(subpopulation_abundance_sub, n_p)
#   
#   # table of p_profile levels
#   levels_p_profile = levels(subpopulation_abundance_sub$p_profile)
#   
#   # convert content of column p_profile from the level of p_profile to character p_profile
#   interaction$p_profile = levels_p_profile[interaction$p_profile]
#   
#   # create a vector of subpopulatoin abundance, replicate n_p times for each element
#   subp_abundance <- rep(subpopulation_abundance_sub$abundance, each = n_p)
#   
#   # fill in the vector to interaction table as a new column
#   interaction <- cbind(interaction, subp_abundance)
#   
#   # replace the copy with 0 for rows with subp_abundance = 0
#   for (i in 1:nrow(interaction)){
#     if(interaction$subp_abundance[i] == 0){
#       interaction$copy[i] = 0
#     }
#   }
#   
#   # group data by t, strain, plsmid_id, and sum up plasmid copies
#   interaction_sum <- interaction %>%
#     group_by(t, strain, plasmid_id) %>%
#     summarize(copy_sum = sum(copy*subp_abundance))
#   
#   # for each t, check if all strain_id exist, if no (due to extinction), add the blank rows for the missing strain
#   times = unique(interaction_sum$t)
#   for (i in times){
#     interaction_sum_sub = subset(interaction_sum, t == i)
#     for(j in 1:n_b){
#       if(!(j %in% interaction_sum_sub$strain)){
#         for(k in 1:n_p){
#           new_row <- data.frame(t = i, strain = j, plasmid_id = k, copy_sum = 0)
#           interaction_sum <- rbind(interaction_sum, new_row)
#         }
#       }
#       else{interaction_sum <- interaction_sum}  
#     }
#   }
#   
#   
#   # add characters B & P before strain & plasmid ids
#   interaction_sum$strain <- paste0("B", interaction_sum$strain)
#   interaction_sum$plasmid_id <- paste0("P", interaction_sum$plasmid_id)
#   
#   # convert dada frame into form tibble
#   interaction_sum = as_tibble(interaction_sum)
#   
#   # set title for the animated heatmap
#   title <- sprintf("Bipartite Interaction (w): %s", trt_group)  # generate title
#   
#   # set file name for saving
#   name <- sprintf("uw_bp_network_%s.gif", trt_group)
#   
#   # set variables that will later be used in transition_states()
#   s_length = t_final/t_output + 1
#   
#   # create animated heatmap across time (using package ggplot2)
#   # replace the 0 values in copy_sum by NA, for a distinctive color filling
#   # interaction_sum$copy_sum[interaction_sum$copy_sum == 0] <- NA
#   
#   p <- ggplot(interaction_sum, aes(plasmid_id, strain, fill = ifelse(copy_sum != 0, copy_sum, NA))) + 
#     geom_tile() +
#     scale_fill_gradientn(colours=c("navy", "skyblue", "greenyellow", "yellow", "#fd8d3c", "#bd0026"),
#                          limits = c(0, 2e4),  # Define the range of values, maximum based on carrying capacity
#                          labels = label_number(scale = 1e-3), na.value = "white", name="Links (K)") +
#     coord_fixed() +
#     paper_figs_theme +
#     theme(text = element_text(size = 18),
#           plot.title = element_text(size = 16),
#           plot.subtitle=element_text(size=14),
#           panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(), 
#           panel.background = element_blank(),
#           plot.margin = unit(c(0,0,0,0), "cm")) +
#     guides(fill = guide_colourbar(theme = theme(legend.key.height = unit(15, "lines")))) +
#     # Here comes the gganimate specific bits
#     transition_states(as.integer(t), state_length = s_length, transition_length = c(rep(5, s_length)), wrap = TRUE) +
#     labs(title = title, subtitle = "t: {next_state}", x = 'Plasmids', y = 'Bacteria') 
#   
#   #animated_p <- animate(p, fps=5)
#   # dev.new()
#  animate(p, fps=5, height = 380, width = 380)
#   # anim_save(name, animation = animated_p) # save plot into .gif file; optional
# }


### below are functions designed for 3x4 system's analysis Rmd

## Function that create interaction table, with result as a subset of data frame of a given t and replicate ##
interaction_table2 <- function(result, K){
  # construct vectors for the table of interaction
  vector_exp <- c()
  vector_rep <- c()
  vector_t <- c()
  vector_p_profile <- factor()
  vector_strain <- c()
  vector_plasmid_id <- c()
  vector_copy <- c()
  
  # fill in the table of interaction
  for (i in 1:nrow(result)){
    plasmid_profile <- as.numeric(strsplit(as.character(result$p_profile[i]), " ")[[1]]) # convert plasmid profile from factor (e.g. 01) to numeric vector c(0,1))
    for (j in 1:n_p){
      # create new row
      vector_exp <- c(vector_exp, result$experiment[1])
      vector_rep <- c(vector_rep, result$replicate[1])
      vector_t <- c(vector_t, result$t[i])
      vector_p_profile <- c(vector_p_profile, result$p_profile[i])
      vector_strain <- c(vector_strain, result$strain[i])
      vector_plasmid_id <- c(vector_plasmid_id, j)
      vector_copy <- c(vector_copy, as.numeric(plasmid_profile[j]*result$abundance[i]/K))
    }
  }
  
  # construct the table of interaction
  interaction <- data.frame(experiment = vector_exp, replicate = vector_rep, t = vector_t, strain = vector_strain, p_profile = vector_p_profile, plasmid_id = vector_plasmid_id, copy = vector_copy)
  return(interaction)
}


# Function to convert bipartite matrix to long format data frame (for plottting using ggplot2)
matrix_to_long_df <- function(mat, name) {
  tb <- melt(mat)
  colnames(tb) <- c("Row", "Col", "Value")
  tb$Matrix <- name
  return(tb)
}


# Function for demonstrating population-level bipartite network dynamics of a treatment group (stacked heat map): weighted; mean & variance

w_bp_network_dynamic2 <- function(df, exp, n_layers, K){
  # determine time layers to extract
  times = seq(0, max(df$t), by = max(df$t)/(n_layers-1))
  
  # extract from the treatment group of choice
  df_sub <- subset(df, experiment == exp & t %in% times)
  
  # create an empty interaction table
  interaction <- data.frame(experiment = numeric(0), replicate = numeric(0), t = numeric(0), strain = numeric(0), p_profile = c(), plasmid_id = numeric(0), copy = numeric(0))
  
  # for each replicate at each time, create an interaction table, then combined the tables
  replicates <- unique(df_sub$replicate)
  
  for (i in times){
    for (j in replicates){
      df_sub2 <- subset(df_sub, t == i & replicate == j)
      interaction_table_sub <- interaction_table2(df_sub2, K)
      interaction <- rbind(interaction, interaction_table_sub)
    }
  }
  
  
  # group dependent variable by treatment groups, t, strain and p_profile, then calculate mean & sd of plasmid number in each strain across replicates 
  sum_interaction <- interaction %>%
    group_by(experiment, replicate, t, strain, plasmid_id) %>%
    summarize(copy_sum = sum(copy)) %>%
    ungroup() %>%
    group_by(experiment, t, strain, plasmid_id) %>%
    summarize(n = n(), copy_mean = mean(copy_sum), copy_sd = sd(copy_sum))
  
  # create layers of bipartite interaction matrix
  matrix_list_mean <- list()
  matrix_list_sd <- list()
  
  for (i in 1:length(times)){
    sum_interaction_sub <- data_sub(times[i], sum_interaction)
    # check if all strain_id exist, if no (due to extinction), add the blank rows for the missing strain
    for(j in 1:n_b){
      if(!(j %in% sum_interaction_sub$strain)){
        for(k in 1:n_p){
          new_row <- data.frame(experiment = sum_interaction_sub$experiment[1], t = times[i], strain = j, plasmid_id = k, n = length(replicates), copy_mean = 0, copy_sd = 0)
          sum_interaction_sub <- rbind(sum_interaction_sub, new_row)
        }
      }
      else{sum_interaction_sub <- sum_interaction_sub}  
    }
    
    matrix_list_mean[[i]] <- matrix(sum_interaction_sub$copy_mean, nrow = n_b, ncol = n_p, byrow = TRUE)
    matrix_list_sd[[i]] <- matrix(sum_interaction_sub$copy_sd, nrow = n_b, ncol = n_p, byrow = TRUE)
    
    colnames(matrix_list_mean[[i]]) <- paste0("P", 1:n_p)                             # Column names
    rownames(matrix_list_mean[[i]]) <- paste0("B", 1:n_b)
    colnames(matrix_list_sd[[i]]) <- paste0("P", 1:n_p)                             # Column names
    rownames(matrix_list_sd[[i]]) <- paste0("B", 1:n_b)
  }
  
  # name each matrix according to t
  matrix_names <- paste("t =", times)
  names(matrix_list_mean) <- matrix_names
  names(matrix_list_sd) <- matrix_names
  
  # convert list of matrices into list of long data frames
  long_df_list_mean <- lapply(names(matrix_list_mean), function(name) {
    matrix_to_long_df(matrix_list_mean[[name]], name)
  })
  
  long_df_list_sd <- lapply(names(matrix_list_sd), function(name) {
    matrix_to_long_df(matrix_list_sd[[name]], name)
  })
  
  # combine data frames
  combined_long_dfmean <- bind_rows(long_df_list_mean)
  combined_long_dfsd <- bind_rows(long_df_list_sd)
  
  # make sure matrix name is ordered
  combined_long_dfmean$Matrix <- factor(combined_long_dfmean$Matrix, levels = matrix_names)
  combined_long_dfsd$Matrix <- factor(combined_long_dfsd$Matrix, levels = matrix_names)
  
  # convert value 0 into NA (to make distinct colors)
  combined_long_dfmean$Value[combined_long_dfmean$Value == 0] <- NA
  combined_long_dfsd$Value[combined_long_dfsd$Value == 0] <- NA
  
  list_of_tables <- list(Mean = combined_long_dfmean, SD = combined_long_dfsd)
  
  return(list_of_tables)
}

# Function that produce a network plot from a single pair of nodes
plotweb_singlepair <- function(edge_list) {
  
  plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 10), ylim = c(0, 10), axes = FALSE)
  
  # Draw the higher rectangular box
  rect(xleft = 4, ybottom = 7.7, xright = 6, ytop = 9, col = colsP[as.character(edge_list$from[1])], border = "black")
  text(x = 5, y = 9.2, labels = edge_list$from[1], col = "black", cex = 0.7)
  
  # Draw the lower rectangular box
  rect(xleft = 4, ybottom = 1, xright = 6, ytop = 2.5, col = colsB[as.character(edge_list$strain[1])], border = "black")
  text(x = 5, y = 0.7, labels = edge_list$to[1], col = "black", cex = 0.7)
  
  # Draw the middle rectangular box
  rect(xleft = 4, ybottom = 2.2, xright = 6, ytop = 7.8, col = "grey", border = "black")
}


# Function for converting data subset (of a given t and replicate) into edge list
network_vis <- function(df, exp, time, rep){
  
  # construct table of the edge list
  edge_list <- tibble(from = character(),
                      to = character(),
                      weight = integer(),
                      strain = character())
  
  # extract from the treatment group of choice
  df_sub <- subset(df, experiment == exp & t == time & replicate == rep)
  
  for (i in 1:length(df_sub$abundance)){
    plasmid_profile <- as.numeric(strsplit(as.character(df_sub$p_profile[i]), " ")[[1]]) # convert plasmid profile from factor (e.g. 01) to numeric vector c(0,1))
    for (j in 1:n_p){
      if (plasmid_profile[j] != 0 && df_sub$abundance[i] > 0){
        new_row <- tibble(from = paste0("P", j), to = as.character(df_sub$subpop_id[i]), weight = df_sub$abundance[i], strain = paste0("B", df_sub$strain[i])) 
        edge_list <- rbind(edge_list, new_row)  
      }
    }
  }
  
  # create monolayer object using package infomapecology
  monolayer_object <- create_monolayer_object(edge_list,
                                              bipartite = T, 
                                              directed = T, 
                                              group_names = c('Pasmid', 'Host'))
  
  if(nrow(edge_list) == 1) {
    # message <- paste0("The network has only one subpouplatoin left: ", edge_list$strain[1], " infected with ", edge_list$from[1]," with an abundance of ", edge_list$weight[1], ".")
    # print(message)
    plotweb_singlepair(edge_list)
  } else {
    
    # plot the adjacency matrix (using package bipartite)
    # visweb(monolayer_object$mat) 
    
    # colsB <- rainbow(length(unique(E(monolayer_object$igraph_object)$strain))
    
    # Define colors for each strain & plasmid
    colsB <- c("B1" = "#F8766D", "B2" = "#00BA38", "B3" = "#619CFF")
    colsP <- c("P1" = "#1b9e77", "P2" = "#d95f02", "P3" = "#7570b3", "P4" = "#e7298a")
    
    
    # create a color vector for the higher nodes
    higher_node_names <- colnames(monolayer_object$mat)  # Get the names of higher nodes
    col_high <- colsP[higher_node_names]
    
    # plot the network from adjacency matrix
    plotweb(monolayer_object$mat, col.low = colsB[as.character(E(monolayer_object$igraph_object)$strain)], col.high = col_high) 
    
  }
  return(monolayer_object)
}


#result <- network_vis(df_LB, 7, 2000, 1)

#names(monolayer_object)
#monolayer_object$edge_list # edge list
#monolayer_object$mat # adjacency matrix
#monolayer_object$mode # unipartite(U)/bipartite(B)
#monolayer_object$directed # directed?
#monolayer_object$nodes 
#monolayer_object$igraph_object
#vertex_attr(monolayer_object$igraph_object) # check the group of nodes
#plot(monolayer_object$igraph_object, layout=layout_as_bipartite) # plot the network from igraph

# Function to calculate Shannon diversity index
shannon_diversity <- function(counts) {
  proportions <- counts / sum(counts)
  -sum(proportions * log(proportions))
}

# Function to visualize perturbation impact
pulse_expd_periodic2 <- function(lag, time, perturb, t_p, k){
  if (time < lag) {impact <- 0} else {
    impact <- k * exp(-(1/perturb) * (time %% t_p))
  }
  return(impact)
}

# e.g.
# t = 0:2000
# impact9 = c()
# impact100 = c()
# impact1000 = c()
#for (i in t){
#  element9 <- pulse_expd_periodic2(100, i, 9, 100, 9.0)
#  element100 <- pulse_expd_periodic2(100, i, 100, 100, 9.0)
#  element1000 <- pulse_expd_periodic2(100, i, 1000, 100, 9.0)
#  impact9 <- c(impact9, element9)
#  impact100 <- c(impact100, element100)
#  impact1000 <- c(impact1000, element1000)
#}
# plot(t, impact9, type = "l", col = "blue", lwd=2, ylab = expression(epsilon[i]))
# lines(t, impact100, col="green", lwd=2)
# lines(t, impact1000, col="red", lwd=2)
# legend("topright", legend=c(expression(k[i] == 9), expression(k[i] == 100), expression(k[i] == 1000)), col=c("blue", "green", "red"), lwd=2, lty=1, bg="white")

# Function to add transparency of a color
add_transparency <- function(color, alpha = 0.5) {
  rgb_vals <- col2rgb(color) / 255
  rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha = alpha)
}