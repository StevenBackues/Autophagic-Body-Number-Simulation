require(tidyverse)
require(dgof)

# Gnarly split violin functions - https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
# https://stackoverflow.com/a/45614547
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
# ************************************************************************************************************
raw_exp_data <- read_csv("Experimental data.csv") %>% drop_na()
char_time <- as.character(Sys.time())
# ************************************************************************************************************

# ensures file structure is correct for saving output
check_folders <- function(char_time) {
  
  stripped_time <- gsub(" ", "_", char_time, fixed = TRUE)
  
  if (!file.exists("resim_output")) {
    dir.create(file.path(getwd(), "resim_output"))
  }
  
  path <- file.path(getwd(), "resim_output")
  if (!file.exists(file.path(path, stripped_time))) {
    dir.create(file.path(path, stripped_time))
  }
  
  path <- file.path(path, stripped_time)
  if (!file.exists(file.path(path, "plots"))) {
    dir.create(file.path(path, "plots"))
  }
  
  path <- file.path(path, "plots")
  if(!file.exists(file.path(path, "QQ_plots"))) {
    dir.create(file.path(path, "QQ_plots"))
  }
}

# general analysis of which body mean fits best
run_ks_tests <- function(df, char_time) {
  # to prevent downstream errors
  check_folders(char_time)
  
  #initiate ks result vectors
  ks_list <- list(length(unique(df$input_body_mean)))
  
  #perform ks test using each input body number mean
  for (i in 1:length(unique(df$input_body_mean))) {
    active_mean <- as.integer((i + (min(unique(df$input_body_mean)) - 1)))
    active_sim_data <- subset(df, input_body_mean == active_mean)
    ks_result <- ks.test(raw_exp_data$exp, ecdf(active_sim_data$sliced_bodies))
    ks_list[[i]] <- data.frame(mean = active_mean, D = ks_result[[1]])
  }
  
  ks_df <- do.call("bind_rows", ks_list)
  
  # plot D stat vs mean
  ggplot(ks_df, aes(x = as.factor(mean), y = D)) +
    geom_point() +
    geom_text(aes(y = D + 0.002, label = round(D, digits = 3)), size = 2, vjust = 0) +
    theme_bw() +
    xlab("Simulation Input Mean") +
    ylab("K-S Test D Statistic") +
    labs(title = "K-S Test Results") +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(colour = "grey90", linetype = "dashed"))
  
  # save in sensible location
  ggsave("KS_Plot.png",
         path = paste("resim_output/", gsub(" ", "_", char_time, fixed = TRUE), "/plots", sep = ""),
         width = 8, 
         height = 8, 
         unit = "cm", 
         dpi = 300)
  
  # return a vector containing the means with lowest D values in ascending order
  ks_df2 <- ks_df %>%
    arrange(D)
  return(ks_df2$mean)
}

plot_split_violins <- function(df, ks_means, char_time) {
  # check folder structure
  check_folders(char_time)
  
  # initialize a dataframe for later
  by_mean <- data_frame()
  
  # arrange the experimental data for split violin plots
  exp_data <- raw_exp_data %>%
    select(exp) %>%
    gather(set, num_body) %>%
    mutate(method = "Experimental")
  
  # find out how many different means are present, use up to 5
  if (length(ks_means) < 5) {
    num_means <- length(ks_means)
  } else {
    num_means <- 5
  }
  
  # split resim data into list of dataframes, show at most 5 of the lowest KS valued
  filtered_df <- df %>% 
    filter(input_body_mean %in% ks_means[1:num_means])
  by_mean_list <- split(filtered_df, filtered_df$input_body_mean)
  
  # loop through at most 5 lowest KS valued means to set up data frames for split violin plotting
  for (i in 1:num_means) {
    X <- by_mean_list[[i]] %>%
      select(input_body_mean, sliced_bodies) %>%
      mutate(method = "Simulated")
    names(X) <- c("set", "num_body", "method")
    exp_data$set <- X$set[i]
    X <- bind_rows(X, exp_data)
    by_mean <- bind_rows(by_mean, X)
  }
  
  # generate plots
  ggplot(by_mean, aes(x = factor(set), y = num_body, fill = method))  + 
    geom_split_violin(adjust = .7, size = 0.25) +
    scale_fill_brewer(palette="Accent", labels = c("Experimental\n(Mean)", "Simulated\n(Mean)")) +
    stat_summary(fun.y="mean", geom="point", shape=95, size=4, position=position_dodge(0.3)) +
    theme_bw() +
    labs(title = "Split Violin Plots For Lowest K-S Means", fill = "Type") +
    xlab("Input Mean") +
    ylab("Bodies/Vacuole") +
    guides(fill=guide_legend(title = NULL)) +
    theme(panel.grid.major.x = element_blank(),
          axis.text = element_text(size = 4),
          axis.title = element_text(size = 6),
          plot.title = element_text(size = 6),
          legend.text = element_text(size = 4))
  
  # save in sensible location
  ggsave("Split_Violins.png",
         path = paste("resim_output/", gsub(" ", "_", char_time, fixed = TRUE), "/plots", sep = ""),
         width = 12, 
         height = 6, 
         unit = "cm", 
         dpi = 300)
}

# iterable function for plot_QQs
plot_QQ <- function(df) {
  
  # plots a single QQ plot for a given simulation input mean, saves it as png
  qqplot(x = raw_exp_data$exp,
         y = df$sliced_bodies,
         main = paste("Simulated Mean:", df$input_body_mean[1]),
         xlab = "Experimental",
         ylab = "Simulated")
  abline(0, 1)
}

plot_QQs <- function(df, char_time) {
  # just to be safe / for testing
  check_folders(char_time)
  
  # split up the data frame for serarate analysis
  unique_means <- split(df, df$input_body_mean)
  
  # open the png writing function
  png(paste("resim_output/", gsub(" ", "_", char_time), "/plots/QQ_plots/QQ_plot_%d.png", sep= ""),
      width = 400,
      height = 400)
  
  # generate a QQ plot for each unique mean
  map(unique_means, plot_QQ)
  
  dev.off()
}

plot_mean_differences <- function(df, char_time) {
  # just to be safe / for testing
  check_folders(char_time)
  
  # prepare data for plotting
  exp_mean <- mean(raw_exp_data$exp, na.rm = TRUE)
  
  means_df <- df %>%
    group_by(input_body_mean) %>%
    summarize(mean = mean(sliced_bodies, na.rm = TRUE)) %>%
    mutate(abs_diff = abs(exp_mean - mean))
  
  # plot, plot, plot it up
  ggplot(means_df, aes(x = factor(input_body_mean), y = abs_diff)) +
    geom_point()
}

plot_means <- function(df, char_time) {
  
  exp_mean <- mean(raw_exp_data$exp, na.rm = TRUE)
  
  # just to be safe / for testing
  check_folders(char_time)
  
  #prepare data for plotting
  means_df <- df %>%
    group_by(input_body_mean) %>%
    summarize(mean = mean(sliced_bodies, na.rm = TRUE))
  
  # compare experimental and simulated observed means
  ggplot(means_df, aes(x = factor(input_body_mean), y = mean)) +
    geom_hline(yintercept = exp_mean, colour = "grey55") +
    geom_point() +
    annotate("text", x = 1.8, y = exp_mean + .05, label = "Experimental", size = 2) +
    theme_bw() +
    labs(title = "Simulated Means") +
    xlab("Input Mean") +
    ylab("Mean Bodies / Slice")

  ggsave("Means_Plot.png",
         path = paste("resim_output/", gsub(" ", "_", char_time, fixed = TRUE), "/plots", sep = ""),
         width = 8, 
         height = 8, 
         unit = "cm", 
         dpi = 300)
}

# replacement for the large write.table calls in older script
write_data <- function(df, char_time) {
  # just to be safe / for testing
  check_folders(char_time)
  
  # to make the write function more readable
  df_name <- deparse(substitute(df))
  file_name <- paste(df_name, ".csv", sep = "")
  file_path <- file.path("resim_output", gsub(" ", "_", char_time, fixed = TRUE), file_name)
  
  write_csv(df, path = file_path)
}
