library(tidyr)
library(mirt)
source("R/theta2arclengrm.R")

grm_4_icc_plots <- function(grm, item, titlestr, theta_interval, Qvec, Qvec_al, labsize = 12,
                            xlab = "x-axis", plotTitle = TRUE, lgdpos = "bottom") {
  surp_base <- length(unique(grm@Data$data[, item]))
  x_lab <- paste0("Test Information (", length(unique(grm@Data$data[, item])), "-bits)")
  p1 <- group_fit_plot(
    model = grm, 
    theta_interval = theta_interval, 
    item_to_plot = item, 
    grouped_probs = NULL, 
    surprisal = TRUE,
    surp_base = surp_base,
    arc_length = FALSE,
    x_lab = TeX("$\\theta$"),
    color_plot = TRUE
  ) + 
    geom_vline(xintercept = Qvec, color="black", linetype = "dashed") +
    geom_hline(yintercept = 0,  color="black", linetype = "dashed") +
    ggtitle(paste0(titlestr, item)) +
    ylim(c(0, 4)) +
    theme_bw() +
    theme(axis.title=element_text(size=labsize,face="bold"),
          legend.title = element_blank(),
          legend.position = lgdpos,
          axis.text  = element_text(size = 10),
          legend.text = element_text(size = 10))
  if (xlab == "") p1 <- p1 + theme(axis.title.x = element_blank())
  if (plotTitle) {
    p1 <- p1 + theme(plot.title = element_text(size = 16, hjust = 0.5))
  } else {
    p1 <- p1 + theme(plot.title = element_blank())
  }
  
  p3 <- group_fit_plot(
    model = grm, 
    theta_interval = theta_interval, 
    item_to_plot = item, 
    grouped_probs = NULL, 
    surprisal = TRUE,
    surp_base = surp_base,
    arc_length = TRUE,
    x_lab = x_lab,
    color_plot = TRUE
  ) + 
    geom_vline(xintercept = Qvec_al, color="black", linetype = "dashed") +
    geom_hline(yintercept = 0,  color="black", linetype = "dashed") +
    ggtitle(paste0(titlestr, item)) +
    ylim(c(0, 4)) +
    theme_bw() +
    theme(axis.title=element_text(size=labsize,face="bold"),
          legend.title = element_blank(),
          legend.position = lgdpos,
          axis.text  = element_text(size = 10),
          legend.text = element_text(size = 10))
  if (xlab == "") p3 <- p3 + theme(axis.title.x = element_blank())
  if (plotTitle) {
    p3 <- p3 + theme(plot.title = element_text(size = 16, hjust = 0.5))
  } else {
    p3 <- p3 + theme(plot.title = element_blank())
  }
  
  p2 <- group_fit_plot(
    model = grm, 
    theta_interval = theta_interval, 
    item_to_plot = item, 
    grouped_probs = NULL, 
    surprisal = FALSE,
    surp_base = surp_base,
    arc_length = FALSE,
    x_lab = TeX("$\\theta$"),
    color_plot = TRUE
  ) + 
    geom_vline(xintercept = Qvec, color="black", linetype = "dashed") +
    geom_hline(yintercept = 0.5,  color="black", linetype = "dashed") +
    ggtitle(paste0(titlestr, item)) +
    theme_bw() +
    theme(axis.title=element_text(size=labsize,face="bold"),
          legend.title = element_blank(),
          legend.position = lgdpos,
          axis.text  = element_text(size = 10),
          legend.text = element_text(size = 10))
  if (xlab == "") p2 <- p2 + theme(axis.title.x = element_blank())
  if (plotTitle) {
    p2 <- p2 + theme(plot.title = element_text(size = 16, hjust = 0.5))
  } else {
    p2 <- p2 + theme(plot.title = element_blank())
  }
  
  p4 <- group_fit_plot(
    model = grm, 
    theta_interval = theta_interval, 
    item_to_plot = item, 
    grouped_probs = NULL, 
    surprisal = FALSE,
    surp_base = surp_base,
    arc_length = TRUE,
    x_lab = x_lab,
    color_plot = TRUE
  ) + 
    geom_vline(xintercept = Qvec_al, color="black", linetype = "dashed") +
    geom_hline(yintercept = 0.5,  color="black", linetype = "dashed") +
    ggtitle(paste0(titlestr, item)) +
    theme_bw()+
    theme(axis.title=element_text(size=labsize,face="bold"),
          legend.title = element_blank(),
          legend.position = lgdpos,
          axis.text  = element_text(size = 10),
          legend.text = element_text(size = 10))
  if (xlab == "") p4 <- p4 + theme(axis.title.x = element_blank())
  if (plotTitle) {
    p4 <- p4 + theme(plot.title = element_text(size = 16, hjust = 0.5))
  } else {
    p4 <- p4 + theme(plot.title = element_blank())
  }
  
  return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4))
}

latent_group_probabilities <- function(
    model,
    data, 
    theta_scores, 
    groups = 10, 
    arc_length_scores = NULL
) {
  sorted_indices <- order(theta_scores)
  theta_scores <- theta_scores[sorted_indices]
  data <- data[sorted_indices, ]
  group_factor <- as.factor((seq_along(theta_scores) - 1) %/% (length(theta_scores) / groups))
  grouped_theta <- split(theta_scores, group_factor)
  grouped_data <- split(data, group_factor)
  
  if (!is.null(arc_length_scores)) {
    arc_length_scores <- arc_length_scores[sorted_indices]
    grouped_entropy <- split(arc_length_scores, group_factor)
    group_averages <- sapply(grouped_entropy, function(group) mean(group))
  } else {
    group_averages <- sapply(grouped_theta, function(group) mean(group))
  }
  
  data_probs <- list()
  model_probs <- list()
  for (item in seq_len(ncol(data))) {
    extr <- extract.item(model, item)
    item_model_probs <- lapply(grouped_theta, function(group_thetas) {
      probtrace(extr, group_thetas) |> 
        colMeans()
    })
    model_probs[[item]] <- do.call(rbind, item_model_probs)
    item_responses <- 0:(ncol(model_probs[[item]]) - 1)
    
    item_data_probs <- lapply(grouped_data, function(group) {
      table(factor(group[, item], levels = item_responses)) / nrow(group)
    })
    data_probs[[item]] <- do.call(rbind, item_data_probs)
  }
  
  return(list(data_probs = data_probs, model_probs = model_probs, group_averages = group_averages))
}


plot_bin_diff <- function(
    curve_scores,
    curve_probabilities, 
    item_group_means, 
    item_model_probs, 
    item_data_probs, 
    surprisal = TRUE,
    surp_base = 2,
    x_lab = TeX("$\\theta$"),
    y_lab = ifelse(surprisal, paste0("Surprisal (", surp_base, "bits)"), "Proportion/Probability"), 
    color_plot = TRUE, 
    x_max = NULL
) {
  if (surprisal) {
    curve_probabilities <- -log(curve_probabilities, surp_base)
    item_model_probs <- -log(item_model_probs, surp_base)
    item_data_probs <- -log(item_data_probs, surp_base)
  }
  
  ncat <- ncol(curve_probabilities)
  x_max <- ifelse(is.null(x_max), max(curve_scores), x_max)
  curve_probabilities <- cbind(curve_scores, curve_probabilities)
  curve_probabilities <- as.data.frame(curve_probabilities)
  colnames(curve_probabilities) <- c("score", 1:ncat)
  curve_probabilities <- pivot_longer(curve_probabilities, 1:ncat + 1, names_to = "category", values_to = "prob")
  
  item_model_probs <- as.data.frame(item_model_probs)
  colnames(item_model_probs) <- 1:ncat
  
  eo_df <- item_model_probs |>
    cbind(item_group_means) |>
    pivot_longer(cols = 1:ncat, names_to = "category", values_to = "e") |>
    cbind(o = as.vector(t(item_data_probs))) |>
    pivot_longer(cols = c("o", "e"), names_to = "oe", values_to = "prob") |>
    mutate(category = as.factor(category), oe = as.factor(oe))
  
  plot <- ggplot(eo_df) +
    # theme_bw() +
    geom_point(aes(x = item_group_means, y = prob, color = category, shape = oe),
               size = 3
    ) +
    geom_line(
      data = curve_probabilities,
      mapping = aes(x = score, y = prob, color = category),
      linewidth = 1
    ) +
    scale_shape_manual(
      values = c(19, 1),
      name = "Bin value", labels = c("Expected", "Observed")
    ) +
    xlab(x_lab) +
    ylab(y_lab) +
    scale_x_continuous(
      limits = c(NA, x_max),
      expand = expansion(mult = c(0, 0), add = c(0.0001, 0.0001))
    )
  if (color_plot) {
    plot <- plot + scale_color_discrete(name = "Score", labels = 0:10)
  } else {
    plot <- plot + scale_color_grey(name = "Score", start = 0.7, end = 0, labels = 0:10)
  }
  if (!surprisal) {
    plot <- plot +
      scale_y_continuous(
        limits = c(0, 1),
        expand = expansion(mult = c(0, 0), add = c(0.05, 0.05))
      )
  }
  if (surprisal) {
    plot <- plot +
      scale_y_continuous(
        limits = c(0, min(6, max(eo_df$prob))),
        expand = expansion(mult = c(0, 0), add = c(0.05, 0.05))
      )
  }
  plot
}

plot_curve <- function(
    curve_scores,
    curve_probabilities, 
    surprisal = TRUE,
    surp_base = 2,
    x_lab = TeX("$\\theta$"),
    y_lab = ifelse(surprisal, paste0("Surprisal (", surp_base, "bits)"), "Proportion/Probability"), 
    color_plot = TRUE, 
    x_max = NULL
) {
  if (surprisal) {
    curve_probabilities <- -log(curve_probabilities, surp_base)
  }
  
  ncat <- ncol(curve_probabilities)
  x_max <- ifelse(is.null(x_max), max(curve_scores), x_max)
  curve_probabilities <- cbind(curve_scores, curve_probabilities)
  curve_probabilities <- as.data.frame(curve_probabilities)
  colnames(curve_probabilities) <- c("score", 1:ncat)
  curve_probabilities <- pivot_longer(curve_probabilities, 1:ncat + 1, names_to = "category", values_to = "prob")
  
  
  plot <- ggplot(curve_probabilities, mapping = aes(x = score, y = prob, color = category)) +
    # theme_bw() +
    geom_line(linewidth = 1) +
    xlab(x_lab) +
    ylab(y_lab) +
    scale_x_continuous(
      limits = c(NA, x_max),
      expand = expansion(mult = c(0, 0), add = c(0.0001, 0.0001))
    )
  if (color_plot) {
    plot <- plot + scale_color_discrete(name = "Score", labels = 0:10)
  } else {
    plot <- plot + scale_color_grey(name = "Score", start = 0.7, end = 0, labels = 0:10)
  }
  if (!surprisal) {
    plot <- plot +
      scale_y_continuous(
        limits = c(0, 1),
        expand = expansion(mult = c(0, 0), add = c(0.05, 0.05))
      )
  }
  if (surprisal) {
    plot <- plot +
      scale_y_continuous(
        limits = c(0, min(6, max(curve_probabilities$prob))),
        expand = expansion(mult = c(0, 0), add = c(0.05, 0.05))
      )
  }
  plot
}

group_fit_plot <- function(
    model, 
    theta_interval, 
    item_to_plot,
    grouped_probs = NULL,
    plot_model_probs = TRUE,
    surprisal = FALSE,
    surp_base = 4,
    arc_length = TRUE,
    x_lab = TeX("$\\theta$"),
    y_lab = ifelse(surprisal, paste0("Surprisal (", surp_base, "bits)"), "Proportion/Probability"),
    color_plot = TRUE
) {
  # Calculate total categories
  total_cat <- sum(apply(model@Data$data, 2, max) + 1)
  
  if (arc_length) {
    scores <- theta2arclengrm(model, theta_interval, total_cat, surp_base = surp_base)[, 1]
  } else {
    scores <- theta_interval
  }
  
  extr <- extract.item(model, item_to_plot)
  
  if (is.null(grouped_probs)) {
    plot <- plot_curve(
      curve_scores = scores,
      curve_probabilities = probtrace(extr, theta_interval), 
      x_lab = x_lab,
      y_lab = y_lab,
      surprisal = surprisal,
      surp_base = surp_base,
      color_plot = color_plot,
      x_max = NULL
    )
  }
  else {
    if (plot_model_probs) {
      plot <- plot_bin_diff(
        curve_scores = scores,
        curve_probabilities = probtrace(extr, theta_interval), 
        item_group_means = grouped_probs$group_averages, 
        item_model_probs = grouped_probs$model_probs[[item_to_plot]],
        item_data_probs = grouped_probs$data_probs[[item_to_plot]],
        x_lab = x_lab,
        y_lab = y_lab,
        surprisal = surprisal,
        surp_base = surp_base,
        color_plot = color_plot,
        x_max = NULL
      )
    } else {
      plot <- plot_obs_probs(
        curve_scores = scores,
        curve_probabilities = probtrace(extr, theta_interval), 
        item_group_means = grouped_probs$group_averages, 
        item_data_probs = grouped_probs$data_probs[[item_to_plot]],
        x_lab = x_lab,
        y_lab = y_lab,
        surprisal = surprisal,
        surp_base = surp_base,
        color_plot = color_plot,
        x_max = NULL
      )
    }
  }
  
  return(plot)
}

plot_obs_probs <- function(
    curve_scores,
    curve_probabilities, 
    item_group_means, 
    item_data_probs, 
    surprisal = TRUE,
    surp_base = 2,
    x_lab = TeX("$\\theta$"),
    y_lab = ifelse(surprisal, paste0("Surprisal (", surp_base, "bits)"), "Proportion/Probability"), 
    color_plot = TRUE, 
    x_max = NULL
) {
  if (surprisal) {
    curve_probabilities <- -log(curve_probabilities, surp_base)
    item_data_probs <- -log(item_data_probs, surp_base)
  }
  
  ncat <- ncol(curve_probabilities)
  x_max <- ifelse(is.null(x_max), max(curve_scores), x_max)
  curve_probabilities <- cbind(curve_scores, curve_probabilities)
  curve_probabilities <- as.data.frame(curve_probabilities)
  colnames(curve_probabilities) <- c("score", 1:ncat)
  curve_probabilities <- pivot_longer(curve_probabilities, 1:ncat + 1, names_to = "category", values_to = "prob")
  
  item_data_probs <- as.data.frame(item_data_probs)
  colnames(item_data_probs) <- 1:ncat
  
  data_prob_df <- item_data_probs |>
    cbind(item_group_means) |>
    pivot_longer(cols = 1:ncat, names_to = "category", values_to = "prob")
  
  plot <- ggplot(data_prob_df) +
    # theme_bw() +
    geom_point(aes(x = item_group_means, y = prob, color = category),
               size = 3
    ) +
    geom_line(
      data = curve_probabilities,
      mapping = aes(x = score, y = prob, color = category),
      linewidth = 1
    ) +
    scale_shape_manual(
      values = c(19, 1),
      name = "Bin value", labels = c("Expected", "Observed")
    ) +
    xlab(x_lab) +
    ylab(y_lab) +
    scale_x_continuous(
      limits = c(NA, x_max),
      expand = expansion(mult = c(0, 0), add = c(0.0001, 0.0001))
    )
  if (color_plot) {
    plot <- plot + scale_color_discrete(name = "Score", labels = 0:10)
  } else {
    plot <- plot + scale_color_grey(name = "Score", start = 0.7, end = 0, labels = 0:10)
  }
  if (!surprisal) {
    plot <- plot +
      scale_y_continuous(
        limits = c(0, 1),
        expand = expansion(mult = c(0, 0), add = c(0.05, 0.05))
      )
  }
  if (surprisal) {
    plot <- plot +
      scale_y_continuous(
        limits = c(0, min(6, max(curve_probabilities$prob))),
        expand = expansion(mult = c(0, 0), add = c(0.05, 0.05))
      )
  }
  plot
}
