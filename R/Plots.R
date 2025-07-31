#' @import patchwork
#' @import ggplot2
#' @import data.table
NULL

#' Plot OM and MP runs with faceted panels
#'
#' Combines Operating Model (OM) performance plot with Management Procedure (MP)
#' panels in a single layout.
#'
#' @param om_dt Data table containing OM performance data
#' @param mp_dt Data table containing MP performance data
#' @param select Metric to plot (default: "mean(C)")
#' @param y_label Y-axis label for OM plot (default: "Catch")
#' @param y_labelmp Y-axis label for MP plot (default: "Mean Catch")
#' @param buffer_filter Vector of buffer types to filter (default: c("cpue3", "cpue36", "cpue367"))
#' @param ncols Number of columns for MP facets (NULL for automatic)
#'
#' @return A patchwork combined plot object
#' @export
plotOMrunsDT <- function(om_dt, mp_dt,
                         select = "mean(C)",
                         y_label = "Catch",
                         y_labelmp = "Mean Catch",
                         buffer_filter = c("cpue3", "cpue36", "cpue367"),
                         ncols = NULL) {
  p_om <- plot_omperf(om_dt, select = select) + labs(y = y_label)
  mp_filtered <- get_tuned_mps(mp_dt, select_name = select, buffer_filter = buffer_filter)
  p_mp <- plot_mp_panels(mp_filtered, metric_label = y_labelmp, ncols = ncols)
  combined <- p_om / p_mp + plot_layout(heights = c(1, 2))
  return(combined)
}

#' Extract profile data from performance data table
#'
#' Filters and processes performance data for profiling analysis, extracting
#' target values from MP names.
#'
#' @param dt Data table containing performance data
#' @param group_pattern Regex pattern to identify MP groups (e.g., "cpue3_buffer")
#' @param target_pattern Regex pattern to identify targets in MP names (default: "target_")
#'
#' @return Processed data table with target values and MP groups extracted
#' @export
get_profile <- function(dt, group_pattern = "cpue3_buffer", target_pattern = "target_") {
  dt[, year := as.integer(year)]
  dt[, iter := as.integer(iter)]
  dt[, data := as.numeric(data)]
  
  prof_dt <- dt[grepl(group_pattern, mp) & grepl(target_pattern, mp)]
  prof_dt[, target := as.integer(sub(".*_target_", "", mp))]
  prof_dt[, mp_group := sub(".*_([a-z0-9]+_buffer)_.*", "\\1", mp)]
  
  return(prof_dt)
}

#' Plot a single profile metric
#'
#' Creates a line plot for a specific performance metric across target values.
#'
#' @param profile_summary Summary data table containing profile metrics
#' @param metric_name Name of the metric to plot (default: "P(Green)")
#' @param y_label Label for y-axis (default: "green")
#'
#' @return ggplot object
#' @export
plot_profile_metric <- function(profile_summary, metric_name = "P(Green)", y_label = "green") {
  ggplot(profile_summary[name == metric_name], aes(x = target, y = value)) +
    geom_line(color = "steelblue", linewidth = 1) +
    geom_point(color = "darkblue", size = 2) +
    labs(x = "CPUE Target", y = y_label, title = paste0(y_label, " vs CPUE Target")) +
    theme_minimal()
}

#' Plot a set of profile metrics
#'
#' Creates a composite plot showing three key metrics (P(Green), P(SB>SB[limit]), IAC(C))
#' across target values with automatic y-axis scaling.
#'
#' @param profile_dt Data table containing profile data
#' @param title Plot title
#' @param outpath Optional path to save plot (default: NULL returns plot object)
#' @param width Plot width in pixels when saving (default: 3000)
#'
#' @return patchwork plot object or saves to file if outpath provided
#' @export
plot_profile_set <- function(profile_dt, title, outpath = NULL, width = 3000) {
  summary_dt <- profile_dt[, .(value = mean(data, na.rm = TRUE)), by = .(target, name)]
  
  ymaxs <- summary_dt[name %in% c("P(Green)", "P(SB>SB[limit])", "IAC(C)"),
                      .(max_val = max(value, na.rm = TRUE)),
                      by = name]
  ymaxs[, ylim := ifelse(name == "IAC(C)", ceiling(max_val / 10) * 2 * 10, 1)]
  limits <- setNames(ymaxs$ylim, ymaxs$name)
  
  plot_profile_metric_scaled <- function(metric_name, y_label) {
    ggplot(summary_dt[name == metric_name], aes(x = target, y = value)) +
      geom_line(color = "steelblue", linewidth = 1) +
      geom_point(color = "darkblue", size = 2) +
      labs(x = "CPUE Target", y = y_label, title = NULL) +
      ylim(0, limits[[metric_name]]) +
      theme_minimal()
  }
  
  a <- plot_profile_metric_scaled("P(Green)", "P(Kobe = Green)")
  b <- plot_profile_metric_scaled("P(SB>SB[limit])", "P(SB > SB[limit])")
  c <- plot_profile_metric_scaled("IAC(C)", "Inter-Annual Catch Variability")
  
  p_combined <- a | b | c +
    plot_layout(nrow = 1) +
    plot_annotation(title = title)
  
  if (!is.null(outpath)) {
    pubpng(outpath, p_combined, width = width)
  } else {
    return(p_combined)
  }
}

#' Plot stacked metrics for selected MPs
#'
#' Creates a vertical stack of plots showing multiple metrics for selected
#' management procedures with uncertainty ribbons.
#'
#' @param dt Data table containing performance data
#' @param mp_list List of MPs to include
#' @param metric_names Vector of metric names (default: c("mean(C)", "SB", "F"))
#' @param y_labels Vector of y-axis labels (default: c("Catch (tonnes)", "SSB (tonnes)", "F"))
#' @param titles Vector of plot titles (default: c("Catches", "Spawning Biomass", "Fishing Mortality"))
#' @param probs Quantiles to plot (default: c(0.1, 0.25, 0.5, 0.75, 0.9))
#'
#' @return patchwork plot object
#' @export
plot_selected_metrics_stack <- function(dt, mp_list,
                                        metric_names = c("mean(C)", "SB", "F"),
                                        y_labels = c("Catch (tonnes)", "SSB (tonnes)", "F"),
                                        titles = c("Catches", "Spawning Biomass", "Fishing Mortality"),
                                        probs = c(0.1, 0.25, 0.5, 0.75, 0.9)) {
  dt[, year := as.integer(year)]
  dt[, iter := as.integer(iter)]
  dt[, data := as.numeric(data)]
  dt <- dt[mp %in% mp_list & name %in% metric_names]
  dt[, mp_label := factor(mp, levels = mp_list)]
  
  plot_metric <- function(metric, y_label, title) {
    dsub <- dt[name == metric]
    summary_dt <- dsub[, as.list(quantile(data, probs = probs, na.rm = TRUE)), by = .(year, mp_label)]
    qnames <- paste0("q", 100 * probs)
    setnames(summary_dt, old = paste0(100 * probs, "%"), new = qnames)
    
    ggplot(summary_dt, aes(x = year, fill = mp_label, color = mp_label)) +
      geom_ribbon(aes(ymin = get(qnames[1]), ymax = get(qnames[5])), alpha = 0.1, color = NA) +
      geom_ribbon(aes(ymin = get(qnames[2]), ymax = get(qnames[4])), alpha = 0.3, color = NA) +
      geom_line(aes(y = get(qnames[3])), linewidth = 0.5, alpha = 0.75) +
      geom_vline(xintercept = 2034, linetype = "dashed", color = "darkgrey") +
      geom_vline(xintercept = 2042, linetype = "dashed", color = "darkgrey") +
      labs(title = title, y = y_label, x = "Year", color = "MP", fill = "MP") +
      theme_minimal()
  }
  
  plots <- Map(plot_metric, metric_names, y_labels, titles)
  wrap_plots(plots, ncol = 1, guides = "collect") & theme(legend.position = "right")
}

#' Plot stacked metrics across target values
#'
#' Creates a vertical stack of plots showing multiple metrics across different
#' target values with uncertainty ribbons.
#'
#' @param dt Data table containing performance data
#' @param metric_names Vector of metric names (default: c("mean(C)", "SB", "F"))
#' @param y_labels Vector of y-axis labels (default: c("Catch (tonnes)", "SSB", "F"))
#' @param titles Vector of plot titles (default: c("Catches", "Spawning Biomass", "Fishing Mortality"))
#' @param probs Quantiles to plot (default: c(0.1, 0.25, 0.5, 0.75, 0.9))
#'
#' @return patchwork plot object
#' @export
plot_profile_targets_stack <- function(dt,
                                       metric_names = c("mean(C)", "SB", "F"),
                                       y_labels = c("Catch (tonnes)", "SSB", "F"),
                                       titles = c("Catches", "Spawning Biomass", "Fishing Mortality"),
                                       probs = c(0.1, 0.25, 0.5, 0.75, 0.9)) {
  dt[, year := as.integer(year)]
  dt[, iter := as.integer(iter)]
  dt[, data := as.numeric(data)]
  dt <- dt[name %in% metric_names]
  dt[, target_label := factor(target, levels = sort(unique(target)))]
  
  plot_metric <- function(metric, y_label, title) {
    dsub <- dt[name == metric]
    summary_dt <- dsub[, as.list(quantile(data, probs = probs, na.rm = TRUE)), by = .(year, target_label)]
    qnames <- paste0("q", 100 * probs)
    setnames(summary_dt, old = paste0(100 * probs, "%"), new = qnames)
    
    ggplot(summary_dt, aes(x = year, fill = target_label, color = target_label)) +
      geom_ribbon(aes(ymin = get(qnames[1]), ymax = get(qnames[5])), alpha = 0.1, color = NA) +
      geom_ribbon(aes(ymin = get(qnames[2]), ymax = get(qnames[4])), alpha = 0.3, color = NA) +
      geom_line(aes(y = get(qnames[3])), linewidth = 0.7, alpha = 0.8) +
      geom_vline(xintercept = 2034, linetype = "dashed", color = "darkgrey") +
      geom_vline(xintercept = 2042, linetype = "dashed", color = "darkgrey") +
      labs(title = title, y = y_label, x = "Year", color = "Target", fill = "Target") +
      theme_minimal()
  }
  
  plots <- Map(plot_metric, metric_names, y_labels, titles)
  wrap_plots(plots, ncol = 1, guides = "collect") & theme(legend.position = "right")
}

#' Plot summary of OM and MP performance
#'
#' Creates a composite plot showing OM trajectories, Kobe plot, and short/long-term
#' performance metrics.
#'
#' @param om Operating model object
#' @param runs MP runs object
#'
#' @return patchwork plot object
#' @export
plotSummary <- function(om, runs) {
  (plot(window(om, start = 2010), runs) + guides(colour = FALSE)) +
    kobeMPs(performance(runs)) +
    (plotBPs(performance(runs)[year == "short"],
             statistics = c("green", "SBMSY", "C", "IACC")) + 
       guides(fill = FALSE) + ggtitle("2024-2028")) +
    (plotBPs(performance(runs)[year == "long"],
             statistics = c("green", "SBMSY", "C", "IACC")) + 
       guides(fill = FALSE) + ggtitle("2034-2042")) +
    plot_layout(design = "AAA\nBCD", guides = "collect")
}

#' Plot Operating Model performance
#'
#' Creates a plot showing OM performance metrics with uncertainty ribbons.
#'
#' @param dt Data table containing OM performance data
#' @param select Metric(s) to plot (default: "mean(C)")
#' @param probs Quantiles to plot (default: c(0.1, 0.25, 0.5, 0.75, 0.9))
#'
#' @return ggplot object
#' @export
plot_omperf <- function(dt, select = c("mean(C)"), probs = c(0.1, 0.25, 0.5, 0.75, 0.9)) {
  dt[, year := as.integer(year)]
  dt[, iter := as.integer(iter)]
  dt[, data := as.numeric(data)]
  dt_subset <- dt[name %in% select]
  summary_dt <- dt_subset[, as.list(quantile(data, probs = probs, na.rm = TRUE)), by = .(year, name)]
  setnames(summary_dt, old = paste0(100 * probs, "%"), new = paste0("q", 100 * probs))
  ggplot(summary_dt, aes(x = year)) +
    geom_ribbon(aes(ymin = q10, ymax = q90, fill = name), alpha = 0.2) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = name), alpha = 0.4) +
    geom_line(aes(y = q50, color = name), linewidth = 1) +
    labs(y = select, title = "Operating Model") +
    scale_fill_discrete(name = NULL) +
    scale_color_discrete(name = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

#' Plot MP performance panels
#'
#' Creates faceted plots showing MP performance across different tuning levels
#' and buffer types.
#'
#' @param dt Data table containing MP performance data
#' @param metric_label Y-axis label (default: "mean(C)")
#' @param probs Quantiles to plot (default: c(0.1, 0.25, 0.5, 0.75, 0.9))
#' @param ncols Number of facet columns (NULL for automatic)
#'
#' @return ggplot object
#' @export
plot_mp_panels <- function(dt, metric_label = "mean(C)",
                           probs = c(0.1, 0.25, 0.5, 0.75, 0.9),
                           ncols = NULL) {
  summary_dt <- dt[, as.list(quantile(data, probs = probs, na.rm = TRUE)),
                   by = .(year, tune_level, buffer_type)]
  setnames(summary_dt, old = paste0(100 * probs, "%"), new = paste0("q", 100 * probs))
  p <- ggplot(summary_dt, aes(x = year)) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = "salmon", alpha = 0.2) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "salmon", alpha = 0.4) +
    geom_line(aes(y = q50), color = "black", linewidth = 0.8) +
    geom_vline(xintercept = 2034, linetype = "dashed", color = "darkgrey", linewidth = 0.8) +
    geom_vline(xintercept = 2042, linetype = "dashed", color = "darkgrey", linewidth = 0.8) +
    labs(title = "Management Procedures", y = metric_label, x = "Year") +
    theme_minimal()
  
  if (is.null(ncols)) {
    p <- p + facet_grid(rows = vars(tune_level), cols = vars(buffer_type))
  } else {
    p <- p + facet_grid(rows = vars(tune_level), cols = vars(buffer_type), labeller = "label_value") +
      theme(strip.text = element_text(size = 10))
  }
  return(p)
}

#' Plot stacked C, SSB, and F metrics
#'
#' Creates a vertical stack of plots showing catch, spawning biomass, and fishing
#' mortality metrics with optional OM comparison.
#'
#' @param dt Data table containing MP performance data
#' @param om_dt Optional data table containing OM performance data (default: NULL)
#' @param metric_names Vector of metric names (default: c("mean(C)", "SB", "F"))
#' @param y_labels Vector of y-axis labels (default: c("Catch (t)", "SSB (t)", "F"))
#' @param titles Vector of plot titles (default: c("Catches", "Spawning Biomass", "Fishing Mortality"))
#' @param probs Quantiles to plot (default: c(0.1, 0.25, 0.5, 0.75, 0.9))
#'
#' @return patchwork plot object
#' @export
plot_C_SSB_F_stack <- function(dt, om_dt = NULL,
                               metric_names = c("mean(C)", "SB", "F"),
                               y_labels = c("Catch (t)", "SSB (t)", "F"),
                               titles = c("Catches", "Spawning Biomass", "Fishing Mortality"),
                               probs = c(0.1, 0.25, 0.5, 0.75, 0.9)) {
  dt[, year := as.integer(year)]
  dt[, iter := as.integer(iter)]
  dt[, data := as.numeric(data)]
  dt <- dt[name %in% metric_names]
  dt[, mp_label := factor(mp, levels = unique(mp))]
  
  if (!is.null(om_dt)) {
    om_dt[, year := as.integer(year)]
    om_dt[, iter := as.integer(iter)]
    om_dt[, data := as.numeric(data)]
    om_dt <- om_dt[name %in% metric_names]
    om_dt[, mp_label := "OM"]
    
    last_om_year <- max(om_dt$year)
    first_mp_year <- min(dt$year)
    if (last_om_year + 1 < first_mp_year) {
      message("⚠️ Gap between OM and MP years: ", last_om_year, " -> ", first_mp_year)
    }
    
    dt <- rbind(om_dt, dt, fill = TRUE)
  }
  
  plot_metric <- function(metric, y_label, title) {
    dsub <- dt[name == metric]
    summary_dt <- dsub[, as.list(quantile(data, probs = probs, na.rm = TRUE)), by = .(year, mp_label)]
    qnames <- paste0("q", 100 * probs)
    setnames(summary_dt, old = paste0(100 * probs, "%"), new = qnames)
    
    p <- ggplot(summary_dt, aes(x = year, fill = mp_label, color = mp_label)) +
      geom_ribbon(aes(ymin = get(qnames[1]), ymax = get(qnames[5])), alpha = 0.1, color = NA) +
      geom_ribbon(aes(ymin = get(qnames[2]), ymax = get(qnames[4])), alpha = 0.3, color = NA) +
      geom_line(aes(y = get(qnames[3])), linewidth = 0.7, alpha = 0.8) +
      geom_vline(xintercept = 2034, linetype = "dashed", color = "darkgrey") +
      geom_vline(xintercept = 2042, linetype = "dashed", color = "darkgrey") +
      labs(title = title, y = y_label, x = "Year", color = "Scenario", fill = "Scenario") +
      theme_minimal()
    
    if (metric == "F") {
      ymax <- max(summary_dt$q90, na.rm = TRUE)
      buffer <- 0.05 * ymax
      p <- p + coord_cartesian(ylim = c(0, ymax + buffer))
    }
    return(p)
  }
  
  plots <- Map(plot_metric, metric_names, y_labels, titles)
  wrap_plots(plots, ncol = 1, guides = "collect") & theme(legend.position = "right")
}