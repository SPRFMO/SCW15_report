extract_tracking_array <- function(fl_mse_obj, metrics = NULL) {
  # Raw array and dimnames
  arr <- tracking(fl_mse_obj)@.Data
  dn <- dimnames(tracking(fl_mse_obj))
  
  # Get all metrics if none specified
  if (is.null(metrics)) metrics <- dn$metric
  
  # Indices for requested metrics
  metric_idx <- which(dn$metric %in% metrics)
  year_vals <- as.integer(dn$year)
  iter_vals <- as.integer(dn$iter)
  
  if (length(metric_idx) == 0) stop("None of the requested metrics are present.")
  
  # Create long-form list of data.frames
  out <- lapply(metric_idx, function(i) {
    mname <- dn$metric[i]
    slice <- arr[i, , 1, 1, 1, ]
    dimnames(slice) <- list(year = year_vals, iteration = iter_vals)
    df <- as.data.frame(as.table(slice))
    names(df)[3] <- mname
    df
  })
  
  # Merge and return
  merged <- Reduce(function(x, y) merge(x, y, by = c("year", "iteration")), out)
  merged <- merged[order(merged$year, merged$iteration), ]
  tibble::as_tibble(merged)
}

# Example usage for all metrics:
# df_all <- extract_tracking_array(bufflows[[1]])
# summary(df_all)
# # Example usage:
# df <- extract_tracking_array(run)
# names(df_all)
# df_all |> dplyr::filter(SB.om<1e5) |> ggplot(aes(x=(SB.om),y=(tac.hcr))) + geom_point(size=.8,alpha=.3) +
#   ggthemes::theme_few()
# summary(tracking(run))
