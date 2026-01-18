get_outcome_status_at_specific_time <- function(df, status_field, time_field, specific_time) {
  specific_time <- as.numeric(specific_time)
  if ((specific_time <= min(as.numeric(df[,time_field]), na.rm = TRUE)) | (specific_time >= max(as.numeric(df[,time_field]), na.rm = TRUE))) {
    output <- list(outcome = "Unsuccessful", message = "The specific time must be more than the minimum time and less than the maximum time.", new_data = df)
  } else {
    current_time <- df[,time_field]
    current_status <- df[, status_field]
    if (! is.factor(current_status)) {
      factor_levels <- sort(unique(current_status))
      factor_levels[! is.na(factor_levels)]
      current_status <- factor(current_status, levels = factor_levels)
    }
    new_time <- rep(specific_time, nrow(df))
    new_status <- rep(NA, nrow(df))
    # If the time is less than or equal to the specific time, time and status is not changed
    new_time[! is.na(current_time) & (current_time <= specific_time)] <- current_time[! is.na(current_time) & (current_time <= specific_time)]
    new_status[! is.na(current_status) & (current_time <= specific_time)] <- as.character(current_status[! is.na(current_status) & (current_time <= specific_time)])
    # If the time is more than the specific time, the time is the specific time and whether it is event or 'no event', the new status is 'no event'
    new_time[! is.na(current_time) & (current_time > specific_time)] <- specific_time
    new_status[! is.na(current_status) & (current_time > specific_time)] <- as.character(levels(current_status)[1])
    new_status <- factor(new_status, levels = levels(current_status))
    df <- cbind.data.frame(df, df[,c(time_field, status_field)])
    colnames(df)[(ncol(df)-1):ncol(df)] <- paste0("unmodified_", c(time_field, status_field))
    df[,time_field] <- new_time
    df[,status_field] <- new_status
    output <- list(outcome = "Successful", message = paste0("The time field and status field are replaced with new time and status. The unmodified fields are available as new columns: ", paste0("'unmodified_", c(time_field, status_field), "'", collapse = " and "), "."), new_data = df)
  }
}
