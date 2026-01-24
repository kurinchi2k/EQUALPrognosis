prepare_datasets <- function(df, simulations, outcome_name, outcome_type, outcome_time = NULL, verbose = TRUE) {
  {
    # Convert this to dataframe
    df <- data.frame(df, check.names = FALSE)
    # Remove the subjects without the outcome
    df <- df[! is.na(df[,outcome_name]),]
    # If outcome_type is time-to-event, then time to event is a minimum requirement. Therefore, remove subjects without time
    if (outcome_type == "time-to-event") {
      df <- df[! is.na(df[,outcome_time]),]
    }
    # Assign an unique id for each subject
    df$PTE_id <- paste0(formatC(1:nrow(df), width = 6, flag = 0))
    # Scale quantitative outcomes
    if (outcome_type == "quantitative") {
      df[,outcome_name] <- scale(df[,outcome_name], center = FALSE, scale = TRUE)
    }
  }
  # Sampling
  {
    simulation_number <- 1
    working_seeds <- 0
    df_training_list <- list()
    df_only_validation_list <- list()
    df_all_subjects_list <- list()
    if (verbose == TRUE) {cat("\nWorking samples\n")}
    repeat {
      # split_factor identifies the rows to include at random, noting that the same subject may be chosen more than one to achieve the desired sample size
      split_factor <- sample(1:nrow(df), size = nrow(df), replace = TRUE, prob = NULL)
      # Include the rows in split_factor: this results in database df_training of same sample size as original and includes a random sample in that database
      df_training <- df[split_factor,]
      # Check whether there are two levels of outcomes in the training data set but this applies only if the outcome is not quantitative
      if (outcome_type == "quantitative") {
        two_levels <- TRUE
      } else {
        two_levels <- (nlevels(as.factor(as.character(df_training[,outcome_name]))) == 2)
      }
      # If there are two levels, also obtain the samples not included in df_training, which is called 'df_validation' and all subjects which is the original database
      if (two_levels == TRUE) {
        df_only_validation <- df[-split_factor,]
        df_all_subjects <- df
        # Additional modifications are required to combine sparse levels in ordinal factors
        each_column <- lapply(1:ncol(df_training), function(y) {
          # If it is ordinal factor, it combines the missing level to the next level (except for the last level where it combines with the previous level) if not, no change is made
          # If it is not ordinal, there is no change
          # For ordinal factors
          if (is.ordered(df_training[,y]) == TRUE) {
            # Get the levels
            summary <- rbind.data.frame(
              table(df_training[,y]),
              table(df_only_validation[,y])
            )
            # Loop through each level
            for (i in 1:nlevels(df_training[,y])) {
              # Assign a variable j which is the next level unless it is the last level in which case it is the previous level
              j <- ifelse(i < nlevels(df_training[,y]), i+1, i-1)
              # If there are no subjects at this level in the training data set but there are such subjects in the validation data set, assign the next level (if level is not the last level and to the previous level if the level is the last level)
              if (summary[1,i] == 0) {
                if (summary[2,i] > 0) {
                  df_only_validation[((! is.na(df_only_validation[,y])) &
                                        (df_only_validation[,y] == levels(df_training[,y])[i]))
                                     ,y] <- levels(df_training[,y])[j]
                  # As these out-of-sample patients with excluded levels of ordinal data in training data set are present in the original dataset, these must also be changed as for df_only_validation
                  df_all_subjects[((! is.na(df_all_subjects[,y])) &
                                     (df_all_subjects[,y] == levels(df_training[,y])[i]))
                                  ,y] <- levels(df_training[,y])[j]
                }
              }
            }
            # Get the changed variable
            output <- list(df_only_validation_column = data.frame(df_only_validation[,y]),
                           df_all_subjects_column = data.frame(df_all_subjects[,y]))
          } else {
            # Get the unchanged variable
            output <- list(df_only_validation_column = data.frame(df_only_validation[,y]),
                           df_all_subjects_column = data.frame(df_all_subjects[,y]))
          }
          return(output)
        })
        # Combine the columns to update the df_only_validation and df_all_subjects
        df_only_validation <- do.call(cbind.data.frame, lapply(each_column, "[[",1))
        df_all_subjects <- do.call(cbind.data.frame, lapply(each_column, "[[",2))
        # Rename the columns to the original column names
        colnames(df_only_validation) <- colnames(df_training)
        colnames(df_all_subjects) <- colnames(df_training)
        # Add to the collection
        df_training_list <- c(df_training_list, list(df_training))
        df_only_validation_list <- c(df_only_validation_list, list(df_only_validation))
        df_all_subjects_list <- c(df_all_subjects_list, list(df_all_subjects))
        # Add to working seeds
        working_seeds <- working_seeds + 1
      }
      if (verbose == TRUE) {cat(paste0(working_seeds, "..."))}
      # Exit loop if we have the desired simulations
      if (working_seeds >= simulations) {break}
      # Go to the next simulation
      simulation_number <- simulation_number + 1
    }
  }
  output <- list(df_training_list = df_training_list, df_only_validation_list = df_only_validation_list,
                 df_all_subjects_list = df_all_subjects_list)
}
