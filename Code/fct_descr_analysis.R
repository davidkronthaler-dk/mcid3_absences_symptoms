################################################################################
# Title: "fct_descr_analysis.R"                                                #
# Description: Set of utility functions needed for the analysis                #
################################################################################

## Function for significance coding
##------------------------------------------------------------------------------
signif_code <- Vectorize(function(p) {
  if (p <= 0.001) {
    "***"
  } else if (p <= 0.01) {
    "**"
  } else if (p <= 0.05) {
    "*"
  } else {
    ""
  }
}, vectorize.args = "p")

## Function to check whether a minimum is unique
##------------------------------------------------------------------------------
is_unique_min <- function(vec) {
  min_value <- min(vec)
  sum(vec == min_value) == 1
}

## Function to link respiratory viruses to absences
##------------------------------------------------------------------------------
link_absences_to_pathogens <- function(absences, range_days) {
  # 1. Extract all positive test results possibly linked to an absence by student id
  # 2. Check whether date of positive test result is "range_days" or less days before absence 
  # or after return from absence
  absences$test_id <- NA
  absences$pathogen <- NA
  absences$pathogen <- factor(absences$pathogen, levels = pathgs)
  
  for (i in 1:nrow(absences)) {
    t <- viral_loads[viral_loads$student_id == absences$student_id[i],
                     c("date", "pathogen", "test_id", "ct")]
    if (nrow(t) == 0) {
      absences$pathogen[i] <- NA
    } else {
      x <- abs(as.numeric(difftime(as.Date(t$date),
                                   as.Date(absences$date_absence[i]),
                                   units = "days")))
      y <- abs(as.numeric(difftime(as.Date(t$date),
                                   as.Date(absences$date_return[i]),
                                   units = "days")))
      if (min(c(x, y)) <= range_days) {
        miner <- apply(cbind(x, y), 1, min)
        if (is_unique_min(miner)) {
          absences$pathogen[i] <- t[which.min(miner), "pathogen"]
          absences$test_id[i] <- t[which.min(miner), "test_id"]
        } else {
          # For absences that cannot be uniquely linked to a pathogen, 
          # select the pathogen with the lowest CT value
          t <- t %>% arrange(ct)
          absences$pathogen[i] <- t[which.min(miner), "pathogen"]
          absences$test_id[i] <- t[which.min(miner), "test_id"]
        }
      } else {
        absences$pathogen[i] <- NA
      }
    }
  }
  return(absences)
}

## Function for normalization
##------------------------------------------------------------------------------
normalize <- function(x) {
  if (max(x) == min(x)) {
    return(rep(0, length(x))) 
  } else {
    return((x - min(x)) / (max(x) - min(x)))
  }
}
