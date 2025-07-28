##------------------------------------------------------------------------------
## Title: "fct_descr_analysis.R"                                                
## Description: Utility functions                
##-------------------------------------------------------------------------------

## Significance coding
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

## Unique minimum
##------------------------------------------------------------------------------
is_unique_min <- function(vec) {
  min_value <- min(vec)
  sum(vec == min_value) == 1
}

## Normalization to [0,1]
##------------------------------------------------------------------------------
normalize <- function(x) {
  if (max(x) == min(x)) {
    return(rep(0, length(x))) 
  } else {
    return((x - min(x)) / (max(x) - min(x)))
  }
}
