#' From: https://jeromyanglim.tumblr.com/post/50228877196/round-numbers-in-data-frame-that-contains-non

round_df <- function(x, digits) {
    # round all numeric variables
    # x: data frame 
    # digits: number of digits to round
    numeric_columns <- sapply(x, mode) == 'numeric'
    x[numeric_columns] <-  round(x[numeric_columns], digits)
    x
}