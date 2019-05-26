#' Calculate mean allele richness at all sampling levels
#'
#' Calculates the mean allele richness across all markers for a sample at all
#' sampling levels. Missing data is not allowed. Based on the formula presented
#' in https://www.genetics.org/content/157/1/389
#'
#' @param pop a data frame with individuals in columns and markers in rows,
#' there must be atleast two individuals
#' @param allele_coding the coding used for indicating the different alleles
#' @param num_cores the number of cores to use with mclapply(), must be 1 on
#' windows, default is detectCores()
#'
#' @importFrom magrittr %>%
#' @importFrom parallel detectCores mclapply
#' @importFrom scrime rowTables
#' 
#' @return a table of expected allele richness for each marker at each 
#' subsampling level with markers in rows and sampling levels in columns
#'
#' @export

allele_richness <- function (
  pop, allele_coding = 1:2, num_cores = detectCores()
) {
  # the total number of alleles observed at each marker
  n <- ncol(pop)
  # probs contains the probability of not observing allele i at each 
  # sub-sampling level (n - k) for each possible count of allele i with
  # allele count in rows and k in columns
  #
  # for each subsampling level
  probs <- mclapply(0:(n - 1), function (k) {
    # a vector for containng the probs of not observing allele i at each count
    # level at each subsampling (n - k) level
    inter <- rep(0, n)
    # if n - k <= 1 then the prob of not observing allele i is 0 at all count
    # levels, the smaller k is compared to n the more levels will have probs of
    # not observing allel i greater than 0
    if (n - k > 1) {
      # probs of not observing allele i are linear decreasing, therefore the top
      # half and bottom half are 1 - mirrors, we can use this fact to skip a lot
      # of computation
      temp <- lapply(1:floor((n - k) / 2), function (n_i) {
        (n - n_i - k) / (n - k)
      }) %>% do.call(c, .)
      # concatenating the calced probs with their 1 - mirror, if n - k is odd
      # the middle value will equal 0.5 which we do not need to mirror
      temp <- c(temp, rev(1 - temp[which(temp != 0.5)]))
      inter[1:length(temp)] <- temp
      inter
    } else {
      inter
    }
  }, mc.cores = num_cores) %>% do.call(cbind, .)
  # creates a data frame containg the count of each allele for each marker
  marker_allele_count <- rowTables(pop, allele_coding)
  # we calcuate the mean allele richness across all markers at each subsampling
  # level (n - k) by calculating the product of not observing each allele at
  # each sub-sampling level then taking the sum of these for each marker and 
  # then taking the mean across all markers
  #
  # for each marker
  mclapply(1:nrow(marker_allele_count), function (marker) {
    (1 - lapply(1:length(marker_allele_counts[marker, ]), function (allele) {
      # for each allele, calc the probability of not observing the allele at
      # each sub-sampling level
      cumprod(probs[marker_allele_count[[marker, allele]], ])
    # rbind the probabilities for each allele at each sub-smapling level,
    # subtract from one to turn them into probabilities of observing the allele,
    # and sum the alleles together
    }) %>% do.call(rbind, .)) %>% colSums()
  # return a table withexpected allele ricness in rows and sub-sampling levels
  # in columns
  }, mc.cores = num_cores) %>% do.call(rbind, .)
}
