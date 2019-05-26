#' Calculate mean allele richness at all sampling levels
#'
#' Calculates the mean allele richness across all markers for a sample at all
#' sampling levels. Missing data is not allowed. Based on the formula presented
#' in https://www.genetics.org/content/157/1/389
#'
#' @param pop a data frame with individuals in columns and markers in rows,
#' there must be atleast two individuals
#' @param allele_coding the coding used for indicating the different alleles
#' @param num_cores the number of cores to use, must be 1 on windows, can use 
#'  qdetectCores() of the parallel package on linux
#'
#' @importFrom magrittr %>%
#' @importFrom parallel mclapply
#' @importFrom scrime rowTables
#' 
#' @return a table of expected allele richness for each marker at each 
#' subsampling level with markers in rows and sampling levels in columns
#'
#' @export

allele_richness <- function (pop, allele_coding = 1:2, num_cores = 1) {
  # the total number of alleles observed at each marker
  n <- ncol(pop)
  # probs contains the probability of not observing allele i at each 
  # sub-sampling level (n - k) for each possible count of allele i
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
  # creates a data frame containg the counts of each allele for each marker
  markers <- rowTables(pop, coding)
  # we calcuate the mean allele richness across all markers at each subsampling
  # level (n - k) by calculating the product of not observing each allele at
  # each sub-sampling level then taking the sum of these for each marker and 
  # then taking the mean across all markers
  #
  # k here is equivalent to that above, column 1 of probs contains probs of not
  # observing allele i when k is equal to zero
  mclapply(1:n, function (k) {
    # for each marker
    lapply(1:nrow(markers), function (marker) {
      # for each allele, calc the product of the probabilities of not observing
      # it at all sub-sampling levels up to (n - k), then sum the allel products
      (1 - lapply(1:length(markers[marker, ]), function (allele) {
        prod(probs[markers[[marker, allele]], 1:k])
      }) %>% unlist()) %>% sum()
    # take the mean acoss all markers
    }) %>% unlist()
  }, mc.cores = num_cores) %>% do.call(cbind, .)
}
