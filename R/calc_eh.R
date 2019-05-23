#' Calculate expected heterozygosity
#'
#' Calculates expected heterozygosity of markers given a matrix where rows are
#' markers, and markers are haploid, encoded as 0 and 2
#'
#' @param genotypes a matrix of haploid markers as rows with alleles encoded as
#' 0 and 2
#'
#' @importFrom magrittr %>%
#' 
#' @return a vector of expected heterozygosities of the markers
#'
#' @export

calc_eh <- function (genotypes) {
  p <- rowSums(genotypes == 0)
  q <- rowSums(genotypes == 2)
  n <- (p + q)
  p_f <- p / n
  q_f <- 1 - p_f
  (n / (n - 1)) * (1 - (p_f ^ 2 + q_f ^ 2))
}