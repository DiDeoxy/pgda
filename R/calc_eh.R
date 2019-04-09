#' Calculate expected heterozygosity
#'
#' Calculates expected heterozygosity of markers given a matrix where rows are
#' markers, and markers are haploid, encoded as 0 and 2
#'
#' @param genotypes a matrix of haploid markers as rows with alleles encoded as
#' 0 and 2
#'
#' @return a vector of expected heterozygosities of the markers
#'
#' @export

calc_eh <- function (genotypes) {
  apply(genotypes, 1, function (snp) {
    n <- 2
    p <- sum(snp == 0) / sum(snp == 0 | 2)
    q <- sum(snp == 2) / sum(snp == 0 | 2)
    (n / (n - 1)) * (1 - (p^2 + q^2))
  })
}