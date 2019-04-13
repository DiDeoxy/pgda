#' Calculate expected heterozygosity
#'
#' Calculates expected heterozygosity of markers given a matrix where rows are
#' markers, and markers are haploid, encoded as 0 and 2
#'
#' @param genotypes a matrix of haploid markers as rows with alleles encoded as
#' 0 and 2
#'
#' @importFrom parallel detectCores mclapply
#' 
#' @return a vector of expected heterozygosities of the markers
#'
#' @export

calc_eh <- function (genotypes) {
  print(genotypes)
  mclapply(1:nrow(genotypes), function (row) {
    n <- sum(genotypes[row, ] == 0 | 2)
    p <- sum(genotypes[row, ] == 0) / n
    q <- sum(genotypes[row, ] == 2) / n
    (n / (n - 1)) * (1 - (p^2 + q^2))
  }, mc.cores = detectCores())
}