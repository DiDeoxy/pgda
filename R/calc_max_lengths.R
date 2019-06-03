#' Calc genome and homoeologous set lengths from chromosomes
#'
#' @param chrom_lengths a list of lengths by chrom in numeric order
#'
#' @return a list containing the lengths of the different genome and 
#' homoeolog sets
#'
#' @export

max_lengths <- function (chrom_lengths) {
  c(
    A = max(chrom_lengths[seq(1, 19, 3)]),
    B = max(chrom_lengths[seq(2, 20, 3)]),
    D = max(chrom_lengths[seq(3, 21, 3)]),
    one = max(chrom_lengths[c(1, 2, 3)]),
    two = max(chrom_lengths[c(4, 5, 6)]),
    three = max(chrom_lengths[c(7, 8, 9)]),
    four = max(chrom_lengths[c(10, 11, 12)]),
    five = max(chrom_lengths[c(13, 14, 15)]),
    six = max(chrom_lengths[c(16, 17, 18)]),
    seven = max(chrom_lengths[c(19, 20, 21)])
  )
}