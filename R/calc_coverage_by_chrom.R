#' Calculate mean allele richness at all sampling levels
#'
#' Calculates the mean allele richness across all markers for a sample at all
#' sampling levels. Based on the formula presented in 
#' https://www.genetics.org/content/157/1/389
#'
#' @param snp_chrom a factor vector inidicating the chrom of each marker
#' @param snp_pos a vector with the position of each snp on its chrom
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble add_row
#' 
#' @return a table of expected allele richness for each marker at each 
#' subsampling level with markers in rows and sampling levels in columns
#'
#' @export

coverage_by_chrom <- function (snp_chrom, snp_pos) {
  half_mean_dist <- (
    by(snp_pos, snp_chrom, max) %>% as.list() %>% unlist() %>% sum()
  ) / (2 * length(snp_chrom))

  by(snp_pos, snp_chrom, function (chrom_pos) {
    intervals <- tibble(start = double(), end = double())
    cur_interval <- list(start = double(), end = double())
    for (i in 1:length(chrom_pos)) {
      if (i == 1) {
        cur_interval$start <- max(0, chrom_pos[i] - half_mean_dist)
        cur_interval$end <- chrom_pos[i] + half_mean_dist
      } else if (i > 1 && i < length(chrom_pos)) {
        if ((chrom_pos[i] - half_mean_dist) <= cur_interval$end) {
          cur_interval$end <- chrom_pos[i] + half_mean_dist
        } else {
          intervals <- intervals %>%
            add_row(start = cur_interval$start, end = cur_interval$end)
          cur_interval$start <- chrom_pos[i] - half_mean_dist
          cur_interval$end <- chrom_pos[i] + half_mean_dist
        }
      } else {
        if ((chrom_pos[i] - half_mean_dist) <= cur_interval$end) {
          cur_interval$end <- chrom_pos[i] + half_mean_dist
          intervals <- intervals %>%
            add_row(start = cur_interval$start, end = cur_interval$end)
        } else {
          intervals <- intervals %>%
            add_row(start = cur_interval$start, end = cur_interval$end)
          cur_interval$start <- chrom_pos[i] - half_mean_dist
          cur_interval$end <- chrom_pos[i] + half_mean_dist
          intervals <- intervals %>%
            add_row(start = cur_interval$start, end = cur_interval$end)
        }
      }
    }
    (intervals[, 2] - intervals[, 1]) %>% sum()
  }) %>% as.list() %>% unlist()
}