#' Calculate the length of each chrom, optionally subtract the min position from
#' it
#'
#' @param snp_chrom a factor vector inidicating the chrom of each marker
#' @param snp_pos a vector with the position of each snp on its chrom
#' 
#' @return a vector of span distances
#'
#' @export

span_by_chrom <- function (snp_chrom, snp_pos, diff = FALSE) {
  if (diff) {
    by(snp_pos, snp_chrom,
      function(chrom_pos) {
        max(chrom_pos) - min(chrom_pos)
      }
    ) %>% as.list() %>% unlist()
  } else {
    by(snp_pos, snp_chrom, max) %>% as.list() %>% unlist()
  }
}