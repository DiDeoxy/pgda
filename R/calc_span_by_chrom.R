#' Calculate the difference between the max and min positions on a chromosome
#'
#' @param snp_chrom a factor vector inidicating the chrom of each marker
#' @param snp_pos a vector with the position of each snp on its chrom
#' 
#' @return a vector of span distances
#'
#' @export

span_by_chrom <- function (snp_chrom, snp_pos) {
  by(snp_pos, snp_chrom,
    function(chrom_pos) {
      max(chrom_pos) - min(chrom_pos)
    }
  ) %>% as.list() %>% unlist()
}

