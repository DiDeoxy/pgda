#' Calculate LNG
#'
#' Calculates the length of each chromsome in each genome, the number of
#' markers and the gaps between neighbouring markers
#'
#' @param snp_data the parsed snp data of the gds object
#'
#' @return a list contianing the above data organised by genome

calc_lng <- function(snp_data, scale) {
  # number of snps and mean distances between the genome
  lng <- by(snp_data, snp_data$chrom,
    function (chrom_data) {
      list(
        leng = max(chrom_data$pos / scale) - min(chrom_data$pos / scale),
        num = length(chrom_data$pos / scale),
        gaps = diff(chrom_data$pos / scale)
      )
    }
  )
  ret <- list(
    A = list(leng = vector(), num = vector(), gaps = vector()),
    B = list(leng = vector(), num = vector(), gaps = vector()),
    D = list(leng = vector(), num = vector(), gaps = vector())
  )
  for (i in 1:length(lng)) {
    if (i %in% seq(1, 19, 3)) {
      ret$A$leng <- c(ret$A$leng, lng[[i]]$leng)
      ret$A$num <- c(ret$A$num, lng[[i]]$num)
      ret$A$gaps <- c(ret$A$gaps, lng[[i]]$gaps)
    } else if (i %in% seq(2, 20, 3)) {
      ret$B$leng <- c(ret$B$leng, lng[[i]]$leng)
      ret$B$num <- c(ret$B$num, lng[[i]]$num)
      ret$B$gaps <- c(ret$B$gaps, lng[[i]]$gaps)
    } else {
      ret$D$leng <- c(ret$D$leng, lng[[i]]$leng)
      ret$D$num <- c(ret$D$num, lng[[i]]$num)
      ret$D$gaps <- c(ret$D$gaps, lng[[i]]$gaps)
    }
  }
  ret
}