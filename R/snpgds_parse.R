#' Parse a gds file
#'
#' Loads the data from a gds object into a form amenable to manipulation
#' using R code
#'
#' @importFrom dplyr tibble
#' @importFrom magrittr %>%
#' @importFrom SNPRelate snpgdsClose snpgdsOpen
#' @importFrom gdsfmt index.gdsn read.gdsn
#'
#' @param gds_file the path to the target gds file
#'
#' @return a list of lists containing the data of the gds object
#'
#' @export

snpgds_parse <- function(gds_file) {
  gds <- snpgdsOpen(gds_file)

  # make a tibble of the snp data
  snp <- tibble(
    id = read.gdsn(index.gdsn(gds, "snp.id")) %>% as.character(),
    chrom = read.gdsn(index.gdsn(gds, "snp.chromosome")),
    pos = read.gdsn(index.gdsn(gds, "snp.position"))
  )

  chrom_lengths <- by(snp$pos, snp$chrom, max)

  # convert chroms values from integer to names
  chroms <- outer(as.character(1:7), c("A", "B", "D"), paste, sep = "") %>% 
    t() %>% as.vector()
  for (i in 1:21) {
    snp$chrom[which(snp$chrom == i)] <- chroms[i]
  }

  # extract the genotypes
  genotypes <- read.gdsn(index.gdsn(gds, "genotype"))

  # extract the categorical information
  sample_annot <- read.gdsn(index.gdsn(gds, "samp_annot"))
  sample_id <- as.character(read.gdsn(index.gdsn(gds, "sample.id")))
  snpgdsClose(gds)

  return(
    list(
      snp = snp, genotypes = genotypes,
      chrom_lengths = chrom_lengths,
      sample = list(id = sample_id, annot = sample_annot)
    )
  )
}