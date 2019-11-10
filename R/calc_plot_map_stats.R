#' Calculate and plot map stats
#'
#' Calcluates various map statistics creating plots and tables
#'
#' @importFrom dplyr tibble
#' @importFrom readr write_csv
#' @importFrom stringr str_c
#'
#' @param phys_gds the file path to the phys gds object
#' @param gen_gds the file path to the gen gds object
#' @param out_name the base name to use for the csv
#' @param plot_title the title of the combined plots
#' @param y_lim the limit for the y axis of the plots
#'
#' @return Outputs plots and tables of the data
#'
#' @export

calc_plot_map_stats <- function (
  phys_gds, gen_gds, out_name, plot_title, y_lim
) {
  phys_data <- snpgds_parse(phys_gds)
  gen_data <- snpgds_parse(gen_gds)

  # calc stats
  phys_lng <- calc_lng(phys_data$snp, 1e6)
  gen_lng <- calc_lng(gen_data$snp, 100)
  genome_ld <- calc_ld_stats(phys_gds, phys_data$snp)

  # plot the ld
  plot_gaps_nbs_ld(
    phys_lng, gen_lng, genome_ld, out_name, plot_title, y_lim
  )

  # calc chrom coverage
  phys_cov <- coverage_by_chrom(phys_data$snp$chrom, phys_data$snp$pos) / 1e6
  gen_cov <- coverage_by_chrom(gen_data$snp$chrom, gen_data$snp$pos) / 100

  # find the maf and mr
  maf_mr <- calc_maf_mr(phys_data)

  # find the min length of the top percentile of gaps
  phys_tp <- quantile(
    c(phys_lng$A$gaps, phys_lng$B$gaps, phys_lng$D$gaps),
    prob = 0.99, na.rm = TRUE
  )
  gen_tp <- quantile(
    c(gen_lng$A$gaps, gen_lng$B$gaps, gen_lng$D$gaps),
    prob = 0.99, na.rm = TRUE
  )

  map_stats <- tibble(
    "Genome" = c("A", "B", "D", "All"),
    "MAF" = c(
      mean(maf_mr$A$maf),
      mean(maf_mr$B$maf),
      mean(maf_mr$D$maf),
      mean(c(maf_mr$A$maf, maf_mr$B$maf, maf_mr$D$maf))
    ),
    "MR" = c(
      mean(maf_mr$A$mr),
      mean(maf_mr$B$mr),
      mean(maf_mr$D$mr),
      mean(c(maf_mr$A$mr, maf_mr$B$mr, maf_mr$D$mr))
    ),
    "Span (Gb)" = c(
      sum(phys_lng$A$leng) / 1000,
      sum(phys_lng$B$leng) / 1000,
      sum(phys_lng$D$leng) / 1000,
      sum(phys_lng$A$leng, phys_lng$B$leng, phys_lng$D$leng) / 1000
    ),
    "% Coverage (Pseudo-chrom.)" = c(
      phys_cov[seq(1, 19, 3)] %>% sum() / sum(phys_lng$A$leng),
      phys_cov[seq(2, 20, 3)] %>% sum() / sum(phys_lng$B$leng),
      phys_cov[seq(3, 21, 3)] %>% sum() / sum(phys_lng$D$leng),
      phys_cov %>% sum() /
        sum(phys_lng$A$leng, phys_lng$B$leng, phys_lng$D$leng)
    ) * 100,
    "Span (cM)" = c(
      sum(gen_lng$A$leng),
      sum(gen_lng$B$leng),
      sum(gen_lng$D$leng),
      sum(gen_lng$A$leng, gen_lng$B$leng, gen_lng$D$leng)
    ),
    "% Coverage (Genetic)" = c(
      gen_cov[seq(1, 19, 3)] %>% sum() / sum(gen_lng$A$leng),
      gen_cov[seq(2, 20, 3)] %>% sum() / sum(gen_lng$B$leng),
      gen_cov[seq(3, 21, 3)] %>% sum() / sum(gen_lng$D$leng),
      gen_cov %>% sum() / sum(gen_lng$A$leng, gen_lng$B$leng, gen_lng$D$leng)
    ) * 100,
    "Num. SNPs" = c(
      sum(phys_lng$A$num),
      sum(phys_lng$B$num),
      sum(phys_lng$D$num),
      sum(phys_lng$A$num, phys_lng$B$num, phys_lng$D$num)
    ),
    "Mean Gap Size (Mb)" = c(
      mean(phys_lng$A$gaps),
      mean(phys_lng$B$gaps),
      mean(phys_lng$D$gaps),
      mean(c(phys_lng$A$gaps, phys_lng$B$gaps, phys_lng$D$gaps))
    ),
    "Mean Gap Size (cM)" = c(
      mean(gen_lng$A$gaps),
      mean(gen_lng$B$gaps),
      mean(gen_lng$D$gaps),
      mean(c(gen_lng$A$gaps, gen_lng$B$gaps, gen_lng$D$gaps))
    ),
    "Num. Top 1% Gaps (Pseudo-chrom.)" = c(
      sum(phys_lng$A$gaps >= phys_tp),
      sum(phys_lng$B$gaps >= phys_tp),
      sum(phys_lng$D$gaps >= phys_tp),
      sum(
        c(
          phys_lng$A$gaps >= phys_tp,
          phys_lng$B$gaps >= phys_tp,
          phys_lng$D$gaps >= phys_tp
        )
      )
    ),
    "Num. Top 1% Gaps (Genetic)" = c(
      sum(gen_lng$A$gaps >= gen_tp),
      sum(gen_lng$B$gaps >= gen_tp),
      sum(gen_lng$D$gaps >= gen_tp),
      sum(
        c(
          gen_lng$A$gaps >= gen_tp,
          gen_lng$B$gaps >= gen_tp,
          gen_lng$D$gaps >= gen_tp
        )
      )
    ),
    "Min Length Top 1% Gap (Mb)" = c(
      min(phys_lng$A$gaps[which(phys_lng$A$gaps >= phys_tp)]),
      min(phys_lng$B$gaps[which(phys_lng$B$gaps >= phys_tp)]),
      min(phys_lng$D$gaps[which(phys_lng$D$gaps >= phys_tp)]),
      phys_tp
    ),
    "Min Length Top 1% Gap (cM)" = c(
      min(gen_lng$A$gaps[which(gen_lng$A$gaps >= gen_tp)]),
      min(gen_lng$B$gaps[which(gen_lng$B$gaps >= gen_tp)]),
      min(gen_lng$D$gaps[which(gen_lng$D$gaps >= gen_tp)]),
      gen_tp
    ),
    "Max Length Gap (Mb)" = c(
      max(phys_lng$A$gaps),
      max(phys_lng$B$gaps),
      max(phys_lng$D$gaps),
      max(c(phys_lng$A$gaps, phys_lng$B$gaps, phys_lng$D$gaps))
    ),
    "Max Length Gap (cM)" = c(
      max(gen_lng$A$gaps),
      max(gen_lng$B$gaps),
      max(gen_lng$D$gaps),
      max(c(gen_lng$A$gaps, gen_lng$B$gaps, gen_lng$D$gaps))
    ),
    "Mean Pairwise LD" = c(
      mean(genome_ld$A$pw, na.rm = TRUE),
      mean(genome_ld$B$pw, na.rm = TRUE),
      mean(genome_ld$D$pw, na.rm = TRUE),
      mean(c(genome_ld$A$pw, genome_ld$B$pw, genome_ld$D$pw), na.rm = TRUE)
    ),
    "Mean Neighbour LD" = c(
      mean(genome_ld$A$nbs, na.rm = TRUE),
      mean(genome_ld$B$nbs, na.rm = TRUE),
      mean(genome_ld$D$nbs, na.rm = TRUE),
      mean(c(genome_ld$A$nbs, genome_ld$B$nbs, genome_ld$D$nbs), na.rm = TRUE)
    )
  ) %>% round_df(., 2) %>% t()
  write_csv(
    map_stats, file.path(
      map_stats_and_plots, str_c(out_name, ".csv")
    )
  )
}