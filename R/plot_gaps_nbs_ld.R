#' Output figures based on map stats
#'
#' Outputs a figure of gap distances and neighbouring ld by genome and overall
#'
#' @param phys_lng data supplied by the calc_lng function for phys map
#' @param gen_lng data supplied by the calc_lng function for gen map
#' @param genome_ld data supplied by the calc_cld_stats function
#' @param plot_file_name the name for the output file
#' @param plot_title the title of the plot
#' @param y_lim the y-axis limit of the plot
#'
#' @importFrom dplyr tibble
#' @importFrom GGally ggmatrix
#' @importFrom ggplot2 aes element_text geom_freqpoly ggplot scale_colour_manual
#' @importFrom ggplot2 scale_x_log10 scale_y_log10 theme xlim ylim xlab ylab ggtitle
#' @importFrom magrittr %>%
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stringr str_c
#'
#' @return None

plot_gaps_nbs_ld <- function(
  phys_lng, gen_lng, genome_ld, out_name, plot_title, y_lim
) {
  # histograms and boxplots depicting the distribution of gaps on each genome
  gaps_ld <- tibble(
    genome = factor(
      c(rep("A", length(phys_lng$A$gaps)),
        rep("B", length(phys_lng$B$gaps)),
        rep("D", length(phys_lng$D$gaps)),
        rep("All", length(c(phys_lng$A$gaps, phys_lng$B$gaps, phys_lng$D$gaps)))
      ),
      levels = c("A", "B", "D", "All")
    ),
    phys_gaps = c(
      phys_lng$A$gaps,
      phys_lng$B$gaps,
      phys_lng$D$gaps,
      phys_lng$A$gaps, phys_lng$B$gaps, phys_lng$D$gaps
    ),
    gen_gaps = c(
      gen_lng$A$gaps,
      gen_lng$B$gaps,
      gen_lng$D$gaps,
      gen_lng$A$gaps, gen_lng$B$gaps, gen_lng$D$gaps
    ),
    ld = c(
      genome_ld$A$nbs,
      genome_ld$B$nbs,
      genome_ld$D$nbs,
      genome_ld$A$nbs, genome_ld$B$nbs, genome_ld$D$nbs
    )
  )

  # log_gaps <- cbind(gaps$gaps, log10(gaps$gaps))
  # print(log_gaps[is.infinite(log10(gaps$gaps)), ])
  plots <- list()
  plots[[1]] <- gaps_ld %>%
    ggplot() +
    geom_freqpoly(aes(phys_gaps * 1e6, colour = genome), size = 0.3) +
    scale_colour_manual(values = brewer.pal(4, "Dark2")) +
    # causes some values to be removed
    scale_x_log10(breaks = c(1, 1e2, 1e4, 1e6, 1e8), limits = c(NA, 1e8)) +
    scale_y_log10(limits = c(NA, y_lim))
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  plots[[2]] <- gaps_ld %>%
    ggplot() +
    geom_freqpoly(aes(gen_gaps, colour = genome), size = 0.3) +
    scale_colour_manual(values = brewer.pal(4, "Dark2")) +
    # causes some values to be removed
    scale_x_log10(breaks = c(0.01, 0.1, 1, 10), limits = c(NA, 30)) +
    scale_y_log10(limits = c(NA, y_lim))
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  plots[[3]] <- gaps_ld %>%
    ggplot() +
    geom_freqpoly(aes(ld, colour = genome), size = 0.3) +
    scale_colour_manual(values = brewer.pal(4, "Dark2")) +
    # causes some values to be removed
    xlim(0, 1) +
    scale_y_log10(limits = c(1, y_lim))
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  # turn plot list into ggmatrix
  plots_matrix <- ggmatrix(
    plots, nrow = 1, ncol = 3, yAxisLabels = "Num Markers",
    xAxisLabels = c(
      "Neighbour Distance in bp",
      "Neighbour Distance in cM",
      "Neighbour Abs. Comp. LD"
    ),
    # title = plot_title,
    legend = c(1, 2)
  )

  # plot the matrix
  png(
    file.path(
      map_stats_and_plots, str_c(out_name, ".png")
    ),
    family = "Times New Roman", width = 240, height = 120, pointsize = 10,
    units = "mm", res = 300
  )
  print(plots[[3]] +
    theme(
      legend.position = "bottom",
      text = element_text(size = 8, lineheight = 0.1)
    ) +
    xlab("Nieghbouring Marker LD") +
    ylab("Num Markers") +
    ggtitle("Hisograms of Neighbouring Marker LD by Genome")
  )
  # print(plots_matrix +
  #   theme(legend.position = "bottom",
  #     text = element_text(size = 8, lineheight = 0.1)
  #   )
  # )
  dev.off()
}