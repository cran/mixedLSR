#' Heatmap Plot of the mixedLSR Coefficient Matrices
#'
#' @param a A coefficient matrix from mixed_lsr model.
#' @param abs A boolean for taking the absolute value of the coefficient matrix.
#'
#' @importFrom ggplot2 ggplot aes geom_tile theme element_blank scale_fill_viridis_c facet_grid vars label_both
#'
#' @return A ggplot2 heatmap of the coefficient matrix, separated by subgroup.
#' @export
#'
#' @examples
#' simulate <- simulate_lsr()
#' plot_lsr(simulate$a)
plot_lsr <- function(a, abs = TRUE){
  a_plot <- NULL
  for(i in seq(1,length(a))){
    if(!is.null(a[[i]])){
      a_i <- a[[i]]
      a_i_p <- as.data.frame(as.table(a_i))
      a_i_p$Subgroup <- paste(i)
      if(abs){a_i_p$Freq <- abs(a_i_p$Freq)}
      a_plot <- rbind(a_plot,a_i_p)
    }
  }
  ggplot2::ggplot(a_plot, ggplot2::aes_string(x = "Var2", y = "Var1", fill = "Freq")) +
    ggplot2::geom_tile() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.text.x=ggplot2::element_blank(),
          axis.ticks.x=ggplot2::element_blank(),
          axis.text.y=ggplot2::element_blank(),
          axis.ticks.y=ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          legend.position = "top") +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::scale_fill_viridis_c(name="Coefficient") +
    ggplot2::facet_grid(rows = "Subgroup", labeller = ggplot2::label_both)
}

