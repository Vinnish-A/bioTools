
# plot --------------------------------------------------------------------


library(ggrepel)


## volcano ----

part_non = function(data_, x_, y_, label_, color_ = 'black', alpha_ = 0.75) {

  data_point_ = data_

  x_ = sym(x_)
  y_ = sym(y_)

  list(
    geom_point(data = data_point_, aes(!!x_, !!y_), color = color_, alpha = alpha_)
  )

}

part_up = function(data_, x_, y_, label_, color_ = '#cc0303', max_ = 12, alpha_ = 0.75) {

  data_point_ = data_
  data_label_ = data_[!is.na(data_[[label_]]), ]

  nudge_x_ = max_ - data_label_[[x_]]
  x_ = sym(x_)
  y_ = sym(y_)

  list(
    geom_point(data = data_point_, aes(!!x_, !!y_), color = color_, alpha = alpha_),
    geom_text_repel(
      data = data_label_, aes(!!x_, !!y_, label = !!label_),
      nudge_x = nudge_x_, color = color_,
      max.overlaps = Inf, hjust = 1, vjust = 0,
      direction = 'y',
      segment.linetype = 2, segment.size = 0.1, segment.alpha = 0.5,
      min.segment.length = 0
    )
  )

}

part_down = function(data_, x_, y_, label_, color_ = '#026401', min_ = -12, alpha_ = 0.75) {

  data_point_ = data_
  data_label_ = data_[!is.na(data_[[label_]]), ]

  nudge_x_ = min_ - data_label_[[x_]]
  x_ = sym(x_)
  y_ = sym(y_)
  label_ = sym(label_)

  list(
    geom_point(data = data_point_, aes(!!x_, !!y_), color = color_, alpha = alpha_),
    geom_text_repel(
      data = data_label_, aes(!!x_, !!y_, label = !!label_),
      nudge_x = nudge_x_, color = color_,
      max.overlaps = Inf, hjust = 1, vjust = 0,
      direction = 'y',
      segment.linetype = 2, segment.size = 0.1, segment.alpha = 0.5,
      min.segment.length = 0
    )
  )

}


#' Create a Volcano Plot
#'
#' Combines up-, down-, and non-significant layers into a complete volcano plot.
#'
#' @param df_ A data frame containing `x_`, `y_`, `label_`, and a grouping column (e.g., `colorBy_`).
#' @param x_ Column name for the x-axis variable. Default is "logFC".
#' @param y_ Column name for the y-axis variable. Default is "log10padj".
#' @param label_ Column name for gene labels. Default is "label".
#' @param colorBy_ Column name used for splitting data into "up", "down", and "non" groups. Default is "direction".
#' @param colors_ Named vector of colors for "up", "down", and "non" categories.
#' @param alpha_ Point transparency. Default is 0.75.
#' @param range_ A numeric vector of length 2 specifying the x-axis range. Default is c(-12, 12).
#'
#' @return A ggplot object representing the volcano plot.
#'
#' @export
plot_volcano = function(
    df_, x_ = 'logFC', y_ = 'log10padj', label_ = 'label',
    colorBy_ = 'direction', colors_ = color_myth,
    alpha_ = 0.75, range_ = c(-12, 12)
) {

  lst_dp_ = split(df_, df_[[colorBy_]])

  ggplot() +
    part_non(lst_dp_$non, x_, y_, label_, colors_[[2]], alpha_) +
    part_down(lst_dp_$down, x_, y_, label_, colors_[[1]], range_[[1]], alpha_) +
    part_up(lst_dp_$up, x_, y_, label_, colors_[[3]], range_[[2]], alpha_) +
    xlab('log2(Fold Change)') +
    ylab('-log1o(Adj Pvalue)') +
    theme_bw()

}

## kegg ----

#' Visualize KEGG Pathway with Differential Gene Expression Data
#'
#' This function overlays gene expression data onto a KEGG pathway graph.
#'
#' @param deg_ A data frame with gene expression results (e.g., log fold change and adjusted p-values).
#' @param id_ The KEGG pathway ID (e.g., "hsa00010").
#' @param logFC_ Column name for log fold change in `deg_`. Default is "log2FoldChange".
#' @param padj_ Column name for adjusted p-values in `deg_`. Default is "padj".
#' @param color_ A vector of colors for the gradient scale. Default is `color_myth`.
#'
#' @return A ggraph object visualizing the KEGG pathway.
#'
#' @export
plot_kegg = function(deg_, id_, logFC_ = 'log2FoldChange', padj_ = 'padj', color_ = color_myth) {

  library(ggraph)
  library(ggkegg)

  graph_ = pathway(id_, 'tmp') |>
    activate(nodes) |>
    mutate(converted_name=convert_id(substr(id_, 1, 3))) |>
    left_join(deg_, by = c('converted_name' = 'symbol'))

  ggraph(graph_, layout="manual", x = x, y = y)+
    geom_node_rect(aes(fill = !!sym(logFC_), filter = !is.na(!!sym(padj_)) & !!sym(padj_)<0.05))+
    # ggfx::with_outer_glow(geom_node_rect(aes(fill=log2FoldChange, filter=!is.na(padj) & padj<0.05)), colour="yellow", expand=2)+
    overlay_raw_map(id_, transparent_colors = c("#cccccc","#FFFFFF","#BFBFFF","#BFFFBF"))+
    scale_fill_gradient2(low=color_[[1]], mid=color_[[2]], high=color_[[3]]) +
    theme_void()

}

## heatmap ----

#' Plot Heatmap with Gene Set and Class Annotations
#'
#' Generates a heatmap for selected genes with hierarchical clustering and class labels.
#'
#' @param df_ A data frame containing expression data with `symbol` as row names.
#' @param table_class_ A named vector mapping sample names to class labels.
#' @param geneset_ A named list of gene sets to visualize.
#'
#' @return A ggplot object representing the heatmap.
#'
#' @import ggplot2
#' @import ggh4x
#' @importFrom ggtree ggtree
#'
#' @export
plot_heatmap = function(df_, table_class_, geneset_) {

  name_ = names(geneset_)[[1]]
  geneset_ = geneset_[[name_]]

  dp_pre_ = df_ |>
    filter(symbol %in% geneset_) |>
    column_to_rownames('symbol') |>
    t() |>
    as_tibble(rownames = 'sample') |>
    mutate(across(where(is.numeric), ~ as.numeric(scale(log2(.x + 1))))) |>
    column_to_rownames('sample') |>
    t() |>
    as_tibble(rownames = 'symbol') |>
    drop_na()

  pic_tree_ = dp_pre_ |>
    column_to_rownames('symbol') |>
    dist() |>
    hclust() |>
    ggtree(color = 'white')

  order_ = pic_tree_$data |>
    filter(isTip) |>
    arrange(y) |>
    pull(label)

  pic_heatmap_ = dp_pre_|>
    pivot_longer(-symbol, names_to = 'sample', values_to = 'expr') |>
    mutate(sample = factor(sample, levels = names(table_class_)),
           class = table_class_[sample],
           symbol = factor(symbol, levels = order_)) |>
    ggplot(aes(interaction(sample, class), symbol, fill = expr)) +
    # ggplot(aes(sample, symbol, fill = expr)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("green", "black", "red"),
      values = scales::rescale(c(-2, 0, 2)),
      limits = c(-2, 2),
      oob = scales::squish
    ) +
    xlab(NULL) +
    ylab(NULL) +
    labs(title = name_) +
    guides(x = "axis_nested") +
    theme_classic() +
    theme(
      axis.text.x = element_text(color = 'white'),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      ggh4x.axis.nesttext.x = element_text(color = "black"),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )

  return(pic_heatmap_)

}

#' Plot Complex Heatmap with Class Splits
#'
#' Uses the ComplexHeatmap package to generate a detailed heatmap with class splits.
#'
#' @param df_ A data frame containing expression data.
#' @param table_class_ A named vector of class annotations.
#' @param geneset_ A named list of gene sets.
#'
#' @return A Heatmap object from the ComplexHeatmap package.
#'
#' @export
plot_chtmap = function(df_, table_class_, geneset_) {

  library(ComplexHeatmap)
  library(circlize)

  title_ = names(geneset_)[[1]]
  geneset_ = geneset_[[title_]]

  dp_pre_ = df_ |>
    select(symbol, all_of(names(table_class_))) |>
    filter(symbol %in% geneset_) |>
    tTidy('symbol', 'sample') |>
    mutate(across(where(is.numeric), ~ as.numeric(scale(log2(.x + 1))))) |>
    tTidy('sample', 'symbol') |>
    drop_na() |>
    column_to_rownames('symbol') |>
    as.matrix()

  pic_tree_ = dp_pre_ |>
    dist() |>
    hclust() |>
    ggtree(color = 'white')

  order_ = pic_tree_$data |>
    filter(isTip) |>
    arrange(y) |>
    pull(label)

  dp_ = dp_pre_[order_, ]

  col_fun = colorRamp2(c(-1, 0, 1), c("green", "black", "red"))

  res_ = Heatmap(
    dp_,
    name = "expr",
    col = col_fun,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    row_names_side = "left",
    column_split = table_class_[colnames(dp_pre_)]
  )

  return(res_)

}

## gsea ----

#' Plot Gene Set Enrichment Analysis (GSEA) Results
#'
#' Generates a GSEA plot for a specific term.
#'
#' @param gsea_ A GSEA result object from the `clusterProfiler` package.
#' @param term_ The term to visualize.
#' @param sizef_ A scaling factor for the line thickness. Default is 2.
#' @param color_ A vector of colors for the gradient scale.
#'
#' @return A ggplot object representing the GSEA plot.
#'
#' @export
plot_gsea = function(gsea_, term_, sizef_ = 2, color_ = color_monika) {

  geneSetID_ = which(gsea_@result$ID == term_)

  ES_ = gsea_@result$enrichmentScore[geneSetID_]
  FDR_ = gsea_@result$pvalue[geneSetID_]
  tag_ = sprintf("ES=%.3f, False Discovery Rate=%.3f", ES_, FDR_)

  gsdata_ = enrichplot:::gsInfo(gsea_, geneSetID_) |> as_tibble()

  p1_ = gsdata_ |>
    ggplot(aes(x = x)) +
    geom_segment(aes(xend = x, y = 0, yend = -runningScore, color = x), linewidth = 0.1 * sizef_, data = subset(gsdata_, position == 1)) +
    theme_void() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0.1), breaks = NULL) +
    scale_color_gradient(low = color_[[1]], high = color_[[2]]) +
    xlab(term_) +
    ylab(NULL) +
    theme(
      legend.position = "none",
      plot.margin = margin(t = -0.1, b = 0, unit = "cm"),
      axis.title.x = element_text(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.line.x = element_blank()
    ) +
    annotation_custom(
      grob = textGrob(
        tag_,
        x = unit(1, "npc"), y = unit(0.05, "npc"),
        hjust = 1, vjust = 0,
        gp = gpar(col = "black", fontsize = 10, family = "italic", fontface = "italic")
      )
    )

  p1_

}
