
# plot --------------------------------------------------------------------

#' Register a Function Across Multiple Environments
#'
#' @export
register = function(fun_, ...) {

  name_full_ = as.character(substitute(fun_))
  name_last_ = strsplit(name_full_, '_')[[1]]
  name_last_ = name_last_[[length(name_last_)]]

  fun_ = get(name_full_, envir = parent.frame())
  envs_ = list(...)

  for (i_ in seq_along(envs_)) assign(name_last_, fun_, envir = envs_[[i_]])

}

#' theme_dropx
#'
#' @export
theme_dropx = function() {

  theme(
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

}

#' theme_dropy
#'
#' @export
theme_dropy = function() {

  theme(
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
  )

}

## universal ----

plot_stacked = function(df_, group_, subgroup_, n_ = 'n', thres_ = 0.1, attention_ = NULL, summarise_ = F) {

  if (summarise_) {

    dp_bar_ = tibble(group = df_[[group_]],
                     subgroup = df_[[subgroup_]]) |>
      group_by(group, subgroup) |>
      summarise(n = n())

  } else {

    dp_bar_ = tibble(group = df_[[group_]],
                     subgroup = df_[[subgroup_]],
                     n = df_[[n_]])

  }

  dp_bar_ = dp_bar_ |>
    arrange(group, subgroup) |>
    group_by(group) |>
    mutate(ratio = signif(n/sum(n), 3),
           label = ifelse(ratio > thres_, paste0(ratio * 100, "%"), NA),
           postion_top = rev(cumsum(rev(ratio))),
           postion_bottom = postion_top - ratio)

  if (!is.null(attention_)) {
    dp_bar_ = dp_bar_ |>
      mutate(label = ifelse(subgroup %in% attention_, label, NA))
  }

  nGroup_ = length(unique(dp_bar_[['subgroup']]))
  colors_ = get_color(nGroup_)

  dp_bar_ |>
    ggplot() +
    geom_bar(aes(group, n, fill = subgroup), color = "#f3f4f4", position = "fill", stat = "identity", size = 1) +
    geom_text(aes(group, (postion_top + postion_bottom)/2, label = label), size = 4, color = "white", fontface = 'bold') +
    scale_y_continuous(labels = paste0(100 * seq(0, 1, 0.25), "%")) +
    scale_fill_manual(values = colors_) +
    xlab(NULL) +
    ylab(NULL) +
    labs(fill = NULL) +
    theme_classic() +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.75)
    )

}

## volcano ----

#' part_non
#'
#' @keywords internal
part_non = function(data_, x_, y_, label_, color_ = 'black', alpha_ = 0.75) {

  data_point_ = data_

  x_ = sym(x_)
  y_ = sym(y_)

  list(
    geom_point(data = data_point_, aes(!!x_, !!y_), color = color_, alpha = alpha_)
  )

}

#' part_up
#'
#' @keywords internal
part_up = function(data_, x_, y_, label_, color_ = '#cc0303', max_ = 12, alpha_ = 0.75) {

  data_point_ = data_
  data_label_ = data_[!is.na(data_[[label_]]), ]

  nudge_x_ = max_ - data_label_[[x_]]
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

#' part_down
#'
#' @keywords internal
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
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' df = data.frame(
#'  logFC = c(-3, 1.5, 2, -4, 0),
#'  log10padj = c(1, 2, 3, 4, 0.5),
#'  label = c("Gene1", "Gene2", NA, "Gene4", NA),
#'  direction = c("down", "up", "up", "down", "non")
#' )
#' plot_volcano(df)
#'
#' @export
plot_volcano = function(
    df_, x_ = 'logFC', y_ = 'log10padj', label_ = 'label',
    colorBy_ = 'direction', colors_ = get_color(3, 'htmap'),
    alpha_ = 0.75, range_ = c(-12, 12)
) {

  lst_dp_ = split(df_, df_[[colorBy_]])

  ggplot() +
    part_non(lst_dp_$non, x_, y_, label_, colors_[[2]], alpha_) +
    part_down(lst_dp_$down, x_, y_, label_, colors_[[1]], range_[[1]], alpha_) +
    part_up(lst_dp_$up, x_, y_, label_, colors_[[3]], range_[[2]], alpha_) +
    xlab('log2(Fold Change)') +
    ylab('-log10(Adj Pvalue)') +
    theme_bw()

}

## enrich ----

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
plot_gsea = function(gsea_, term_, sizef_ = 2, color_ = get_color(2, 'ddlc')) {

  checkReliance('enrichplot')

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

#' Visualize KEGG Pathway with Differential Gene Expression Data
#'
#' This function overlays gene expression data onto a KEGG pathway graph.
#'
#' @param deg_ A data frame with gene expression results (e.g., log fold change and adjusted p-values).
#' @param id_ The KEGG pathway ID (e.g., "hsa00010").
#' @param logFC_ Column name for log fold change in `deg_`. Default is "log2FoldChange".
#' @param padj_ Column name for adjusted p-values in `deg_`. Default is "padj".
#' @param color_ A vector of colors for the gradient scale. Default is `get_color(3, 'htmap')`.
#'
#' @return A ggraph object visualizing the KEGG pathway.
#'
#' @export
plot_pathway = function(deg_, id_, logFC_ = 'log2FoldChange', padj_ = 'padj', color_ = get_color(3, 'htmap')) {

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

#' Plot GO Terms
#'
#' This function creates a bar plot for the top GO terms based on their -log10(q-value).
#' It supports grouping the data by a specified feature, such as the ontology type.
#'
#' @param df_ A dataframe, result attr of a GO obj.
#' @param title_ A character string specifying the title of the plot.
#' @param nTerm_ An integer specifying the number of top terms to display (default is 5).
#' @param colorBy_ A character string indicating which feature to group and color the bars by. Defaults to 'ONTOLOGY'.
#' @return A `ggplot` object displaying the bar plot for the top GO terms, grouped and colored by the specified feature.
#'
#' @examples
#' plot_GO(res_GO@result, 'Test')
#'
#' @export
plot_GO = function(df_, title_, nTerm_ = 5, colorBy_ = 'ONTOLOGY') {

  facet_ = as.formula(paste('~', colorBy_))

  colorBy_ = sym(colorBy_)

  dp_ = df_ |>
    group_by(!!colorBy_) |>
    slice_min(qvalue, n = nTerm_, with_ties = F) |>
    mutate(Description = factor(Description, levels = rev(Description))) |>
    ungroup()

  dp_ |>
    ggplot() +
    geom_col(aes(-log10(qvalue), Description, fill = !!colorBy_)) +
    geom_text(aes(0, Description, label = capitalize(Description)), hjust = 0, nudge_x = 0.1, size = 4) +
    geom_vline(xintercept = -log10(0.05), linetype = 2, alpha = 0.5) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = get_color(4, 'morandi')) +
    xlab('-Log10(FDR)') +
    ylab(NULL) +
    labs(title = paste(title_, '| GO Terms')) +
    facet_wrap(facet_, ncol = 1, scales = 'free_y', strip.position = 'left') +
    theme_classic() +
    theme(legend.position = 'none',
          axis.ticks.y = element_blank(),
          panel.grid.major = element_line(color = "#fafafa"),
          panel.grid.minor = element_line(color = "#fafafa"),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0.5))

}

#' Plot KEGG Terms
#'
#' This function creates a bar plot for the top KEGG pathway terms
#' based on their -log10(q-value) for a given KEGG result object.
#'
#' @param df_ A dataframe, result attr of a KEGG obj.
#' @param title_ A character string specifying the title of the plot.
#' @param nTerm_ An integer specifying the number of top terms to display (default is 5).
#' @return A `ggplot` object displaying the bar plot for the top KEGG pathway terms.
#'
#' @examples
#' plot_KEGG(res_KEGG@result, 'Test')
#'
#' @export
plot_KEGG = function(df_, title_, nTerm_ = 5) {

  dp_ = df_ |>
    slice_min(qvalue, n = nTerm_, with_ties = F) |>
    mutate(ONTOLOGY = 'KEGG',
           Description = factor(Description, levels = rev(Description)))

  dp_ |>
    ggplot() +
    geom_col(aes(-log10(qvalue), Description, fill = ONTOLOGY)) +
    geom_text(aes(0, Description, label = capitalize(Description)), hjust = 0, nudge_x = 0.1, size = 4) +
    geom_vline(xintercept = -log10(0.05), linetype = 2, alpha = 0.5) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = get_color(4, 'morandi')) +
    xlab('-Log10(FDR)') +
    ylab(NULL) +
    labs(title = paste(title_, '| KEGG Terms')) +
    theme_classic() +
    theme(legend.position = 'none',
          axis.ticks.y = element_blank(),
          panel.grid.major = element_line(color = "#fafafa"),
          panel.grid.minor = element_line(color = "#fafafa"),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0.5))

}

#' Plot Gene Ontology (GO) and KEGG Pathway Terms Together
#'
#' This function combines the GO and KEGG results and creates a bar plot for the top terms based on their -log10(q-value).
#'
#' @param go_ A GO result object containing the GO data.
#' @param kegg_ A KEGG result object containing the KEGG pathway data.
#' @param title_ A character string specifying the title of the plot.
#' @param nTerm_ An integer specifying the number of top terms to display (default is 5).
#' @param avoid_ A character string to exclude specific terms from the plot. Default is an empty string, meaning no terms are excluded.
#' @return A `ggplot` object displaying the combined bar plot for the top GO and KEGG terms.
#'
#' @examples
#' plot_GOKEGG(res_GO, res_KEGG, 'Test')
#'
#' @export
plot_GOKEGG = function(go_, kegg_, title_, nTerm_ = 5, avoid_ = '') {

  dp_ = bind_rows(
    go_@result |>
      select(Description, qvalue, ONTOLOGY),
    kegg_@result |>
      mutate(ONTOLOGY = 'KEGG') |>
      select(Description, qvalue, ONTOLOGY)
  )

  if (avoid_ != '') {

    dp_ = dp_ |>
      filter(!map(avoid_, ~ str_detect(Description, .x)) |> pmap(`|`) |> unlist())

  }

  dp_ = dp_ |>
    group_by(ONTOLOGY) |>
    slice_min(qvalue, n = nTerm_, with_ties = F) |>
    mutate(ONTOLOGY = factor(ONTOLOGY, levels = c('BP', 'CC', 'MF', 'KEGG')),
           Description = factor(Description, levels = rev(Description))) |>
    ungroup()

  dp_ |>
    ggplot() +
    geom_col(aes(-log10(qvalue), Description, fill = ONTOLOGY)) +
    geom_vline(xintercept = -log10(0.05), linetype = 2, alpha = 0.5) +
    geom_text(aes(0, Description, label = capitalize(Description)), hjust = 0, nudge_x = 0.1, size = 4) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = get_color(4, 'morandi')) +
    xlab('-Log10(FDR)') +
    ylab(NULL) +
    labs(title = title_) +
    facet_wrap(~ ONTOLOGY, ncol = 1, scales = 'free_y', strip.position = 'left') +
    theme_classic() +
    theme(legend.position = 'none',
          axis.ticks.y = element_blank(),
          panel.grid.major = element_line(color = "#fafafa"),
          panel.grid.minor = element_line(color = "#fafafa"),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0.5))

}

### register ----

plot_enrich = new.env()

register(plot_gsea, plot_enrich)
register(plot_pathway, plot_enrich)
register(plot_GO, plot_enrich)
register(plot_KEGG, plot_enrich)
register(plot_GOKEGG, plot_enrich)

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
plot_heatmap = function(df_, table_class_, geneset_, scale_ = T, color_ = get_color(3), range_ = NULL) {

  checkReliance('ggtree')

  if (scale_) {
    handleRow_ = \(vec__) as.numeric(scale(vec__))
  } else {
    handleRow_ = \(vec__) as.numeric(vec__)
  }

  name_ = names(geneset_)[[1]]
  geneset_ = geneset_[[name_]]

  dp_pre_ = df_ |>
    filter(symbol %in% geneset_) |>
    tidyT('sample', 'symbol') |>
    mutate(across(where(is.numeric), ~ as.numeric(scale(log2(.x + 1))))) |>
    tidyT('symbol', 'sample') |>
    drop_na()

  pic_tree_ = dp_pre_ |>
    column_to_rownames('symbol') |>
    dist() |>
    hclust() |>
    ggtree::ggtree(color = 'white')

  order_ = pic_tree_$data |>
    filter(isTip) |>
    arrange(y) |>
    pull(label)

  dp_ = dp_pre_|>
    pivot_longer(-symbol, names_to = 'sample', values_to = 'expr') |>
    mutate(sample = factor(sample, levels = names(table_class_)),
           class = table_class_[sample],
           symbol = factor(symbol, levels = order_))

  limits = range(range_) %||% dp_[[expr]]

  p_ = dp_ |>
    ggplot(aes(interaction(sample, class), symbol, fill = expr)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = color_,
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

  return(p_)

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

  checkReliance('ComplexHeatmap', 'circlize')

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

  col_fun = colorRamp2(c(-1, 0, 1), get_color(3))

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

## sc ----

plot_sc = new.env()

#' Plot UMAP visualization
#'
#' Creates a UMAP plot with points colored by the specified grouping variable.
#'
#' @param df_ Data frame. Contains UMAP coordinates and grouping information.
#' @param x_ Character. The column name for the UMAP x-axis.
#' @param y_ Character. The column name for the UMAP y-axis.
#' @param colorBy_ Character. The column name for the grouping variable used for coloring points.
#' @param label_ Logic. Whether to add label.
#' @param border_ Logic. Whether to add border.
#' @param discrete_ Logic. Whether to use discrete colors.
#' @param arrowSize_ Numeric. Relative size of arrow.
#'
#' @return A ggplot object displaying the UMAP visualization.
#'
#' @importFrom grid arrow
#'
#' @examples
#' plot_umap(df_, x_ = "UMAP_1", y_ = "UMAP_2", colorBy_ = "CellType")
#'
#' @export
plot_umap = function(df_, x_ = 'umap_1', y_ = 'umap_2', colorBy_ = 'seurat_clusters', label_ = T, border_ = T, discrete_ = T, arrowSize_ = 0.5) {

  params_chr_ = c(x_, y_, colorBy_)

  x_ = sym(x_)
  y_ = sym(y_)
  colorBy_ = sym(colorBy_)

  axis_ = guide_axis_truncated(
    trunc_lower = unit(0, "npc"),
    trunc_upper = unit(8 * arrowSize_, "cm")
  )

  p_ = df_ |>
    ggplot()

  if (border_) {
    p_ = p_ +
      ggunchull::stat_unchull(aes(!!x_, !!y_, color = !!colorBy_, fill = !!colorBy_), alpha = 0.2, size = 1, lty = 2, qval = 0.8, delta = 1)
  }

  p_ = p_ |>
    geom_point(aes(!!x_, !!y_, color = !!colorBy_), size = 1) +
    guides(fill = "none", x = axis_, y = axis_) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    labs(color = "", title = paste("UMAP by", colorBy_)) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.key = element_blank(),
      aspect.ratio = 1,
      legend.position = "bottom",
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_text(hjust = 0.05, face = "italic"),
      axis.line.x = element_line(arrow = arrow(type = "open", length = unit(arrowSize_, "cm"))),
      axis.line.y = element_line(arrow = arrow(type = "open", length = unit(arrowSize_, "cm")))
    )


  if (label_) {

    dp_label_ = df_ |>
      group_by(!!colorBy_) |>
      summarise(x = mean(!!x_), y = mean(!!y_))

    p_ = p_ +
      geom_shadowtext(aes(x, y, label = !!colorBy_), data = dp_label_, color = 'black', bg.color = 'white', bg.r = 0.2)

  }

  if (discrete_) {

    nColor_ = length(unique(df_[[params_chr_[[3]]]]))
    color_ = get_color(nColor_)

    p_ = p_ +
      scale_fill_manual(values = color_) +
      scale_color_manual(values = color_)

  }

  return(p_)

}

#' Plot Single-Cell Differential Expression (DEG) Results with GO Enrichment
#'
#' This function creates a combined visualization of single-cell differential expression (DEG) results.
#' It includes a heatmap of log2 fold changes, a dot plot of cluster annotations, and a bar plot for Gene Ontology (GO) enrichment results.
#'
#' @param seu_ A Seurat object containing the single-cell RNA-seq data.
#' @param marker_all_ A data frame with all marker gene data containing columns `gene`, `cluster`, and `avg_log2FC`.
#' @param num_each_ An integer specifying the number of genes to display per cluster (default is 15).
#' @param ident_ A vector of cluster identities to use in the analysis (default is the cluster identities from the Seurat object).
#' @param breaks_ A numeric vector specifying the breaks for the color scale in the heatmap (default is c(-5, 0, 5)).
#' @param organism_ A character string specifying the organism for GO enrichment analysis. Options are `'hsa'` for human (default) or `'mmu'` for mouse.
#' @return A combined plot consisting of:
#'   - A heatmap showing the average log2 fold change for genes across clusters.
#'   - A dot plot showing the cluster annotations.
#'   - A bar plot showing the top 5 GO terms enriched in each cluster.
#'
#' @examples
#' plot_sc_deg(seurat_object, marker_data, num_each_ = 20, ident_ = c('Cluster1', 'Cluster2'), organism_ = 'hsa')
#'
#' @export
plot_sc_deg = function(seu_, marker_all_, num_each_ = 15, ident_ = Idents(seu_), breaks_ = c(-5, 0, 5), organism_ = 'hsa') {

  ident_ = levels(Idents(seu_))

  deg_sliced_ = marker_all_ |>
    distinct(gene, .keep_all = T) |>
    filter(avg_log2FC > 0) |>
    group_by(cluster) |>
    slice_max(avg_log2FC, n = 6*num_each_) |>
    slice_sample(n = num_each_) |>
    ungroup() |>
    arrange(cluster) |>
    pull(gene, cluster)

  dp_logfc_ = map(
    ident_, \(cluster_) {
      marker_each_ = FindMarkers(seu_, ident.1 = cluster_, only.pos = F, features = deg_sliced_)
      marker_each_$cluster = cluster_
      return(as_tibble(marker_each_, rownames = "gene"))
    }
  ) |> bind_rows() |>
    arrange(cluster) |>
    mutate(cluster = factor(cluster),
           gene = factor(gene, levels = deg_sliced_)) |>
    select(cluster, avg_log2FC, gene)

  patch_ = expand.grid(cluster = ident_, gene = deg_sliced_) |>
    as_tibble() |>
    mutate(V = paste(cluster, gene, sep = '_')) |>
    filter(!V %in% with(dp_logfc_, paste(cluster, gene, sep = '_'))) |>
    select(-V) |>
    mutate(avg_log2FC = runif(length(row_number()), -0.75, 0.25))

  dp_logfc_ = bind_rows(dp_logfc_, patch_)

  if (organism_ == 'hsa') {
    enrich_ = compareCluster(gene ~ cluster, data = marker_all_, fun = 'enrichGO', OrgDb = 'org.Hs.eg.db', keyType = 'SYMBOL')
  } else if (organism_ == 'mmu') {
    enrich_ = compareCluster(gene ~ cluster, data = marker_all_, fun = 'enrichGO', OrgDb = 'org.Mm.eg.db', keyType = 'SYMBOL')
  }

  dp_enrich_ = enrich_@compareClusterResult |>
    filter(cluster != "NA") |>
    group_by(cluster) |>
    slice_max(-p.adjust, n = 5, with_ties = F) |>
    ungroup() |>
    mutate(cluster = factor(cluster, levels = rev(ident_)))

  limits_ = min(abs(range(dp_logfc_$avg_log2FC)))
  limits_ = c(-limits_, limits_)

  p1_ = dp_logfc_ |>
    ggplot(aes(x = cluster,
               y = reorder(gene, -as.numeric(cluster)),
               fill = avg_log2FC,
               color = avg_log2FC)) +
    geom_tile() +
    xlab("") +
    ylab("") +
    labs(fill = "", color = "") +
    scale_color_gradient2(low = "#01665e", mid = "white", high = "#8c510a", breaks = breaks_, limits = limits_, oob = scales::squish) +
    scale_fill_gradient2(low = "#01665e", mid = "white", high = "#8c510a", breaks = breaks_, limits = limits_, oob = scales::squish) +
    scale_x_discrete(breaks = NULL) +
    theme_bw() +
    theme(legend.position = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", linewidth = 0.2, fill = NA),
          plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
          legend.key.size = unit(0.5, "cm"),
          legend.key.width = unit(0.5, "cm"))

  ncolor_ = length(ident_)
  colors_ = get_color(ncolor_)

  p2_ = tibble(x = 0, y = ident_) |>
    ggplot(aes(x = y, y = 1, color = factor(y))) +
    geom_point(size = 6, show.legend = F) +
    scale_color_manual(values = colors_) +
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.position = "none",
          panel.spacing = unit(0, "lines"),
          panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "pt"),
          axis.text.x = element_text(angle = 30,
                                     size = 12,
                                     hjust = 1,
                                     vjust = 1,
                                     color = "black"),
          axis.title  = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          axis.text.y = element_blank())

  limits_ = c(0, mean(sort(-log10(dp_enrich_$p.adjust), T)[1:6]))

  p3_ = dp_enrich_ |>
    ggplot(aes(x = reorder(Description, -log10(p.adjust)),
               y = -log10(p.adjust),
               fill = cluster)) +
    geom_bar(position = position_dodge(),
             stat = "identity",
             show.legend = F) +
    scale_x_discrete(position = "top") +
    scale_y_continuous(expand = c(0,0), limits = limits_, oob = scales::squish, labels = scales::number_format(accuracy = 1)) +
    scale_fill_manual(values = rev(colors_))+
    facet_wrap(~ cluster, ncol = 1, scales = "free_y") +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title.y = element_blank(),
          panel.spacing = unit(0, "lines"),
          axis.text = element_text(size = 12),
          panel.border = element_rect(colour = "grey50")) +
    ylab("-Log10(Padj)") +
    coord_flip()

  p1_ + p3_ + p2_ + plot_spacer() + plot_layout(widths = c(2, 1), height = c(12, 1))

}

### register ----

plot_sc = new.env()

register(plot_umap, plot_sc)
register(plot_sc_deg, plot_sc)

