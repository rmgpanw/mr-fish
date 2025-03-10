
# Figure 3 ---------------------------------------------------------------

plot_figure_3 <- function(all_mr_fish_results_df, file_path = "Figure 3.png") {
  plot_replication_forest_plot_figure(
    all_mr_fish_results_df = all_mr_fish_results_df,
    selected_outcomes = c("PBC", "Any ischaemic stroke", "CD", "FEV1"),
    file_path = file_path
  )
}

# Figure 4 ----------------------------------------------------------------

plot_figure_4 <- function(all_mr_fish_results_df, file_path = "Figure 4.png") {
  plot_replication_forest_plot_figure(
    all_mr_fish_results_df = all_mr_fish_results_df,
    selected_outcomes = c("Alzheimer's", "CHD", "T2DM", "AMD"),
    file_path = file_path
  )
}

plot_replication_forest_plot_figure <- function(all_mr_fish_results_df, selected_outcomes, file_path) {
  
  # validate args
  stopifnot(identical(length(selected_outcomes), 4L))
  stopifnot(all(selected_outcomes %in% all_mr_fish_results_df$Outcome))
  stopifnot(rlang::is_string(file_path))
  
  # make list of replication forest plots
  forest_plots <- selected_outcomes %>%
    purrr::set_names() %>%
    purrr::map(
      \(x) plot_replication_forest_plot_single(
        all_mr_fish_results_df,
        OUTCOME = x,
        TEXT_SIZE = 8,
        LEGEND_TEXT_SIZE = 12,
        X_TITLE_TEXT_SIZE = 12,
        PVALUE_THRESHOLD = 0.0048
      )
    )
  
  # manualy adjust axes
  if ("T2DM" %in% names(forest_plots)) {
    forest_plots$T2DM <- forest_plots$T2DM + ggplot2::coord_cartesian(xlim = c(1e-03, 2e06))
  }
  
  if ("CD" %in% names(forest_plots)) {
    forest_plots$CD <- forest_plots$CD + ggplot2::coord_cartesian(xlim = c(1e-07, 1.1e04))
  }
  
  if ("FEV1" %in% names(forest_plots)) {
    forest_plots$FEV1 <- forest_plots$FEV1 + ggplot2::coord_cartesian(xlim = c(-2, 1))
  }
  
  # remove scientific notation from axes (see
  # https://ggplot2.tidyverse.org/articles/articles/faq-axes.html and
  # http://www.sthda.com/english/wiki/ggplot2-axis-scales-and-transformations#log-and-sqrt-transformations
  # )
  
  outcome_type <- all_mr_fish_results_df |>
    dplyr::mutate(binary_outcome = dplyr::case_when(`Outcome units` == "logOR" ~ "binary",
                                                    TRUE ~ "continuous")) |>
    dplyr::select(binary_outcome,
                  Outcome) |>
    dplyr::group_by(binary_outcome) |>
    tidyr::nest() |>
    tibble::deframe() |>
    purrr::map(\(x) x$Outcome)
  
  forest_plots <- forest_plots %>%
    purrr::imap( ~ {
      if (.y %in% outcome_type$binary) {
        .x +
          ggplot2::scale_x_continuous(
            labels = scales::label_number(drop0trailing = TRUE),
            trans = 'log',
            breaks = scales::log_breaks(n = 7)
          )
      }
      else {
        .x
      }
    })
  
  # extract legend and move to bottom (See http://www.sthda.com/english/wiki/ggplot2-legend-easy-steps-to-change-the-position-and-the-appearance-of-a-graph-legend-in-r-software#change-the-legend-position-for-multiple-guides)
  p_legend <- cowplot::get_legend(forest_plots[[1]] +
                                    ggplot2::theme(legend.position = "bottom", legend.box = "horizontal"))
  
  forest_plots <- forest_plots %>%
    purrr::map( ~ .x + ggplot2::theme(legend.position = 'None'))
  
  # plot - egg
  egg_main_plot <- egg::ggarrange(
    plots = forest_plots,
    byrow = TRUE,
    nrow = 2, ncol = 2
  )
  
  # add legend
  cowplot_plot <-
    cowplot::plot_grid(cowplot::plot_grid(egg_main_plot),
                       p_legend,
                       nrow = 2,
                       ncol = 1,
                       rel_heights = c(1, 0.05)) +
    # without this, background is black when saved as .png
    # theme(plot.background = element_rect(fill = "white"))
    cowplot::theme_nothing() + ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white",
                                                                                      colour = NA))
  
  ggplot2::ggsave(filename = file_path,
                  plot = cowplot_plot,
                  width = 32,
                  height = 48,
                  units = "cm",
                  dpi = 300)  
  
  return(file_path)
}

plot_replication_forest_plot_single <- function(all_mr_fish_results_df,
                                                OUTCOME,
                                                TEXT_SIZE = 8,
                                                LEGEND_TEXT_SIZE = 12,
                                                X_TITLE_TEXT_SIZE = 12,
                                                PVALUE_THRESHOLD = 0.0048) {
  
  all_mr_fish_results_df <- remove_optimised_genomewide_crp_mr(all_mr_fish_results_df)
  
  # concordance df
  concordance_df <-
    make_concordance_df(all_mr_fish_results_df, pval_threshold = PVALUE_THRESHOLD)
  
  # get list of proteins with at least one significant result across all analyses
  proteins_with_at_least_one_sig_result <- concordance_df %>% 
    dplyr::filter(Outcome == !!OUTCOME) %>% 
    dplyr::filter(total_sig > 0) %>% 
    dplyr::pull(Protein) %>% 
    unique() %>% 
    sort()
  
  # Prepare for plotting a single outcome
  df <- all_mr_fish_results_df |>
    dplyr::filter(Outcome == !!OUTCOME,
                  Protein %in% c(proteins_with_at_least_one_sig_result,
                                 "CRP",
                                 "CRP (genome-wide)"))
  
  # label analysis groups
  WEIGHT_ORDER <- c("Trans", "Cis", "Genomewide")
  
  df <- df %>%
    dplyr::mutate(group = dplyr::case_when(
      !(Protein %in% c("CRP", "CRP (genome-wide)")) ~ "Trans",
      Protein == "CRP" ~ "Cis",
      Protein == "CRP (genome-wide)" ~ "Genomewide",
      TRUE ~ "Error"
    )) %>% 
    dplyr::mutate(group = factor(group, levels = WEIGHT_ORDER))
  
  stopifnot(identical(0L, nrow(df[df$group == "Error", ])))
  
  df <- df |>
    dplyr::mutate(group = forcats::fct_relabel(group, \(x) dplyr::case_when(
      x == "Trans" ~ "*trans*-for-CRP",
      x == "Cis" ~ "*cis*-for-CRP",
      x == "Genomewide" ~ "Genome-wide"
    )))
  
  # order studies
  STUDY_ORDER <- c(
    "CRP Ligthart",
    "CRP UKBB",
    "pQTL (DECODE)",
    "eQTL (GTEx, whole blood)"
  )
  
  # determine protein order
  # constants
  PROTEIN_ORDER <- df %>%
    # order by CRP 2018 estimates
    dplyr::filter(Weight == "CRP Ligthart") %>% 
    dplyr::arrange(Estimate) %>%
    dplyr::mutate(protein_order = as.numeric(rownames(.))) %>%
    dplyr::select(Protein,
                  protein_order)
  
  df <- df %>%
    # order proteins
    dplyr::left_join(PROTEIN_ORDER, by = "Protein") %>%
    dplyr::arrange(protein_order) %>%
    # bold proteins with >=x concordant results
    dplyr::left_join(
      concordance_df %>%
        dplyr::filter(Outcome == !!OUTCOME) %>%
        dplyr::select(Protein,
                      Outcome,
                      total_sig_dir),
      by = c("Protein", "Outcome")
    ) %>%
    dplyr::mutate(
      Protein = dplyr::case_when(
        total_sig_dir >= 3 ~ paste0("**", Protein, "**"),
        total_sig_dir == 2 ~ paste0("*", Protein, "*"),
        TRUE ~ Protein
      )
    ) %>%
    # relabel and reorder 'study' column
    dplyr::mutate(
      Weight = dplyr::case_when(
        Weight == "eQTL" ~ "eQTL (GTEx, whole blood)",
        TRUE ~ Weight
      )
    ) %>%
    dplyr::mutate(Weight = factor(Weight, levels = rev(STUDY_ORDER)))
  
  stopifnot(identical(0L, sum(is.na(df$Weight))))
  
  # labels; logodds -> odds
  if (identical(unique(df[["Outcome units"]]), "logOR")) {
    LOGODDS <- TRUE
    XLABEL <- paste0(OUTCOME, " (OR)")
  } else {
    LOGODDS <- FALSE
    XLABEL <- paste0 (OUTCOME, " (", unique(df[["Outcome units"]]), ")")
  }
  
  
  remove_eqtl <- TRUE
  if (remove_eqtl) {
    df <- df %>%
      dplyr::filter(Weight != "eQTL (GTEx, whole blood)") %>%
      dplyr::mutate(Weight = forcats::fct_drop(Weight, only = "eQTL (GTEx, whole blood)"))
  }
  
  # Ligthart et al only for genome-wide estimate
  expected_nrow <- as.integer(nrow(df) - 1)
  
  df <- df |>
    dplyr::filter(!(Weight == "CRP UKBB" & (Protein == "*CRP (genome-wide)*" |
                                              Protein == "CRP (genome-wide)")))
  
  stopifnot(identical(expected_nrow, nrow(df)) |
              identical(expected_nrow + 1L, nrow(df)))
  
  # plot
  result <- df %>%
    dplyr::mutate(shape = ifelse(Pvalue < !!PVALUE_THRESHOLD, "Significant", "Non-significant")) %>% 
    dplyr::filter(!is.na(Estimate)) %>% 
    ggforestplot::forestplot(
      estimate = Estimate,
      name = Protein,
      pvalue = Pvalue,
      se = SE,
      psignif = PVALUE_THRESHOLD,
      logodds = LOGODDS,
      colour = Weight,
      xlab = XLABEL
    ) +
    ggplot2::aes(shape = shape) + 
    ggplot2::scale_shape_manual(
      values = c("Significant" = 19, "Non-significant" = 1),  # Filled for significant, hollow for non-significant
      name = "P-value",  # Legend title
      labels = c(stringr::str_glue("Estimate (p ≥ {PVALUE_THRESHOLD})"), 
                 stringr::str_glue("Estimate (p < {PVALUE_THRESHOLD})"))  # Custom legend labels
    ) +
    ggplot2::theme(
      axis.text.x = ggtext::element_markdown(size = TEXT_SIZE),
      axis.text.y = ggtext::element_markdown(size = TEXT_SIZE),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.x = ggtext::element_markdown(size = X_TITLE_TEXT_SIZE),
      axis.title.y = ggtext::element_markdown(size = TEXT_SIZE),
      legend.title = ggplot2::element_blank(),
      legend.text = ggtext::element_markdown(size = LEGEND_TEXT_SIZE)
    ) +
    ggplot2::scale_colour_discrete(drop = FALSE) +
    # group by cis/genomewide/trans
    ggforce::facet_col(facets = ~ group,
                       scales = "free_y",
                       space = "free") +
    ggplot2::theme(strip.text = ggtext::element_markdown(),
                   legend.position = "bottom")
  
  return(result)
}

make_concordance_df <- function (all_mr_fish_results_df, pval_threshold) 
{
  concordance_df <- all_mr_fish_results_df %>%
    dplyr::select(Protein,
                  Outcome,
                  Estimate,
                  Pvalue,
                  Weight) %>% 
    dplyr::mutate(effect_dir = sign(Estimate)) %>%
    dplyr::mutate(sig_result = ifelse(Pvalue < !!pval_threshold, TRUE, FALSE)) %>% 
    dplyr::mutate(Weight = dplyr::case_when(
      Weight == "CRP Ligthart" ~ "crp2018",
      Weight == "CRP UKBB" ~ "crpukb",
      Weight == "pQTL (DECODE)" ~ "pqtl",
      Weight == "eQTL" ~ "eqtl_whole_blood",
    )) %>% 
    tidyr::pivot_wider(names_from = Weight,
                       values_from = c(Estimate,
                                       Pvalue,
                                       effect_dir,
                                       sig_result))
  
  concordance_df[["pvalue_threshold"]] <- pval_threshold
  
  stopifnot(identical(
    0L,
    concordance_df |>
      dplyr::count(Protein, Outcome) |>
      dplyr::filter(n > 1) |>
      nrow()
  ))
  
  concordance_df <- concordance_df |>
    dplyr::mutate(total_sig = rowSums(dplyr::across(
      tidyselect::starts_with("sig_result_")
    ),
    na.rm = TRUE)) %>%
    dplyr::mutate(crp2018_noteql_crpukb = ifelse(effect_dir_crp2018 != effect_dir_crpukb,
                                                 TRUE,
                                                 FALSE)) %>%
    dplyr::mutate(eqtl_noteql_pqtl = ifelse(effect_dir_pqtl != effect_dir_eqtl_whole_blood,
                                            TRUE,
                                            FALSE)) %>%
    # subtract '1' if both crp2018 and crpukb are significant, but in opposite
    # effect directions
    dplyr::mutate(total_sig_dir = total_sig) %>%
    dplyr::mutate(total_sig_dir = dplyr::case_when(((sig_result_crp2018 == TRUE) &
                                                      (sig_result_crpukb == TRUE) &
                                                      (crp2018_noteql_crpukb == TRUE)
    ) ~ total_sig_dir - 1,
    TRUE ~ total_sig_dir)) %>%
    # subtract '1' if both pQTL and eQTL are significant but in opposite directions
    dplyr::mutate(total_sig_dir = dplyr::case_when(((sig_result_pqtl == TRUE) &
                                                      (sig_result_eqtl_whole_blood == TRUE) &
                                                      (eqtl_noteql_pqtl == TRUE)
    ) ~ total_sig_dir - 1,
    TRUE ~ total_sig_dir)) %>%
    dplyr::select(
      total_sig_dir,
      total_sig,
      tidyselect::contains("sig_result"),
      tidyselect::contains("effect_dir"),
      # crp2018_noteql_crpukb,
      tidyselect::everything()
    )
  
  return(concordance_df)
}


# Figure S1 ---------------------------------------------------------------

plot_figure_s1 <- function(all_mr_fish_results_df, file_path = "Figure S1.png") {
  p <- plot_cis_v_gw_manhattan(all_mr_fish_results_df, genomewide_type = "CRP (genome-wide)")
  
  ggplot2::ggsave(
    filename = file_path,
    plot = p,
    width = 28,
    height = 24,
    units = "cm",
    dpi = 300
  )
  
  return(file_path)
}

plot_cis_v_gw_manhattan <- function(all_mr_fish_results_df,
                                    genomewide_type,
                                    phenotype_order = NULL,
                                    cis_gw_pval_threshold = 0.001,
                                    interaction_pval_threshold = 0.001,
                                    point_size = 4,
                                    legend_text_size = 10) {
  gw_v_crp_manhattan_df <- all_mr_fish_results_df %>%
    create_cis_v_gw_table(
      cis_gw_pval_threshold = cis_gw_pval_threshold,
      genomewide_type = genomewide_type,
      interaction_pval_threshold = interaction_pval_threshold,
      prettify = FALSE
    )
  
  # reshape and add back required cols
  gw_v_crp_manhattan_df <- gw_v_crp_manhattan_df %>%
    tidyr::pivot_longer(
      cols = tidyselect::ends_with(c("_Cis", "_Genomewide")) ,
      names_to = c(".value", "Protein"),
      names_pattern = "(.*)_([Cis|Genomewide])"
    ) %>%
    
    # remove concordant/<blank> values in 'Concordance' col
    dplyr::mutate(Concordance = dplyr::case_when(
      Concordance %in% c("", "Concordant") ~ "Not discordant",
      TRUE ~ Concordance
    )) %>%
    
    dplyr::mutate(Concordance = factor(
      Concordance,
      levels = c(
        'Discordant: both significant, opposite effect direction',
        'Discordant: CRP significant',
        'Discordant: genome-wide significant',
        'Discordant: both significant, same effect direction',
        'Not discordant'
      )
    )) %>%
    
    # re-add `cis_vs_gw_pval_raw_alpha`,
    dplyr::mutate(
      cis_vs_gw_pval_raw_alpha = dplyr::case_when(
        cis_vs_gw_pval_raw < !!interaction_pval_threshold ~ 1,
        TRUE ~ 0.5
      ),
      
      # cis vs gw
      cis_vs_gw_log10pval = log10(cis_vs_gw_pval_raw) * -1
    ) %>%
    dplyr::mutate(
      cis_vs_gw_log10pval = dplyr::case_when(
        cis_vs_gw_log10pval > 15 ~ 15,
        # limit size of extremely significant pvals
        cis_vs_gw_log10pval < -15 ~ -15,
        TRUE ~ cis_vs_gw_log10pval
      )
    )
  
  # plot
  p <- gw_v_crp_manhattan_df %>%
    plot_interaction_pvalues_manhattan(
      phenotype_order = phenotype_order,
      geom_text_size = 3,
      xlab_text_size = 10,
      ylab_text_size = 10,
      xlab_title_size = 12,
      ylab_title_size = 12,
      xhjust = 0.95,
      xvjust = 0.2,
      point_size = 3,
      ggrepel_labels = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(0.001),
      linetype = "dotted",
      colour = "purple",
      linewidth = 1,
      alpha = 0.6
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      colour = "dark green",
      linewidth = 1,
      alpha = 0.6
    ) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())
  
  # colour by concordance category
  
  # first reset geom_point layer
  p$layers[[2]] <- NULL
  
  p +
    ggplot2::geom_point(
      ggplot2::aes(colour = Concordance, shape = Concordance),
      size = point_size,
      show.legend = TRUE,
      inherit.aes = TRUE
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = legend_text_size)
    ) +
    ggplot2::scale_shape_manual(values = c(15, 16, 17, 18, 20)) +
    ggplot2::guides(shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
}

plot_interaction_pvalues_manhattan <- function (df,
                                                phenotype_order = NULL,
                                                point_size = NULL,
                                                geom_text_size = 3,
                                                xlab_text_size = 8,
                                                ylab_text_size = 8,
                                                xlab_title_size = 10,
                                                ylab_title_size = 10,
                                                xhjust = 0.95,
                                                xvjust = 0.2,
                                                ggrepel_labels = FALSE)
{
  df_by_phenotype <- df %>% dplyr::group_by(Outcome) %>%
    dplyr::slice(1L)
  if (!is.null(phenotype_order)) {
    if (all(!phenotype_order %in% df_by_phenotype[["Outcome"]])) {
      stop("all values in phenotype_order are not present in df_by_phenotype[['Outcome']]")
    }
    if (any(!phenotype_order %in% df_by_phenotype[["Outcome"]])) {
      warning("phenotype_order contains values that are not present in df[[x]]")
    }
    df_by_phenotype <- df_by_phenotype %>%
      dplyr::filter(.data[["Outcome"]] %in% phenotype_order)
    df_by_phenotype[["Outcome"]] <- factor(df_by_phenotype[["Outcome"]], levels = phenotype_order)
  }
  my_plot <- df_by_phenotype %>% ggplot2::ggplot(ggplot2::aes(x = Outcome, y = cis_vs_gw_log10pval, alpha = cis_vs_gw_pval_raw_alpha)) +
    ggplot2::scale_y_continuous(breaks = c(0, 5, 10, 15)) +
    if (ggrepel_labels) {
      ggrepel::geom_label_repel(
        data = subset(df_by_phenotype, cis_vs_gw_pval_raw_alpha == 1),
        ggplot2::aes(Outcome, cis_vs_gw_log10pval, label = Outcome),
        size = geom_text_size,
        colour = "white",
        fill = "#0066FF",
        fontface = "bold"
      )
    }
  else {
    ggplot2::geom_label(
      data = subset(df_by_phenotype, cis_vs_gw_pval_raw_alpha == 1),
      ggplot2::aes(Outcome, cis_vs_gw_log10pval + 1, label = Outcome),
      size = geom_text_size,
      colour = "black",
      fill = "white"
    )
  }
  my_plot <- my_plot + ggplot2::scale_alpha_continuous(guide = "none") +
    ggplot2::labs(x = "Outcome", y = "-log<sub>10</sub>(interaction p-value)") +
    ggplot2::theme_classic()
  my_plot <- my_plot + ggplot2::theme(
    axis.text.x = ggtext::element_markdown(
      angle = 90,
      size = xlab_text_size,
      hjust = xhjust,
      vjust = xvjust
    ),
    axis.text.y = ggtext::element_markdown(angle = 0, size = ylab_text_size),
    axis.title.x = ggtext::element_markdown(size = xlab_title_size),
    axis.title.y = ggtext::element_markdown(size = ylab_title_size)
  )
  if (!is.null(point_size)) {
    my_plot <- my_plot + ggplot2::geom_point(size = point_size, show.legend = FALSE)
  }
  else {
    my_plot <- my_plot + ggplot2::geom_point()
  } + return(my_plot)
}


create_cis_v_gw_table <- function(df,
                                  genomewide_type,
                                  cis_gw_pval_threshold = 0.001,
                                  interaction_pval_threshold = 0.001,
                                  log_odds_to_odds = TRUE,
                                  prettify = TRUE) {
  
  rlang::arg_match(genomewide_type,
                   c("CRP (genome-wide)", "CRP (genome-wide, optimised)"))
  
  result <- df %>%
    
    # CRP analyses only, Cis and Genomewide
    dplyr::filter(Weight == "CRP Ligthart",
                  Protein %in% c("CRP", genomewide_type)) %>% 
    
    dplyr::mutate(Protein = dplyr::case_when(Protein == "CRP" ~ "Cis",
                                             TRUE ~ "Genomewide")) %>% 
    
    # adjust significant_p threshold
    dplyr::mutate(significant_p = dplyr::case_when(Pvalue < !!cis_gw_pval_threshold ~ TRUE,
                                                   TRUE ~ FALSE)) %>% 
    dplyr::mutate(significant_p_nominal = dplyr::case_when(Pvalue < 0.05 ~ TRUE, 
                                                           TRUE ~ FALSE)) %>% 
    
    dplyr::select(
      Weight,
      Protein,
      `Phenotype category`,
      Outcome,
      Estimate,
      significant_p,
      significant_p_nominal,
      Pvalue,
      SE,
      unit = `Outcome units`,
      Estimator,
      `No. variants`,
      `F statistic`
    ) %>%
    tidyr::pivot_wider(
      names_from = "Protein",
      values_from = c("Estimate",
                      "significant_p",
                      "significant_p_nominal",
                      "Phenotype category",
                      "Pvalue",
                      "SE",
                      "unit",
                      "Estimator",
                      "No. variants",
                      "F statistic")
    )
  
  # Calculate interaction p values
  result <- result %>% 
    dplyr::mutate(cis_vs_gw_pval_raw = calc_interaction_pvalues(effect1 = Estimate_Cis,
                                                                effect2 = Estimate_Genomewide,
                                                                se1 = SE_Cis,
                                                                se2 = SE_Genomewide,
                                                                pvals_only = TRUE))
  
  # indicate which are significant
  result <- result %>%
    dplyr::mutate(
      sig_cis_vs_gw_pval_raw = dplyr::case_when(
        cis_vs_gw_pval_raw < interaction_pval_threshold ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    
    dplyr::mutate(# scientific notation
      cis_vs_gw_pval_sci_not = .scientific_format(cis_vs_gw_pval_raw)) %>%
    
    # bold if significant
    dplyr::mutate(
      cis_vs_gw_pval_sci_not = dplyr::case_when(
        sig_cis_vs_gw_pval_raw == TRUE ~ paste0("**", cis_vs_gw_pval_sci_not, "**"),
        TRUE ~ cis_vs_gw_pval_sci_not
      )
    ) %>%
    
    # mutate concordance column
    dplyr::mutate(
      Genomewide = dplyr::case_when(
        significant_p_Genomewide == TRUE &
          sign(Estimate_Genomewide) == 1 ~ "+",
        significant_p_Genomewide == TRUE &
          sign(Estimate_Genomewide) == -1 ~ "-",
        TRUE ~ ""
      ),
      `Cis-CRP` = dplyr::case_when(
        significant_p_Cis == TRUE & sign(Estimate_Cis) == 1 ~ "+",
        significant_p_Cis == TRUE & sign(Estimate_Cis) == -1 ~ "-",
        TRUE ~ ""
      )
    ) %>%
    
    ## concordant if both significant in same direction and non-significant interaction p-value
    dplyr::mutate(Concordance = dplyr::case_when(
      ((sign(Estimate_Cis) == sign(Estimate_Genomewide)) &
         (sig_cis_vs_gw_pval_raw == FALSE) &
         (significant_p_nominal_Cis &
            significant_p_nominal_Genomewide)
      ) ~ "Concordant",
      
      ## Discordant (opposite effect dir) if both significant in opposite directions, with significant interaction pval
      ((sign(Estimate_Cis) != sign(Estimate_Genomewide)) &
         (sig_cis_vs_gw_pval_raw == TRUE) &
         (significant_p_nominal_Cis &
            significant_p_nominal_Genomewide)
      ) ~ "Discordant: both significant, opposite effect direction",
      TRUE ~ ""
    )) %>%
    
    ## Discordant (same effect direction) - significant interaction p value, both significant
    dplyr::mutate(Concordance = dplyr::case_when(
      ((sign(Estimate_Cis) == sign(Estimate_Genomewide)) &
         (sig_cis_vs_gw_pval_raw == TRUE) &
         (significant_p_nominal_Cis &
            significant_p_nominal_Genomewide)
      ) ~ "Discordant: both significant, same effect direction",
      
      ## only CRP significant
      (sig_cis_vs_gw_pval_raw == TRUE) &
        (significant_p_nominal_Cis &
           !significant_p_nominal_Genomewide) ~ "Discordant: CRP significant",
      
      ## only genome-wide significant
      (sig_cis_vs_gw_pval_raw == TRUE) &
        (!significant_p_nominal_Cis &
           significant_p_nominal_Genomewide) ~ "Discordant: genome-wide significant",
      TRUE ~ Concordance
    ))
  
  # relabel Genomewide and cis-CRP from "+/-" to Estimate (CI)
  # convert log(OR) to OR for binary traits
  
  if (log_odds_to_odds) {
    result <- result %>%
      dplyr::mutate(
        Genomewide = dplyr::case_when(
          unit_Cis == "logOR" ~ paste0(
            .format(exp(Estimate_Genomewide)),
            # "\\\n(",
            " (",
            .format(exp((
              Estimate_Genomewide - 2 * SE_Genomewide
            ))),
            ", ",
            .format(exp((
              Estimate_Genomewide + 2 * SE_Genomewide
            ))),
            ")"
            # ,
            # model_Genomewide
          ),
          TRUE ~ paste0(
            .format(Estimate_Genomewide),
            # "\\\n(",
            " (",
            .format((Estimate_Genomewide - 2 * SE_Genomewide)),
            ", ",
            .format((Estimate_Genomewide + 2 * SE_Genomewide)),
            ")"
            # ,
            # model_Genomewide
          )
        ),
        `Cis-CRP` = dplyr::case_when(
          unit_Cis == "logOR" ~ paste0(
            .format(exp(Estimate_Cis)),
            # "\\\n(",
            " (",
            .format(exp((
              Estimate_Cis - 2 * SE_Cis
            ))),
            ", ",
            .format(exp((
              Estimate_Cis + 2 * SE_Cis
            ))),
            ")"
            # ,
            # model_Cis
          ),
          TRUE ~ paste0(
            .format(Estimate_Cis),
            # "\\\n(",
            " (",
            .format((Estimate_Cis - 2 * SE_Cis)),
            ", ",
            .format((Estimate_Cis + 2 * SE_Cis)),
            ")"
            # ,
            # model_Cis
          )
        )
      ) %>%
      dplyr::mutate(across(
        tidyselect::starts_with("unit_"),
        ~ dplyr::case_when(stringr::str_detect(.x, "logOR") ~ "OR",
                           TRUE ~ .x)
      ))
  } else {
    result <- result %>%
      dplyr::mutate(
        Genomewide = paste0(
          .format(Estimate_Genomewide),
          # "\\\n(",
          " (",
          .format((Estimate_Genomewide - 2 * SE_Genomewide)),
          ", ",
          .format((Estimate_Genomewide + 2 * SE_Genomewide)),
          ")"
          # ,
          # model_Genomewide
        ),
        `Cis-CRP` = paste0(
          .format(Estimate_Cis),
          # "\\\n(",
          " (",
          .format((Estimate_Cis - 2 * SE_Cis)),
          ", ",
          .format((Estimate_Cis + 2 * SE_Cis)),
          ")"
          # ,
          # model_Cis
        )
      )
  }
  
  # # bold (p<0.001, `cis_gw_pval_threshold`)/italicise (p<0.05) significant cis/GW results
  result <- result %>%
    dplyr::mutate(
      Genomewide = dplyr::case_when(
        (Pvalue_Genomewide < 0.05) &
          (Pvalue_Genomewide >= !!cis_gw_pval_threshold) ~ paste0("_",
                                                                  Genomewide,
                                                                  "_"),
        Pvalue_Genomewide < !!cis_gw_pval_threshold ~ paste0("**",
                                                             Genomewide,
                                                             "**"),
        TRUE ~ Genomewide
      ),
      `Cis-CRP` = dplyr::case_when(
        (Pvalue_Cis < 0.05) &
          (Pvalue_Cis >= !!cis_gw_pval_threshold) ~ paste0("_",
                                                           `Cis-CRP`,
                                                           "_"),
        Pvalue_Cis < !!cis_gw_pval_threshold ~ paste0("**",
                                                      `Cis-CRP`,
                                                      "**"),
        TRUE ~ `Cis-CRP`
      )
    ) %>%
    
    # remove NA results
    dplyr::mutate(
      Genomewide = dplyr::case_when(is.na(Estimate_Genomewide) ~ "–",
                                    TRUE ~ Genomewide),
      `Cis-CRP` = dplyr::case_when(is.na(Estimate_Cis) ~ "–",
                                   TRUE ~ `Cis-CRP`)
    )
  
  return(result) 
}

#' Calculate interaction p-values
#'
#' @param effect1 numeric vector
#' @param effect2 numeric vector
#' @param se1 numeric vector
#' @param se2 numeric vector
#' @param pvals_only logical. Indicate whether to return only interaction
#'   p-values, or a dataframe also including the input betas and standard
#'   errors. Default = \code{TRUE}
#'
#' @return df or numeric vector of interaction p-values (determined by
#'   \code{pvals_only} argument)
#' @export
calc_interaction_pvalues <- function(effect1,
                                     effect2,
                                     se1,
                                     se2,
                                     pvals_only = TRUE) {
  # difference
  difference <- effect1 - effect2
  
  # SE of the difference
  se_diff <- sqrt(se1^2 + se2^2)
  
  # test
  pvalue <- stats::pnorm(abs(difference / se_diff), lower.tail = F) * 2
  
  # results
  res <- data.frame(point = difference, se = se_diff, pval = pvalue)
  
  if (pvals_only) {
    return(res$pval)
  } else if (pvals_only == FALSE) {
    return(res)
  }
}

.scientific_format <- function(x, 
                               digits = 1,
                               min_value = "<1.0 X 10^-16^") {
  dplyr::case_when(
    !(x %in% c(0, 1)) ~ formatC(x,
                                format = "e",
                                digits = digits) %>%
      stringr::str_replace(pattern = "e-0",
                           replacement = "e-") %>%
      stringr::str_replace(pattern = "e-",
                           replacement = " x 10^-") %>%
      paste0("^"),
    x == 0 ~ min_value,
    x == 1 ~ "1"
  )
}

.format <- function(x, digits = 2) {
  round(x, digits = digits) %>% 
    format(scientific = FALSE,
           trim = TRUE)
}

# Figure S2 ---------------------------------------------------------------

plot_figure_s2 <- function(all_mr_fish_results_df,
                           druggable_proteins,
                           file_path = "Figure S2.png") {
  
  result <- plot_mr_heatmap(all_mr_fish_results_df,
                            Weight = "CRP Ligthart",
                            proteins_to_bold = druggable_proteins)
  
  ggplot2::ggsave(
    file_path,
    plot = result,
    width = 13,
    height = 13,
    units = "cm",
    dpi = 600
  )
  
  return(file_path)
}


# Figure S3 ---------------------------------------------------------------

plot_figure_s3 <- function(all_mr_fish_results_df,
                           druggable_proteins,
                           file_path = "Figure S3.png") {
  p_list <- list(`pQTL (deCODE)` = plot_mr_heatmap(all_mr_fish_results_df,
                                                   Weight = "pQTL (DECODE)",
                                                   proteins_to_bold = druggable_proteins),
                 `eQTL (GTEx, whole blood)` = plot_mr_heatmap(all_mr_fish_results_df,
                                                              Weight = "eQTL",
                                                              proteins_to_bold = druggable_proteins))
  
  p_list <- p_list %>% 
    purrr::imap(~ .x + ggplot2::labs(subtitle = .y))
  
  p_list$`pQTL (deCODE)` <- p_list$`pQTL (deCODE)` + ggplot2::theme(legend.position = 'None')
  
  # combine
  cowplot_plot <-
    cowplot::plot_grid(plotlist = p_list,
                       ncol = 1,
                       rel_heights = c(1, 2))
  
  ggplot2::ggsave(file_path, 
                  plot = cowplot_plot, 
                  width = 13, 
                  height = 14,
                  units = "cm",
                  dpi = 600)
  
  return(file_path)
}

plot_mr_heatmap <- function(all_mr_fish_results_df,
                            proteins_to_bold,
                            Weight,
                            lwd = 0.3,
                            linetype = 1,
                            pval_threshold = 0.05) {
  
  stopifnot(Weight %in% unique(all_mr_fish_results_df$Weight))
  
  df <- all_mr_fish_results_df |>
    dplyr::filter(Outcome != "CRP") |>
    dplyr::filter(stringr::str_detect(Protein, "genome", negate = TRUE)) |>
    dplyr::filter(Weight == !!Weight)
  
  # fill empty Protein-Outcome pairs
  all_commbinations <- tidyr::expand_grid(Protein = unique(df$Protein),
                                          Outcome = unique(df$Outcome))
  
  combinations_to_add <- all_commbinations |>
    dplyr::anti_join(df |>
                       dplyr::distinct(Protein, Outcome),
                     by = dplyr::join_by(Protein, Outcome)) |>
    dplyr::left_join(df |>
                       dplyr::distinct(Outcome,
                                       `Phenotype category`),
                     by = "Outcome")
  
  df <- df |>
    dplyr::bind_rows(combinations_to_add)
  
  # set up protein and outcome ordering for plot
  phenotype_order <- df |>
    dplyr::distinct(`Phenotype category`,
                    Outcome) |>
    dplyr::filter(Outcome != "CRP") |>
    dplyr::arrange(`Phenotype category`,
                   Outcome) |>
    dplyr::pull(Outcome)
  
  
  protein_order <- df |>
    dplyr::distinct(Protein) |>
    dplyr::filter(Protein != "CRP") |>
    dplyr::filter(stringr::str_detect(Protein, "genome", negate = TRUE)) |>
    dplyr::pull(Protein) |>
    sort() |>
    rev()
  
  protein_order <- c(protein_order, "CRP")
  
  # Highlight selected proteins with bold
  stopifnot(is.character(proteins_to_bold))
  
  df <- df %>%
    dplyr::mutate(Protein = ifelse(
      Protein %in% proteins_to_bold,
      yes = paste0("**", Protein, "**"),
      no = Protein
    ))
  
  protein_order <- protein_order %>%
    purrr::map_chr(\(x) ifelse(
      x %in% proteins_to_bold,
      yes = paste0("**", x, "**"),
      no = x
    ))
  
  result <- df %>%
    # relabel categories with long names
    dplyr::mutate(
      `Phenotype category` = dplyr::case_when(
        `Phenotype category` == "Eye and adnexa" ~ "Eyes",
        `Phenotype category` == "Genitourinary system" ~ "Renal",
        TRUE ~ `Phenotype category`
      )) |>
    dplyr::mutate(log10pval_dir = -log10(Pvalue) * sign(Estimate)) |>
    dplyr::mutate(log10pval_dir = dplyr::case_when(abs(log10pval_dir) > 15 ~ 15 * sign(log10pval_dir),
                                                   TRUE ~ log10pval_dir)) |>
    dplyr::mutate(starred = dplyr::case_when(Pvalue < !!pval_threshold ~ "*",
                                             TRUE ~ "")) |>
    tidyr::complete(Outcome, Protein,
                    fill = list(log10pval_dir = NA, 
                                starred = NA)) |>
    plot_neglog10dir_heatmap(
      x = "Outcome",
      y = "Protein",
      x_order = phenotype_order,
      y_order = protein_order,
      fill = "log10pval_dir",
      lwd = lwd,
      linetype = linetype,
      pval_threshold = 0.05,
      xlab_text_size = 4,
      ylab_text_size = 4,
      xlab_text_angle = 90,
      ylab_text_angle = 0,
      xtitle_text_size = 6,
      ytitle_text_size = 6,
      legend_title_size = 6,
      legend_text_size = 5,
      label = "starred"
    ) +
    ggplot2::facet_grid( ~ .data[["Phenotype category"]], scales = "free", space = "free") +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(size = 3.5),
      legend.key.height = ggplot2::unit(0.4, 'cm'),
      legend.key.width = ggplot2::unit(0.8, 'cm'),
      legend.position = "bottom", 
      plot.background = ggplot2::element_rect(fill = "white"),
      axis.title.x =  ggplot2::element_blank(),
      axis.title.y =  ggplot2::element_blank()
    ) 
  
  return(result)
}

plot_neglog10dir_heatmap <- function(df,
                                     x = "Outcome",
                                     y = "Protein",
                                     fill = "log10pval_dir",
                                     pval_threshold = NULL,
                                     xlab = "Phenotype",
                                     ylab = "Protein",
                                     legend_lab = "Direction x -log<sub>10</sub>(p-value)",
                                     xlab_text_angle = 90,
                                     ylab_text_angle = 0,
                                     label = "starred",
                                     label_size = TRUE,
                                     ...)
{
  if (!is.null(pval_threshold)) {
    df <- df %>% dplyr::mutate(heatmap_fill = dplyr::case_when(dplyr::between(
      .data[[fill]], log10(pval_threshold), -log10(pval_threshold)
    ) ~
      0, TRUE ~ .data[[fill]]))
  }
  plot_heatmap(
    df = df,
    x = x,
    y = y,
    fill = "heatmap_fill",
    xlab = xlab,
    ylab = ylab,
    legend_lab = legend_lab,
    xlab_text_angle = xlab_text_angle,
    ylab_text_angle = ylab_text_angle,
    label = label,
    label_size = label_size,
    ...
  ) + ggplot2::scale_fill_distiller(palette = "RdBu", limits = c(-15, 15))
}

plot_heatmap <- function(df,
                         x,
                         y,
                         fill,
                         label = NULL,
                         label_size = TRUE,
                         x_order = NULL,
                         y_order = NULL,
                         xlab = "x title",
                         ylab = "y title",
                         legend_lab = "legend title",
                         xlab_text_angle = 90,
                         ylab_text_angle = 0,
                         xlab_text_size = NULL,
                         ylab_text_size = NULL,
                         xhjust = 0.95,
                         xvjust = 0.2,
                         yhjust = NULL,
                         yvjust = NULL,
                         xtitle_text_size = NULL,
                         ytitle_text_size = NULL,
                         legend_title_size = NULL,
                         legend_text_size = NULL,
                         geom_tile_color = "white",
                         lwd = 0.5,
                         linetype = 1)
{
  if (!is.null(x_order)) {
    if (all(!x_order %in% df[[x]])) {
      stop("all values in x_order are not present in df[[x]]")
    }
    if (any(!x_order %in% df[[x]])) {
      warning("x_order contains values that are not present in df[[x]]")
    }
    df <- df %>% dplyr::filter(.data[[x]] %in% x_order)
    df[[x]] <- factor(df[[x]], levels = x_order)
  }
  if (!is.null(y_order)) {
    if (all(!y_order %in% df[[y]])) {
      stop("all values in y_order are not present in df[[x]]")
    }
    if (any(!y_order %in% df[[y]])) {
      warning("y_order contains values that are not present in df[[y]]")
    }
    df <- df %>% dplyr::filter(.data[[y]] %in% y_order)
    df[[y]] <- factor(df[[y]], levels = y_order)
  }
  heatmap <- df %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])) +
    ggplot2::geom_tile(color = geom_tile_color,
                       lwd = lwd,
                       linetype = linetype) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) + ggplot2::labs(fill = legend_lab) +
    ggplot2::theme(
      axis.text.x = ggtext::element_markdown(
        angle = xlab_text_angle,
        size = xlab_text_size,
        hjust = xhjust,
        vjust = xvjust
      ),
      axis.text.y = ggtext::element_markdown(
        angle = ylab_text_angle,
        size = ylab_text_size,
        hjust = yhjust,
        vjust = yvjust
      ),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.x = ggtext::element_markdown(size = xtitle_text_size),
      axis.title.y = ggtext::element_markdown(size = ytitle_text_size),
      legend.title = ggtext::element_markdown(size = legend_title_size),
      legend.text = ggtext::element_markdown(size = legend_text_size),
    ) + ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0))
  
  if (!is.null(label)) {
    heatmap <- heatmap + 
      ggplot2::geom_text(ggplot2::aes(label = .data[[label]]), size = label_size)
  }
  return(heatmap)
}


# Figure S4 ---------------------------------------------------------------

plot_figure_s4 <- function(all_mr_fish_results_df,
                           druggable_proteins,
                           file_path = "Figure S4.png") {
  # plot
  replication_heatmap <- plot_replication_heatmap(
    all_mr_fish_results_df = all_mr_fish_results_df,
    druggable_proteins = druggable_proteins
  )
  
  # save plot
  ggplot2::ggsave(
    filename = file_path,
    plot = replication_heatmap,
    width = 45,
    height = 25,
    units = "cm",
    dpi = 600
  )
  
  return(file_path)
}

plot_replication_heatmap <- function(all_mr_fish_results_df,
                                     replication_pval = 0.0048,
                                     n_sig_replicates = 3,
                                     TEXT_SIZE = 10,
                                     lwd = 0.8,
                                     linetype = 1,
                                     druggable_proteins = NULL,
                                     facet_by_phenotype_category = TRUE) {
  # validate args
  if (!is.null(druggable_proteins)) {
    stopifnot(is.character(druggable_proteins))
  }
  
  # make concordance df with a 'concordance' indicator column
  required_cols <- c("Outcome", "Protein", "Estimate", "Pvalue", "Weight")
  
  concordance_df <- all_mr_fish_results_df %>%
    dplyr::select(tidyselect::all_of(required_cols)) %>%
    dplyr::mutate(effect_dir = sign(Estimate)) %>%
    dplyr::mutate(sig_result = ifelse(Pvalue < !!replication_pval, TRUE, FALSE)) %>%
    dplyr::mutate(
      Weight = dplyr::case_when(
        Weight == "CRP Ligthart" ~ "crp2018",
        Weight == "CRP UKBB" ~ "crpukb",
        Weight == "eQTL" ~ "eqtl_whole_blood",
        Weight == "pQTL (DECODE)" ~ "pqtl"
      )
    ) |>
    tidyr::pivot_wider(
      names_from = Weight,
      values_from = c(Estimate, Pvalue, effect_dir, sig_result)
    ) |>
    dplyr::filter(Outcome != "CRP",
                  stringr::str_detect(Protein, "genome", negate = TRUE))
  
  stopifnot(identical(
    0L,
    concordance_df |>
      dplyr::count(Outcome, Protein) |>
      dplyr::filter(n > 1) |>
      nrow()
  ))
  
  concordance_df[["pvalue_threshold"]] <- replication_pval
  
  # mutate total_sig_dir col. NOTE: only count crp_ukb result if effect direction
  # is concordant with crp2018
  concordance_df <- concordance_df %>%
    dplyr::mutate(total_sig = rowSums(dplyr::across(tidyselect::all_of(
      c(
        "sig_result_crp2018",
        "sig_result_crpukb",
        "sig_result_pqtl",
        "sig_result_eqtl_whole_blood"
      )
    )), na.rm = TRUE)) %>%
    dplyr::mutate(
      crp2018_noteql_crpukb = ifelse(effect_dir_crp2018 != effect_dir_crpukb, TRUE, FALSE),
      pqtl_noteql_eqtl = ifelse(effect_dir_pqtl != effect_dir_eqtl_whole_blood, TRUE, FALSE)
    ) %>%
    # subtract '1' if both crp2018 and crpukb are significant, but in opposite
    # effect directions
    dplyr::mutate(total_sig_dir = dplyr::case_when(((sig_result_crp2018 == TRUE) &
                                                      (sig_result_crpukb == TRUE) &
                                                      (crp2018_noteql_crpukb == TRUE)
    ) ~ total_sig - 1, TRUE ~ total_sig)) %>%
    
    # select required cols
    dplyr::select(
      total_sig_dir,
      total_sig,
      tidyselect::contains("sig_result"),
      tidyselect::contains("effect_dir"),
      tidyselect::everything()
    )
  
  # summarise for each protein, number of outcomes with significant results
  concordance_df_n_outcomes_sig_per_protein <-
    concordance_df %>%
    dplyr::group_by(Protein) %>%
    dplyr::summarize(n_sig_outcomes = sum(total_sig_dir >= !!n_sig_replicates, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    # order
    dplyr::arrange(n_sig_outcomes)
  
  # summarise for each protein, number of proteins with significant results
  concordance_df_n_proteins_sig_per_outcome <-
    concordance_df %>%
    dplyr::group_by(Outcome) %>%
    dplyr::summarize(n_sig_proteins = sum(total_sig_dir >= !!n_sig_replicates, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    # order
    dplyr::arrange(n_sig_proteins)
  
  # Order by marginal counts
  protein_order <-  concordance_df_n_outcomes_sig_per_protein$Protein %>%
    subset(., !. == "CRP") %>%
    c("CRP", .)
  
  outcome_order <- concordance_df_n_proteins_sig_per_outcome$Outcome
  
  # bold selected proteins
  if (!is.null(druggable_proteins)) {
    concordance_df <- concordance_df %>%
      dplyr::mutate(Protein = ifelse(
        Protein %in% !!druggable_proteins,
        yes = paste0("**", Protein, "**"),
        no = Protein
      ))
    
    concordance_df_n_outcomes_sig_per_protein <- concordance_df_n_outcomes_sig_per_protein %>%
      dplyr::mutate(Protein = ifelse(
        Protein %in% !!druggable_proteins,
        yes = paste0("**", Protein, "**"),
        no = Protein
      ))
    
    protein_order <- protein_order %>%
      purrr::map_chr( ~ ifelse(
        .x %in% druggable_proteins,
        yes = paste0("**", .x, "**"),
        no = .x
      ))
  }
  
  # plot concordance heatmap
  p <- concordance_df %>%
    # add phenotype_category column
    dplyr::left_join(
      all_mr_fish_results_df %>%
        dplyr::rename(phenotype_category = `Phenotype category`) %>%
        dplyr::distinct(Outcome, phenotype_category),
      by = "Outcome"
    ) %>%
    
    # relabel categories with long names
    dplyr::mutate(
      phenotype_category = dplyr::case_when(
        phenotype_category == "Eye and adnexa" ~ "Eyes",
        phenotype_category == "Genitourinary system" ~ "Renal",
        TRUE ~ phenotype_category
      )
    ) %>%
    # set NAs to zero
    dplyr::mutate(total_sig_dir = ifelse(is.na(total_sig_dir), 0, total_sig_dir)) %>%
    # make factor - so can work with scale_fill_manual
    dplyr::mutate(total_sig_dir = as.factor(total_sig_dir)) %>%
    plot_heatmap(
      x = "Outcome",
      y = "Protein",
      y_order = protein_order,
      x_order = outcome_order,
      fill = "total_sig_dir",
      xlab = "**Phenotype**",
      ylab = "**Gene**",
      legend_lab = "N significant\n MR replicates",
      xlab_text_angle = 90,
      ylab_text_angle = 0,
      xhjust = 0.95,
      xvjust = 0.2,
      xlab_text_size = TEXT_SIZE,
      ylab_text_size = TEXT_SIZE,
      xtitle_text_size = TEXT_SIZE,
      ytitle_text_size = TEXT_SIZE,
      legend_title_size = TEXT_SIZE,
      legend_text_size = TEXT_SIZE,
      # plotly_plot = TRUE,
      lwd = lwd,
      linetype = linetype
    ) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(size = TEXT_SIZE),
      legend.key.height = ggplot2::unit(2, 'cm'),
      legend.key.width = ggplot2::unit(1, 'cm'),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    ) +
    
    ggplot2::scale_fill_manual(
      values = c("grey97", "yellow", "orange", "red", "darkred"),
      na.translate = FALSE
    )
  
  if (facet_by_phenotype_category) {
    # facet by phenotype category
    p <- p +
      ggplot2::facet_grid( ~ .data[["phenotype_category"]], scales = "free", space = "free")
  }
  
  # get p without legend
  p_minus_legend <-
    p + ggplot2::theme(legend.position = 'None')
  
  p_legend <- cowplot::get_legend(p)
  
  # add marginal boxplot
  # https://stackoverflow.com/questions/38856309/r-grid-arrange-marginal-plots-to-ggplot2-heatmap-geom-tile/46415093
  gg_cols <-
    concordance_df_n_outcomes_sig_per_protein %>%
    # order proteins
    dplyr::mutate(Protein = factor(Protein, levels = protein_order)) %>%
    ggplot2::ggplot(ggplot2::aes(x = Protein, y = n_sig_outcomes)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::ylab(paste0(
      "N outcomes with >",
      n_sig_replicates - 1,
      "\nsignificant MR replicates"
    )) +
    ggplot2::coord_flip() +
    ggplot2::theme_classic() +
    
    ggplot2::scale_y_continuous(expand = c(0, 0),
                                breaks = seq(
                                  from = 0,
                                  to = max(
                                    concordance_df_n_outcomes_sig_per_protein$n_sig_outcomes,
                                    na.rm = TRUE
                                  ),
                                  by = 2
                                )) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = TEXT_SIZE),
      axis.title.x = ggplot2::element_text(size = TEXT_SIZE)
    )
  
  gg_rows <-
    concordance_df_n_proteins_sig_per_outcome %>%
    # add phenotype_category column
    dplyr::left_join(
      all_mr_fish_results_df %>%
        dplyr::rename(phenotype_category = `Phenotype category`) %>%
        dplyr::distinct(Outcome, phenotype_category),
      by = "Outcome"
    ) %>%
    
    # relabel categories with long names
    dplyr::mutate(
      phenotype_category = dplyr::case_when(
        phenotype_category == "Eye and adnexa" ~ "Eyes",
        phenotype_category == "Genitourinary system" ~ "Renal",
        TRUE ~ phenotype_category
      )
    ) %>%
    # order outcomes
    dplyr::mutate(Outcome = factor(Outcome, levels = outcome_order)) %>%
    # plot
    ggplot2::ggplot(ggplot2::aes(x = Outcome, y = n_sig_proteins)) +
    ggplot2::geom_bar(stat = "identity",
                      position = "dodge",
                      width = 0.9) +
    ggplot2::ylab(paste0(
      "N proteins >",
      n_sig_replicates - 1,
      "\nsignificant MR replicates"
    )) +
    ggplot2::theme_classic() +
    
    ggplot2::scale_y_continuous(expand = c(0, 0),
                                breaks = seq(
                                  from = 0,
                                  to = max(
                                    concordance_df_n_proteins_sig_per_outcome$n_sig_proteins,
                                    na.rm = TRUE
                                  ),
                                  by = 2
                                )) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.length = ggplot2::unit(0, "pt"),
      axis.title.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = TEXT_SIZE),
      axis.title.y = ggplot2::element_text(size = TEXT_SIZE)
    ) +
    
    ggplot2::facet_grid( ~ .data[["phenotype_category"]], scales = "free", space = "free") +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_blank()
    )
  
  gg_empty <-
    concordance_df %>%
    ggplot2::ggplot(ggplot2::aes(x = Protein, y = Outcome)) +
    ggplot2::geom_blank() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      line = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    )
  
  # plot - egg
  egg_main_plot <-
    egg::ggarrange(
      # plots
      gg_rows,
      gg_empty,
      # p_legend, does not work with egg package
      p_minus_legend,
      gg_cols,
      # arrange
      nrow = 2,
      ncol = 2,
      widths = c(4, 1),
      heights = c(1, 4)
    )
  
  # add legend
  cowplot_plot <-
    cowplot::plot_grid(cowplot::plot_grid(egg_main_plot),
                       p_legend,
                       rel_widths = c(1, 0.1)) +
    
    cowplot::theme_nothing() +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  
  return(cowplot_plot)
}

# Figure A1 ---------------------------------------------------------------

plot_figure_a1 <- function(all_mr_fish_results_df, file_path = "Figure A1.png") {
  
  crp_pqtl_eqtl_df <- create_crp_pqtl_eqtl_df(all_mr_fish_results_df)
  
  crp_pqtl_eqtl_df_summarised <- c(1, 0.5, 0.1, 0.05, 0.01, 0.001, 0.0001) %>%
    purrr::set_names() %>%
    purrr::map(
      \(x) compare_effect_direction_trans_v_cis_weighted_cis_mr(crp_pqtl_eqtl_df = crp_pqtl_eqtl_df, pval_threshold = x)
    )
  
  figure_a1_ggplot <- c("eqtl", "pqtl") %>%
    purrr::set_names() %>%
    purrr::map(
      \(x) summarise_summary_crp_pqtl_eqtl_df(summary_crp_pqtl_eqtl_df = crp_pqtl_eqtl_df_summarised, pqtl_eqtl = x)
    ) %>%
    dplyr::bind_rows(.id = "eqtl_pqtl") %>%
    dplyr::mutate(eqtl_pqtl = dplyr::case_when(eqtl_pqtl == "eqtl" ~ "eQTL", eqtl_pqtl == "pqtl" ~ "pQTL")) %>%
    forest_plot_summarise_summary_crp_pqtl_eqtl_df() + ggplot2::facet_grid(~ eqtl_pqtl)
  
  ggplot2::ggsave(
    filename = file_path,
    plot = figure_a1_ggplot,
    width = 25,
    height = 10,
    units = "cm",
    dpi = 600
  )
  
  return(file_path)
}

forest_plot_summarise_summary_crp_pqtl_eqtl_df <-
  function(summarised_summary_crp_pqtl_eqtl_df,
           ylab_title = "Agreement (%)",
           nudge_x = 0.3) {
    summarised_summary_crp_pqtl_eqtl_df %>%
      dplyr::filter(pval_threshold %in% c("1", "0.1", "0.01", "0.001", "0.0001")) %>%
      ggplot2::ggplot(ggplot2::aes(x = pval_threshold, y = pct_expected)) +
      ggplot2::geom_pointrange(
        ggplot2::aes(ymin = ymin, ymax = ymax),
        colour = "black",
        alpha = 0.9,
        size = 0.5
      ) +
      
      ggplot2::xlab("P value threshold") +
      ggplot2::ylab(ylab_title) +
      ggplot2::geom_hline(
        yintercept = 50,
        linetype = "dashed",
        colour = "darkred"
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 100)) +
      ggplot2::scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100)) +
      ggplot2::theme_bw() +
      ggplot2::geom_text(ggplot2::aes(label = paste0("(", total_n, ")")), nudge_x = nudge_x)
  }

summarise_summary_crp_pqtl_eqtl_df <-
  function(summary_crp_pqtl_eqtl_df, pqtl_eqtl) {
    result <- summary_crp_pqtl_eqtl_df %>%
      purrr::map(\(x) x[[pqtl_eqtl]] %>%
                   dplyr::count(got_expected_result) %>%
                   dplyr::mutate(pct = n / sum(.$n))) %>%
      dplyr::bind_rows(.id = "pval_threshold") %>%
      dplyr::group_by(pval_threshold) %>%
      dplyr::mutate(total_n = sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(got_expected_result == "Yes") %>%
      dplyr::mutate(pval_threshold = as.numeric(pval_threshold)) %>%
      dplyr::arrange(pval_threshold) %>%
      dplyr::select(pval_threshold = pval_threshold,
                    n_expected_result = n,
                    total_n,
      ) %>%
      # format pval threshold, keeping order from largest to smallest
      dplyr::arrange(dplyr::desc(pval_threshold)) %>%
      tibble::rownames_to_column(var = "pval_order") %>%
      dplyr::mutate(pval_order = as.integer(pval_order)) %>%
      dplyr::mutate(
        pval_threshold = format(pval_threshold, scientific = FALSE) %>%
          stringr::str_remove_all("\\.?0+$")
      ) %>%
      dplyr::mutate(pval_threshold = factor(pval_order, levels = pval_order, labels = pval_threshold))
    
    # calculate pct_expected, confidence intervals and p values (see
    # http://www.sthda.com/english/wiki/one-proportion-z-test-in-r )
    result <- result %>%
      dplyr::mutate(binom_test = purrr::map2(n_expected_result, total_n, ~ broom::tidy(
        binom.test(
          x = .x,
          n = .y,
          p = 0.5,
          alternative = "two.sided"
        )
      ))) %>%
      dplyr::mutate(
        ymin = purrr::map_dbl(binom_test, ~ .x$conf.low) * 100,
        ymax = purrr::map_dbl(binom_test, ~ .x$conf.high * 100),
        p.value = purrr::map_dbl(binom_test, ~ .x$p.value),
        pct_expected = purrr::map_dbl(binom_test, ~ .x$estimate * 100)
      ) %>%
      dplyr::mutate(p.value = as.numeric(p.value)) %>%
      dplyr::mutate(significance = dplyr::case_when(
        (p.value < 0.05) & (p.value >= 0.01) ~ "*",
        (p.value < 0.01) & (p.value >= 0.001) ~ "**",
        (p.value < 0.001) ~ "***",
        TRUE ~ ""
      ))
    
    return(result)
  }

compare_effect_direction_trans_v_cis_weighted_cis_mr <- function(crp_pqtl_eqtl_df, pval_threshold = 0.05) {
  # Get causal effects of trans-for-CRP proteins on CRP
  crp_outcome_mr_df <- crp_pqtl_eqtl_df |>
    dplyr::filter(Outcome == "CRP") |>
    dplyr::select(-Estimate_crp, -Pvalue_crp, -Outcome) |>
    dplyr::rename_with(\(x) paste0(x, "_cause_crp"), .cols = -Protein)
  
  # Reformat crp_pqtl_eqtl_df
  crp_pqtl_eqtl_df <- crp_pqtl_eqtl_df |>
    dplyr::filter(Outcome != "CRP") |>
    dplyr::full_join(crp_outcome_mr_df, by = "Protein")
  
  # Append col with expected effect dir for trans-weighted (CRP-weighted) cis-MR
  crp_pqtl_eqtl_df <- crp_pqtl_eqtl_df |>
    dplyr::mutate(
      eqtl_expected_effect_dir_crp = sign(Estimate_eqtl) * sign(Estimate_eqtl_cause_crp),
      pqtl_expected_effect_dir_crp = sign(Estimate_pqtl) * sign(Estimate_pqtl_cause_crp)
    )
  
  # GTEx eQTL result -------------------------------------------------------------
  
  eqtl_result <- crp_pqtl_eqtl_df |>
    dplyr::select(
      Outcome,
      Protein,
      Estimate_crp,
      Pvalue_crp,
      Estimate_eqtl,
      Pvalue_eqtl,
      Estimate_eqtl_cause_crp,
      Pvalue_eqtl_cause_crp,
      eqtl_expected_effect_dir_crp
    ) |>
    dplyr::filter(!is.na(eqtl_expected_effect_dir_crp)) |>
    
    # filter for significant causal associations (CRPwt-D, QTLwt-D, QTL-CRP)
    dplyr::filter(
      (Pvalue_crp < !!pval_threshold) &
        (Pvalue_eqtl < !!pval_threshold) &
        (Pvalue_eqtl_cause_crp < !!pval_threshold)
    ) |>
    # append column - got expected effect dir for CRP-weighted result?
    dplyr::mutate(got_expected_result = ifelse(
      eqtl_expected_effect_dir_crp == sign(Estimate_crp),
      yes = "Yes",
      no = "No"
    ))
  
  # DECODE pqtl result -------------------------------------------------------------
  
  pqtl_result <- crp_pqtl_eqtl_df |>
    dplyr::select(
      Outcome,
      Protein,
      Estimate_crp,
      Pvalue_crp,
      Estimate_pqtl,
      Pvalue_pqtl,
      Estimate_pqtl_cause_crp,
      Pvalue_pqtl_cause_crp,
      pqtl_expected_effect_dir_crp
    ) |>
    dplyr::filter(!is.na(pqtl_expected_effect_dir_crp)) |>
    
    # filter for significant causal associations (CRPwt-D, QTLwt-D, QTL-CRP)
    dplyr::filter(
      (Pvalue_crp < !!pval_threshold) &
        (Pvalue_pqtl < !!pval_threshold) &
        (Pvalue_pqtl_cause_crp < !!pval_threshold)
    ) |>
    # append column - got expected effect dir for CRP-weighted result?
    dplyr::mutate(got_expected_result = ifelse(
      pqtl_expected_effect_dir_crp == sign(Estimate_crp),
      yes = "Yes",
      no = "No"
    ))
  
  # Return results as list ---------------------------------------
  
  result <- list(eqtl = eqtl_result, pqtl = pqtl_result)
  
  # Append pval threshold
  result <- result |>
    purrr::map(\(x) x |>
                 dplyr::mutate(sig_threshold = !!pval_threshold))
  
  return(result)
}

create_crp_pqtl_eqtl_df <- function(all_mr_fish_results_df) {
  all_mr_fish_results_df |>
    dplyr::filter(Weight %in% c("CRP Ligthart", "eQTL", "pQTL (DECODE)")) |>
    dplyr::mutate(
      Weight = dplyr::case_when(
        Weight == "CRP Ligthart" ~ "crp",
        Weight == "eQTL" ~ "eqtl",
        Weight == "pQTL (DECODE)" ~ "pqtl"
      )
    ) |>
    dplyr::select(Weight, Protein, Outcome, Estimate, Pvalue) |>
    
    tidyr::pivot_wider(names_from = Weight,
                       values_from = c(Estimate, Pvalue))
}


# sTable 1 ----------------------------------------------------------------

create_table_s1 <- function(all_mr_fish_results_df) {
  create_stable_cis_v_gw(
    all_mr_fish_results_df,
    genomewide_type = "CRP (genome-wide)",
    return_flextable = TRUE
  )
}


# sTable2 -----------------------------------------------------------------

create_table_s2 <- function(all_mr_fish_results_df) {
  create_stable_cis_v_gw(
    all_mr_fish_results_df,
    genomewide_type = "CRP (genome-wide, optimised)",
    return_flextable = TRUE
  )
}

create_stable_cis_v_gw <- function(all_mr_fish_results_df, genomewide_type, return_flextable = TRUE) {
  result <- create_cis_v_gw_table(all_mr_fish_results_df, genomewide_type) |>
    dplyr::select(-Weight, -`Phenotype category_Cis`) |>
    dplyr::rename(Category = `Phenotype category_Genomewide`) |>
    dplyr::mutate(phenotype_units = as.character(stringr::str_glue("{Outcome} ({unit_Cis})"))) |>
    dplyr::select(
      Category,
      `Phenotype (units)` = phenotype_units,
      `Genome-wide` = Genomewide,
      `*cis*-CRP` = `Cis-CRP`,
      `Interaction p value` = cis_vs_gw_pval_sci_not,
      Concordance,
      `Genome-wide Estimator` = Estimator_Genomewide,
      `*cis*-CRP Estimator` = Estimator_Cis,
      `Genome-wide no. variants` = `No. variants_Genomewide`,
      `*cis*-CRP no. variants` = `No. variants_Cis`,
      `Genome-wide F statistic` = `F statistic_Genomewide`,
      `*cis*-CRP F statistic` = `F statistic_Cis`
    ) |>
    dplyr::arrange(Category, `Phenotype (units)`) |>
    # dplyr::rename_with(\(x) paste0("**", x, "**")) |>
    flextable::as_grouped_data(groups = "Category")
  
  if (return_flextable) {
    result <- result |>
      flextable::flextable() |>
      ftExtra::colformat_md(part = "all") |>
      flextable::bold(part = "header") |>
      flextable::add_footer_lines(
        "Estimates in italics are nominally significant (p<0.05) and those in bold are significant after adjusting for multiple testing."
      ) |>
      flextable::autofit()
  }
  
  return(result)
}


# sTable3 and sTable 4 -----------------------------------------------------------------

create_table_s3_and_s4 <- function(all_mr_fish_results_df,
                                   return_flextable = TRUE) {
  message("Retrieving Ensembl IDs...")
  genes <- lookup_gene_symbols(gene_symbols = unique(all_mr_fish_results_df$Protein))
  
  message("Querying OpenTargets...")
  
  ot_data_all <- list(
    known_drugs = ot_knowndrugs(ensgIds = na.omit(genes$id)),
    tractability = ot_target_tractability(ensgIds = na.omit(genes$id))
  )
  
  gene_symbols_ids <- genes %>%
    dplyr::select(approvedSymbol = display_name, ensgId = id)
  
  selected_druggable_categories <- c(
    "approved_drug",
    "advanced_clinical",
    "phase_1_clinical",
    "high_quality_ligand",
    "high_quality_pocket",
    "druggable_family",
    "uni_prot_loc_high_conf",
    "go_cc_high_conf"
  )
  
  tractability <- ot_data_all$tractability %>%
    dplyr::mutate(druggable_target_score = rowSums(dplyr::across(
      tidyselect::all_of(selected_druggable_categories)
    ), na.rm = TRUE))
  
  tractable_targets <- tractability %>%
    dplyr::filter(druggable_target_score > 0) %>%
    dplyr::pull(approvedSymbol) %>%
    unique()
  
  known_drugs <- ot_data_all$known_drugs %>%
    dplyr::left_join(gene_symbols_ids, by = "ensgId") %>%
    dplyr::relocate(approvedSymbol, .after = ensgId)
  
  if ("drug.name" %in% names(known_drugs)) {
    existing_drug <- known_drugs %>%
      dplyr::filter(!is.na(drug.name)) %>%
      dplyr::pull(approvedSymbol) %>%
      unique()
  } else {
    existing_drug <- character()
  }
  
  druggable_or_drugged <- list(
    existing_drug = sort(existing_drug),
    tractable_targets = sort(tractable_targets),
    union = sort(union(existing_drug, tractable_targets))
  )
  
  tab_druggable_genes <- ot_data_all$known_drugs %>%
    dplyr::full_join(gene_symbols_ids, by = "ensgId") %>%
    dplyr::relocate(approvedSymbol, .after = ensgId) |>
    # dplyr::filter(!is.na(drug.name)) %>%
    dplyr::group_by(approvedSymbol, drug.name) %>%
    dplyr::mutate(phase = dplyr::case_when(is.na(phase) ~ -1L, TRUE ~ phase)) %>%
    dplyr::mutate(max_phase = max(phase, na.rm = TRUE)) %>%
    dplyr::slice(which.max(phase)) %>%
    dplyr::ungroup() %>%
    dplyr::select(approvedSymbol,
                  drug.name,
                  max_phase,
                  status,
                  drugType,
                  mechanismOfAction) %>%
    dplyr::arrange(approvedSymbol,
                   drug.name,
                   max_phase,
                   status,
                   drugType,
                   mechanismOfAction) %>%
    dplyr::bind_rows(tibble::tibble(
      approvedSymbol = subset(
        druggable_or_drugged$tractable_targets,
        !druggable_or_drugged$tractable_targets %in% druggable_or_drugged$existing_drug
      )
    )) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), \(x) ifelse(is.na(x), "-", x))) |>
    dplyr::arrange(approvedSymbol) %>%
    dplyr::mutate(max_phase = dplyr::case_when(max_phase == -1L ~ "-", TRUE ~ as.character(max_phase))) |>
    dplyr::distinct() |>
    dplyr::filter(approvedSymbol %in% druggable_or_drugged$union)
  
  tab_druggable_genes_compounds <- ot_data_all$known_drugs %>%
    dplyr::left_join(gene_symbols_ids, by = "ensgId") |>
    dplyr::filter(approvedSymbol %in% gene_symbols_ids$approvedSymbol) %>%
    dplyr::filter(!is.na(disease.name)) %>%
    dplyr::group_by(drug.name) %>%
    dplyr::summarise(target_disease = paste(sort(unique(disease.name)), sep = "", collapse = "; ")) %>%
    dplyr::rename(`Compound name` = drug.name, `Target disease` = target_disease)
  
  result <- list(tab_druggable_genes = tab_druggable_genes,
                 tab_druggable_genes_compounds = tab_druggable_genes_compounds)
  
  if (return_flextable) {
    result$tab_druggable_genes <- result$tab_druggable_genes |>
      flextable::as_grouped_data(groups = "approvedSymbol") |>
      flextable::flextable() |>
      flextable::bold(j = 1) |>
      flextable::autofit()
    
    result$tab_druggable_genes_compounds <- result$tab_druggable_genes_compounds |>
      flextable::flextable() |>
      flextable::theme_zebra() |>
      flextable::autofit()
  }
  
  return(result)
}

# Miscellaneous -------------------------------------------------------------------

remove_optimised_genomewide_crp_mr <- function(all_mr_fish_results_df) {
  all_mr_fish_results_df |>
    dplyr::filter(Protein != "CRP (genome-wide, optimised)")
}

#' Get drugs known to target a given EnsemblId product
#'
#' @param ensgIds Character vector of EnsemblIDs
#' @param size Maximum number of rows to return per EnsemblID queried
#'
#' @return A tibble
#' @export
#' @examples
#' 
#' ot_knowndrugs("ENSG00000157764")
#' 
ot_knowndrugs <- function(ensgIds,
                          size = 1000) {
  
  ot_knowndrugs_single_possibly <- purrr::possibly(ot_knowndrugs_single,
                                                   otherwise = tibble::tibble(ensgId = NA))
  
  ensgIds |>
    purrr::set_names() |>
    purrr::map(\(x) ot_knowndrugs_single_possibly(ensgId = x, size = size), .progress = TRUE) |>
    dplyr::bind_rows(.id = "ensgId")
}

#' Get tractability data for EnsemblIDs
#'
#' Returns a tibble of available [tractability
#' data](https://platform-docs.opentargets.org/target/tractability) for one or
#' more EnsemblIDs.
#'
#' Each queried EnsemblID will return 4 rows of data, one for each [assessment
#' modality](https://platform-docs.opentargets.org/target/tractability#assessments):
#'
#' - 'SM' = Small molecule
#' - 'AB' = Antibody
#' - 'PR' = PROTAC
#' - 'OC' = Other modalities (*see [example](https://platform.opentargets.org/target/ENSG00000169083) on OpenTargets web browser*)
#'
#' @inheritParams ot_knowndrugs
#'
#' @return A tibble
#' @export
#' @examples
#' # get tractability data for 2 targets (note ENSG00000162711 does not have tractability data)
#' result <- ot_target_tractability(c("ENSG00000157764", "ENSG00000122482", "ENSG00000162711"))
#'
#' # see if these are 'druggable', as per the pipeline by Finan et al
#' result |>
#'   dplyr::filter(modality == "SM") |>
#'   dplyr::select(ensgId, druggable_family)
ot_target_tractability <- function(ensgIds) {
  ensgIds |>
    purrr::set_names() |>
    purrr::map(\(x) ot_target_tractability_single(ensgId = x), .progress = TRUE) |>
    dplyr::bind_rows() |>
    dplyr::rename("ensgId" = "id")
}

ot_target_tractability_single <- function(ensgId) {
  result <- submit_ot_query(
    variables = list("ensgId" = ensgId),
    query = "
      query tractability(
        $ensgId: String!
      ) {
        target(ensemblId: $ensgId){
          id
          approvedSymbol
          biotype
          tractability {
            label
            modality
            value
          }
        }
      }
      "
  )
  
  x <- result$data$target
  
  result <- x[1:3] |>
    tibble::as_tibble()
  
  if (!rlang::is_empty(x$tractability)) {
    tractability <- x$tractability |>
      dplyr::bind_rows() |>
      tidyr::pivot_wider(names_from = "label", 
                         values_from = "value") |>
      janitor::clean_names()
    
    result <- result |>
      dplyr::bind_cols(tractability)
  }
  
  result
}

ot_knowndrugs_single <- function(ensgId,
                                 size = 1000) {
  result <- submit_ot_query(
    variables = list("ensgId" = ensgId, "size" = size),
    query = "
      query KnownDrugsQuery(
        $ensgId: String!
        $cursor: String
        $freeTextQuery: String
        $size: Int = 10
      ) {
        target(ensemblId: $ensgId) {
          id
          knownDrugs(cursor: $cursor, freeTextQuery: $freeTextQuery, size: $size) {
            count
            cursor
            rows {
              phase
              status
              urls {
                name
                url
              }
              disease {
                id
                name
              }
              drug {
                id
                name
                mechanismsOfAction {
                  rows {
                    actionType
                    targets {
                      id
                    }
                  }
                }
              }
              drugType
              mechanismOfAction
            }
          }
        }
      }
    
      "
  )
  
  result$data$target$knownDrugs$rows |>
    purrr::map(\(x) {
      if (is.null(x$status)) {
        x$status <- "Unknown"
      }
      
      x
    }) |>
    purrr::map(as.data.frame) |>
    dplyr::bind_rows() |>
    tidyr::unite(col = "urls.name", 
                 tidyselect::starts_with("urls.name"), 
                 sep = "|", 
                 remove = TRUE, 
                 na.rm = TRUE) |>
    tidyr::unite(col = "urls.url", 
                 tidyselect::starts_with("urls.url"), 
                 sep = "|", 
                 remove = TRUE, 
                 na.rm = TRUE) |>
    tidyr::unite(col = "drug.mechanismsOfAction", 
                 tidyselect::starts_with("drug.mechanismsOfAction"), 
                 sep = "|", 
                 remove = TRUE, 
                 na.rm = TRUE) |>
    tibble::as_tibble()
}


ot_associatedDiseases <- function(ensemblId) {
  submit_ot_query(
    variables = list("ensemblId" = ensemblId),
    query = "
  query target($ensemblId: String!){
    target(ensemblId: $ensemblId){
      id
      approvedSymbol
      biotype
      associatedDiseases {
        count
        rows {
          disease {
            id
            name
          }
          datasourceScores {
            id
            score
          }
        }
      }
      geneticConstraint {
        constraintType
        exp
        obs
        score
        oe
        oeLower
        oeUpper
      }
      tractability {
        label
        modality
        value
      }
    }
  }
"
  )
}

#' Base function for submitting queries to OpenaTargets GraphQL APi
#'
#' See  [OpenTargets GraphQL API
#' documentation](https://platform-docs.opentargets.org/data-access/graphql-api),
#' the [API playground](https://platform.opentargets.org/api) and the online
#' browser, which includes option to download results or retrieve the API query
#' that will retrieve these
#' e.g.[ENSG00000169083](https://platform.opentargets.org/target/ENSG00000169083)
#'
#' @param query A query string
#' @param variables A named list of variables, corresponding to `query`
#'
#' @return A response in JSON format.
#' @noRd
submit_ot_query <- function(variables,
                            query) {
  # Set base URL of GraphQL API endpoint
  base_url <- "https://api.platform.opentargets.org/api/v4/graphql"
  
  # Construct POST request body object with query string and variables
  post_body <- list(query = query, variables = variables)
  
  # Perform POST request
  httr::POST(url = base_url,
             body = post_body,
             encode = 'json') |>
    httr::content()
}

#' Request info for a gene symbol
#'
#' Queries Ensembl REST API (build [38](https://rest.ensembl.org/) or
#' [37](https://grch37.rest.ensembl.org/)). Includes gene coordinates,
#' description and biotype (e.g. whether protein coding or not).
#'
#' @param gene_symbols A character vector of gene symbols e.g. "ESPN"
#' @param build Either 'GRCh38' (default) or 'GRCh37'.
#'
#' @return A tibble.
#' @export
#' @examples
#'  lookup_gene_symbols(c("ARMS2", "CFH"))
lookup_gene_symbols <- function(gene_symbols,
                                build = "GRCh38") {
  gene_symbols |>
    purrr::map(\(x) lookup_gene_symbol_single(gene_symbol = x,
                                              build = build), 
               .progress = TRUE) |>
    dplyr::bind_rows()
}

lookup_gene_symbol_single <- function(gene_symbol,
                                      build) {
  req <- req_build_ensembl(
    build = build,
    ext = paste0("lookup/symbol/homo_sapiens/", gene_symbol, "?")
  )
  
  tryCatch(
    expr = req |>
      httr2::req_perform() |>
      httr2::resp_body_json() |>
      tibble::as_tibble(),
    error = function(e)
      tibble::tibble(display_name = gene_symbol,
                     assembly_name = build)
  )
}

req_build_ensembl <- function(build,
                              ext) {
  # validate args
  match.arg(build,
            c("GRCh38", "GRCh37"))
  
  server <- switch(build,
                   GRCh38 = "https://rest.ensembl.org/",
                   GRCh37 = "https://grch37.rest.ensembl.org/")
  
  # define REST query to get the gene ID from the gene name
  req <- httr2::request(paste0(server, ext))
  
  req |>
    httr2::req_headers("Accept" = "application/json")
}