library(targets)
library(tarchetypes)
library(magrittr)

# functions for recreating final manuscript figures
source(file.path("code", "utilities.R"))

list(
  tar_target(
    SUPPLEMENTARY_XLSX,
    "mr-fish_supplementary_tables.xlsx",
    format = "file"
  ),
  
  tar_target(supplementary_tables,
             SUPPLEMENTARY_XLSX |>
               readxl::excel_sheets() |>
               purrr::set_names() |>
               purrr::map(\(x) readxl::read_excel(SUPPLEMENTARY_XLSX,
                                                  sheet = x, 
                                                  skip = 1))),
  
  tar_target(druggable_proteins,
             supplementary_tables$sTable3 |>
               dplyr::distinct(approvedSymbol) |>
               dplyr::filter(!is.na(approvedSymbol)) |>
               dplyr::pull(approvedSymbol)),
  
  # tar_target(FIGURE_1,
  #            FALSE),
  # 
  # tar_target(FIGURE_2,
  #            FALSE),
  
  tar_target(FIGURE_3,
             plot_figure_3(all_mr_fish_results_df = supplementary_tables$sTable6,
                           file_path = "output/Figure 3.png"),
             format = "file"),
  
  tar_target(FIGURE_4,
             plot_figure_4(all_mr_fish_results_df = supplementary_tables$sTable6,
                           file_path = "output/Figure 4.png"),
             format = "file"),
  
  tar_target(FIGURE_S1,
             plot_figure_s1(supplementary_tables$sTable6,
                            file_path = "output/Figure S1.png"),
             format = "file"),
  
  tar_target(FIGURE_S2,
             plot_figure_s2(supplementary_tables$sTable6, 
                            druggable_proteins, 
                            file_path = "output/Figure S2.png"),
             format = "file"),
  
  tar_target(
    FIGURE_S3,
    plot_figure_s3(supplementary_tables$sTable6, 
                   druggable_proteins, 
                   file_path = "output/Figure S3.png"),
    format = "file"
  ), 
  
  tar_target(FIGURE_S4,
             plot_figure_s4(supplementary_tables$sTable6, 
                            druggable_proteins, 
                            file_path = "output/Figure S4.png"),
             format = "file"
  ), 
  
  tar_target(
    TABLE_S1,
    flextable::save_as_html(
      values = list(sTable1 = create_table_s1(supplementary_tables$sTable6)),
      path = "output/sTable1.html"
    ),
    format = "file"
  ), 
  
  tar_target(
    TABLE_S2,
    flextable::save_as_html(
      values = list(sTable2 = create_table_s2(supplementary_tables$sTable6)),
      path = "output/sTable2.html"
    ),
    format = "file"
  ), 
  
  tar_target(
    table_s3_and_s4,
    create_table_s3_and_s4(supplementary_tables$sTable6)
  ),
  
  tar_target(
    TABLE_S3,
    flextable::save_as_html(
      values = list(sTable3 = table_s3_and_s4$tab_druggable_genes),
      path = "output/sTable3.html"
    ),
    format = "file"
  ), 
  
  tar_target(
    TABLE_S4,
    flextable::save_as_html(
      values = list(sTable4 = table_s3_and_s4$tab_druggable_genes_compounds),
      path = "output/sTable4.html"
    ),
    format = "file"
  ), 
  
  tar_target(FIGURE_A1,
             plot_figure_a1(all_mr_fish_results_df = supplementary_tables$sTable6,
                            file_path = "output/Figure A1.png"),
             format = "file")
)