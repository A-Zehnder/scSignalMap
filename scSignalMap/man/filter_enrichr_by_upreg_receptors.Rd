#' Filter Enrichr results to keep only pathways containing upregulated receptors
#'
#' Takes the full Enrichr results and the upregulated receptors (after intersection with LR pairs)
#' and keeps only those pathways that contain at least one of the receptor gene symbols.
#' It adds a column "Matching_Receptors" with the matching receptor(s) and overwrites/saves
#' the filtered enrichr results.
#'
#' @param enrichr_results Dataframe returned by find_enriched_pathways (combined across databases)
#' @param upreg_receptors_filtered_and_compared Dataframe from intersect_upreg_receptors_with_lr_interactions
#' @param directory Directory where the filtered file should be saved
#' @param sender_clean Cleaned sender cell type name (for filename)
#' @param receiver_clean Cleaned receiver cell type name (for filename)
#' @param suffix Suffix for the output file, default "_enrichr_results_filtered.csv"
#' @return The filtered enrichr_results dataframe (also saved to disk)
#' @export
filter_enrichr_by_upreg_receptors = function(enrichr_results, upreg_receptors_filtered_and_compared) {

  # Safety checks
  if (is.null(enrichr_results) || nrow(enrichr_results) == 0) {
    message("No enrichr results to filter.")
    return(NULL)
  }
  if (is.null(upreg_receptors_filtered_and_compared) || nrow(upreg_receptors_filtered_and_compared) == 0) {
    message("No upregulated receptors – skipping enrichr filtering.")
    return(NULL)
  }

  # Standardise column names (in case Adjusted.P.value vs Adjusted_P_value)
  if ("Adjusted.P.value" %in% colnames(enrichr_results)) {
    enrichr_results = enrichr_results %>%
      dplyr::rename(Adjusted_P_value = Adjusted.P.value)
  }
  if ("Combined.Score" %in% colnames(enrichr_results)) {
    enrichr_results = enrichr_results %>%
      dplyr::rename(Combined_Score = Combined.Score)
  }

  # Unique receptor symbols (case-insensitive matching)
  receptor_symbols = tolower(unique(upreg_receptors_filtered_and_compared$gene_symbol))

  # Split Genes column into list and lowercase for matching
  enrichr_results = enrichr_results %>%
    dplyr::mutate(Genes_list = lapply(strsplit(Genes, ";"), function(x) trimws(tolower(x))))

  # Find matching receptors (original case)
  enrichr_results = enrichr_results %>%
    dplyr::mutate(Matching_Receptors = sapply(Genes_list, function(pathway_genes) {
      matches = intersect(pathway_genes, receptor_symbols)
      if (length(matches) == 0) return("")
      # Map back to original case
      orig = upreg_receptors_filtered_and_compared$gene_symbol[
        tolower(upreg_receptors_filtered_and_compared$gene_symbol) %in% matches
      ]
      paste(unique(orig), collapse = ";")
    })) %>%
    dplyr::select(-Genes_list)

  # Keep only pathways with at least one matching receptor
  enrichr_filtered = enrichr_results %>%
    dplyr::filter(Matching_Receptors != "" & !is.na(Matching_Receptors))
  message("Filtered to ", nrow(enrichr_filtered), " pathways containing upregulated receptors.")
  return(enrichr_filtered)
}
