#' Export scSignalMap results in Neo4J-friendly format
#'
#' Takes the full MapInteractions output and creates three clean CSV files:
#' 1. senders_ligands: One row per sender-ligand with expression stats
#' 2. ligands_receptor_pairs: One row per valid L-R pair with stats
#' 3. receptors_receivers: One row per receptor-receiver with expression stats
#'
#' These are ideal for importing into Neo4j as nodes and relationships.
#'
#' @param interactions Dataframe from MapInteractions() (full pairs_data)
#' @param output_dir Directory to save the three CSV files (default: "Neo4J/")
#' @param prefix Prefix for filenames (e.g., "Q1_norm_test")
#' @return Invisibly returns list of file paths created
#' @export
export_for_neo4j = function(
  interactions,
  output_dir = "Neo4J/",
  prefix = "scSignalMap"
) {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # 1. Sender - Ligand
  sender_ligands = interactions %>%
    dplyr::select(
      Sender,
      Ligand_Symbol,
      Ligand_Cells_Exp,
      Ligand_Avg_Exp,
      Ligand_secreted,
      Ligand_Counts,
      Ligand_gte_3,
      Ligand_gte_10
    ) #%>%
    #dplyr::distinct()

  sender_file = file.path(output_dir, paste0(prefix, "_senders_ligands.csv"))
  write.csv(sender_ligands, sender_file, row.names = FALSE)
  message("Saved senders to ligands: ", sender_file)

  # 2. Ligand - Receptor pairs
  lr_pairs = interactions %>%
    dplyr::select(
      Ligand_Symbol,
      Ligand_Counts,
      Ligand_gte_3,
      Ligand_gte_10,
      Receptor_Symbol,
      Receptor_Counts,
      Receptor_gte_3,
      Receptor_gte_10
    ) #%>%
    #dplyr::distinct()

  lr_file = file.path(output_dir, paste0(prefix, "_ligands_receptor_pairs.csv"))
  write.csv(lr_pairs, lr_file, row.names = FALSE)
  message("Saved L-R pairs: ", lr_file)

  # 3. Receptor - Receiver
  receiver_receptors = interactions %>%
    dplyr::select(
      Receptor_Symbol,
      Receptor_Counts,
      Receptor_gte_3,
      Receptor_gte_10,
      Receptor_Cluster_Marker,
      Receiver,
      Receptor_Avg_Exp
    ) #%>%
    #dplyr::distinct()

  receiver_file = file.path(output_dir, paste0(prefix, "_receptors_receivers.csv"))
  write.csv(receiver_receptors, receiver_file, row.names = FALSE)
  message("Saved receptors to receivers: ", receiver_file)

  invisible(list(
    sender_ligands = sender_file,
    ligand_receptor_pairs = lr_file,
    receiver_receptors = receiver_file
  ))
}
