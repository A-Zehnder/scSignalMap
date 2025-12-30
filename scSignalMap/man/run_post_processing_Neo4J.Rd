#' Run Post-Processing & Neo4J Export
#'
#' Takes the output from run_scSignalMap_core() and saves all CSV files,
#' exports core Neo4J tables, and generates load scripts.
#'
#' @param all_results List returned by run_full_scSignalMap_pipeline()
#' @param neo4j_prefix Prefix for core Neo4J CSV filenames (default: "scSignalMap")
#' @param upload_to_drive Logical; upload Neo4J CSVs to Google Drive and generate cloud script
#' @param drive_folder_name Name of Google Drive folder (default: auto-generated with date)
#' @param generate_local_script Logical; create local file:/// load script
#' @param generate_cloud_script Logical; create HTTPS cloud load script (auto-enabled if upload_to_drive=TRUE)
#' @export
run_post_processing_Neo4J = function(
    all_results,
    neo4j_prefix = "scSignalMap",
    upload_to_drive = FALSE,
    drive_folder_name = NULL,
    generate_local_script = TRUE,
    generate_cloud_script = NULL
) {
  if (is.null(generate_cloud_script)) {
    generate_cloud_script = upload_to_drive
  }
  
  # Ensure directories exist
  dir.create("enrichr_results", showWarnings = FALSE, recursive = TRUE)
  dir.create("Neo4J", showWarnings = FALSE, recursive = TRUE)
  
  # Grab the global LR_interactions (same for all pairs)
  LR_interactions = all_results[[1]]$LR_interactions
  
  message("Saving per-pair CSV files and filtered enrichr results...")
  
  for (pair_name in names(all_results)) {
    res = all_results[[pair_name]]
    sender_clean   = res$sender_clean
    receiver_clean = res$receiver_clean
    
    dir_pair = file.path("enrichr_results", paste0(sender_clean, "_", receiver_clean))
    dir.create(dir_pair, showWarnings = FALSE, recursive = TRUE)
    
    # Save individual results
    write.csv(res$LR_interactions, file.path(dir_pair, paste0(sender_clean, "_", receiver_clean, "_LR_interactions2.csv")), row.names = FALSE)
    write.csv(res$de_cond_celltype, file.path(dir_pair, paste0(sender_clean, "_", receiver_clean, "_de_cond_celltype2.csv")), row.names = FALSE)
    write.csv(res$upreg_receptors, file.path(dir_pair, paste0(sender_clean, "_", receiver_clean, "_upreg_receptors2.csv")), row.names = FALSE)
    write.csv(res$interactions_filtered, file.path(dir_pair, paste0(sender_clean, "_", receiver_clean, "_interactions_filtered2.csv")), row.names = FALSE)
    write.csv(res$upreg_receptors_filtered_and_compared, file.path(dir_pair, paste0(sender_clean, "_", receiver_clean, "_upreg_receptors_filtered_and_compared2.csv")), row.names = FALSE)
    write.csv(res$enrichr_results, file.path(dir_pair, paste0(sender_clean, "_", receiver_clean, "_enrichr_results2.csv")), row.names = FALSE)
    
    # Save filtered enrichr for Neo4J (only if non-empty)
    if (!is.null(res$enrichr_filtered) && nrow(res$enrichr_filtered) > 0) {
      neo4j_file = file.path("Neo4J",
                             paste0(sender_clean, "_", receiver_clean, "_enrichr_results_DATABASE2.csv"))
      write.csv(res$enrichr_filtered, neo4j_file, row.names = FALSE)
      message("Saved ", nrow(res$enrichr_filtered),
              " signaling-relevant pathways for ", pair_name)
    } else {
      message("No signaling-relevant pathways for ", pair_name, " (skipped)")
    }
  }
  
  # Export the three core Neo4J tables (from full LR_interactions)
  export_for_neo4j(
    interactions = LR_interactions,
    output_dir = "Neo4J/",
    prefix = neo4j_prefix
  )
  
  # Generate scripts
  if (generate_local_script) {
    generate_neo4j_local_load_script(
      neo4j_dir = "Neo4J/",
      output_file = "Neo4J/load_scSignalMap_local.cypher"
    )
  }
  
  if (upload_to_drive) {
    if (!requireNamespace("googledrive", quietly = TRUE)) {
      stop("Package 'googledrive' needed for upload_to_drive = TRUE.")
    }
    library(googledrive)
    
    message("\n=== Uploading Neo4J files to Google Drive ===")
    drive_auth()
    
    if (is.null(drive_folder_name)) {
      drive_folder_name = paste0("scSignalMap_Neo4j_", format(Sys.Date(), "%Y%m%d"))
    }
    
    drive_folder = try(drive_get(drive_folder_name), silent = TRUE)
    if (inherits(drive_folder, "try-error") || nrow(drive_folder) == 0) {
      drive_folder = drive_mkdir(drive_folder_name)
      message("Created folder: ", drive_folder_name)
    } else {
      drive_folder = drive_folder[1, ]
    }
    
    local_files = list.files("Neo4J/", pattern = "\\.csv$", full.names = TRUE)
    if (length(local_files) == 0) {
      warning("No CSV files in Neo4J/ to upload.")
    } else {
      message("Uploading ", length(local_files), " files...")
      uploaded = lapply(local_files, function(f) drive_upload(media = f, path = drive_folder, overwrite = TRUE))
      uploaded = do.call(rbind, uploaded)
      
      drive_share(uploaded, role = "reader", type = "anyone")
      
      file_urls = setNames(
        paste0("https://drive.google.com/uc?id=", uploaded$id, "&export=download"),
        uploaded$name
      )
      
      message("Direct download URLs:")
      print(file_urls)
      
      if (generate_cloud_script) {
        generate_neo4j_cloud_load_script(
          file_urls = file_urls,
          output_file = "Neo4J/load_scSignalMap_cloud.cypher"
        )
        message("Cloud load script generated: Neo4J/load_scSignalMap_cloud.cypher")
      }
    }
  }
  
  message("Post-processing and Neo4J export complete.")
  invisible(all_results)
}                   
