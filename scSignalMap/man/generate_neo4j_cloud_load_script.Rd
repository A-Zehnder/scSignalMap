#' Generate Neo4j LOAD CSV script for cloud instances (Sandbox/Aura) using HTTPS URLs
#'
#' This version is for Neo4j cloud deployments (Sandbox, AuraDB) where local file:/// paths are not supported.
#' Instead, provide direct public HTTPS URLs to your CSV files hosted on Google Drive (or GitHub, Dropbox, etc.).
#'
#' @param file_urls A named vector where names are the original filenames (as expected by the script)
#'                  and values are the direct HTTPS URLs (e.g., Google Drive download links).
#'                  Example: c("my_senders_ligands.csv" = "https://drive.google.com/uc?id=abc123&export=download", ...)
#' @param output_file Path to save the generated .cypher script (default: "Neo4J/load_scSignalMap_cloud.cypher")
#' @return Invisibly returns the path to the generated script
#' @export
generate_neo4j_cloud_load_script = function(
  file_urls,
  output_file = "Neo4J/load_scSignalMap_cloud.cypher"
) {
  # Validate that required core files are provided
  required_core = c(
    grep("_senders_ligands\\.csv$", names(file_urls), value = TRUE),
    grep("_ligands_receptor_pairs\\.csv$", names(file_urls), value = TRUE),
    grep("_receptors_receivers\\.csv$", names(file_urls), value = TRUE)
  )
  if (length(required_core) < 3) {
    stop("Missing one or more core files in file_urls. Need: *_senders_ligands.csv, *_ligands_receptor_pairs.csv, *_receptors_receivers.csv")
  }
  
  sender_url = file_urls[grep("_senders_ligands\\.csv$", names(file_urls))]
  lr_url     = file_urls[grep("_ligands_receptor_pairs\\.csv$", names(file_urls))]
  receiver_url = file_urls[grep("_receptors_receivers\\.csv$", names(file_urls))]
  
  # Pathway files (optional)
  pathway_urls = file_urls[grep("_enrichr_results_DATABASE2\\.csv$", names(file_urls))]
  
  cypher = c(
    "// ================================================",
    "// Auto-generated Neo4j import script for scSignalMap (CLOUD VERSION)",
    "// Generated on: ", Sys.Date(),
    "// Uses direct HTTPS URLs (e.g., Google Drive direct download links)",
    "// Run this directly in Neo4j Browser on your Sandbox/Aura instance",
    "// ================================================\n",
    "// Indexes",
    "CREATE INDEX IF NOT EXISTS FOR (s:Sender) ON (s.name);",
    "CREATE INDEX IF NOT EXISTS FOR (l:Ligand_Symbol) ON (l.name);",
    "CREATE INDEX IF NOT EXISTS FOR (p:Receptor_Symbol) ON (p.name);",
    "CREATE INDEX IF NOT EXISTS FOR (c:Receiver) ON (c.name);",
    "CREATE INDEX IF NOT EXISTS FOR (z:Signal_Pathway_Info) ON (z.name);\n",
    "// 1. Sender to Ligand_Symbol",
    paste0('LOAD CSV WITH HEADERS FROM "', sender_url, '" AS row'),
    "MERGE (s:Sender {name: row.Sender})",
    "MERGE (l:Ligand_Symbol {name: row.Ligand_Symbol})",
    " ON CREATE SET l += {",
    " counts: toInteger(coalesce(row.Ligand_Counts, 0)),",
    " l_exp_lvl_3: toFloat(coalesce(row.Ligand_gte_3, 0.0)),",
    " l_exp_lvl_10: toFloat(coalesce(row.Ligand_gte_10, 0.0)),",
    " Ligand_secreted: toBoolean(coalesce(row.Ligand_secreted, false))",
    " }",
    "MERGE (s)-[r:sender2ligand]->(l)",
    " ON CREATE SET r.Ligand_Avg_Exp = toFloat(coalesce(row.Ligand_Avg_Exp, 0.0));\n",
    "// 2. Ligand_Symbol to Receptor_Symbol",
    paste0('LOAD CSV WITH HEADERS FROM "', lr_url, '" AS row'),
    "MERGE (l:Ligand_Symbol {name: row.Ligand_Symbol})",
    "MERGE (p:Receptor_Symbol {name: row.Receptor_Symbol})",
    " ON CREATE SET p += {",
    " counts: toInteger(coalesce(row.Receptor_Counts, 0)),",
    " r_exp_lvl_3: toFloat(coalesce(row.Receptor_gte_3, 0.0)),",
    " r_exp_lvl_10: toFloat(coalesce(row.Receptor_gte_10, 0.0))",
    " }",
    "MERGE (l)-[:ligand2receptorsymbol]->(p);\n",
    "// 3. Receptor_Symbol to Receiver",
    paste0('LOAD CSV WITH HEADERS FROM "', receiver_url, '" AS row'),
    "MERGE (p:Receptor_Symbol {name: row.Receptor_Symbol})",
    " ON CREATE SET p.Receptor_Cluster_Marker = toBoolean(coalesce(row.Receptor_Cluster_Marker, false))",
    "MERGE (c:Receiver {name: row.Receiver})",
    "MERGE (p)-[r:receptor2receivecluster]->(c)",
    " ON CREATE SET r.Receptor_Avg_Exp = toFloat(coalesce(row.Receptor_Avg_Exp, 0.0));"
  )
  
  if (length(pathway_urls) > 0) {
    unwinds = character()
    for (i in seq_along(pathway_urls)) {
      f_name = names(pathway_urls)[i]
      url = pathway_urls[i]
      # Extract receiver from filename (same logic as original)
      receiver = sub(".*_", "", sub("_enrichr_results_DATABASE2\\.csv$", "", f_name))
      receiver = gsub("[.////]", "", receiver)
      unwinds = c(unwinds, paste0(' {url: "', url, '", receiver: "', receiver, '"}'))
    }
    cypher = c(cypher,
      "\n// 4. Receiver to Signal_Pathway_Info (filtered pathways)",
      "UNWIND [",
      paste(unwinds, collapse = ",\n"),
      "] AS item",
      'LOAD CSV WITH HEADERS FROM item.url AS row',
      "WITH item.receiver AS recvName, row",
      "WHERE row.Term IS NOT NULL",
      "MERGE (c:Receiver {name: recvName})",
      "MERGE (z:Signal_Pathway_Info {name: row.Term})",
      " ON CREATE SET z += {",
      " Database: coalesce(row.database, 'Unknown'),",
      " Pathway_Adjusted_P_Value: toFloat(coalesce(row.Adjusted_P_value, 999)),",
      " Combined_Score: toFloat(coalesce(row.Combined_Score, 0.0)),",
      " Matching_Receptors: split(coalesce(row.Matching_Receptors, ''), ';')",
      " }",
      "MERGE (c)-[r:receiver2signalpathway]->(z)",
      " ON CREATE SET r.Addl_Linked_Genes = split(coalesce(row.Genes, ''), ';');"
    )
  }
  
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  writeLines(cypher, output_file)
  message("Cloud Neo4j load script generated: ", output_file)
  message("Run this script directly in your Neo4j Sandbox/Aura Browser.")
  invisible(output_file)
}
