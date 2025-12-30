#' Generate Neo4j LOAD CSV script using local files from scSignalMap pipeline
#'
#' Automatically finds:
#' - The three core files in Neo4J/ (*_senders_ligands.csv, etc.)
#' - All filtered enrichr files in Neo4J/ (*_enrichr_results_DATABASE2.csv)
#' and generates a complete Cypher script using file:/// paths.
#'
#' @param neo4j_dir Directory with all CSV files (default: "Neo4J/")
#' @param output_file Where to save the .cypher script (default: "Neo4J/load_scSignalMap.cypher")
#' @return Invisibly returns path to generated script
#' @export
generate_neo4j_local_load_script = function(
  neo4j_dir = "Neo4J/",
  output_file = file.path(neo4j_dir, "load_scSignalMap_local.cypher")
) {
  # List all CSV files
  all_files = list.files(neo4j_dir, pattern = "\\.csv$", full.names = FALSE)

  # Identify core files
  sender_file = all_files[grep("_senders_ligands\\.csv$", all_files)]
  lr_file     = all_files[grep("_ligands_receptor_pairs\\.csv$", all_files)]
  receiver_file = all_files[grep("_receptors_receivers\\.csv$", all_files)]
  pathway_files = all_files[grep("_enrichr_results_DATABASE2\\.csv$", all_files)]

  if (length(sender_file) == 0 || length(lr_file) == 0 || length(receiver_file) == 0) {
    stop("Missing one or more core files. Found:\n", paste(all_files, collapse = "\n"))
  }

  cypher = c(
    "// ================================================",
    "// Auto-generated Neo4j import script for scSignalMap",
    "// Generated on: ", Sys.Date(),
    "// All files are in one folder - copy everything to Neo4j import/",
    "// ================================================\n",
    "// Indexes",
    "CREATE INDEX IF NOT EXISTS FOR (s:Sender) ON (s.name);",
    "CREATE INDEX IF NOT EXISTS FOR (l:Ligand_Symbol) ON (l.name);",
    "CREATE INDEX IF NOT EXISTS FOR (p:Receptor_Symbol) ON (p.name);",
    "CREATE INDEX IF NOT EXISTS FOR (c:Receiver) ON (c.name);",
    "CREATE INDEX IF NOT EXISTS FOR (z:Signal_Pathway_Info) ON (z.name);\n",
    "// 1. Sender → Ligand_Symbol",
    paste0('LOAD CSV WITH HEADERS FROM "file:///', sender_file, '" AS row'),
    "MERGE (s:Sender {name: row.Sender})",
    "MERGE (l:Ligand_Symbol {name: row.Ligand_Symbol})",
    "  ON CREATE SET l += {",
    "    counts: toInteger(coalesce(row.Ligand_Counts, 0)),",
    "    l_exp_lvl_3: toFloat(coalesce(row.Ligand_gte_3, 0.0)),",
    "    l_exp_lvl_10: toFloat(coalesce(row.Ligand_gte_10, 0.0)),",
    "    Ligand_secreted: toBoolean(coalesce(row.Ligand_secreted, false))",
    "  }",
    "MERGE (s)-[r:sender2ligand]->(l)",
    "  ON CREATE SET r.Ligand_Avg_Exp = toFloat(coalesce(row.Ligand_Avg_Exp, 0.0));\n",
    "// 2. Ligand_Symbol to Receptor_Symbol",
    paste0('LOAD CSV WITH HEADERS FROM "file:///', lr_file, '" AS row'),
    "MERGE (l:Ligand_Symbol {name: row.Ligand_Symbol})",
    "MERGE (p:Receptor_Symbol {name: row.Receptor_Symbol})",
    "  ON CREATE SET p += {",
    "    counts: toInteger(coalesce(row.Receptor_Counts, 0)),",
    "    r_exp_lvl_3: toFloat(coalesce(row.Receptor_gte_3, 0.0)),",
    "    r_exp_lvl_10: toFloat(coalesce(row.Receptor_gte_10, 0.0))",
    "  }",
    "MERGE (l)-[:ligand2receptorsymbol]->(p);\n",
    "// 3. Receptor_Symbol to Receiver",
    paste0('LOAD CSV WITH HEADERS FROM "file:///', receiver_file, '" AS row'),
    "MERGE (p:Receptor_Symbol {name: row.Receptor_Symbol})",
    "  ON CREATE SET p.Receptor_Cluster_Marker = toBoolean(coalesce(row.Receptor_Cluster_Marker, false))",
    "MERGE (c:Receiver {name: row.Receiver})",
    "MERGE (p)-[r:receptor2receivecluster]->(c)",
    "  ON CREATE SET r.Receptor_Avg_Exp = toFloat(coalesce(row.Receptor_Avg_Exp, 0.0));"
  )

  if (length(pathway_files) > 0) {
    unwinds = character()
    for (f in pathway_files) {
      # Extract receiver from filename: "Sender_Receiver_enrichr_results_DATABASE2.csv" → Receiver
      receiver = sub(".*_", "", sub("_enrichr_results_DATABASE2\\.csv$", "", f))
      receiver = gsub("[.////]", "", receiver)
      unwinds = c(unwinds, paste0('  {file: "', f, '", receiver: "', receiver, '"}'))
    }

    cypher = c(cypher,
      "\n// 4. Receiver to Signal_Pathway_Info (filtered pathways with linked receptors)",
      "UNWIND [",
      paste(unwinds, collapse = ",\n"),
      "] AS item",
      'LOAD CSV WITH HEADERS FROM "file:///" + item.file AS row',
      "WITH item.receiver AS recvName, row",
      "WHERE row.Term IS NOT NULL",
      "MERGE (c:Receiver {name: recvName})",
      "MERGE (z:Signal_Pathway_Info {name: row.Term})",
      "  ON CREATE SET z += {",
      "    Database: coalesce(row.database, 'Unknown'),",
      "    Pathway_Adjusted_P_Value: toFloat(coalesce(row.Adjusted_P_value, 999)),",
      "    Combined_Score: toFloat(coalesce(row.Combined_Score, 0.0)),",
      "    Matching_Receptors: split(coalesce(row.Matching_Receptors, ''), ';')",
      "  }",
      "MERGE (c)-[r:receiver2signalpathway]->(z)",
      "  ON CREATE SET r.Addl_Linked_Genes = split(coalesce(row.Genes, ''), ';');"
    )
  }

  # Save script
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  writeLines(cypher, output_file)
  message("Neo4j load script generated: ", output_file)
  message("   1. Copy all files from 'Neo4J/' into your Neo4j import/ directory")
  message("   2. Run the script in Neo4j Browser")

  invisible(output_file)
}
