\name{run_post_processing_Neo4J}
\alias{run_post_processing_Neo4J}
\title{Run Post-Processing and Neo4J Export}

\description{
  Performs all post-processing steps after the core scSignalMap analysis:
  \itemize{
    \item Saves detailed per-pair CSV files (DE genes, interactions, enrichment results)
    \item Saves filtered pathway enrichment files for Neo4j import
    \item Exports the three core relationship tables using \code{\link{export_for_neo4j}}
    \item Generates a local Neo4j import script (optional)
    \item Optionally uploads all Neo4j CSVs to Google Drive, makes them public,
          and generates a cloud-compatible import script
  }
}

\usage{
run_post_processing_Neo4J(
  all_results,
  neo4j_prefix = "scSignalMap",
  upload_to_drive = FALSE,
  drive_folder_name = NULL,
  generate_local_script = TRUE,
  generate_cloud_script = NULL
)
}

\arguments{
  \item{all_results}{
    Named list returned by \code{\link{run_full_scSignalMap_pipeline}} (or \code{run_scSignalMap_core}).
    Each element corresponds to a sender-receiver pair and contains all intermediate results.
  }
  \item{neo4j_prefix}{
    Character string used as prefix for the three core Neo4j CSV filenames
    (senders_ligands, ligands_receptor_pairs, receptors_receivers).
    Default: \code{"scSignalMap"}.
  }
  \item{upload_to_drive}{
    Logical; if \code{TRUE}, uploads all CSV files in the \code{Neo4J/} directory
    to Google Drive, sets them to "Anyone with the link" sharing, and prints direct download URLs.
    Requires the \code{googledrive} package and authentication.
    Default: \code{FALSE}.
  }
  \item{drive_folder_name}{
    Name of the Google Drive folder to create or use. If \code{NULL} (default),
    a folder named \code{"scSignalMap_Neo4j_YYYYMMDD"} is created.
  }
  \item{generate_local_script}{
    Logical; if \code{TRUE}, generates a local \code{file:///}-based Cypher script
    (\code{load_scSignalMap_local.cypher}). Default: \code{TRUE}.
  }
  \item{generate_cloud_script}{
    Logical; if \code{TRUE}, generates an HTTPS-based Cypher script for cloud instances
    (\code{load_scSignalMap_cloud.cypher}). Automatically set to \code{TRUE} when
    \code{upload_to_drive = TRUE}. Default: \code{NULL} (inherits from \code{upload_to_drive}).
  }
}

\details{
  This function is the final step in the scSignalMap → Neo4j workflow. It organizes and saves
  all intermediate results for inspection and prepares the data for import into Neo4j.

  Files are saved in two directories:
  \itemize{
    \item \code{enrichr_results/<Sender>_<Receiver>/} – per-pair detailed CSVs
    \item \code{Neo4J/} – core relationship tables + filtered pathway files
  }
}

\value{
  Invisibly returns the input \code{all_results} list (for piping if desired).
}

\examples{
\donttest{
# After running the core analysis
results <- run_full_scSignalMap_pipeline(
  seurat_obj = my_seurat,
  celltype_column = "celltype",
  sender_celltypes = c("Macrophages"),
  receiver_celltypes = c("Tumor"),
  cond_column = "condition",
  cond_name1 = "Treatment",
  cond_name2 = "Control"
)

# Basic export (local only)
run_post_processing_Neo4J(
  all_results = results,
  neo4j_prefix = "my_analysis_20251230",
  generate_local_script = TRUE
)

# Full export with Google Drive upload
run_post_processing_Neo4J(
  all_results = results,
  neo4j_prefix = "my_analysis_20251230",
  upload_to_drive = TRUE,
  drive_folder_name = "MyProject_Neo4j_Results"
)
}
}

\seealso{
  \code{\link{run_full_scSignalMap_pipeline}},
  \code{\link{export_for_neo4j}},
  \code{\link{generate_neo4j_local_load_script}},
  \code{\link{generate_neo4j_cloud_load_script}}
}

\author{
  Andrew Zehnder
}

\keyword{Neo4j}
\keyword{export}
\keyword{post-processing}
