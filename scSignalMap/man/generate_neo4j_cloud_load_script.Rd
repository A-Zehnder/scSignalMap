\name{generate_neo4j_cloud_load_script}
\alias{generate_neo4j_cloud_load_script}
\title{Generate Neo4j LOAD CSV script for cloud instances using HTTPS URLs}

\description{
  Generates a complete Cypher script for importing scSignalMap data into Neo4j cloud deployments
  (e.g., Neo4j Sandbox or AuraDB), where local \code{file:///} paths are not supported.

  Instead of local file paths, this function uses direct public HTTPS URLs pointing to CSV files
  hosted on Google Drive, GitHub releases, Dropbox, or any other publicly accessible location.
}

\usage{
generate_neo4j_cloud_load_script(file_urls,
                                 output_file = "Neo4J/load_scSignalMap_cloud.cypher")
}

\arguments{
  \item{file_urls}{
    A named character vector where the names are the original CSV filenames
    (exactly as produced by \code{\link{export_for_neo4j}} and the filtered enrichr files)
    and the values are direct downloadable HTTPS URLs.

    Required core files (names must match):
    \itemize{
      \item \code{*_senders_ligands.csv}
      \item \code{*_ligands_receptor_pairs.csv}
      \item \code{*_receptors_receivers.csv}
    }

    Optional pathway files (automatically detected):
    \itemize{
      \item Any files ending in \code{_enrichr_results_DATABASE2.csv}
    }

    Example:
    \preformatted{
    file_urls = c(
      "mydata_senders_ligands.csv" = "https://drive.google.com/uc?id=abc123&export=download",
      "mydata_ligands_receptor_pairs.csv" = "https://drive.google.com/uc?id=def456&export=download",
      "mydata_receptors_receivers.csv" = "https://drive.google.com/uc?id=ghi789&export=download",
      "Tcells_Tumor_enrichr_results_DATABASE2.csv" = "https://drive.google.com/uc?id=jkl012&export=download"
    )
    }
  }
  \item{output_file}{
    Path and filename for the generated Cypher script.
    Default: \code{"Neo4J/load_scSignalMap_cloud.cypher"}.
  }
}

\details{
  The generated script includes:
  \itemize{
    \item Index creation for fast lookups
    \item Loading and merging of sender to ligand relationships
    \item Loading and merging of ligand to receptor relationships
    \item Loading and merging of receptor to receiver relationships
    \item Optional loading of filtered pathway enrichment results linked to receivers
  }

  The script can be copied and run directly in the Neo4j Browser of a cloud instance.
}

\value{
  Invisibly returns the path to the generated Cypher script file.
}

\examples{
\donttest{
# Example URLs (replace with your actual public links)
cloud_urls = c(
  "scSignalMap_GB126_senders_ligands.csv" =
    "https://drive.google.com/uc?id=ABC123&export=download",
  "scSignalMap_GB126_ligands_receptor_pairs.csv" =
    "https://drive.google.com/uc?id=DEF456&export=download",
  "scSignalMap_GB126_receptors_receivers.csv" =
    "https://drive.google.com/uc?id=GHI789&export=download",
  "Macrophages_Tumor_enrichr_results_DATABASE2.csv" =
    "https://drive.google.com/uc?id=JKL012&export=download"
)

generate_neo4j_cloud_load_script(
  file_urls = cloud_urls,
  output_file = "Neo4J/my_analysis_cloud.cypher"
)
}
}

\seealso{
  \code{\link{generate_neo4j_local_load_script}},
  \code{\link{export_for_neo4j}},
  \code{\link{run_post_processing_Neo4J}}
}

\author{
  Andrew Zehnder
}

\keyword{Neo4j}
\keyword{import}
\keyword{cloud}
