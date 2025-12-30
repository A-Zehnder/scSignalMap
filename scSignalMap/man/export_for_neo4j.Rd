\name{export_for_neo4j}
\alias{export_for_neo4j}
\title{Export scSignalMap results in Neo4J-friendly format}

\description{
  Takes the full \code{MapInteractions} output and creates three clean CSV files:
  \enumerate{
    \item \code{senders_ligands}: One row per sender-ligand with expression stats
    \item \code{ligands_receptor_pairs}: One row per valid L-R pair with stats
    \item \code{receptors_receivers}: One row per receptor-receiver with expression stats
  }

  These files are ideal for importing into Neo4j as nodes and relationships.
}

\usage{
export_for_neo4j(interactions, output_dir = "Neo4J/", prefix = "scSignalMap")
}

\arguments{
  \item{interactions}{
    Data frame returned by \code{\link{MapInteractions}} (the full \code{pairs_data}).
  }
  \item{output_dir}{
    Directory in which to save the three CSV files. Default is \code{"Neo4J/"}.
  }
  \item{prefix}{
    Character string used as a prefix for the output filenames (e.g., \code{"Q1_norm_test"}).
    Default is \code{"scSignalMap"}.
  }
}

\value{
  Invisibly returns a named list containing the full paths to the three created files:
  \item{sender_ligands}{Path to the sender-ligand CSV}
  \item{ligand_receptor_pairs}{Path to the ligand-receptor pairs CSV}
  \item{receiver_receptors}{Path to the receptor-receiver CSV}
}

\examples{
\donttest{
export_for_neo4j(interactions, output_dir = "Neo4J/", prefix = "scSignalMap")
}
}

\author{
  Andrew Zehnder
}

\keyword{IO}
\keyword{Neo4j}
