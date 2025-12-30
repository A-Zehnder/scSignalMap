\name{filter_enrichr_by_upreg_receptors}
\alias{filter_enrichr_by_upreg_receptors}
\title{Filter Enrichr results to pathways containing upregulated receptors}

\description{
  Takes the full Enrichr results (combined across databases) and a set of upregulated receptors
  (after intersection with ligand-receptor pairs) and retains only those pathways that contain
  at least one of the upregulated receptor gene symbols.

  A new column \code{Matching_Receptors} is added listing the matching receptor(s) (original case)
  for each pathway.
}

\usage{
filter_enrichr_by_upreg_receptors(enrichr_results, upreg_receptors_filtered_and_compared)
}

\arguments{
  \item{enrichr_results}{
    Data frame returned by \code{\link{find_enriched_pathways}} (the combined enrichment results
    across all databases).
  }
  \item{upreg_receptors_filtered_and_compared}{
    Data frame returned by \code{\link{intersect_upreg_receptors_with_lr_interactions}},
    containing upregulated receptors that overlap with known ligand-receptor interactions.
  }
}

\details{
  The function performs case-insensitive matching on gene symbols and maps matches back to their
  original capitalization. It standardizes column names (\code{Adjusted.P.value} → \code{Adjusted_P_value},
  \code{Combined.Score} → \code{Combined_Score}) for consistency.

  Pathways without any matching upregulated receptor are removed. A message is printed indicating
  how many pathways remain after filtering.
}

\value{
  A filtered data frame containing only enriched pathways that include at least one upregulated
  receptor, with an additional column \code{Matching_Receptors} listing the matching receptor
  symbols (separated by ";"). Returns \code{NULL} if there are no enrichr results or no upregulated
  receptors to filter on.
}

\examples{
\donttest{
# Assuming you have run the full pipeline
enrichr_all   = results[[1]]$enrichr_results
upreg_int     = results[[1]]$upreg_receptors_filtered_and_compared

filtered_pathways = filter_enrichr_by_upreg_receptors(
  enrichr_results = enrichr_all,
  upreg_receptors_filtered_and_compared = upreg_int
)

head(filtered_pathways)
}
}

\seealso{
  \code{\link{find_enriched_pathways}},
  \code{\link{intersect_upreg_receptors_with_lr_interactions}}
}

\author{
  Andrew Zehnder
}

\keyword{enrichment}
\keyword{pathway}
\keyword{receptor}
