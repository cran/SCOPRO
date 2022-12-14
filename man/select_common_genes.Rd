% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/markers_functions.R
\name{select_common_genes}
\alias{select_common_genes}
\title{select_common_genes}
\usage{
select_common_genes(
  SCOPRO_output,
  marker_stages,
  selected_stages,
  name_vivo,
  cluster_vitro,
  name_vitro
)
}
\arguments{
\item{SCOPRO_output}{output given by function \emph{SCOPRO}}

\item{marker_stages}{Second element of the list given as output by the function \emph{select_top_markers}}

\item{selected_stages}{In vivo stages for which the markers where computed with the function \emph{select_top_markers}}

\item{name_vivo}{name of the in vivo stage on which SCOPRO is run}

\item{cluster_vitro}{cluster for in vitro dataset}

\item{name_vitro}{name of the in vitro stage for which we want to know the conserved markers with the \emph{name_vivo} stage}
}
\value{
Character vector with the names of the conserved markers of \emph{name_vivo} stage in the \emph{name_vitro} stage
}
\description{
select_common_genes
}
\author{
Gabriele Lubatti \email{gabriele.lubatti@helmholtz-muenchen.de}
}
