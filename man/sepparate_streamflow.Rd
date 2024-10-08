% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Events_classification.R
\name{sepparate_streamflow}
\alias{sepparate_streamflow}
\title{Time series formatting (P and Q)}
\usage{
sepparate_streamflow(q_obs, p_obs, par_beta = 0.925)
}
\arguments{
\item{q_obs}{zoo time series of observed streamflow for a defined period}

\item{p_obs}{zoo time series of observed streamflow for a defined period}

\item{par_beta}{numeric representing the filter parameter. Default value 
is 0.925 as recommended by Arnold and Allen (1999)}
}
\value{
zoo object with streamflow ('Q_tot', 'Q_bf', 'Q_dr', 'Q_dr_bool') 
and precipitation data ('P_obs', 'P_bool')
}
\description{
Time series formatting (P and Q)
}
