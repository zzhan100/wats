\name{wats}
\alias{wats}
\title{Window Adjusted Time Series}
\usage{
wats()
}
\description{
Window adjusted methods to deal with time series.
}
\examples{
data(toy)
eg <- wats(DV=dv0, IV_list=ivs0, L=3, Total_period=3, r=1.2, BW=2, alpha=0.5) ;eg
}
