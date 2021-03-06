% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/module_CHM_G.R
\name{v_chmFoldChangeLog10}
\alias{v_chmFoldChangeLog10}
\title{Calculate log10 fold-change value for module: CHM_G/T (both Gene and Transcript)}
\usage{
v_chmFoldChangeLog10(
  outputFolderPath,
  log10Threshold,
  grpName_fc,
  TPM2RHKG,
  filterNoTPM
)
}
\description{
(Internal) Helper function, used to calculate log10 fold-change value for module CHM_G/T (both Gene and Transcript), should be called within the respective main function (v_chmSignaturePanel for Gene and v_chmTranscript for Transcript).
}
\keyword{internal}
