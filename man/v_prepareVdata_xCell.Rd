% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/module_xCell.R
\name{v_prepareVdata_xCell}
\alias{v_prepareVdata_xCell}
\title{Prepare data for module: xCell}
\usage{
v_prepareVdata_xCell(
  doGE = FALSE,
  doLocalAnalysis = FALSE,
  colSliceOrder = grpName
)
}
\arguments{
\item{doGE}{logic, whether to prepare Gene Expression (GE) data, if FALSE, no GE data will be available for downstream v_chmXcell function; if TRUE, please make sure GE data files are properly named, see Details for more information about file naming.}

\item{doLocalAnalysis}{logic, whether to perform local Cell Type Enrichment Analysis (based on xCell), if FALSE, no local analysis results will be available for downstream v_chmXcell function.}

\item{colSliceOrder}{character vector, the order (from left to right) of groups to be shown on the output plots and files, should be set to the same value when called from different functions within the same module (e.g., v_prepareVdata_xCell and v_chmXcell). By default will use the internal character vector "grpName" and can be left unchanged. Groups not included in the "colSliceOrder" will also be excluded from the output plots and files (in most cases), and from certain analysis process (depending on the situation).}
}
\value{
list, because R CMD check discourages assignments to the global environment within functions, user needs to run the function with explicitly assigning the return value to a global variable named "prepareVdata_xCell_returnList", which will be a list containing the required variables for downstream analyses.
}
\description{
The workflow of vigilante is highly module-based. To ensure a successful and smooth run, vigilante needs to prepare input data before continuing.
}
\details{
Oftentimes input data files generated by upstream tools came with diverse naming conventions. It might be easy for the user to recognize those files, but not for vigilante if there is no recognizable patterns.

To make input data files clear to vigilante, it would be nice to have them named something like "studyID_sampleID_(other descriptions).file extension". Here "studyID" is the name of the study or project, and it will be used in multiple naming situations (such as on the plot, or in the output file names), so it is recommended to be concise and meaningful.

For module xCell, currently supported input data files are listed below, please contact the author if you want to add more files to the supported list:
Gene Expression (GE): *cDNA_genes.sf for Salmon
}
