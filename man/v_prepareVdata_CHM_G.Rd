% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/module_CHM_G.R
\name{v_prepareVdata_CHM_G}
\alias{v_prepareVdata_CHM_G}
\title{Prepare data for module: CHM_G (Gene)}
\usage{
v_prepareVdata_CHM_G(doGE = FALSE)
}
\arguments{
\item{doGE}{logic, whether to prepare Gene Expression (GE) data, if FALSE, no GE data will be available for downstream v_chmSignaturePanel function; if TRUE, please make sure GE data files are properly named, see Details for more information about file naming.}
}
\value{
list, because R CMD check discourages assignments to the global environment within functions, user needs to run the function with explicitly assigning the return value to a global variable named "prepareVdata_CHM_G_returnList", which will be a list containing the required variables for downstream analyses.
}
\description{
The workflow of vigilante is highly module-based. To ensure a successful and smooth run, vigilante needs to prepare input data before continuing.
}
\details{
Oftentimes input data files generated by upstream tools came with diverse naming conventions. It might be easy for the user to recognize those files, but not for vigilante if there is no recognizable patterns.

To make input data files clear to vigilante, it would be nice to have them named something like "studyID_sampleID_(other descriptions).file extension". Here "studyID" is the name of the study or project, and it will be used in multiple naming situations (such as on the plot, or in the output file names), so it is recommended to be concise and meaningful.

For module CHM_G, currently supported input data files are listed below, please contact the author if you want to add more files to the supported list:
Gene Expression (GE): *cDNA_genes.sf for Salmon
}
