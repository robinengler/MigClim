\name{MigClim.validate}
\alias{MigClim.validate}
\title{Validation of genetic clusters migration simulation results.}
\description{Compare the output of a genetic clusters migration simulation to an observed genetic clusters distribution.}
\usage{MigClim.validate (validateFile="Validation.txt", nrPoints=0, simFile="out1.asc",
  nrClusters=4)}
\arguments{
  \item{validateFile}{The name of the file containing the observed distribution. See details below for the required file format. A full file name (including file extension) is expected.}
  \item{nrPoints}{The number of points in the observed distribution in the 'validateFile' file.}
  \item{simFile}{The name of an output file from the genetic cluster migration simulation. Here, a full file name is expected, i.e., including an iteration number and the '.asc' extension, for example 'out5.asc'.}
  \item{nrClusters}{The number of genetic clusters.}
}
\details{
This function calculates a matching score between the genetic clusters distribtion predicted by a simulated scenario and an observed population genetic structure. The observed distribution (i.e., the validation file) should be in the following format:
\preformatted{   N  X        Y         C
   3  24.62192	45.60239  1
  18  19.23068	48.88922  1
  23  24.63611	46.59972  2
  25  23.59028	46.70118  3
  30  23.32314	42.19968  4
   :    :        :       :
}
The first line in the file is a header line, and each next line should contain four values, being an identification number (not used here), an X-coordinate, a Y-coordinate, and finally the number of the genetic cluster the point belongs to. There should be 'nrPoints' lines in this file.}
\value{A score (between 0 and 1) indicating the fit of the distribution in the output file. This is given as an array with two elements. The first element gives the "total" score, i.e., measured over all points regardless of which cluster each point belongs to. The second element give the "average" score, where the fit is calculated over each cluster separately, and then averaged over all cluster scores (so small clusters get the same weight as large clusters).}
\seealso{MigClim.genClust ()}
