\name{MigClim.genClust}
\alias{MigClim.genClust}
\title{Run a simulation of migration of genetic clusters.}
\description{Simulate the migration of genetic clusters. Centers of origin of the genetic clusters are picked randomly or defined by the user as the initial distribution. The simulation makes the genetic clusters migrate until the last time step for which data files are provided. Unlike the 'migrate' function in this package, this genetic clusters migration simulation assumes large time scales (e.g., 1000 years per step).}
\usage{MigClim.genClust (hsMap="hsMap", barrier="barrier", nrClusters=4,
  nrIterations=1, threshold=445, outFile="out", initFile="")}
\arguments{
  \item{hsMap}{The 'base' name of the raster files that contain the habitat suitability maps for each time step in ASCII grid format. Iteration numbers (1,2,3,...) and the file extension '.asc' are automatically added to this 'base' name to get the file name for the habitat suitability map for each next iteration. For example, if the habitat suitability raster files are named 'hsMap1.asc', 'hsMap2.asc', etc., the value of this argument should be 'hsMap'. Habitat suitability maps indicate the suitability of each cell to be colonized as a value between 0 (fully unsuitable) and 1000 (fully suitable).}
  \item{barrier}{The 'base' name of the raster files that contain the barriers for each time step in ASCII grid format. Iteration numbers (1,2,3,...) and the file extension '.asc' are automatically added to this 'base' name to get the file name for the barriers for each next iteration. For example, if the barrier raster files are named 'barrier1.asc', 'barrier2.asc', etc., the value of this argument should be 'barrier'. Barrier files indicate whether there is a barrier to migration present (1) or absent (0 or nodata_value) in each cell.}
  \item{nrClusters}{The number of genetic clusters to use.}
  \item{nrIterations}{The number of iterations for which the simulation will be run. This number must match the number of habitat suitability maps and barrier files.}
  \item{threshold}{The threshold value (in [0:1000]) above which a cell is considered suitable.}
  \item{outFile}{The 'base' name of the raster files that will contain the output for each time step in ASCII grid format. Iteration numbers (1,2,3,...) and the file extension '.asc' are automatically added to this 'base' name. For example, is the value of this argument is 'out', the output raster files will be named 'out1.asc', 'out2.asc', etc.}
  \item{initFile}{If an empty string (default value), initial starting points for the genetic clusters are generated at random, and then saved as a raster file with iteration number 0 (e.g., 'out0.asc'). Otherwise, the initial distribution is read from a file with the name as given for this argument. The file name is assumed to be the full name (including the file extension), and to be a raster file in ASCII grid format.}
}
\details{
'nrClusters' origins of the genetic clusters represented by suitable pixels are randomly picked as the inital state. The remaining suitable pixels are assigned to one of these clusters using a nearest neighbor rule. Then, for each following time-step (e.g. every thousand years) up to the present, any suitable pixel in any timeframe t is colonized by the genetic cluster from the closest suitable pixel from timeframe t-1. Alternatively, a user defined distribution of the genetic clusters may be provided in ASCII grid format, with 0 (or nodata_value) as unsuitable, and a value of 1 to 'nrClusters' attributed to each suitable pixel. The habitat suitability maps should be provided in ASCII grid format with value from 0 (totally unsuitable) to 1000 (fully suitable), as typical ouput from BIOMOD. The function output will be written in ASCII grid format.

It is assumed that the input files are in the current working directory. If not, either set the working directory using 'setwd', or include the full pathname in the base name of the input files. The output files are also written to the current working directory, unless a full pathname is included.}
\value{No value is returned, but the distribution of the genetic clusters for each iteration is written to a file in raster ASCII grid format.}
\references{Espindola A., Pellissier L., Maiorano L., Hordijk W., Guisan A., Alvarez N. Niche hindcasting across 24 millennia predicts the current intra-specific species genetic structure. Ecology Letters, in revision, 2012.}
\seealso{Migclim.validate ()}
\examples{
\dontrun{
### Some example data files can be downloaded from the following web page:
###   http://www.unil.ch/ecospat/page89413.html
###
### Run the example as follows (set the current working directory to the
###   folder where the example data files are located):
###
MigClim.genClust (hsMap="hsMap", barrier="barrier", nrClusters=4,
  nrIterations=5, threshold=445, outFile="out", initFile="")
MigClim.validate (validateFile="Validation.txt", nrPoints=79,
  simFile="out5.asc", nrClusters=4)
}}
