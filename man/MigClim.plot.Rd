\name{MigClim.plot}
\alias{MigClim.plot}
\title{Plots ascii raster ouputs of the MigClim.migrate function}
\description{Plots ascii raster ouputs of the MigClim.migrate function with an adequate color scale and saves the result as a JPEG or PNG file (or does simply display it in the R console).}
\usage{MigClim.plot(asciiFile, outDir="", fileFormat="jpeg", fullOutput=FALSE)}
\arguments{
  \item{asciiFile}{Name of the main ascii raster ouput "_raster.asc" obtained when running a simulation using MigClim.migrate. If the ascii file is not located in the working directory, then the file name should include its path as well. File names should end in "_raster.asc" or "_raster".}
  \item{outDir}{Optional ouput directory. When this value is set to its default ("") then the ouput files are saved in the same directory as the ascii grids that are being plotted.}
  \item{fileFormat}{Format in which the ouput graphics are saved. Can be either "jpeg", "png" or "inR". If "inR" is passed as parameter value, the plotted image is not saved to disk but simply displayed in the R console.}
  \item{fullOutput}{If "TRUE" then all the intermediate ouputs for the given simulation are also plotted and saved to disk.}
}
\value{The function does not return any value but saves "jpeg" or "png" image files to disk, or displays a graphic in the R consile (fileFormat="inR").}
\seealso{MigClim.migrate ()}
