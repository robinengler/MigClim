\name{MigClim.migrate}
\alias{MigClim.migrate}
\title{Initialize and run the MigClim simulation}
\description{Initialize the main MigClim simulation by passing the relevant
  parameter values, which will be written to a parameter file. The
  simulation is then run with these parameter values.}
\usage{MigClim.migrate (iniDist="InitialDist", hsMap="HSmap", rcThreshold=0, 
  envChgSteps=1, dispSteps=1, dispKernel=c(1.0,1.0), 
  barrier="", barrierType="strong", 
  iniMatAge=1, propaguleProd=c(1.0),
  lddFreq=0.0, lddMinDist=NULL, lddMaxDist=NULL, 
  simulName="MigClimTest", replicateNb=1, overWrite=FALSE, 
  testMode=FALSE, fullOutput=FALSE, keepTempFiles=FALSE)}
\arguments{
  \item{iniDist}{The initial distribution of the species. This can be given either a string indicating the name of a raster file (see 'Details' for supported formats) or as a data frame object (see 'Details' for how to structure your data frame). Please note that the inputs for 'iniDist', 'hsMap' and 'barrier' (optional) must always be given in the same format. Note that the values of the species' initial distribution layer must be binary and integer numbers: 1 (species is present) or 0 (species is absent).}
  \item{hsMap}{The habitat suitability values. This can be given as a string indicating the 'base name' of the raster files that contain the habitat suitability maps. Iteration numbers (1,2,3,...) are automatically added to this 'base name' to get the file name for the habitat suitability map for each successive environmental change iteration (see the 'Details' section for supported formats). Alternatively, the habitat suitability information can also be given as a data frame object, where each column indicates a successive habitat suitability map (see the 'Details' section for further information on how this data frame must be structured). Note that the values of the habitat suitability layers must be integer numbers in the range 0 to 1000.}
  \item{rcThreshold}{The reclassification threshold: an integer value between 0 and 1000; default=0). If 'rcThreshold > 0', then the continuous values of the habitat suitability maps (in the range 0:1000) will be reclassified according to 'rcThreshold'. Values of habitat suitability < 'rcThreshold' are reclassified to '0' (unsuitable habitat) and values >= 'rcThreshold' are reclassified to '1000' (fully suitable habitat).  In the case where 'rcThreshold=0', the habitat suitability values are not reclassified, and are instead considered as habitat 'invasibility', modulating the probability of an unoccupied cell to become colonized (probabilities are computed as 'habitat suitability / 1000').}
  \item{envChgSteps}{The number of environmental change steps to perform. At each environmental change step the habitat suitability values are updated with the values of the corresponding habitat suitability map (and therefore the number of environmental change steps must match the number of habitat suitability maps available).}
  \item{dispSteps}{The number of dispersal steps to perform within each environmental change step. For instance, if one wants to simulate dispersal to occur once a year, and the habitat suitability maps represent 5 years intervals, then 'dispSteps' should be set to 5.}
  \item{dispKernel}{The dispersal kernel. A vector of dispersal probabilities (values in the range 0.0 to 1.0) giving the conditional probability for a source cell to colonize an empty cell given the distance between both cells. The distance unit is the 'pixel', with the first value in the vector representing the probability for a source cell to colonize a directly adjacent cell. See also the MigClim user guide (available by typing 'MigClim.userGuide' in R) for more details on this parameter.}
  \item{barrier}{The name of the raster file that contains barrier information or a single column data frame (or vector) containing this information. If an empty string is given (default value), no barrier information is used. The values of the barrier layer must integer numbers and binary: either 1 (indicating that the cell is a barrier) or 0 (indicating that the cell is not a barrier).}
  \item{barrierType}{The barrier type to use. Values can be either 'strong' (default value) or 'weak'. Not relevant if barrier information is not used. 'weak' barriers will allow dispersal to proceed through two diagonally adjacent barrier pixels, 'strong' barriers won't. See the MigClim user guide (type 'MigClim.userGuide()' in R) for detailed explanations of the difference between these two barrier types.}
  \item{iniMatAge}{The initial maturity age of newly colonized cells. Newly colonized cells younger than this age cannot produce propagules and hence cannot colonize other cells. When newly colonized cells reach an age equal to 'iniMatAge', then their probability to produce propagules is set to the first value indicated in the 'propaguleProd' vector. The time unit that measures cell 'age' is a dispersal step, which usually should be equal to a year.}
  \item{propaguleProd}{The propagule production probability as a function of cell 'age'. A vector where each successive value indicates the propagule production probability of a cell that has reached its 'iniMatAge' age. The first value of the vector corresponds to the cells having an age equal to 'iniMatAge' and successive values correspond to an increase in 1 unit of age from the 'iniMatAge'. When the probability of propagule production reaches 1 (full maturity age), then it is no longer needed to indicate this value in the 'propaguleProd' vector as it will be considered to be 1 from then on. The length of the 'propaguleProd' vector is thus equal to 'full maturity age'-'iniMatAge' (or a length of 1 if both ages are equal). Propagule production probabilities must be given in the range 0.0 to 1.0.}
  \item{lddFreq}{The long-distance dispersal frequency, i.e., the probability for an occupied cell with full propagule production potential to generate a long distance dispersal event. If set to 0.0 (default), no long-distance dispersal is performed. Value should be given in the range 0.0 and 1.0.}
  \item{lddMinDist}{The minimum distance for long-distance dispersal (not used if 'lddFreq=0'). This value must be larger than the regular dispersal distance given by the length of 'dispKernel'.}
  \item{lddMaxDist}{The maximum distance for long-distance dispersal (not used if 'lddFreq=0'). This value must be >= 'lddMinDist'.}
  \item{simulName}{The 'base name' to be used for the different outputs produced by the MigClim simulation. Three different types of outputs are produced by the 'MigClim.migrate()' function: ascii grid files named 'simulName'+'_raster.asc' that contains the final state of the simulation, 'simulName'+'_stats.txt' files that contain the simulation's outputs after each dispersal step, and 'simulName'+'_summary.txt' files that contain a single-line summary of the entire simulation.}
  \item{replicateNb}{Number of times a simulation should be replicated. The final outputs include all the outputs from individual runs as well as the average of all runs.}
  \item{overWrite}{If 'TRUE' then any existing file with the same name as an ouput of the MigClim.migrate function will be mercilessly overwritten. If 'FALSE' then the function will stop if any output file does already exist.}
  \item{testMode}{If 'TRUE' then the MigClim.migrate function will check all the provided input data but will not run the actual simulation. Useful for testing your data before running several successive simulations or simulations that might take a long time.}
  \item{fullOutput}{If 'TRUE', the current state of the simulation is written to an ASCII raster file after each dispersal step (allowing to reconstruct the dispersal process at each step). If 'FALSE' (default), only the final state of the simulation is written to an ASCII grid file.}
  \item{keepTempFiles}{If 'FALSE' (default), then any '.asc' file created from a conversion process in the function will be deleted when the simulation completes. If you wish to keep these files then set the value of this parameter to 'TRUE'.}
}
\details{The input data for initial distribution ('iniDist'), habitat suitability ('hsMap'), and (optionally) barriers ('barrier') can be provided as either a string giving the name of a raster file (the name should be given relative to the working directory) or as a data frame object. For a given simulation, all these inputs must be given in the same format.
Option 1: Input data provided as raster files. In this case, the input must be a string that contains the name of the raster files relative to the working directory. The following raster formats are supported: (i) ascii grid (files must have a '.asc' extension), (ii) R rasterLayer (see 'raster' package), (iii) ESRI GRID, (iv) GeoTIFF (files must have a '.tif' extension). Note that all input grids need to have exactly the same pixel size and the same extent (i.e. the same number of rows and columns).

The standard ASCII grid Raster format looks as follows (actual values depend on file content):
\preformatted{ncols         100
nrows         100
xllcorner     0
yllcorner     0
cellsize      10
nodata_value -9999
0 0 1 0 0 1 1 1 0 ...
1 0 1 0 1 0 1 0 1 ...
:
:
}

Option 2: Input data provided as data frames. It is important to note that, unlike when inputing data in a raster format, inputing it as data frames requires to provide the actual data frame object to the function, not a string indicating its name. 'iniDist' shoud have 3 colums (in this order): X coordinate, Y coordinate and initial distribution of the species. Only values of 0 or 1 are accpted: 1 = species present, 0 = species absent.
\preformatted{Xcoordinate  Ycoordinate IniDist
520000       200100      1
520050       200100      1
520100       200100      0
...
}
          'hsMap' shoud have a number of colums equal to 'envChgSteps'. Each column giving the suitability of the habitat for a given environmental change step (in chronological order). Values of habitat suitability must be in the range [0:1000].
\preformatted{
HSmap1  HSmap2  HSmap3  ...
789     800     230     ...
560     720     150     ...
120     540      86     ...
...
}

'barrier' is optional and must have only one column (it can also be a vector). It must contain only values of either 1 (pixel is a barrier feature) or 0 (pixel is not a barrier feature).

The function output(s) will be written in ascii GRID format (with .asc extension).
}
\value{The number of environmental change steps performed. The function also writes the following outputs into the current working directory: an ASCII grid raster file named 'simulName'+'_raster.asc' that contains the final state of the simulation, a 'simulName'+'_stats.txt' file that contains the simulation's outputs after each dispersal event, and a 'simulName'+'_summary.txt' file that contains a single-line summary of the entire simulation. If fullOutput=TRUE then an ASCII raster file containing the state of the simulation at the end of each dispersal step is also saved as output with the following name structure: 'simulName' + '_step_' + dispersal step code + '.asc'. The output file 'simulName' + '_stats.txt' contains summary statistics for each individual dispersal step. The output file 'simulName' + '_summary.txt' contains summary statistics over the entire simulation (including the running time in seconds).}
\references{Engler R., Hordijk W. and Guisan A. The MigClim R package - seamless integration of dispersal constraints into projections of species distribution models. Ecography, in review.}
\seealso{MigClim.plot(), MigClim.userGuide()}
\examples{
  
  data(MigClim.testData)

  ### Run MigClim with a data frame type input.
  n<-MigClim.migrate (iniDist=MigClim.testData[,1:3],
     hsMap=MigClim.testData[,4:8], rcThreshold=500, 
     envChgSteps=5, dispSteps=5, dispKernel=c(1.0,0.4,0.16,0.06,0.03),
     barrier=MigClim.testData[,9], barrierType="strong", 
     iniMatAge=1, propaguleProd=c(0.01,0.08,0.5,0.92),
     lddFreq=0.1, lddMinDist=6, lddMaxDist=15, 
     simulName="MigClimTest", replicateNb=1, overWrite=TRUE, 
     testMode=FALSE, fullOutput=FALSE, keepTempFiles=FALSE)
  
\dontrun{  ### Convert data from a dataframe to ASCII grid format
  library(SDMTools)
  dataframe2asc(MigClim.testData[,c(2,1,3:9)])
  ### and then run MigClim with "raster" type input.
  n<-MigClim.migrate (iniDist="InitialDist",
     hsMap="HSmap", rcThreshold=500, 
     envChgSteps=5, dispSteps=5, dispKernel=c(1.0,0.4,0.16,0.06,0.03), 
     barrier="Barrier", barrierType="strong", 
     iniMatAge=1, propaguleProd=c(0.01,0.08,0.5,0.92), 
     lddFreq=0.1, lddMinDist=6, lddMaxDist=15, 
     simulName="MigClimTest", replicateNb=1, overWrite=TRUE, 
     testMode=FALSE, fullOutput=FALSE, keepTempFiles=FALSE)}}

