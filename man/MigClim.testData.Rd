\name{MigClim.testData}
\alias{MigClim.testData}
\title{Test data for the MigClim dispersal simulation tool}
\description{Data frame that contains data for simulating the dispersal of a species under environmental change conditions. The data frame contains information on the species' initial distribution in the landscape, how the suitability of its habitat evolves over time, and some information on which areas in the landscape represent barriers to dispersal for the species. Each row of the data frame contains the information for a given 50 x 50m pixel in the landscape}
\usage{data(MigClim.testData)}
\format{MigClim.testData is a data frame with the following columns
Xcoord: The X coordinate of the given pixel.
Ycoord: The Y coordinate of the given pixel.
InitialDist: The species's initial distribution. 1 = species presence, 0 = species absence.
HSmap1: Habitat suitability in the range 0:1000 for the first time step.
HSmap2: Habitat suitability in the range 0:1000 for the second time step.
HSmap3: ...
...
Barrier: Indicate whether a given pixel is a barrier feature or not (barrier features impede dispersal accross them). 1 = pixel is a barrier, 0 = pixel is not a barrier.}

