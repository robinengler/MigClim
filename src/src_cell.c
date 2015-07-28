/*
** src_cell.c: Function for finding potential souce cells for a given sink cell.
**
** Wim Hordijk    Last modified: 03 October 2011
**
** This C code is based on the original Visual Basic code of Robin Engler.
*/

#include "migclim.h"


/*
** mcSrcCell: This function will search, for a given input location, for a
**            suitable "source" pixel that can lead to the colonization of
**            the input location (sink pixel).
**
** Parameters:
**   - i:        The row number of the sink pixel.
**   - j:        The column number of the sink pixel.
**   - curState: The matrix that contains the current state of the cellular
**               automaton.
**   - pxlAge:   A matrix giving the "age" of each colonized pixel.
**   - loopID:   Indicates in which "loop" the simulation currently is. The
**               'loopID' enables to retrieve the 'environmental change' loop
**               and the 'dispersal' loop that the simulation is currently in.
**   - habSuit:  The habitat suitability of the current pixel.
**   - barriers: A pointer to the barriers matrix.
**
** Returns:
**   If a suitable source cell was found: true.
**   Otherwise:                           false.
*/

bool mcSrcCell (int i, int j, int **curState, int **pxlAge, int loopID,
		int habSuit, int **barriers)
{
  int    k, l, realDist, pxlSizeFactor;
  double probCol, rnd;
  bool   sourceFound;

  /*
  ** For now let's set these paramters to fixed values. Later we can implement
  ** them as variables.
  */
  pxlSizeFactor = 1;
  sourceFound = false;
        
  /*
  ** Search for a potential source cell. i and j are the coordinates of the
  ** sink cell. k and l are the coordinates of the potential source cell.
  */
  for (k = i - dispDist; k <= i + dispDist; k++)
  {
    for (l = j - dispDist; l <= j + dispDist; l++)
    {
      /*
      ** 1. Test of basic conditions to see if a pixel could be a potential
      **    source cell:
      **    - The pixel must be within the limits of the matrix's extent.
      **    - The pixel must be colonized, but not during the current loop.
      **    - The pixel must have reached its age of "initial maturity"
      **      (otherwise it cannot produce seeds).
      */
      if ((k >= 0) && (k < nrRows) && (l >= 0) && (l < nrCols))
      {
	if ((curState[k][l] > 0) &&
	    (curState[k][l] != loopID))
	{
	  if (pxlAge[k][l] >= iniMatAge)
	  {
	    /*
	    ** 2. Compute the distance between sink and (potential) source pixel
	    **    and check if it is <= maximum dispersal distance. The distance
	    **    is computed in pixel units.
	    ** realDist = Fix(Sqr((K - I) ^ 2 + (L - J) ^ 2) + 0.5)
	    */
	    realDist = (int)round (sqrt ((k-i)*(k-i) + (l-j)*(l-j)));
	    if ((realDist > 0) && (realDist <= dispDist))
	    {
	      /*
	      ** 3. Compute the probability of colonization of the sink pixel.
	      **    This probability depends on several factors:
	      **    - Disance between source and sink cells.
	      **    - Age of the source cell.
	      **    - "Invasability" of the sink cell.
	      */
	      if (pxlAge[k][l] >= fullMatAge)
	      {
		probCol = dispKernel[realDist-1] * pxlSizeFactor *
			    (habSuit / 1000.0);
	      }
	      else
	      {
		probCol = dispKernel[realDist-1] * pxlSizeFactor *
			    propaguleProd[pxlAge[k][l] - iniMatAge] *
			    (habSuit / 1000.0);
	      }

	      rnd = UNIF01;
	      if (rnd < probCol || probCol == 1.0)
	      {
		/*
		** When we reach this stage, the last thing we need to check for
		** is whether there is a "barrier" obstacle between the source
		** and sink pixel. We check this last as it requires significant
		** computing time.
		*/
		if (useBarrier)
		{
		  if (!mcIntersectsBarrier (i, j, k, l, barriers))
		  {
		    sourceFound = true;
		    goto End_of_Routine;
		  }
		}
		else
		{
		  sourceFound = true;
		  goto End_of_Routine;
		}
	      }
	    }
	  }
	}
      }
    }
  }
    
 End_of_Routine:
  /*
  ** Return the result.
  */
  return (sourceFound);
}


/*
** EoF: src_cell.c
*/
