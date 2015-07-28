/*
** genclust.c: Implementation of the main genetic clusters migration
**             simulation.
**
** Wim Hordijk   Last modified: 25 January 2012
*/

#include "migclim.h"


/*
** Global variables.
*/
int    nrRows, nrCols, noData;
double xllCorner, yllCorner, cellSize;


/*
** genClust: The core of the genetic cluster migration method. Select random
**           starting points (or read them from file) and perform the migration
**           steps.
**
** Parameters:
**   - nrow:         The number of rows in the data files.
**   - ncol:         The number of columns in the data files.
**   - ncls:         The number of genetic clusters.
**   - niter:        The number of iterations to perform (i.e., number of
**                   data files).
**   - thrs:         The suitability threshold.
**   - suitBaseName: The base name of the suitability data files.
**   - barrBaseName: The base name of the barrier data files.
**   - outBaseName:  The base name of the output data files.
**   - initFile:     The name of the file to read the initial cluster
**                   distribution (starting points) from. If empty, it will
**                   be generated at random and written to a file.
*/

void genClust (int *nrow, int *ncol, int *ncls, int *niter, int *thrs,
	       char **suitBaseName, char **barrBaseName, char **outBaseName,
	       char **initFile)
{
  int    i, j, x, y, **curState, **prevState, **suitability, **barrier, iter,
         nrClusters, nrIterations, threshold;
  double minDist, dist;
  bool   done;
  char   fileName[128];

  /*
  ** Initialize the variables.
  */
  nrRows = *nrow;
  nrCols = *ncol;
  nrClusters = *ncls;
  nrIterations = *niter;
  threshold = *thrs;
  curState = NULL;
  prevState = NULL;
  suitability = NULL;
  barrier = NULL;
  srand (time(NULL));
  
  /*
  ** Allocate the necessary memory.
  */
  curState = (int **)malloc (nrRows * sizeof (int *));
  prevState = (int **)malloc (nrRows * sizeof (int *));
  suitability = (int **)malloc (nrRows * sizeof (int *));
  barrier = (int **)malloc (nrRows * sizeof (int *));
  for (i = 0; i < nrRows; i++)
  {
    curState[i] = (int *)malloc (nrCols * sizeof (int));
    prevState[i] = (int *)malloc (nrCols * sizeof (int));
    suitability[i] = (int *)malloc (nrCols * sizeof (int));
    barrier[i] = (int *)malloc (nrCols * sizeof (int));
  }

  /*
  ** Read the first suitability and barrier data files and initialize the
  ** state matrices.
  */
  sprintf (fileName, "%s1.asc", *suitBaseName);
  if (readMat (fileName, suitability) == -1)
  {
    goto End_of_Routine;
  }
  sprintf (fileName, "%s1.asc", *barrBaseName);
  if (readMat (fileName, barrier) == -1)
  {
    goto End_of_Routine;
  }
  for (i = 0; i < nrRows; i++)
  {
    for (j = 0; j < nrCols; j++)
    {
      if (suitability[i][j] == noData)
      {
	curState[i][j] = prevState[i][j] = noData;
      }
      else
      {
	curState[i][j] = prevState[i][j] = 0;
      }
    }
  }
  
  /*
  ** Get an initial distribution for the genetic clusters. If a file name
  ** is given, read the initial distribution from file. Otherwise create
  ** it at random.
  */
  if (strlen (*initFile) > 1)
  {
    /*
    ** Read the initial state from file.
    */
    if (readMat (*initFile, prevState) == -1)
    {
      goto End_of_Routine;
    }
  }
  else
  {
    /*
    ** Generate starting points at random.
    */
    for (i = 1; i <= nrClusters; i++)
    {
      done = false;
      while (!done)
      {
	x = rand () % nrRows;
	y = rand () % nrCols;
	if ((prevState[x][y] == 0) && (suitability[x][y] >= threshold) &&
	    (barrier[x][y] != 1))
	{
	  prevState[x][y] = i;
	  done = true;
	}
      }
    }
    /*
    ** Save the initial state matrix.
    */
    sprintf (fileName, "%s0.asc", *outBaseName);
    if (writeMat (fileName, prevState) == -1)
    {
      goto End_of_Routine;
    }
  }

  /*
  ** Iterate the migration simulation.
  */
  Rprintf ("Running the genetic clusters migration simulation.\n");
  for (iter = 1; iter <= nrIterations; iter++)
  {
    Rprintf ("  %d...\n", iter);
    /*
    ** Read the next suitability and barrier data files.
    */
    if (iter > 1)
    {
      sprintf (fileName, "%s%d.asc", *suitBaseName, iter);
      if (readMat (fileName, suitability) == -1)
      {
	goto End_of_Routine;
      }
      sprintf (fileName, "%s%d.asc", *barrBaseName, iter);
      if (readMat (fileName, barrier) == -1)
      {
	goto End_of_Routine;
      }
    }

    /*
    ** For each suitable cell in the current state, find the nearest occupied
    ** cell in the previous state and copy the cluster number to the current
    ** cell.
    */
    for (i = 0; i < nrRows; i++)
    {
      for (j = 0; j < nrCols; j++)
      {
	if (suitability[i][j] == noData)
	{
	  curState[i][j] = noData;
	}
	else
	{
	  curState[i][j] = 0;
	}
	/*
	** Check if the current cell is suitable.
	*/
	if ((curState[i][j] == 0) && (suitability[i][j] >= threshold) &&
	    (barrier[i][j] != 1))
	{
	  if (prevState[i][j] > 0)
	  {
	    /*
	    ** If the current cell was already occupied, keep the cluster
	    ** number.
	    */
	    curState[i][j] = prevState[i][j];
	  }
	  else
	  {
	    /*
	    ** Find the closest occupied cell in the previous state and copy
	    ** its cluster number to the current cell.
	    */
	    minDist = nrRows + nrCols;
	    for (x = 0; x < nrRows; x++)
	    {
	      for (y = 0; y < nrCols; y++)
	      {
		if (prevState[x][y] > 0)
		{
		  dist = sqrt ((x-i)*(x-i) + (y-j)*(y-j));
		  if (dist < minDist)
		  {
		    curState[i][j] = prevState[x][y];
		    minDist = dist;
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    /*
    ** Write the current state matrix to file.
    */
    sprintf (fileName, "%s%d.asc", *outBaseName, iter);
    if (writeMat (fileName, curState) == -1)
    {
      goto End_of_Routine;
    }

    /*
    ** Replace the previous state with the current state.
    */
    for (i = 0; i < nrRows; i++)
    {
      for (j = 0; j < nrCols; j++)
      {
	prevState[i][j] = curState[i][j];
      }
    }
  }
  Rprintf ("done.\n");
  
 End_of_Routine:
  /*
  ** Free the allocated memory.
  */
  if (curState != NULL)
  {
    for (i = 0; i < nrRows; i++)
    {
      free (curState[i]);
    }
    free (curState);
  }
  if (prevState != NULL)
  {
    for (i = 0; i < nrRows; i++)
    {
      free (prevState[i]);
    }
    free (prevState);
  }
  if (suitability != NULL)
  {
    for (i = 0; i < nrRows; i++)
    {
      free (suitability[i]);
    }
    free (suitability);
  }
  if (barrier != NULL)
  {
    for (i = 0; i < nrRows; i++)
    {
      free (barrier[i]);
    }
    free (barrier);
  }
}    


/*
** EoF: genclust.c
*/
