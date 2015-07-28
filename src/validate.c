/*
** validate.c: Validate a genetic clusters simulation result by comparing it
**             with an observed cluster distribution and calculating a matching
**             score.
**
** Wim Hordijk   Last modified: 25 January 2012
*/

#include "migclim.h"


/*
** validate: Compare a genetic clusters simulation result with an observed
**           cluster and calculate a matching score. This score is the fraction
**           of observed points that are closest to their own cluster in the
**           simulated result (after the re-labeling of clusters).
**
** Parameters:
**   - obsFileName: The file name with the observed cluster data.
**   - npts:        The number of points (locations) in the data file.
**   - simFileName: The file name with the simulated cluster data.
**   - ncls:        The number of genetic groups.
**   - bestScore:   A pointer to a double array with at least 2 elements,
**                  which will contains the scores (total and average).
*/

void validate (char **obsFileName, int *npts, char **simFileName, int *ncls,
	       double *bestScore)
{
  int     i, j, x, y, n, cl, *obsCluster[2], **simCluster, *match, obs, sim,
          row, col, *points, nrPoints, nrClusters;
  float   c0, c1;
  double *coord[2], dist, minDist, *score, totScore, totScoreMax, avgScore,
          avgScoreMax;
  char    line[1024];
  FILE   *fp;
  
  /*
  ** Initialize the variables.
  */
  nrPoints = *npts;
  nrClusters = *ncls;
  score = NULL;
  points = NULL;
  fp = NULL;
  coord[0] = NULL;
  coord[1] = NULL;
  obsCluster[0] = NULL;
  obsCluster[1] = NULL;
  simCluster = NULL;
  match = NULL;
  bestScore[0] = -1;
  bestScore[1] = -1;
  
  /*
  ** Allocate the necessary memory.
  */
  coord[0] = (double *)malloc (nrPoints * sizeof (double));
  coord[1] = (double *)malloc (nrPoints * sizeof (double));
  obsCluster[0] = (int *)malloc (nrPoints * sizeof (int));
  obsCluster[1] = (int *)malloc (nrPoints * sizeof (int));
  simCluster = (int **)malloc (nrRows * sizeof (int *));
  for (i = 0; i < nrRows; i++)
  {
    simCluster[i] = (int *)malloc (nrCols * sizeof (int));
  }
  match = (int *)malloc ((nrClusters+1) * sizeof (int));
  points = (int *)malloc ((nrClusters+1) * sizeof (int));
  score = (double *)malloc ((nrClusters+1) * sizeof (double));
  points[0] = nrPoints;
  for (i = 1; i <= nrClusters; i++)
  {
    points[i] = 0;
  }
  
  /*
  ** Read the observed cluster data.
  */
  if ((fp = fopen(*obsFileName, "r")) == NULL)
  {
    Rprintf ("Can't open data file %s.\n", obsFileName);
    goto End_of_Routine;
  }
  if (fgets (line, 1024, fp) == NULL)
  {
    Rprintf ("No data in file %s.\n", obsFileName);
    goto End_of_Routine;
  }
  for (i = 0; i < nrPoints; i++)
  {
    if (fgets (line, 1024, fp) == NULL)
    {
      Rprintf ("Invalid number of data points in file %s.\n", obsFileName);
      goto End_of_Routine;
    }
    if (sscanf (line, "%d %g %g %d", &n, &c0, &c1, &cl) != 4)
    {
      Rprintf ("Invalid data format in file %s.\n", obsFileName);
      goto End_of_Routine;
    }
    coord[0][i] = c0;
    coord[1][i] = c1;
    obsCluster[0][i] = cl;
    points[cl]++;
  }
  fclose (fp);
  fp = NULL;

  /*
  ** Read the simulated cluster data.
  */
  if (readMat (*simFileName, simCluster) == -1)
  {
    goto End_of_Routine;
  }

  /*
  ** For each observed cluster point, find the nearest simulated cluster point.
  */
  for (i = 0; i < nrPoints; i++)
  {
    minDist = (nrRows + nrCols) * cellSize;
    for (row = 0; row < nrRows; row++)
    {
      for (col = 0; col < nrCols; col++)
      {
	if (simCluster[row][col] > 0)
	{
	  x = col;
	  y = (nrRows-1) - row;
	  c0 = xllCorner + (x*cellSize) + (cellSize/2.0);
	  c1 = yllCorner + (y*cellSize) + (cellSize/2.0);
	  dist = sqrt ((coord[0][i]-c0)*(coord[0][i]-c0) +
		       (coord[1][i]-c1)*(coord[1][i]-c1));
	  if (dist < minDist)
	  {
	    obsCluster[1][i] = simCluster[row][col];
	    minDist = dist;
	  }
	}
      }
    }
  }

  /*
  ** Create random matches and calculate scores.
  */
  for (i = 1; i <= nrClusters; i++)
  {
    match[i] = i;
  }
  totScoreMax = 0.0;
  avgScoreMax = 0.0;
  for (i = 1; i <= 1000; i++)
  {
    for (j = nrClusters; j > 1; j--)
    {
      x = (rand () % j) + 1;
      y = match[j];
      match[j] = match[x];
      match[x] = y;
    }
    for (j = 0; j <= nrClusters; j++)
    {
      score[j] = 0.0;
    }
    for (j = 0; j <= nrPoints; j++)
    {
      obs = obsCluster[0][j];
      sim = obsCluster[1][j];
      if (sim == match[obs])
      {
	score[0]++;
	score[obs]++;
      }
    }
    for (j = 0; j <= nrClusters; j++)
    {
      score[j] /= points[j];
    }
    totScore = score[0];
    avgScore = 0.0;
    for (j = 1; j <= nrClusters; j++)
    {
      avgScore += (score[j]/nrClusters);
    }
    if (totScore > totScoreMax)
    {
      totScoreMax = totScore;
    }
    if (avgScore > avgScoreMax)
    {
      avgScoreMax = avgScore;
    }
  }
  /*
  ** Set the best score that was found.
  */
  bestScore[0] = totScoreMax;
  bestScore[1] = avgScoreMax;
  
 End_of_Routine:
  /*
  ** Free the allocated memory and close the file, if necessary.
  */
  if (coord[0] != NULL)
  {
    free (coord[0]);
  }
  if (coord[1] != NULL)
  {
    free (coord[1]);
  }
  if (obsCluster[0] != NULL)
  {
    free (obsCluster[0]);
  }
  if (obsCluster[1] != NULL)
  {
    free (obsCluster[1]);
  }
  if (simCluster != NULL)
  {
    for (i = 0; i < nrRows; i++)
    {
      free (simCluster[i]);
    }
    free (simCluster);
  }
  if (match != NULL)
  {
    free (match);
  }
  if (points != NULL)
  {
    free (points);
  }
  if (score != NULL)
  {
    free (score);
  }
  if (fp != NULL)
  {
    fclose (fp);
  }
}


/*
** EoF: validate.c
*/
