/*
** barriers.c: Functions for performing the barriers related steps.
** Wim Hordijk & Robin Engler   Last modified: 11 May 2012 (RE)
*/

#include "migclim.h"


/*
** mcFilterMatrix: Filter the input matrix (inMatrix) using the filter matrix (filterMatrix).
**                 The different filter possibilities are the following:
**                  -> replace any value < 0 by 0 (this removes NoData values)
**                  -> replace any value of 1 in filterMatrix by 0
**                  -> replace any value of NoData (-9999) in filterMatrix by NoData.
**                    
** Parameters:
**   -> **inMatrix: A pointer to the input matrix.
**   -> **filterMatrix: A pointer to the barrier matrix.
**   ->   filterNoData: if true, removes any value < 0 in inMatrix.
**   ->   filterOnes: if true, replace any value of 1 in filterMatrix by 0.
**   ->   insertNoData: if true, replace any value of NoData (-9999) in filterMatrix by NoData.
*/

void mcFilterMatrix(int **inMatrix, int **filterMatrix, bool filterNoData, bool filterOnes, bool insertNoData)
{
  int i, j;

  /* Set any value < 0 to 0 (removes NoData) */
  if(filterNoData){
    for (i = 0; i < nrRows; i++){
    for (j = 0; j < nrCols; j++){
      if (inMatrix[i][j] < 0) inMatrix[i][j] = 0;
    }
    }
  }

  /* Filter the input matrix for values of 1. */
  if(filterOnes){
    for (i = 0; i < nrRows; i++){
    for (j = 0; j < nrCols; j++){
      if (filterMatrix[i][j] == 1) inMatrix[i][j] = 0;
    }
    }
  }
  
  /* Add NoData where it is present in filterMatrix. */
  if(insertNoData){
    for (i = 0; i < nrRows; i++){
    for (j = 0; j < nrCols; j++){
      if (filterMatrix[i][j] == -9999) inMatrix[i][j] = -9999;
    }
    }
  }
  
}




/*
** mcIntersectsBarrier: Check whether there is a barrier between the source
**                      and sink pixels.
** Parameters:
**   - snkX:     The x-coordinate of the sink pixel.
**   - snkY:     The y-coordinate of the sink pixel.
**   - srcX:     The x-coordinate of the source pixel.
**   - srcY:     The y-coordinate of the source pixel.
**   - barriers: A pointer to the barriers matrix.
**
** Returns:
**   If there is a barrier: True.
**   Otherwise:             False.
*/

bool mcIntersectsBarrier (int snkX, int snkY, int srcX, int srcY,
			              int **barriers)
{
  int  dstX, dstY, i, pxlX, pxlY, distMax, barCounter;
  bool barFound;

  barFound = false;
  
  /*
  ** Calculate the distance in both dimensions between the source and sink
  ** pixels and take the largest of the two.
  */
  dstX = srcX - snkX;
  dstY = srcY - snkY;
  if (abs (dstX) >= abs (dstY))
  {
    distMax = abs (dstX);
  }
  else
  {
    distMax = abs (dstY);
  }

  /*
  ** Check the possible paths from source to sink and see if there is a path
  ** without barriers.
  */
  if (barrierType == WEAK_BARRIER)
  {
    /*
    ** Weak barrier: If there is at least one free path we're good.
    **
    ** BARRIER MIDDLE
    */
    barFound = false;
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX + (1.0 * i / distMax * dstX));
      pxlY = (int)round (snkY + (1.0 * i / distMax * dstY));
      if (barriers[pxlX][pxlY] == 1)
      {
	barFound = true;
	break;
      }
    }
    if (!barFound)
    {
      goto End_of_Routine;
    }
    /*
    ** BARRIER TOP_LEFT
    */
    barFound = false;
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX - 0.49 + (1.0 * i / distMax * dstX));
      pxlY = (int)round (snkY - 0.49 + (1.0 * i / distMax * dstY));
      if (barriers[pxlX][pxlY] == 1)
      {
	barFound = true;
	break;
      }
    }
    if (!barFound)
    {
      goto End_of_Routine;
    }
    /*
    ** BARRIER TOP_RIGHT
    */
    barFound = false;
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX + 0.49 + (1.0 * i / distMax * dstX));
      pxlY = (int)round (snkY - 0.49 + (1.0 * i / distMax * dstY));
      if (barriers[pxlX][pxlY] == 1)
      {
	barFound = true;
	break;
      }
    }
    if (!barFound)
    {
      goto End_of_Routine;
    }
    /*
    ** Barrier DOWN_LEFT
    */
    barFound = false;
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX - 0.49 + (1.0 * i / distMax * dstX));
      pxlY = (int)round (snkY + 0.49 + (1.0 * i / distMax * dstY));
      if (barriers[pxlX][pxlY] == 1)
      {
	barFound = true;
	break;
      }
    }
    if (!barFound)
    {
      goto End_of_Routine;
    }
    /*
    ** Barrier DOWN_RIGHT
    */
    barFound = false;
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX + 0.49 + (1.0 * i / distMax * dstX));
      pxlY = (int)round (snkY + 0.49 + (1.0 * i / distMax * dstY));
      if (barriers[pxlX][pxlY] == 1)
      {
	barFound = true;
	break;
      }
    }
    if (!barFound)
    {
      goto End_of_Routine;
    }
  }
  else if (barrierType == STRONG_BARRIER)
  {
    /*
    ** Strong barrier: If more than one way is blocked by a barrier then
    **                 colonization fails.
    */
    barCounter = 0;
    barFound = false;
    /*
    ** BARRIER MIDDLE
    */
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX + (1.0 * i / distMax * dstX));
      pxlY = (int)round (snkY + (1.0 * i / distMax * dstY));
      if (barriers[pxlX][pxlY] == 1)
      {
	barCounter++;
	break;
      }
    }
    /*
    ** BARRIER TOP_LEFT
    */
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX - 0.49 + (((i-1.0) / distMax * dstX) +
					((1.0 / distMax * dstX) / 2.0)));
      pxlY = (int)round (snkY - 0.49 + (((i-1.0) / distMax * dstY) +
					((1.0 / distMax * dstY) / 2.0)));
      if (barriers[pxlX][pxlY] == 1)
      {
	barCounter++;
	break;
      }
    }
    if (barCounter > 1)
    {
      barFound = true;
      goto End_of_Routine;
    }
    /*    
    ** BARRIER TOP_RIGHT
    */
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX + 0.49 + (((i-1.0) / distMax * dstX) +
					((1.0 / distMax * dstX) / 2.0)));
      pxlY = (int)round (snkY - 0.49 + (((i-1.0) / distMax * dstY) +
					((1.0 / distMax * dstY) / 2.0)));
      if (barriers[pxlX][pxlY] == 1)
      {
	barCounter++;
	break;
      }
    }
    if (barCounter > 1)
    {
      barFound = true;
      goto End_of_Routine;
    }
    /*
    ** BARRIER DOWN_LEFT
    */
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX - 0.49 + (((i-1.0) / distMax * dstX) +
					((1.0 / distMax * dstX) / 2.0)));
      pxlY = (int)round (snkY + 0.49 + (((i-1.0) / distMax * dstY) +
					((1.0 / distMax * dstY) / 2.0)));
      if (barriers[pxlX][pxlY] == 1)
      {
	barCounter++;
	break;
      }
    }
    if (barCounter > 1)
    {
      barFound = true;
      goto End_of_Routine;
    }
    /*
    ** BARRIER DOWN_RIGHT
    */
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX + 0.49 + (((i-1.0) / distMax * dstX) +
					((1.0 / distMax * dstX) / 2.0)));
      pxlY = (int)round (snkY + 0.49 + (((i-1.0) / distMax * dstY) +
					((1.0 / distMax * dstY) / 2.0)));
      if (barriers[pxlX][pxlY] == 1)
      {
	barCounter++;
	break;
      }
    }
    if (barCounter > 1)
    {
      barFound = true;
      goto End_of_Routine;
    }
  }        

 End_of_Routine:
  /*
  ** Return the result.
  */
  return (barFound);
}
    

/*
** EoF: barriers.c
*/
