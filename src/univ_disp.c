/*
** univ_disp.c: Function for performing the universal dispersal count.
** Wim Hordijk & Robin Engler:   Last modified: 11 May 2012
*/

#include "migclim.h"


/*
** mcUnivDispCnt: This function return the number of pixels in the habitat
**                suitability matrix that are suitable (and that would thus
**                become colonized in the case of unlimited dispersal).   
** Parameters:
**   - habSuit: A pointer to the matrix that contains the current
**              habitat suitability.
** Returns:
**   The number of suitable and non-barrier pixels.
*/

int mcUnivDispCnt (int **habSuit)
{
  int i, j, count;

  /* Count the number of suitable and non-barrier pixels. */
  count = 0;
  for (i = 0; i < nrRows; i++){
    for (j = 0; j < nrCols; j++){
      if (habSuit[i][j] > 0) count++;
    }
  }
	
  /* Return the result. */
  return (count);
}




/*
** updateNoDispMat: This function updates the "NoDispersal_Matrix" with the
**                  habitat suitability values that are contained in the
**                  current "HS_Matrix".
**
** Parameters:
**   - hsMat:       A pointer to the habitat suitability matrix.
**   - niDispMat:   A pointer to the no-dispersal matrix.
**   - noDispCount: A pointer to the no-dispersal count variable (its value
**                  will be updated!).
*/

void updateNoDispMat (int **hsMat, int **noDispMat, int *noDispCount)
{                   
  int i, j;

  if (*noDispCount > 0){
    for(i = 0; i < nrRows; i++){
      for(j = 0; j < nrCols; j++){
	    if((noDispMat[i][j] == 1) && (hsMat[i][j] == 0)){
	      noDispMat[i][j] = 0;
	      (*noDispCount)--;
	    }
      }
    }
  }
}


 
/*
** EoF: univ_disp.c
*/
