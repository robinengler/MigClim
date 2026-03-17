
#include "migclim.h"


int mcUnivDispCnt (int **habSuit){
  /****************************************************************************
  ** Count the number of cells in the habitat suitability matrix that are
  ** suitable, and that would thus become colonized in the case of unlimited
  ** dispersal.
  **
  ** Parameters:
  **   - habSuit: A pointer to the matrix that contains the current
  **              habitat suitability.
  ** Returns:
  **   >- the number of suitable and non-barrier pixels.
  ** *********************************************************************** */
  int i, j, count;

  /* Count the number of suitable and non-barrier cells. */
  count = 0;
  for (i = 0; i < nrRows; i++){
    for (j = 0; j < nrCols; j++){
      if (habSuit[i][j] > 0) count++;
    }
  }
  return (count);
}


void updateNoDispMat (int **hsMat, int **noDispMat, int *noDispCount){
  /****************************************************************************
  ** Update the "no dispersal matrix" with the  habitat suitability values
  ** that are contained in the current habitat suitability matrix.
  **
  ** Parameters:
  **   - hsMat      : pointer to the habitat suitability matrix.
  **   - niDispMat  : pointer to the no-dispersal matrix.
  **   - noDispCount: pointer to the no-dispersal count variable (its value
  **                  will be updated!).
  ** *********************************************************************** */
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
