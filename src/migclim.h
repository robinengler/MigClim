/*
** migclim.h: Header file for the MigClim methods.
** Wim Hordijk   Last modified: 11 May 2012 (RE)
*/

#ifndef _MIGCLIM_H_
#define _MIGCLIM_H_

/*
** Include files.
*/
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <R.h>


/*
** Defines.
**
** UNIF01:         Draw a uniform random number in [0;1]. Note that 'random'
**                 does not work on Windows %-/  so we use 'rand' instead.
** WEAK_BARRIER:   Weak barrier type.
** STRONG_BARRIER: Strong barrier type.
*/
#define UNIF01         ((double)rand () / RAND_MAX)
#define WEAK_BARRIER   1
#define STRONG_BARRIER 2


/*
** Global variables (we just use many global var's here to avoid passing too
** many arguments all the time).
*/

extern int     nrRows, nrCols, envChgSteps, dispSteps, dispDist, iniMatAge,
               fullMatAge, rcThreshold, barrierType, lddMinDist, lddMaxDist, 
               noData, replicateNb;
extern double *dispKernel, *propaguleProd, lddFreq, xllCorner, yllCorner,
               cellSize;
extern char    iniDist[128], hsMap[128], simulName[128], barrier[128];
extern bool    useBarrier, fullOutput;


/*
** Function prototypes.
*/
void mcMigrate           (char **paramFile, int *nrFiles);
bool mcSrcCell           (int i, int j, int **curState, int **pxlAge,
			  int loopID, int habSuit, int **barriers);
int  mcUnivDispCnt       (int **habSuit);
void updateNoDispMat     (int **hsMat, int **noDispMat, int *noDispCount);
void mcFilterMatrix      (int **inMatrix, int **filterMatrix, bool filterNoData, bool filterOnes, bool insertNoData);
bool mcIntersectsBarrier (int snkX, int snkY, int srcX, int srcY, int **barriers);
int  mcInit              (char *paramFile);
int  readMat             (char *fName, int **mat);
int  writeMat            (char *fName, int **mat);
void genClust            (int *nrow, int *ncol, int *ncls, int *niter, int *thrs, char **suitBaseName,
                          char **barrBaseName, char **outBaseName, char **initFile);
void validate            (char **obsFileName, int *npts, char **simFileName, int *ncls, double *bestScore);



#endif  /* _MIGCLIM_H_ */

/*
** EoF: migclim.h
*/
