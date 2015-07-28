/*
** migrate.c: Implementation of the main MigClim function.
** Wim Hordijk & Robin Engler:  Last modified: 11 May 2012 (RE)
*/

#include "migclim.h"


/*
** Global variables, so we won't have to pass too many arguments all the time.
*/
int     nrRows, nrCols, envChgSteps, dispSteps, dispDist, iniMatAge,
        fullMatAge, rcThreshold, barrierType, lddMinDist, lddMaxDist,
        replicateNb;
double *dispKernel, *propaguleProd, lddFreq;
char    iniDist[128], hsMap[128], simulName[128], barrier[128];
bool    useBarrier, fullOutput;
typedef struct _pixel
{
  int row, col;
} pixel;
pixel rndPixel;


/*
** Function prototypes.
*/
void mcRandomPixel   (pixel *pix);
bool mcSinkCellCheck (pixel pix, int **curState, int **habSuit);


/*
** mcMigrate: The core of the MigClim method. Perform the main migration steps.
**            Parameter values are read from a file.
**
** Parameters:
**   - paramFile: The name of the parameter file.
**   - nrFiles:   A pointer to an integer to contain the number of output
**                files created. A value of -1 is returned if an error occurred.
*/

void mcMigrate (char **paramFile, int *nrFiles)
{
  int     i, j, RepLoop, envChgStep, dispStep, loopID, simulTime;  
  bool    advOutput, habIsSuitable, cellInDispDist, tempResilience;
  char    fileName[128], simulName2[128];
  FILE   *fp=NULL, *fp2=NULL;
  double  lddSeedProb;
  time_t  startTime;
  /*
  ** These variables are not (yet) used.
  **
  ** int vegResTime, seedBankResTime, *seedBank, tSinceDecol, tWhenCol;
  */
  
  /*
  ** Pixel counter variables. These variables allow us to record interesting
  ** values to be printed into the output file:
  **
  ** ==> Many of these aren't used yet. <==
  **
  **   - nrColonized:           The number of pixels in "Colonized" state at
  **                            any given time.
  **   - nrAbsent:              The number of pixels in "Absent" state at any
  **                            given time.
  **   - nrDecolonized:         The number of pixels in "Decolonized" state at
  **                            any given time.
  **   - nrTmpResilient:        The number of pixels that in "Temporary
  **                            Resilience" state at any given time.
  **   - nrVegResilient:        The number of pixels that are in "Vegetative
  **                            Resilience" state at any given time.
  **   - nrSeedBank:            The number of pixels that are in Seed Bank
  **                            Resilience state at any given time.
  **   - nrStepColonized:       The number of pixels which changed to
  **                            "Colonized" state during a given step.
  **   - nrStepDecolonized:     The number of pixels which changed to
  **                            "Decolonized" state during a given step.
  **   - nrStepTmpResilient:    The number of pixels that which changed to
  **                            "Temporary Resilience" state during a given
  **                            step.
  **   - nrStepVegResilient:    The number of pixels that are which changed to
  **                            "Vegetative Resilience" state during a given
  **                            step.
  **   - nrStepSeedBank:        The number of pixels that are which changed to
  **                            "Seed Bank Resilience" state during a given
  **                            step.
  **   - nrStepLostUnsuit:      The number of pixels that were lost due to
  **                            unsuitable habitat conditons during a given
  **                            dispersal step.
  **   - nrStepLostExtreme:     The number of pixels that were lost due to
  **                            "Extreme Event" during a given dispersal step.
  **   - nrStepLDDSuccess:      The number of LDD events that were successful
  **                            (can be any of a "normal" LDD or a River or
  **                            Road LDD).
  **   - nrStepVegResRecover:   The number of "Vegetative Resilience" recovery
  **                            events that occured within a given dispersal
  **                            step.
  **   - nrStepSeedBankRecover: The number of "Seed Bank Resilience" recovery
  **                            events that occured within a given dispersal
  **                            step.
  **   - nrTotColonized:        \
  **   - nrTotLostUnsuit:        |
  **   - nrTotLostExtreme:       |   Same as above but over the total
  **   - nrTotLostNatural:        >  simulation instead of one step.
  **   - nrTotLDDSuccess:        |
  **   - nrTotVegResRecover:     |
  **   - nrTotSeedBankRecover:  /
  **   - nrInitial:             The number of initial pixels that are occupied
  **                            by the species.
  **   - nrNoDispersal:         The number of pixels that would be colonized at
  **                            the end of the simulation under the
  **                            "no-dispersal" hypothesis.
  **   - nrUnivDispersal:       The number of pixels that would be colonized
  **                            at the end of the simulation under the
  **                            "unlimited-dispersal" hypothesis.
  */
  int nrColonized, nrAbsent, nrStepColonized, nrStepDecolonized, nrStepLDDSuccess, 
      nrTotColonized, nrTotLDDSuccess, nrInitial, nrNoDispersal, nrUnivDispersal, nrTotDecolonized;
      
  /* As of yet unused count variables:
  ** int nrTotVegResRecover, nrTotSeedBankRecover, nrStepVegResilient,
  **     nrStepSeedBank, nrVegResilient, nrStepSeedBankRecover,
  **     nrStepVegResRecover, nrSeedBank, nrDecolonized; */

  
  /* Matrices:
  **   - currentState:   Values in [-32768;32767]. NoData values are represented by -9999
  **   - habSuitability: Values in [0;1000].
  **   - barriers:       Values in [0;255].
  **   - pixelAge:       Values in [0;255].
  **   - noDispersal:    Values in [0;255].
  */
  int **currentState, **habSuitability, **barriers, **pixelAge, **noDispersal;

  
  /* Initialize the variables. */
  advOutput = false;
  currentState = NULL;
  habSuitability = NULL;
  barriers = NULL;
  pixelAge = NULL;
  noDispersal = NULL;
  propaguleProd = NULL;
  dispKernel = NULL;
  if(mcInit(*paramFile) == -1){ /* Reads the "_param.txt" file */
    *nrFiles = -1;
    goto End_of_Routine;
  }
  
  /* We'll use 'srand' and 'rand' here, as 'srandom' and 'random' do not work on Windows :-( */
  srand (time (NULL));

  /* Allocate the necessary memory. */
  currentState = (int **)malloc (nrRows * sizeof (int *));
  habSuitability = (int **)malloc (nrRows * sizeof (int *));
  barriers = (int **)malloc (nrRows * sizeof (int *));
  pixelAge = (int **)malloc (nrRows * sizeof (int *));
  noDispersal = (int **)malloc (nrRows * sizeof (int *));
  for(i = 0; i < nrRows; i++){
    currentState[i] = (int *)malloc (nrCols * sizeof (int));
    habSuitability[i] = (int *)malloc (nrCols * sizeof (int));
    barriers[i] = (int *)malloc (nrCols * sizeof (int));
    pixelAge[i] = (int *)malloc (nrCols * sizeof (int));
    noDispersal[i] = (int *)malloc (nrCols * sizeof (int));
  }
  
  /* Replicate the simulation replicateNb times. If replicateNb > 1 then the
  ** simulation's output names are "simulName1", "simulName2", etc... */
  for(RepLoop = 1; RepLoop <= replicateNb; RepLoop++){
    
	/* Remember the current time */
    startTime = time(NULL);

    /* If replicateNb > 1 then we need to change the simulation name */	  
    if(replicateNb == 1){
      strcpy(simulName2, simulName);
    }
    else if(replicateNb > 1){
      sprintf(simulName2, "%s%d", simulName, RepLoop);           /* sprinf(): puts a string into variable simulName2 */
    }

    
    /* Load and prepare the data. */
    
    /* Species initial distribution */
    sprintf(fileName, "%s.asc", iniDist);                        /* sprinf(): puts a string into variable fileName */
    if(readMat(fileName, currentState) == -1){
      *nrFiles = -1;
      goto End_of_Routine;
    }
    
    /* Barrier options */
    for(i = 0; i < nrRows; i++){
      for(j = 0; j < nrCols; j++){
	    barriers[i][j] = 0;
      }
    }
    if(useBarrier){
      sprintf(fileName, "%s.asc", barrier);
      if(readMat(fileName, barriers) == -1){
	    *nrFiles = -1;                                           /* if readMat() return -1, an error occured  */
	    goto End_of_Routine;
      }
    } 
    /* Filter the barrier matrix in two ways:
    **  -> reclass any value < 0 as 0 (this is to remove NoData values of -9999).
    **  -> set the cells with NoData in 'currentState' to NoData in 'barriers'
    **     so that the NoData in 'currentState' and 'barriers' are identical */
    mcFilterMatrix(barriers, currentState, true, false, true);
    
    /* Filter the values of current state matrix by the barriers matrix
    ** (when barriers = 1 we set currentState = 0) */
    if(useBarrier) mcFilterMatrix(currentState, barriers, false, true, false);
    
    
    
    /* Time to reach dispersal maturity.
    **
    ** Before a pixel (= a population) can disperse, it needs to reach a certain
    ** user-defined age. The "age" is a matrix that keeps track of the age
    ** of all pixels:
    **   0 = pixel is not colonized.
    **   1 = pixel is colonized since 1 "dispersal event".
    **   2 = pixel is colonized since 2 "dispersal event".
    **   3 = etc...
    ** Fill the "age" to reflect the initial distribution of the species:
    ** where the species is present, pixels get a value of 'FullMaturity', where
    ** the species is absent, pixels get a value of 0. */
    for(i = 0; i < nrRows; i++){
      for(j = 0; j < nrCols; j++){
	    if(currentState[i][j] == 1){
		  pixelAge[i][j] = fullMatAge;
	    }
	    else{
	      pixelAge[i][j] = 0;
	    }
      }
    }
    
    
    /* The "no dispersal" matrix.
    ** This Matrix will keep track of the species distribution under the
    ** "no dispersal" scenario. */  
    for(i = 0; i < nrRows; i++){
      for (j = 0; j < nrCols; j++){
	    noDispersal[i][j] = currentState[i][j];
      }
    }
        
        
    /* Initialize counter variables.
    ** Reset pixel counters to zero before we start the dispersal simulation. */
    nrInitial = 0;
    nrColonized = 0;
    nrAbsent = 0;
    nrNoDispersal = 0;
    nrUnivDispersal = 0;
    nrTotColonized = 0;
    nrTotDecolonized = 0;
    nrTotLDDSuccess = 0;
    tempResilience = true;
    nrStepColonized = 0;
    nrStepDecolonized = 0;
    nrStepLDDSuccess = 0;
    /* Currently unused counters:
    ** nrVegResilient = 0;
    ** nrSeedBank = 0;
    ** nrTotVegResRecover = 0;
    ** nrTotSeedBankRecover = 0; */

    
    /* Count the number of initially colonized pixels (i.e. initial species distribution)
    ** as well as the number of empty (absence) cells. */
    for(i = 0; i < nrRows; i++){
      for(j = 0; j < nrCols; j++){
	    if(currentState[i][j] == 1) nrInitial++;
        if(currentState[i][j] == 0) nrAbsent++;
      }
    }
    nrColonized = nrInitial;
    nrNoDispersal = nrInitial;
    nrUnivDispersal = nrInitial;

    
    /* Write the initial state to the data file. */
    sprintf(fileName, "%s/%s_stats.txt", simulName, simulName2);
    if((fp = fopen (fileName, "w")) != NULL){
      fprintf (fp, "envChgStep\tdispStep\tstepID\tunivDispersal\tNoDispersal\toccupied\tabsent\tstepColonized\tstepDecolonized\tstepLDDsuccess\n");
      fprintf (fp, "0\t0\t1\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", nrUnivDispersal, nrNoDispersal,
	           nrColonized, nrAbsent, nrStepColonized, nrStepDecolonized, nrStepLDDSuccess);
    }
    else{
      *nrFiles = -1;
      Rprintf ("Could not open statistics file for writing.\n");
      goto End_of_Routine;
    }
    
    
    
    /* **************************************************************** */
    /* Simulate plant dispersal and migration (the core of the method). */
    /* **************************************************************** */
    Rprintf("Running MigClim simulation %s.\n", simulName2);
    
    /* Start of environmental change step loop (if simulation is run without change in environment this loop runs only once). */
    for(envChgStep = 1; envChgStep <= envChgSteps; envChgStep++){
	  
      /* Print the current environmental change iteration. */
      Rprintf ("  %d...\n", envChgStep);

      /* Load the habitat suitability layer for the current envChgStep. */
      sprintf (fileName, "%s%d.asc", hsMap, envChgStep);
      if (readMat (fileName, habSuitability) == -1){
	    *nrFiles = -1;
	    goto End_of_Routine;
      }
      
      /* if the user chose a rcThreshold > 0, reclass the habitat suitability into 0 or 1000
      ** if rcThreshold == 0 then suitability values are left unchanged. */
      if(rcThreshold > 0){
	    for(i = 0; i < nrRows; i++){
	      for(j = 0; j < nrCols; j++){
	        if(habSuitability[i][j] < rcThreshold){
	          habSuitability[i][j] = 0;
	        }
	        else{
	          habSuitability[i][j] = 1000;
	        }
	      }
	    }
      }
      
      /* Filter the habitat suitability matrix in three ways:
      **  -> replace any value < 0 by 0 (this removes NoData).
      **  -> set habitat suitability to 0 where barrier = 1.
      **  -> set habitat suitability values to NoData where barrier = NoData. 
      ** and also  */
      mcFilterMatrix(habSuitability, barriers, true, true, true);
      
      
      /* Set the values that will keep track of pixels colonized during the next
      ** climate change loop.
      ** "loopID" is the value that will be given to the pixel colonized
      ** during the current loop. The pixels being decolonized during the
      ** current loop will be given a value of "-loopID". The limits in
      ** terms of number of dispersal events and environmental change
      ** events is the following:
      **  - Without environmental change, the maximum number of dispStep
      **    is 29'500.
      **  - With environmental change, the maximum number of dispStep is
      **    98, the number of environmental change loops is limited to 250.
      ** The coding for the loopID is as follows: envChgStep * 100 + dispStep. */
      if(envChgStep == 0){
        loopID = 1;
      }
      else{
	    loopID = envChgStep * 100;
      }
      
      
      /* "Unlimited" and "no dispersal" scenario pixel count. Here we compute the number of pixels
      ** that would be colonized if dispersal was unlimited or null. This is simply the sum of all
      ** potentially suitable habitats */
      nrUnivDispersal = mcUnivDispCnt(habSuitability);
      updateNoDispMat(habSuitability, noDispersal, &nrNoDispersal);

      /* Reset number of decolonized cells within current dispersal step pixel counter */
	  nrStepDecolonized = 0;
	    
      /* Update for temporarily resilient pixels. */
      for(i = 0; i < nrRows; i++){
	    for(j = 0; j < nrCols; j++){
	      
	      /* Udate non-suitable pixels. If a pixel turned unsuitable, we update its status to "Temporarily Resilient". */
	      if((habSuitability[i][j] == 0) && (currentState[i][j] > 0)){
	        
	        /* If the user selected TemporaryResilience==T, then the pixel is set to "Temporary Resilient" status. */
	        if(tempResilience == true){
	          currentState[i][j] = 29900;
	        }
	        else{
	          /* If not temporary resilience was specified, then the pixel is set to "decolonized" status. */
	          currentState[i][j] = -1 - loopID;
	          pixelAge[i][j] = 0;
	          /* NOTE: Later we can add "Vegetative" and "SeedBank" resilience options at this location. */
	        }
	        
	        /* The number of decolonized cells within current step is increased by one */
	        nrStepDecolonized++;
	      }
	    }
      }

      
      /* *************************************** */
      /* Dispersal event step loop starts here.  */
      /* *************************************** */
      for(dispStep = 1; dispStep <= dispSteps; dispStep++){
	    
	    /* Set the value of "loopID" for the current iteration of the dispersal loop. */
	    loopID++;
          
	    /* Reset pixel counters that count pixels within the current loop. */
	    nrStepColonized = 0;
	    nrStepLDDSuccess = 0;
	    if(dispStep > 1) nrStepDecolonized = 0;

	    /* Currently unused variables:
	    ** nrStepVegResilient = 0;
	    ** nrStepSeedBank = 0;
	    ** nrStepVegResRecover = 0;
	    ** nrStepSeedBankRecover = 0; */
          
	    /* Source cell search: Can the sink pixel be colonized? There are four
	    ** conditions to be met for a sink pixel to become colonized:
	    **   1. Sink pixel is currently suitable and not already colonised.
	    **   2. Sink pixel is within dispersal distance of an already colonised
	    **      and mature pixel.
	    **   3. Source pixel has reached dispersal maturity.
	    **   4. There is no obstacle (barrier) between the pixel to be colonised
	    **      (sink pixel) and the pixel that is already colonised (source
	    **      pixel).
	    **
	    ** Loop through the cellular automaton. */
	    for(i = 0; i < nrRows; i++){
	      for(j = 0; j < nrCols; j++){
	        
		    /* Reset variables. */
	        habIsSuitable = false;
	        cellInDispDist = false;
	        
	        /* 1. Test whether the pixel is a suitable sink (i.e., its habitat
	        **    is suitable, it's unoccupied and is not on a barrier or filter
	        **    pixel). */
	        if((habSuitability[i][j] > 0) && (currentState[i][j] <= 0)) habIsSuitable = true;

	        /* 2. Test whether there is a source cell within the dispersal
	        **    distance. To be more time efficient, this code runs only if
	        **    the answer to the first question is positive. Additionally,
	        **    if there is a source cell within dispersion distance and if
	        **    a barrier was asked for, then we also check that there is no
	        **    barrier between this source cell and the sink cell (this
	        **    verification is carried out in the "SearchSourceCell"
	        **    function). */
	        if(habIsSuitable){
		      /* Now we search if there is a suitable source cell to colonize the sink cell. */
	          if (mcSrcCell (i, j, currentState, pixelAge, loopID, habSuitability[i][j], barriers)) cellInDispDist = true;
	        }
	        
	        /* Update pixel status. */
	        if(habIsSuitable && cellInDispDist){
	          
		      /* Only if the 2 conditions are fullfilled the cell's status is set to colonised. */
	          currentState[i][j] = loopID;
	          nrStepColonized++;
	          
	          /* If the pixel was in seed bank resilience state, then we
	          ** update the corresponding counter. Currently not used.
	          ** if (pixelAge[i][j] == 255) nrStepSeedBank--; */
	          
	          /* Update "age" value. We do this only now because we needed the old "age" value just before 
	          ** to determine whether a pixel was in "Decolonized" or "SeedBank resilience" status. */
	          pixelAge[i][j] = 0;
	        }
	      }
	    }
        
	    /* If the LDD frequence is larger than zero, perform it. */
	    if(lddFreq > 0.0){
	      
		  /* Loop through the entire cellular automaton. */
	      for(i = 0; i < nrRows; i++){
	        for(j = 0; j < nrCols; j++){
	          
		      /* Check if the pixel is a source cell (i.e. is it colonised since at least 1 dispersal Loop) 
	          ** and check if the pixel has reached dispersal maturity. */
	          if((currentState[i][j]) > 0 && (currentState[i][j] != loopID)){
	    	    if(pixelAge[i][j] >= iniMatAge){
	    	  
		    	  /* Set the probability of generating an LDD event. This
	    	      ** probability is weighted by the age of the cell. */
	    	      if(pixelAge[i][j] >= fullMatAge){
	    	        lddSeedProb = lddFreq;
	    	      }
	    	      else{
	    	        lddSeedProb = lddFreq * propaguleProd[pixelAge[i][j] - iniMatAge];
	    	      }
	    	  
	    	      /* Now we can try to generate a LDD event with the calculated probability. */
	    	      if(UNIF01 < lddSeedProb || lddSeedProb == 1.0){
	    	        
		    	    /* Randomly select a pixel within the distance "lddMinDist - lddMaxDist". */
	    	        mcRandomPixel (&rndPixel);
	    	        rndPixel.row = rndPixel.row + i;
	    	        rndPixel.col = rndPixel.col + j;
	    	        
	    	        /* Now we check if this random cell is a suitable sink cell.*/
	    	        if(mcSinkCellCheck (rndPixel, currentState, habSuitability)){
	    	          
		    	      /* if condition is true, the pixel gets colonized.*/
	    	          currentState[rndPixel.row][rndPixel.col] = loopID;
	    	          nrStepColonized++;
	    	          nrStepLDDSuccess++;
	    	          
	    	          /* If the pixel was in seed bank resilience state, then we
	    	          ** update the corresponding counter. Currently not used.
	    	          ** if (pixelAge[rndPixel.row][rndPixel.col] == 255) nrStepSeedBank--;  */
	    	    
	    	          /* Reset pixel age. */
	    	          pixelAge[rndPixel.row][rndPixel.col] = 0;
	    	        }
	    	      }
	    	      
	    	    }
	          }
	        }
	      }
	    }
            
	    /* Update pixel age: At the end of a dispersal loop we want to
	    ** increase the "age" of each colonized pixel.
	    **
	    ** Reminder: pixel "age" structure is as follows:
	    **   0 = Pixel is either "Absent", "Decolonized" or has just been
	    **       "Colonized" during this dispersal step.
	    **   1 to 250 = Pixel is in "Colonized" or "Temporarily Resilient"
	    **       status. The value indicates the number of "dispersal events
	    **       (usually years) since when the pixel was colonized.
	    **   255 = Pixel is in "SeedBank Resilience" state. */
	    for(i = 0; i < nrRows; i++){
	      for(j = 0; j < nrCols; j++){
	        
		    /* If the pixel is in "Colonized" or "Temporarily Resilient" state, update it's age value. */
	        if(currentState[i][j] > 0) pixelAge[i][j] += 1;

	        /* If a pixel is in "Temporarily Resilient" state, we also increase its "currentState" value by 1
	        ** so that the pixels gains 1 year of "Temporarily Resilience" age. */
	        if (currentState[i][j] >= 29900) currentState[i][j] += 1;

	      }
	    }
        
	    /* Update pixel counters. */
	    nrColonized = nrColonized + nrStepColonized - nrStepDecolonized;
	    nrAbsent = nrAbsent - nrStepColonized + nrStepDecolonized;
	    nrTotColonized += nrStepColonized;
	    nrTotDecolonized += nrStepDecolonized;
	    nrTotLDDSuccess += nrStepLDDSuccess;
	    /* Currently unused variables:
	    ** nrTotVegResRecover += nrStepVegResRecover;
	    ** nrTotSeedBankRecover += nrStepSeedBankRecover; */
          
	    /* Write current iteration data to the statistics file. */
	    fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", envChgStep, dispStep, loopID, nrUnivDispersal,
	    	    nrNoDispersal, nrColonized, nrAbsent, nrStepColonized, nrStepDecolonized, nrStepLDDSuccess);
	    	 
	    /* If the user has requested full output, also write the current state matrix to file. */
	    if(fullOutput){
	      sprintf (fileName, "%s/%s_step_%d.asc", simulName, simulName2, loopID);
	      if(writeMat (fileName, currentState) == -1){
	        *nrFiles = -1;
	        goto End_of_Routine;
	      }
	    }
      } /* END OF: dispStep */
    
      
      /* Update temporarily resilient pixels.
      ** Temporarily resilient pixels can be distinguished by:
      **   -> CurrentState_Matrix = 29'900 to 29'999. Increases by 1 at each year.
      **   -> Age_Matrix has a positive value. */
      for (i = 0; i < nrRows; i++){
	    for (j = 0; j < nrCols; j++){
	      if (currentState[i][j] >= 29900){
	        currentState[i][j] = dispSteps - loopID - 1;
	        pixelAge[i][j] = 0;
	      }
	    }
      }
    
    } /* END OF: envChgStep loop */
    Rprintf("All dispersal steps completed. Final output in progress...\n");
  
    
    /* Update currentState matrix for pixels that are suitable but
    ** could not be colonized due to dispersal limitations.
    ** These pixels are assigned a value of 30'000 */
    for(i = 0; i < nrRows; i++){
      for(j = 0; j < nrCols; j++){
	    if((habSuitability[i][j] > 0) && (currentState[i][j] <= 0)) currentState[i][j] = 30000;
      }
    }
  
    /* Write the final state matrix to file. */
    sprintf(fileName, "%s/%s_raster.asc", simulName, simulName2);
    if(writeMat (fileName, currentState) == -1){
      *nrFiles = -1;
      goto End_of_Routine;
    }
  
    /* Write summary output to file. */
    simulTime = time (NULL) - startTime;
    sprintf(fileName, "%s/%s_summary.txt", simulName, simulName2);
    if((fp2 = fopen (fileName, "w")) != NULL){
      fprintf(fp2, "simulName\tiniCount\tnoDispCount\tunivDispCount\toccupiedCount\tabsentCount\ttotColonized\ttotDecolonized\ttotLDDsuccess\trunTime\n");
      fprintf(fp2, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", simulName2, nrInitial, nrNoDispersal, nrUnivDispersal,
	          nrColonized, nrAbsent, nrTotColonized, nrTotDecolonized, nrTotLDDSuccess, simulTime);
      fclose (fp2);
    }
    else{
      *nrFiles = -1;
      Rprintf ("Could not write summary output to file.\n");
      goto End_of_Routine;
    }  
  
    /* Close the data file. */
    if (fp != NULL) fclose (fp);
    fp = NULL;
    
  } /* end of "RepLoop" */
  
  
  /* Set the number of output files created. */
  *nrFiles = envChgSteps;
 
  
 
 End_of_Routine:
  
  /* Close the data file. */
  if(fp != NULL) fclose (fp);

  /* Free the allocated memory. */
  if(currentState != NULL){
    for (i = 0; i < nrRows; i++) free(currentState[i]);
    free(currentState);
  }
  if(habSuitability != NULL){
    for(i = 0; i < nrRows; i++) free(habSuitability[i]);
    free(habSuitability);
  }
  if(barriers != NULL){
    for(i = 0; i < nrRows; i++) free(barriers[i]);
    free(barriers);
  }
  if(pixelAge != NULL){
    for(i = 0; i < nrRows; i++) free(pixelAge[i]);
    free(pixelAge);
  }
  if(noDispersal != NULL){
    for (i = 0; i < nrRows; i++) free(noDispersal[i]);
    free(noDispersal);
  }
  if (dispKernel != NULL) free(dispKernel);
  if (propaguleProd != NULL) free(propaguleProd);

  
  /* If an error occured, display failure message to the user... */
  if(*nrFiles == -1) Rprintf("MigClim simulation aborted...\n");
  
}






/*
** mcRandomPixel: Select a random pixel from a central point (0;0) and within a
**                radius of at least lddMinDist and at most lddMaxDist.
**
** Parameters:
**   - pix:     A pointer to the pixel data structure to put the row and col in.
*/

void mcRandomPixel (pixel *pix)
{
  double rndDist, rndAngle;
  
  /*
  ** Select a random distance between lddMinDist and lddMaxDist, and a random
  ** angle between 0 and 2*pi.
  */
  rndDist = (UNIF01 * (lddMaxDist - lddMinDist)) + lddMinDist;
  rndAngle = UNIF01 * 6.283185;
    
  /*
  ** Convert distance and angle into pixel row and column values.
  */
  pix->row = (int)(rndDist * cos (rndAngle));
  pix->col = (int)(rndDist * sin (rndAngle));
}





/*
** mcSinkCellCheck: Perform a basic check to see whether a given cell fulfills
**                  the conditions to be a "sink" cell (i.e. a cell to be
**                  potentially colonized).
**                  The conditions that are checked are the following:
**                   -> 1. Cell must be within the limits of the cellular automaton.
**                   -> 2. Cell must be empty (i.e. non-occupied).
**                   -> 3. Cell must contain suitable habitat.
**
** Parameters:
**   -> pix: The cell/pixel to consider.
**   -> curState: A pointer to the current state matrix.
**   -> habSuit:  A pointer to the habitat suitability matrix.
**
** Returns:
**   If the cell is suitable: true.
**   Otherwise:               false.
*/

bool mcSinkCellCheck (pixel pix, int **curState, int **habSuit)
{
  bool suitable;
  double rnd;
  suitable = false;

  /* 1. Verify the cell is within the limits of the cellular automaton. */
  if((pix.row < 0) || (pix.row >= nrRows) || (pix.col < 0) || (pix.col >= nrCols)) return(suitable);

  /* 2. Verify the cell is empty. */
  if(curState[pix.row][pix.col] > 0) return(suitable);

  /* 3. Verify the cell contains suitable habitat for the species. */
  if(habSuit[pix.row][pix.col] == 0) return(suitable);            /* check for case when habitat suitability = 0                       */
  rnd = UNIF01 * 1000;                                            /* generate random number: UNIF01 is declared in the header section  */
  if(rnd > (double)habSuit[pix.row][pix.col]) return(suitable);   /* The "(double)" is typecasting the value to a double               */
  
  /* If the function has not exited by now then it means the cell is suitable. */
  suitable = true;
  return (suitable);
}


/*
** EoF: migrate.c
*/
