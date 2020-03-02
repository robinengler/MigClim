/*
####################################################################################################
## migrate.c: implementation of the main MigClim function.
##
## Authors: Robin Engler and Wim Hordijk.
## Last modified: 30 Aug 2019.
####################################################################################################
*/


/* Load header file */
#include "migclim.h"


/* Declare global variables.
** ************************/
int    nrRows, nrCols, envChgSteps, dispSteps, dispDist, iniMatAge, fullMatAge, 
       rcThreshold, barrierType, lddMinDist, lddMaxDist, replicateNb;
float  *dispKernel; 
double *propaguleProd, lddFreq;
char   iniDist[128], hsMap[128], simulName[128], barrier[128], dipersalKernelBasename[128];
bool   useBarrier, fullOutput, singleKernelMode;
typedef struct _pixel{
    int row, col;
} pixel;
pixel rndPixel;


/* Function prototypes.
** *******************/
void mcRandomPixel   (pixel *pix);
bool mcSinkCellCheck (pixel pix, int **curState, int **habSuit);


void mcMigrate (char **paramFile, int *nrFiles){
    /* *******************************************************************************************
    ** The core of the MigClim method that performs species dispersal in the landscape. The 
    ** input parameter values are read from a file that is written to disk by the R code of 
    ** miglicm.
    **
    ** Parameters:
    **  -> paramFile: the name of the parameter file.
    **  -> nrFiles:   a pointer to an integer to contain the number of output files created. 
    **                A value of -1 is returned if an error occurred.
    ** Returns:
    ** ->
    ** ******************************************************************************************/

    int    i, j, RepLoop, envChgStep, dispStep, loopID, simulTime;
    bool   habIsSuitable, cellInDispDist, tempResilience;
    char   fileName[256], simulName2[128];
    FILE   *fp=NULL, *fp2=NULL;
    double lddSeedProb;
    time_t startTime;

    
    // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
    // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
    char tmpString[128];
    int tmpInt;
    FILE *fp3=NULL;
    // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
    // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
    
    
    /* Variables that were not ported to the R version of MigClim.
    ** int vegResTime, seedBankResTime, *seedBank, tSinceDecol, tWhenCol; */

    /* Pixel counter variables. These variables allow us to record interesting
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
        nrTotColonized, nrTotLDDSuccess, nrInitial, nrNoDispersal, nrUnivDispersal, 
        nrTotDecolonized;
        
    
    /* Matrices:
    **   - currentState:   Values in [-32768;32767]. NoData values are represented by -9999
    **   - habSuitability: Values in [0;1000].
    **   - barriers:       Values in [0;255].
    **   - pixelAge:       Values in [0;255].
    **   - noDispersal:    Values in [0;255]. */
    int **currentState, **habSuitability, **barriers, **pixelAge, **noDispersal;


    /* Initialize variables. */
    currentState = NULL;
    habSuitability = NULL;
    barriers = NULL;
    pixelAge = NULL;
    noDispersal = NULL;
    propaguleProd = NULL;
    dispKernel = NULL;
    singleKernelMode = true;

    /* Read the '_param.txt' file */
    if(mcInit(*paramFile) == -1){ 
        *nrFiles = -1;
        goto END_OF_FUNCTION;
    }


    /* Set seed for the random number generator. 
    ** ****************************************
    ** This is no longer done in the C code, as we now set the seed directly in the 
    ** MigClim.migrate() R function. This is possible because we now use the 'unif_rand()' 
    ** function from the Rmath.h library, instead of the regular 'rand()' function of C.
    ** The problem with using 'rand()' is that R also has its own 'rand()' C functions and 
    ** this creates problems when trying to set the seed for the random number generator. 
    ** The seed is set for the original rand() function but when calling rand() the function 
    ** from R is used and the previously set seed is ignored.
    ** An alternative to rand()/srand() would be to use random()/srandom(), which would avoid 
    ** the colliding with R's rand() function, but apparently random() does not work on windows 
    ** machines. */
    if(VERBOSE){
        Rprintf("Verbose mode: Test of R's random number generator:\n");
        Rprintf("Verbose mode: random number [%lf]\n", UNIF01);
        Rprintf("Verbose mode: random number [%lf]\n", UNIF01);
        Rprintf("Verbose mode: random number [%lf]\n", UNIF01);
    }


    /* Allocate memory to matrices.
    ** ***************************/
    currentState = (int **)malloc(nrRows * sizeof (int *));
    habSuitability = (int **)malloc(nrRows * sizeof (int *));
    barriers = (int **)malloc(nrRows * sizeof (int *));
    pixelAge = (int **)malloc(nrRows * sizeof (int *));
    noDispersal = (int **)malloc(nrRows * sizeof (int *));
    for(i = 0; i < nrRows; i++){
        currentState[i] = (int *)malloc(nrCols * sizeof (int));
        habSuitability[i] = (int *)malloc(nrCols * sizeof (int));
        barriers[i] = (int *)malloc(nrCols * sizeof (int));
        pixelAge[i] = (int *)malloc(nrCols * sizeof (int));
        noDispersal[i] = (int *)malloc(nrCols * sizeof (int));
    }


    
    /* *******************************************************************************************
    ** Simulation replication loop starts here.
    ** *******************************************************************************************
    ** Main loop of migclim that repeats the entire simulation process 'replicateNb' times. */
    for(RepLoop = 1; RepLoop <= replicateNb; RepLoop++){

        /* Remember the current time */
        startTime = time(NULL);

        /* If the user selected to run more than 1 replicate of the simulation, add a suffix with
        ** the replicate number to the simulation name. E.g. "testrun" becomes "testrun_1", 
        ** "testrun_2", etc. */
        if(replicateNb == 1){
            strncpy(simulName2, simulName, 128);
        }
        else if(replicateNb > 1){
            snprintf(simulName2, 128, "%s_%d", simulName, RepLoop);
        }


        /* Load and prepare the data.
        ** **************************
        ** Species initial distribution */
        sprintf(fileName, "%s.asc", iniDist);
        if(readMat(fileName, currentState) == -1){
            *nrFiles = -1;
            goto END_OF_FUNCTION;
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
                /* For readMat(), a return value of -1 indicates that an error occurred. */
                *nrFiles = -1;
                goto END_OF_FUNCTION;
            }
        } 

        /* Filter the barrier matrix in two ways:
        **  1. reclass any value < 0 as 0. This is to remove NoData values of -9999.
        **  2. set the cells with NoData in 'currentState' to NoData in 'barriers'
        **     so that the NoData in 'currentState' and 'barriers' are identical */
        mcFilterMatrix(barriers, currentState, true, false, true);

        /* Filter the values of current state matrix by the barriers matrix
        ** (when barriers = 1 we set currentState = 0) */
        if(useBarrier) mcFilterMatrix(currentState, barriers, false, true, false);

        
        /* Time to reach dispersal maturity.
        ** Before a pixel (= a population) can disperse, it needs to reach a certain user-defined
        ** age. pixelAge is a matrix that keeps track of the age of all cells:
        **   0 = cell is not colonized.
        **   1 = cell is colonized since 1 "dispersal event".
        **   2 = cell is colonized since 2 "dispersal event".
        **   3 = etc...
        ** Here we fill the pixelAge matrix to reflect the initial distribution of the species:
        ** where the species is present, cells get a value of 'FullMaturity', where the species 
        ** is absent, cells get a value of 0. */
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


        /* "no dispersal" matrix.
        ** ********************* 
        ** Keeps track of the species distribution under the "no dispersal" scenario. */
        for(i = 0; i < nrRows; i++){
        for (j = 0; j < nrCols; j++){
            noDispersal[i][j] = currentState[i][j];
        }
        }
        
        
        /* Initialize counter variables.
        ** ****************************
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


        /* Count the number of initially colonized pixels, i.e. the initial species distribution, 
         * as well as the number of empty (absence) cells. */
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
        if((fp = fopen(fileName, "w")) != NULL){
            fprintf(fp, "envChgStep\tdispStep\tstepID\tunivDispersal\tNoDispersal\toccupied\tabsent\tstepColonized\tstepDecolonized\tstepLDDsuccess\n");
            fprintf(fp, "0\t0\t1\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", nrUnivDispersal, nrNoDispersal, 
                    nrColonized, nrAbsent, nrStepColonized, nrStepDecolonized, nrStepLDDSuccess);
        }
        else{
            *nrFiles = -1;
            Rprintf ("Could not open statistics file for writing.\n");
            goto END_OF_FUNCTION;
        }

        
        /* ***************************************************************************************
        ** Simulate dispersal and migration - the core of the method.
        ** ***************************************************************************************
        ** Start the loop through all environmental change steps. If the simulation is run 
        ** without change in habitat suitability over time, then this loop runs only once. */
        Rprintf("Running MigClim simulation %s.\n", simulName2);
        for(envChgStep = 1; envChgStep <= envChgSteps; envChgStep++){

            /* Print the current environmental change iteration. */
            Rprintf("  %d...\n", envChgStep);
            
            
            
            
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            // Only load habitat suitability if envChgStep <= 4, because for >= 5 we keep the 
            // same values are 2050.
            if(envChgStep <= 4){
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            
            
            /* Load habitat suitability for the current environmental change step.
            ** ******************************************************************
            ** Read habitat suitability values from raster file. */
            sprintf(fileName, "%s%d.asc", hsMap, envChgStep);
            if(readMat(fileName, habSuitability) == -1){
                *nrFiles = -1;
                goto END_OF_FUNCTION;
            }
            /* If the user passed a reclassification threshold > 0, reclassify the habitat 
            ** suitability values into 0 or 1000. If rcThreshold == 0 habitat suitability 
            ** values are left unchanged. */
            if(rcThreshold > 0){
                for(i = 0; i < nrRows; i++){
                for(j = 0; j < nrCols; j++){
                    if(habSuitability[i][j] < rcThreshold){
                        habSuitability[i][j] = 0;
                    } else habSuitability[i][j] = 1000;
                }
                }
            }
            /* Filter the habitat suitability matrix in three ways:
            **  -> replace any value < 0 by 0 (this removes NoData).
            **  -> set habitat suitability to 0 where barrier = 1.
            **  -> set habitat suitability values to NoData where barrier = NoData. */
            mcFilterMatrix(habSuitability, barriers, true, true, true);
            
            
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            }
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            
            
            
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            /* Load dispersal kernel values for the current environmental change step.
            ** **********************************************************************
            ** If we are running the simulation in multi-kernel mode, i.e. with a custom 
            ** dispersal kernel for each cell, load the kernels for the current environmental
            ** change step. */
            if(singleKernelMode == false){
                /* Read all dispersal kernel matrices into the dispKernel float vector. */
                
                // In fact we only have to relaod the kernels for loops 1 - 4.
                if(envChgStep <= 4){
                    for(i = 0; i < dispDist; i++){
                        if(envChgStep == 1) sprintf(fileName, "%s_2020_%d.asc", dipersalKernelBasename, i+1);
                        if(envChgStep == 2) sprintf(fileName, "%s_2030_%d.asc", dipersalKernelBasename, i+1);
                        if(envChgStep == 3) sprintf(fileName, "%s_2040_%d.asc", dipersalKernelBasename, i+1);
                        if(envChgStep == 4) sprintf(fileName, "%s_2050_%d.asc", dipersalKernelBasename, i+1);
                        if(envChgStep == 5) sprintf(fileName, "%s_2060_%d.asc", dipersalKernelBasename, i+1);
                        if(envChgStep == 6) sprintf(fileName, "%s_2070_%d.asc", dipersalKernelBasename, i+1);
                        if(envChgStep == 7) sprintf(fileName, "%s_2080_%d.asc", dipersalKernelBasename, i+1);
                        if(envChgStep == 8) sprintf(fileName, "%s_2090_%d.asc", dipersalKernelBasename, i+1);
                        if(envChgStep == 9) sprintf(fileName, "%s_2100_%d.asc", dipersalKernelBasename, i+1);
                        
                        /* Important: the pointer we pass to 'readMatrixFloat' is not a pointer to the first
                        ** element of the array, but a pointer to the element 'i * nrRows * nrCols'. This 
                        ** allows to sequentially fill the 'dispKernel' array with the values contained in each 
                        ** of the dispersal kernel rasters. Basically we are simply putting all values in a 
                        ** long vector, and when we want to access the values later, we jump a number of cells 
                        ** equal to the number of cells in a kernel raster. 
                        ** When passing the address of a given cell of an array, C will keep filling the 
                        ** array starting from this position - it's as if this position is the first position
                        ** in the array. */
                        if( readMatrixFloat(fileName, &dispKernel[i * nrRows * nrCols]) != 0 ){
                            Rprintf("ERROR: failed to read dispersal kernel [%s].\n", fileName);
                        }
                    }
                }
            }
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            
            
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            /* Load simple kernel for the current envCHgStep. */
            if(singleKernelMode == true){
                if(envChgStep <= 4){
                    sprintf(fileName, "simple_kernel_%s%d.txt", hsMap, envChgStep);
                    
                    /* Open the input parameter file for reading. */
                    if((fp3 = fopen(fileName, "r")) == NULL){
                        Rprintf ("Can't open parameter file %s\n", fileName);
                        goto END_OF_FUNCTION;
                    }
                    
                    for(i = 0; i < dispDist; i++){
                        if(fscanf(fp3, "%s", &tmpString[0]) != 1){
                            Rprintf("Invalid dispersal kernel values on line %d in parameter file %s.\n", i+1, fileName);
                            goto END_OF_FUNCTION;
                        }
                        dispKernel[i] = strtof(&tmpString[0], NULL);
                        /* Read the carriage return character so that we move the pointer in the file to the 
                        ** start of the next line. */
                        tmpInt = fscanf(fp3, "\n");
                    }
                    fclose(fp3);
                }
            }
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            
            
            
            
            /* Set the values that will keep track of pixels colonized during the next climate 
            ** change loop. 'loopID' is the value that will be given to the pixel colonized during
            ** the current loop. The pixels being decolonized during the current loop will be 
            ** given a value of '-loopID'. The limits in terms of number of dispersal events and 
            ** environmental change events is the following:
            **  - Without environmental change, the maximum number of dispStep is 29'500.
            **  - With environmental change, the maximum number of dispStep is 98, the number 
            **    of environmental change loops is limited to 250.
            **
            ** The coding for the loopID is as follows: envChgStep * 100 + dispStep. */
            if(envChgStep == 0) loopID = 1; else loopID = envChgStep * 100;
            
            
            /* 'Unlimited' and 'no dispersal' scenario pixel count.
            ** ***************************************************
            ** Compute the number of cells that would be colonized if dispersal was unlimited or 
            ** null. This is simply the sum of all potentially suitable habitats */
            nrUnivDispersal = mcUnivDispCnt(habSuitability);
            updateNoDispMat(habSuitability, noDispersal, &nrNoDispersal);

            /* Reset number of decolonized cells within current dispersal step pixel counter */
            nrStepDecolonized = 0;
            
            
            /* Update for temporarily resilient pixels.
            ** ***************************************/
            for(i = 0; i < nrRows; i++){
            for(j = 0; j < nrCols; j++){
                
                /* Update unsuitable pixels. If a pixel turned unsuitable during the last 
                ** environmental change step, we update its status to 'Temporarily Resilient'. */
                if((habSuitability[i][j] == 0) && (currentState[i][j] > 0)){
                    /* If the user enabled the use of temporary resilience status, the cell's 
                    ** status is set to 'Temporary Resilient'. Otherwise the cell's status is 
                    ** set to 'decolonized'. */
                    if(tempResilience == true){
                        currentState[i][j] = 29900;
                    }
                    else{
                        currentState[i][j] = -1 - loopID;
                        pixelAge[i][j] = 0;
                    }
                    
                    /* Increase counter of decolonized cells within the current step. */
                    nrStepDecolonized++;
                }
            }
            }

            
            
            /* ***********************************************************************************
            ** Dispersal event loop starts here.
            ** **********************************************************************************/
            for(dispStep = 1; dispStep <= dispSteps; dispStep++){
                
                /* Set the value of 'loopID' for the current iteration of the dispersal loop. */
                loopID++;
                    
                /* Reset pixel counters that count pixels within the current loop. */
                nrStepColonized = 0;
                nrStepLDDSuccess = 0;
                if(dispStep > 1) nrStepDecolonized = 0;

                
                /* Source cell search:
                ** ******************
                ** Can the sink cell be colonized? There are four conditions to be met for a sink 
                ** cell to become colonized:
                **   1. Sink cell is currently suitable and not already colonized.
                **   2. Sink cell is within dispersal distance of an already colonized
                **      and mature cell.
                **   3. Source pixel has reached dispersal maturity.
                **   4. No obstacle (barrier) between sink and source cells.
                **
                ** Loop through the cell matrix. */
                for(i = 0; i < nrRows; i++){
                for(j = 0; j < nrCols; j++){
                    
                    /* Reset variables. */
                    habIsSuitable = false;
                    cellInDispDist = false;
                    
                    /* Test whether the pixel is a suitable sink cell:
                    **  -> its habitat is suitable.
                    **  -> it is unoccupied.
                    **  -> it is not on a barrier or filter pixel). */
                    if((habSuitability[i][j] > 0) && (currentState[i][j] <= 0)) habIsSuitable = true;

                    /* Search for a source cell within dispersal distance:
                    ** To be more time efficient, this code runs only if the cell's habitat is 
                    ** suitable. If there is a source cell within dispersal distance and if a 
                    ** barrier was provided, then we also check that there is no barrier between 
                    ** the source and sink cell. */
                    if(habIsSuitable){
                        if(mcSrcCell(i, j, currentState, pixelAge, loopID, 
                                     habSuitability[i][j], barriers)) cellInDispDist = true;
                    }
                    
                    /* Update pixel status. 
                    ** Only if both conditions are met is the cell's status set to colonized. */
                    if(habIsSuitable && cellInDispDist){
                        currentState[i][j] = loopID;
                        nrStepColonized++;
                        /* Update the cell's 'age' value. We do this only now because we needed 
                        ** the old 'age' value just before to determine whether a pixel was in 
                        ** 'Decolonized' or 'SeedBank resilience' status. */
                        pixelAge[i][j] = 0;
                    }
                }
                }
                
                /* Long distance dispersal (LDD).
                ** *****************************
                ** If LDD frequency is larger than zero, run LDD events for all suitable cells. */
                if(lddFreq > 0.0){
                    for(i = 0; i < nrRows; i++){
                    for(j = 0; j < nrCols; j++){
                        
                        /* Check whether the pixel is a potential source cell: 
                        **  -> it is colonized since at least 1 dispersal Loop.
                        **  -> it has reached dispersal maturity. */
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
                        
                            /* Try to generate a LDD event with the calculated probability. */
                            if(UNIF01 < lddSeedProb || lddSeedProb == 1.0){
                                
                                /* Randomly select a pixel within distance lddMinDist to 
                                ** lddMaxDist. */
                                mcRandomPixel(&rndPixel);
                                rndPixel.row = rndPixel.row + i;
                                rndPixel.col = rndPixel.col + j;
                                
                                /* Now we check if this random cell is a suitable sink cell.
                                ** If yes, then the cell gets colonized. */
                                if(mcSinkCellCheck (rndPixel, currentState, habSuitability)){
                                    /* Set cell status to colonized */
                                    currentState[rndPixel.row][rndPixel.col] = loopID;
                                    nrStepColonized++;
                                    nrStepLDDSuccess++;
                                    /* Reset pixel age. */
                                    pixelAge[rndPixel.row][rndPixel.col] = 0;
                                }
                            }
                            
                        }
                        }
                    }
                    }
                }
                
                
                /* Update cell age value (i.e. time since colonized).
                ** ************************************************* 
                ** At the end of a dispersal loop we increase the 'age' of each colonized pixel.
                ** Reminder: cell 'age' structure is as follows:
                **   0     = cell is either 'Absent', 'Decolonized' or has just been 'Colonized' 
                **           during current dispersal step.
                **   1-250 = cell is in 'Colonized' or 'Temporarily Resilient' status. The value
                **            indicates the number of 'dispersal events' (usually years) since
                **            the cell was colonized.
                **   255   = cell is in 'SeedBank Resilience' state. */
                for(i = 0; i < nrRows; i++){
                for(j = 0; j < nrCols; j++){
                    /* If the cell is in 'Colonized' or 'Temporarily Resilient' status, update
                    ** it's age value. */
                    if(currentState[i][j] > 0) pixelAge[i][j] += 1;

                    /* If a cell is in 'Temporarily Resilient' state, we also increase its 
                    ** 'currentState' value by 1 so that it gains 1 year of temporarily resilience 
                    ** age. */
                    if (currentState[i][j] >= 29900) currentState[i][j] += 1;
                }
                }
                
                
                /* Update cell counters. 
                ** ********************/
                nrColonized = nrColonized + nrStepColonized - nrStepDecolonized;
                nrAbsent = nrAbsent - nrStepColonized + nrStepDecolonized;
                nrTotColonized += nrStepColonized;
                nrTotDecolonized += nrStepDecolonized;
                nrTotLDDSuccess += nrStepLDDSuccess;

                
                /* Write current iteration data to the statistics file.
                ** ***************************************************/
                fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", envChgStep, dispStep, 
                        loopID, nrUnivDispersal, nrNoDispersal, nrColonized, nrAbsent, 
                        nrStepColonized, nrStepDecolonized, nrStepLDDSuccess);
                    
                /* If the user has requested full output, also write the current state matrix 
                ** to file. */
                if(fullOutput){
                    sprintf(fileName, "%s/%s_step_%d.asc", simulName, simulName2, loopID);
                    if(writeMat(fileName, currentState) == -1){
                        *nrFiles = -1;
                        goto END_OF_FUNCTION;
                    }
                }
            }
            /* ***********************************************************************************
            ** End of dispersal step loop. 
            ** **********************************************************************************/

            
            /* Update temporarily resilient pixels.
            ** ***********************************
            ** Temporarily resilient pixels can be distinguished by:
            **  -> CurrentState_Matrix = 29'900 to 29'999. Increases by 1 at each year.
            **  -> Age_Matrix has a positive value. */
            for (i = 0; i < nrRows; i++){
            for (j = 0; j < nrCols; j++){
                if (currentState[i][j] >= 29900){
                    currentState[i][j] = dispSteps - loopID - 1;
                    pixelAge[i][j] = 0;
                }
            }
            }
            
        }
        /* ***************************************************************************************
        ** End of environmental change step loop. 
        ** **************************************************************************************/
        
        
        
        /* Progress report in R console. */
        Rprintf("All dispersal steps completed. Final output in progress...\n");

        
        /* Update current status matrix.
        ** ****************************
        ** Update currentState matrix for pixels that are suitable but could not be colonized due 
        ** to dispersal limitations. These pixels are assigned a value of 30'000 */
        for(i = 0; i < nrRows; i++){
        for(j = 0; j < nrCols; j++){
            if((habSuitability[i][j] > 0) && (currentState[i][j] <= 0)) currentState[i][j] = 30000;
        }
        }

        /* Write the final state matrix to file. */
        sprintf(fileName, "%s/%s_raster.asc", simulName, simulName2);
        if(writeMat (fileName, currentState) == -1){
            *nrFiles = -1;
            goto END_OF_FUNCTION;
        }

        
        /* Write summary output to file. 
        ** ****************************/
        simulTime = time (NULL) - startTime;
        sprintf(fileName, "%s/%s_summary.txt", simulName, simulName2);
        if((fp2 = fopen (fileName, "w")) != NULL){
            fprintf(fp2, "simulName\tiniCount\tnoDispCount\tunivDispCount\toccupiedCount\tabsentCount\ttotColonized\ttotDecolonized\ttotLDDsuccess\trunTime\n");
            fprintf(fp2, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", simulName2, nrInitial, 
                    nrNoDispersal, nrUnivDispersal, nrColonized, nrAbsent, nrTotColonized, 
                    nrTotDecolonized, nrTotLDDSuccess, simulTime);
            fclose(fp2);
        }
        else{
            *nrFiles = -1;
            Rprintf ("Could not write summary output to file.\n");
            goto END_OF_FUNCTION;
        }

        /* Close the data file. */
        if(fp != NULL) fclose(fp);
        fp = NULL;

    } 
    /* *******************************************************************************************
    ** End of simulation replication loop. 
    ** ******************************************************************************************/
        

    /* Set the number of output files created. */
    *nrFiles = envChgSteps;



 END_OF_FUNCTION:
    /* Close open data file if needed. */
    if(fp != NULL) fclose(fp);

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
    if(dispKernel != NULL) free(dispKernel);
    if(propaguleProd != NULL) free(propaguleProd);

    /* If an error occurred, display failure message to the user... */
    if(*nrFiles == -1) Rprintf("MigClim simulation aborted...\n");
    
}



void mcRandomPixel (pixel *pix){
    /* *******************************************************************************************
    ** Select a random pixel from a central point (0;0), within a radius of at least lddMinDist 
    ** and at most lddMaxDist. Note that 'lddMinDist' and 'lddMaxDist' are defined as global 
    ** variables for the MigClim code and are therefore already available without needing to pass
    ** them to the function. It is OK to have these two values as global variables because their 
    ** value is never modified (they are fixed for a given run of MigClim).
    ** 
    ** Parameters:
    **  -> pix: a pointer to the pixel data structure to put the row and col in.
    ** Returns:
    **  ->
    *********************************************************************************************/

    /* Select a random distance between lddMinDist and lddMaxDist, and a random angle 
    ** between 0 and 2*pi.*/
    double rndDist, rndAngle;
    rndDist = (UNIF01 * (lddMaxDist - lddMinDist)) + lddMinDist;
    rndAngle = UNIF01 * 6.283185;

    // Convert distance and angle into pixel row and column values.
    pix->row = (int)(rndDist * cos(rndAngle));
    pix->col = (int)(rndDist * sin(rndAngle));

}



bool mcSinkCellCheck (pixel pix, int **curState, int **habSuit){
    /* *******************************************************************************************
    ** Perform a basic check to see whether a given cell fulfills the conditions to be a 'sink' 
    ** cell, i.e. a cell that can be potentially colonized. The conditions that are checked are
    ** the following:
    **  1. Cell must be within the limits of the cellular automaton.
    **  2. Cell must be empty, i.e. non-occupied.
    **  3. Cell must contain suitable habitat.
    **
    ** Parameters:
    **  -> pix     : the cell/pixel to consider.
    **  -> curState: a pointer to the current state matrix.
    **  -> habSuit : a pointer to the habitat suitability matrix.
    **
    ** Returns:
    **  -> true : if cell is suitable.
    **  -> false: otherwise
    ** ******************************************************************************************/
    double rnd;

    /* Verify the cell is within the limits of the cellular automaton. */
    if((pix.row < 0) || (pix.row >= nrRows) || (pix.col < 0) || (pix.col >= nrCols)) return(false);

    /* Verify the cell is empty. */
    if(curState[pix.row][pix.col] > 0) return(false);

    /* Verify the cell contains suitable habitat for the species:
    **  -> habitat suitability is > 0.
    **  -> habitat suitability is > random number. Note that UNIF01 is declared in the header 
    **     section and the "(double)" is typecasting the value to a double.
     */
    if(habSuit[pix.row][pix.col] == 0) return(false);
    rnd = UNIF01 * 1000;
    if(rnd > (double)habSuit[pix.row][pix.col]) return(false);
    
    /* If the function has not exited by now then it means the cell is suitable. */
    return (true);
}

