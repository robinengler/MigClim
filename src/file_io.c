/*
####################################################################################################
## file_io.c: functions for file input/output.
####################################################################################################
*/
#include "migclim.h"


int mcInit (char *paramFile){
    /* *******************************************************************************************
    ** Initialize the MigClim model by reading the parameter values from file.
    **
    ** Parameters:
    **  paramFile: The name of the file from which to read the parameter values.
    **
    ** Return value:
    **  '0' if everything went fine, '-1' otherwise.
    ** ******************************************************************************************/
    
    int          i, exitStatus, lineNr, tmpInt;
    unsigned int randomGeneratorSeed;
    long         tmpLong;
    double       tmpDouble;
    char         line[1024], parameterName[256], fileName[256], tmpString[128];
    char         *end;
    FILE         *fp;

    /* Set default value for 'exitStatus' to 1 (error occured) and fp to NULL. When fp is NULL 
    ** we know that the file is closed, when it's non-NULL then the file is open. */
    exitStatus = 1;
    fp = NULL;

    /* Set default parameter values. */
    nrRows = 0;
    nrCols = 0;
    strcpy(iniDist, "");
    strcpy(hsMap, "");
    strcpy(barrier, "");
    useBarrier = false;
    barrierType = STRONG_BARRIER;
    envChgSteps = 0;
    dispSteps = 0;
    dispDist = 0;
    iniMatAge = 0;
    fullMatAge = 0;
    rcThreshold = 0;
    lddMinDist = 0;
    lddMaxDist = 0;
    lddFreq = 0.0;
    fullOutput = false;
    replicateNb = 1;
    tmpInt = 0;
    randomGeneratorSeed = 0;
    //strcpy(simulName, "MigClimTest");


    /* Open the input parameter file for reading. */
    if((fp = fopen(paramFile, "r")) == NULL){
        Rprintf ("Can't open parameter file %s\n", paramFile);
        goto END_OF_FUNCTION;
    }


    /* Read the input parameter file line by line. Here is an example of input parameter file: 
    **   nrRows 408
    **   nrCols 321
    **   iniDist InitialDist
    **   hsMap HSmap
    **   rcThreshold 500
    **   envChgSteps 5
    **   dispSteps 5
    **   dispDist 5
    **   dispKernel 1 0.5 0.25 0.125 0.0625
    **   barrier Barrier
    **   barrierType strong
    **   iniMatAge 1
    **   fullMatAge 2
    **   propaguleProd 1
    **   fullOutput false
    **   replicateNb 1
    **   simulName MigClimTest_1  
    **
    ** Note: the line 'dispKernel' can be replaced by 'dispKernelFile'. */
    
    
    /* Size of migclim input matrices (number of rows and cols):
    ** ********************************************************
    ** All input matrices in migclim have the same number of rows and columns.
    ** nrRows: number of rows. */
    lineNr = 1;
    if((fgets(line, 1024, fp) == NULL) || 
       (sscanf(line,"%s %d\n", parameterName, &nrRows) != 2) || 
       (strcasecmp(parameterName, "nrRows") != 0)){
        Rprintf("ERROR: line %d of file [%s] must contain value for 'nrRows'.\n", 
                lineNr, paramFile);
        goto END_OF_FUNCTION;
    }
    /* 'nrCols': the number of columns. */
    lineNr++;
    if((fgets(line, 1024, fp) == NULL) || 
       (sscanf(line,"%s %d\n", parameterName, &nrCols) != 2) || 
       (strcasecmp(parameterName, "nrCols") != 0)){
        Rprintf("ERROR: line %d of file [%s] must contain value for 'nrCols'.\n", 
                lineNr, paramFile);
        goto END_OF_FUNCTION;
    }

    
    /* Habitat suitability rasters and initial distribution.
    ** **************************************************** 
    ** 'iniDist': string giving the name of the species initial distribution raster file. */
    lineNr++;
    if((fgets(line, 1024, fp) == NULL) || 
       (sscanf(line,"%s %s\n", parameterName, iniDist) != 2) || 
       (strcasecmp(parameterName, "iniDist") != 0)){
        Rprintf("ERROR: line %d of file [%s] must contain value for 'iniDist'.\n", 
                lineNr, paramFile);
        goto END_OF_FUNCTION;
    }
    /* hsMap: a string giving the basename of the habitat suitability raster files. */
    lineNr++;
    if((fgets(line, 1024, fp) == NULL) || 
       (sscanf(line,"%s %s\n", parameterName, hsMap) != 2) || 
       (strcasecmp(parameterName, "hsMap") != 0)){
        Rprintf("ERROR: line %d of file [%s] must contain value for 'hsMap'.\n", 
                lineNr, paramFile);
        goto END_OF_FUNCTION;
    }
    /* rcThreshold: reclassification threshold for values in habitat suitability maps.*/
    lineNr++;
    if((fgets(line, 1024, fp) == NULL) || 
       (sscanf(line,"%s %d\n", parameterName, &rcThreshold) != 2) || 
       (strcasecmp(parameterName, "rcThreshold") != 0)){
        Rprintf("ERROR: line %d of file [%s] must contain value for 'rcThreshold'.\n", 
                lineNr, paramFile);
        Rprintf("ERROR: 'rcThreshold' must be a value in the range [0,1000].\n", lineNr, paramFile);
        goto END_OF_FUNCTION;
    }

    
    /* Environmental change and dispersal steps:
    ** ****************************************
    ** envChgSteps: number of times the habitat suitability maps should be updated. */
    lineNr++;
    if((fgets(line, 1024, fp) == NULL) || 
       (sscanf(line,"%s %d\n", parameterName, &envChgSteps) != 2) || 
       (strcasecmp(parameterName, "envChgSteps") != 0)){
        Rprintf("ERROR: line %d of file [%s] must contain value for 'envChgSteps'.\n", 
                lineNr, paramFile);
        goto END_OF_FUNCTION;
    }
    /* dispSteps: number dispersal events within each environmental change step. */
    lineNr++;
    if((fgets(line, 1024, fp) == NULL) || 
       (sscanf(line,"%s %d\n", parameterName, &dispSteps) != 2) || 
       (strcasecmp(parameterName, "dispSteps") != 0)){
        Rprintf("ERROR: line %d of file [%s] must contain value for 'dispSteps'.\n", 
                lineNr, paramFile);
        goto END_OF_FUNCTION;
    }


    /* Dispersal distance:
    ** ******************
    ** dispDist: maximum dispersal distance (in cell units) for regular dispersal events.*/
    lineNr++;
    if((fgets(line, 1024, fp) == NULL) || 
       (sscanf(line,"%s %d\n", parameterName, &dispDist) != 2) || 
       (strcasecmp(parameterName, "dispDist") != 0)){
        Rprintf("ERROR: line %d of file [%s] must contain value for 'dispDist'.\n", 
                lineNr, paramFile);
        goto END_OF_FUNCTION;
    }


    /* Dispersal kernel (dispKernel):
    ** ***************************** 
    ** Either a list of float values giving the probability of dispersal PDisp associated to each
    ** distance class, or a string giving the basename of the dispersal kernel ascii raster files.
    ** Dispersal kernel raster are files giving the values of PDisp for each pixel individually,
    ** and each raster gives values for a given distance class. */
    lineNr++;
    if(fscanf(fp, "%s", parameterName) != 1){
        Rprintf("ERROR: unable to read line %d of file [%s].\n", lineNr, paramFile);
        goto END_OF_FUNCTION;
    }
    
    /* Case 1: single kernel mode, the user passed a list of float values as kernel, the kernel 
    ** is the same for all cells in the landscape (this is the standard mode). */
    if(strcmp(parameterName, "dispKernel") == 0){
        singleKernelMode = true;
        
        /* Allocate memory so that dispKernel is a vector of size dispDist */
        dispKernel = (float *)malloc(dispDist * sizeof(float));
        for(i = 0; i < dispDist; i++){
            if(fscanf(fp, " %s", &tmpString[0]) != 1){
                Rprintf("Invalid dispersal kernel values on line %d in parameter file %s.\n", 
                        lineNr, paramFile);
                goto END_OF_FUNCTION;
            }
            /* use the 'strtof' (string-to-float) function to convert string to float number. */
            dispKernel[i] = strtof(&tmpString[0], NULL);
        }
        /* Read the carriage return character so that we move the pointer in the file to the 
        ** start of the next line. */
        tmpInt = fscanf(fp, "\n");
    }
    
    /* Case 2: multi-kernel mode. The user passed the basename of the dispersal kernel ascii raster 
    **         files. Read the ascii rasters and store their values in the 'dispKernel' vector. */
    else if( strcmp(parameterName, "dispKernelFile") == 0 ){
        singleKernelMode = false;
        
        /* Allocate memory so that dispKernel is a vector that can contain all values from all 
        ** dispersal kernel rasters. Values are put in the vector one after the other until we 
        ** completed looping through all maps). */
        dispKernel = (float *)malloc(dispDist * nrRows * nrCols * sizeof(float));
        
        /* Read basename of kernel rasters. */
        if (fscanf(fp, " %s\n", dipersalKernelBasename) != 1){
            Rprintf("Invalid dispersal kernel values on line %d in parameter file %s.\n", 
                    lineNr, paramFile);
            goto END_OF_FUNCTION;
        }
    
        /* Read all dispersal kernel matrices into the dispKernel float vector. */
        for(i = 0; i < dispDist; i++){
            sprintf(fileName, "%s%d.asc", dipersalKernelBasename, i+1);
            
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            // TEMPORARY CHANGE IN FILE NAME PATTERN FOR BRYOPHYTE PROJECT.
            sprintf(fileName, "%s_2020_%d.asc", dipersalKernelBasename, i+1);
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            // TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY  --  TEMPORARY
            
            
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
    else{
        Rprintf("ERROR: line %d of parameter file must be 'dispDist' or 'dispKernelFile'.\n", 
                lineNr);
    }


    /* barrier and barrier type (optional values): 
    ** ******************************************
    ** barrier: a string giving the name of the barrier raster file. */
    lineNr++;
    if((fgets(line, 1024, fp) == NULL) || 
       (sscanf(line,"%s %s\n", parameterName, barrier) != 2) || 
       (strcasecmp(parameterName, "barrier") != 0)){
        Rprintf("ERROR: line %d of file [%s] must contain the value for 'barrier'.\n", 
                lineNr, paramFile);
        goto END_OF_FUNCTION;
    }
    if(strcmp(barrier, "NA") == 0) useBarrier = false; else useBarrier = true; 

    /* barrierType: a string giving the name of the barrier ascii raster files. */
    lineNr++;
    if((fgets(line, 1024, fp) == NULL) || 
       (sscanf(line,"%s %s\n", parameterName, tmpString) != 2) || 
       (strcasecmp(parameterName, "barrierType") != 0) || 
       (strlen(tmpString) == 0)){
        Rprintf("ERROR: line %d of file [%s] must contain value for 'barrierType'.\n", 
                lineNr, paramFile);
        goto END_OF_FUNCTION;
    }
    if(useBarrier==true){
        /* If a barrier raster layer was specified, then barrierType must be either 'weak' or 
        ** 'strong'. WEAK_BARRIER and STRONG_BARRIER are simply macros that expand to respectively 
        ** 1 and 2, so the type of barrierType is int and not a char. */
        if(strcmp(tmpString, "weak") == 0){
            barrierType = WEAK_BARRIER;
        }
        else if(strcmp(tmpString, "strong") == 0){
            barrierType = STRONG_BARRIER;
        }
        else{
            Rprintf ("Invalid barrier type on line %d in parameter file %s\n", lineNr, paramFile);
            goto END_OF_FUNCTION;
        }
    }

    

    /* Species propagule production potential.
    ** **************************************
    ** iniMatAge: age at which cells start to produce propagules. Must be a value >= 1. */
    lineNr++;
    if((fgets(line, 1024, fp) == NULL) || 
       (sscanf(line,"%s %d\n", parameterName, &iniMatAge) != 2) || 
       (strcasecmp(parameterName, "iniMatAge") != 0) || 
       (iniMatAge < 1)){
        Rprintf("ERROR: line %d of file [%s] must contain value for 'iniMatAge'.\n", 
                lineNr, paramFile);
        Rprintf("ERROR: 'iniMatAge' must be an integer >= 1.\n", lineNr, paramFile);
        goto END_OF_FUNCTION;
    }

    /* fullMatAge: age at which cells reach full propagule production potential. */
    lineNr++;
    if((fgets(line, 1024, fp) == NULL) || 
       (sscanf(line,"%s %d\n", parameterName, &fullMatAge) != 2) || 
       (strcasecmp(parameterName, "fullMatAge") != 0) || 
       (fullMatAge < iniMatAge)){
        Rprintf("ERROR: line %d of file [%s] must contain value for 'fullMatAge'.\n", 
                lineNr, paramFile);
        Rprintf("ERROR: 'fullMatAge' must be an integer >= 1 and >= 'iniMatAge'.\n", 
                lineNr, paramFile);
        goto END_OF_FUNCTION;
    }

    /* propaguleProd: list of float values giving the probability of propagule production for 
    ** each age category from iniMatAge to fullMatAge. If iniMatAge == fullMatAge then 
    ** propaguleProd has a single value that is equal to 1 (100%). */
    lineNr++;
    if((fscanf(fp, "%s", parameterName) != 1) || (strcasecmp(parameterName, "propaguleProd") != 0)){
        Rprintf("ERROR: line %d of file [%s] must contain values for 'propaguleProd'.\n", 
                lineNr, paramFile);
        goto END_OF_FUNCTION;
    }
    /* Allocate memory so that propaguleProd is a vector of size age */
    if(fullMatAge > iniMatAge){
        propaguleProd = (double *)malloc((fullMatAge - iniMatAge) * sizeof(double));
        for(i = 0; i < (fullMatAge - iniMatAge); i++){
            if(fscanf(fp, " %lf", &tmpDouble) != 1){
                Rprintf("Invalid propaguleProd kernel values on line %d in parameter file [%s].\n", 
                        lineNr, paramFile);
                goto END_OF_FUNCTION;
            }
            propaguleProd[i] = tmpDouble;
        }
    }
    else{
        propaguleProd = NULL;
    }
    /* Read the carriage return character so that we move the pointer in the file to the start 
    ** of the next line. */
    tmpInt = fscanf(fp, "\n");


    
    /* Long Distance Dispersal (LDD).
    ** ***************************** 
    ** lddFreq: frequence of long distance dispersal events. */
    lineNr++;
    if((fgets(line, 1024, fp) == NULL) || 
       (sscanf(line,"%s %lf\n", parameterName, &lddFreq) != 2) || 
       (strcasecmp(parameterName, "lddFreq") != 0) ){
        Rprintf("ERROR: line %d of file [%s] must contain value for 'lddFreq'.\n", 
                lineNr, paramFile);
        Rprintf("ERROR: parameterName[%s], tmpString[%s].\n", parameterName, tmpString);
        goto END_OF_FUNCTION;
    }
    if( (lddFreq < 0.0) || (lddFreq > 1.0) ){
        Rprintf("ERROR: 'lddFreq' must be an float in the range [0,1].\n");
        goto END_OF_FUNCTION;
    }

    /* lddMinDist: minimum distance for long distance dispersal events. */
    lineNr++;
    if((fgets(line, 1024, fp) == NULL) || 
       (sscanf(line,"%s %d\n", parameterName, &lddMinDist) != 2) || 
       (strcasecmp(parameterName, "lddMinDist") != 0)){
        Rprintf("ERROR: line %d of file [%s] must contain value for 'lddMinDist'.\n", 
                lineNr, paramFile);
        goto END_OF_FUNCTION;
    }

    /* lddMaxDist: maximum distance for long distance dispersal events. */
    lineNr++;
    if((fgets(line, 1024, fp) == NULL) || 
       (sscanf(line,"%s %d\n", parameterName, &lddMaxDist) != 2) || 
       (strcasecmp(parameterName, "lddMaxDist") != 0)){
        Rprintf("ERROR: line %d of file [%s] must contain value for 'lddMaxDist'.\n", 
                lineNr, paramFile);
        goto END_OF_FUNCTION;
    }


    /* Other parameters.
    ** ****************
    ** fullOutput: indicates whether migclim should output more (true) or less (false) data. */
    lineNr++;
    if((fgets(line, 1024, fp) == NULL) || 
       (sscanf(line,"%s %s\n", parameterName, tmpString) != 2) || 
       (strcasecmp(parameterName, "fullOutput") != 0)){
        Rprintf("ERROR: line %d of file [%s] must contain value for 'fullOutput'.\n", 
                lineNr, paramFile);
        goto END_OF_FUNCTION;
    }
    if(strcmp(tmpString, "true") == 0) fullOutput = true;
    if(strcmp(tmpString, "false") == 0) fullOutput = false;


    /* replicateNb: number of times the migclim simulation should be replicated. This must be an integer >= 1. */
    lineNr++;
    if((fgets(line, 1024, fp) == NULL) || 
       (sscanf(line,"%s %d\n", parameterName, &replicateNb) != 2) || 
       (strcasecmp(parameterName, "replicateNb") != 0)){
        Rprintf("ERROR: line %d of file [%s] must contain value for 'replicateNb'.\n", 
                lineNr, paramFile);
        goto END_OF_FUNCTION;
    }

    /* simulName: string giving a name for the current migclim run. */
    lineNr++;
    if((fgets(line, 1024, fp) == NULL) || 
       (sscanf(line,"%s %s\n", parameterName, simulName) != 2) || 
       (strcasecmp(parameterName, "simulName") != 0)){
        Rprintf("ERROR: line %d of file [%s] must contain value for 'simulName'.\n", 
                lineNr, paramFile);
        goto END_OF_FUNCTION;
    }

    /* randomGeneratorSeed: either "NA" or an unsigned integer [0,65535] to be used to see the random generator
    **                      in MigClim. If "NA" was passed a value, then the random number generator is seeded
    **                      with the current time. */
    lineNr++;
    if((fgets(line, 1024, fp) == NULL) || 
       (sscanf(line,"%s %s\n", parameterName, &tmpString[0]) != 2) || 
       (strcasecmp(parameterName, "randomGeneratorSeed") != 0) ){
        Rprintf("ERROR: line %d of file [%s] must contain value for 'replicateNb'.\n", 
                lineNr, paramFile);
        goto END_OF_FUNCTION;
    }
    
    /* If 'NA' was passed as randomGeneratorSeed, then we set randomGeneratorSeed = -1, a code 
    ** used to indicate that the random number generator should be seeded with the current time. */
    if(strcasecmp(tmpString, "NA") == 0){
        randomGeneratorSeed = 0;
    }
    else{
        errno = 0;
        tmpLong = strtol(&tmpString[0], &end, 10);
        if(end == &tmpString[0] || 
           ((tmpLong == LONG_MAX || tmpLong == LONG_MIN) && errno == ERANGE)){
            Rprintf("ERROR: failed to convert input string to long integer.\n");
        } 
        if((tmpLong < 1) || (tmpLong > 65535)){
            Rprintf("ERROR: value for 'randomGeneratorSeed' must be either 'NA' or an integer in the range [1,65535].\n");
            goto END_OF_FUNCTION;
        }
        randomGeneratorSeed = (unsigned int)tmpLong;
    }

    
    
    /* Progress report to user.
    ** ***********************/
    if(VERBOSE){
        Rprintf("Input values summary:\n");
        Rprintf("nrRows = [%d]\n", nrRows);
        Rprintf("nrCols = [%d]\n", nrCols);
        Rprintf("iniDist = [%s]\n", iniDist);
        Rprintf("hsMap = [%s]\n", hsMap);
        Rprintf("rcThreshold = [%d]\n", rcThreshold);
        Rprintf("envChgSteps = [%d]\n", envChgSteps);
        Rprintf("dispSteps = [%d]\n", dispSteps);
        Rprintf("dispDist = [%d]\n", dispDist);
        //Rprintf("dispKernel = [%].\n", );
        Rprintf("barrier = [%s]\n", barrier);
        if(barrierType == WEAK_BARRIER) Rprintf("barrierType = [weak]\n");
        if(barrierType == STRONG_BARRIER) Rprintf("barrierType = [strong]\n");
        Rprintf("iniMatAge = [%d]\n", iniMatAge);
        Rprintf("fullMatAge = [%d]\n", fullMatAge);
        //Rprintf("propaguleProd = [%lf].\n", );
        Rprintf("lddFreq = [%lf]\n", lddFreq);
        Rprintf("lddMinDist = [%d]\n", lddMinDist);
        Rprintf("lddMaxDist = [%d]\n", lddMaxDist);
        //Rprintf("fullOutput = [%s].\n", fullOutput);
        Rprintf("replicateNb = [%d]\n", replicateNb);
        Rprintf("simulName = [%s]\n", simulName);
        if(randomGeneratorSeed == -1) Rprintf("randomGeneratorSeed = [NA]\n");
        if(randomGeneratorSeed != -1) Rprintf("randomGeneratorSeed = [%d]\n", randomGeneratorSeed);
    } 

    /* If we reach this point, no error has occurred. Set exit status value to 0. */
    exitStatus = 0;

END_OF_FUNCTION:
    /* Close the file if it is open and return the exit status. */
    if(fp != NULL) fclose(fp);
    return (exitStatus);
}




/*###################################################################################################################### 
### function readMat().
### ******************
### Read a data matrix from an ascii grid file.
### Note: This should eventually be merged with the above "mcReadMatrix" function, but we'll keep it separate for now 
### just to make sure the basic functionality works fine.
###
### Parameters:
###  -> fName: The name of the file to read from.
###  -> mat:   The matrix to put the data in (assumed to be large enough).
###
### Returns: 0 if everything went fine, -1 otherwise.
######################################################################################################################*/
int readMat (char *fName, int **mat)
{
  int   i, j, intVal, status;
  char  line[1024], parameterName[128], dblVal[128];
  FILE *fp;

  status = 0;
  fp = NULL;
  
  /* Open the file for reading. */
  if((fp = fopen(fName, "r")) == NULL){
    status = -1;
    Rprintf ("Can't open data file %s\n", fName);
    goto END_OF_FUNCTION;
  }

  /* Get the meta data. */
  if((fgets (line, 1024, fp) == NULL) || 
      (sscanf (line, "%s %d\n", parameterName, &intVal) != 2) || 
      (strcasecmp (parameterName, "ncols") != 0)){
    status = -1;
    Rprintf ("'ncols' expected in data file %s.\n", fName);
    goto END_OF_FUNCTION;
  }
  
  if(intVal != nrCols){
    status = -1;
    Rprintf ("Invalid number of columns in data file %s\n", fName);
    goto END_OF_FUNCTION;
  }
  
  if((fgets (line, 1024, fp) == NULL) ||
      (sscanf (line, "%s %d\n", parameterName, &intVal) != 2) ||
      (strcasecmp (parameterName, "nrows") != 0)){
    status = -1;
    Rprintf ("'nrows' expected in data file %s\n", fName);
    goto END_OF_FUNCTION;
  }

  if (intVal != nrRows){
    status = -1;
    Rprintf ("Invalid number of rows in data file %s.\n", fName);
    goto END_OF_FUNCTION;
  }

  if((fgets (line, 1024, fp) == NULL) ||
      (sscanf (line, "%s %s\n", parameterName, dblVal) != 2) ||
      (strcasecmp (parameterName, "xllcorner") != 0)){
    status = -1;
    Rprintf ("'xllcorner' expected in data file %s\n", fName);
    goto END_OF_FUNCTION;
  }
  xllCorner = strtod(dblVal, NULL);
  
  if((fgets (line, 1024, fp) == NULL) ||
      (sscanf (line, "%s %s\n", parameterName, dblVal) != 2) ||
      (strcasecmp (parameterName, "yllcorner") != 0)){
    status = -1;
    Rprintf ("'yllcorner' expected in data file %s\n", fName);
    goto END_OF_FUNCTION;
  }
  yllCorner = strtod(dblVal, NULL);
  
  if((fgets (line, 1024, fp) == NULL) ||
      (sscanf (line, "%s %s\n", parameterName, dblVal) != 2) ||
      (strcasecmp (parameterName, "cellsize") != 0)){
    status = -1;
    Rprintf ("'cellsize' expected in data file %s\n", fName);
    goto END_OF_FUNCTION;
  }
  cellSize = strtod(dblVal, NULL);
  
  if((fgets (line, 1024, fp) == NULL) ||
      (sscanf (line, "%s %d\n", parameterName, &noData) != 2) ||
      (strcasecmp (parameterName, "nodata_value") != 0)){
    status = -1;
    Rprintf ("'NODATA_value' expected in data file %s\n", fName);
    goto END_OF_FUNCTION;
  }
  
  /* Read the values into the matrix. */
  for(i = 0; i < nrRows; i++){
    for(j = 0; j < nrCols; j++){

      if(fscanf(fp, "%d", &intVal) != 1){
        status = -1;
        Rprintf ("Invalid value in data file %s\n", fName);
        goto END_OF_FUNCTION;
      }
      mat[i][j] = intVal;
    }
    intVal = fscanf (fp, "\n");
  }



END_OF_FUNCTION:

  /* Close the file and return the function's return status. */
  if (fp != NULL) fclose(fp);
  return (status);
  
}
/*####################################################################################################################*/



/*###################################################################################################################### 
### function readMatrixFloat().
### **************************
### Reads float precision data from an ascii grid file (fName)
 Reads a data matrix from an ascii grid file.
### Note: This should eventually be merged with the above "mcReadMatrix" function, but we'll keep it separate for now 
### just to make sure the basic functionality works fine.
###
### Parameters:
###  -> fName:  A pointer to the name of the file to read from.
###  -> matrix: a pointer to a vector of float values. The vector must be already declared and is assumed to be 
###             large enough to store all values. More specifically, this is a pointer to the location in memory
###             from where we must start writing values into the array.
###
### Returns: 0 if everything went fine, -1 if an error occured.
######################################################################################################################*/
int readMatrixFloat(char *fName, float *matrixToFill)
{
  int   i, j, intValue, exitStatus;
  float floatValue;
  char  line[1024], parameterName[128], dblVal[128];
  FILE *fp;

  /* Set default value for 'exitStatus' to 1 (error occured) and fp to NULL. */
  /* When fp is NULL we know that the file is closed, when it's non-NULL then the file is open. */
  exitStatus = 1;
  fp = NULL;
  

  /* Open input ascii grid file.
  ** **************************
  ** Make sure the file was opened successfully by checking that 'fp' is not NULL. */
  fp = fopen(fName, "r");
  if(fp == NULL){
    Rprintf("Can't open data file %s\n", fName);
    goto END_OF_FUNCTION;
  }


  /* Read ascii grid file metadata.
  ** *****************************
  ** Read the ascii grid's metedata located on top of the file. The metadata has the following structure:
  **     NCOLS 321 
  **     NROWS 408 
  **     XLLCORNER 553437.5 
  **     YLLCORNER 115087.5 
  **     CELLSIZE 100 
  **     NODATA_value -9999  
  **
  ** Note: fgets() returns NULL on error. sscanf() returns the nubmer of elements that were read and assiged
  ** to pointer values specified as input. */

  /* Read the first line of the file and check that it has the structure "NCOLS XXX". Make sure that the number
  ** of columns of the matrix is correct, as all input matrices in MigClim must have the same size.*/
  if((fgets(line,1024,fp) == NULL) || (sscanf(line,"%s %d\n",parameterName,&intValue) != 2) || (strcasecmp(parameterName,"ncols") != 0)){
    Rprintf("'ncols' line expected in header of the input ascii grid file %s.\n", fName);
    goto END_OF_FUNCTION;
  }
  if(intValue != nrCols){
    Rprintf("Invalid number of columns in data file %s\n", fName);
    goto END_OF_FUNCTION;
  }
  
  /* Read the second line of the file. Check it has the structure "NROWS XXX". Check the number of rows is correct */
  if((fgets(line, 1024, fp) == NULL) || (sscanf(line, "%s %d\n", parameterName, &intValue) != 2) || (strcasecmp(parameterName, "nrows") != 0)){
    Rprintf("'nrows' expected in data file %s\n", fName);
    goto END_OF_FUNCTION;
  }
  if (intValue != nrRows){
    Rprintf("Invalid number of rows in data file %s.\n", fName);
    goto END_OF_FUNCTION;
  }

  /* Read 3rd line of file. This line must have the structure 'XLLCORNER XXX'. */
  if((fgets(line, 1024, fp) == NULL) || (sscanf(line, "%s %s\n", parameterName, dblVal) != 2) || (strcasecmp(parameterName, "xllcorner") != 0)){
    Rprintf ("'xllcorner' expected in data file %s\n", fName);
    goto END_OF_FUNCTION;
  }
  xllCorner = strtod(dblVal, NULL);
  
  /* Read 4th line of file. This line must have the structure 'YLLCORNER XXX'. */
  if((fgets(line, 1024, fp) == NULL) || (sscanf(line, "%s %s\n", parameterName, dblVal) != 2) || (strcasecmp(parameterName, "yllcorner") != 0)){
    Rprintf("'yllcorner' expected in data file %s\n", fName);
    goto END_OF_FUNCTION;
  }
  yllCorner = strtod (dblVal, NULL);
  
  /* Read 5th line of file. This line must have the structure 'CELLSIZE XXX'. */
  if((fgets(line, 1024, fp) == NULL) || (sscanf(line, "%s %s\n", parameterName, dblVal) != 2) || (strcasecmp(parameterName, "cellsize") != 0)){
    Rprintf("'cellsize' expected in data file %s\n", fName);
    goto END_OF_FUNCTION;
  }
  cellSize = strtod(dblVal, NULL);
  
  /* Read 6th line of file. This line must have the structure 'NODATA_value XXX'. */
  if((fgets(line, 1024, fp) == NULL) || (sscanf(line, "%s %d\n", parameterName, &noData) != 2) || (strcasecmp(parameterName, "nodata_value") != 0)){
    Rprintf("'NODATA_value' expected in data file %s\n", fName);
    goto END_OF_FUNCTION;
  }
  
  /* Read the ascii grid values into the matrix. */
  for(i = 0; i < nrRows; i++){
    for(j = 0; j < nrCols; j++){

      /* On the input file stream, read the next element - it must be a float - and store it in the matrix to fill. */
      if(fscanf(fp, "%f", &floatValue) != 1){
        Rprintf("Invalid value in data file %s\n", fName);
        goto END_OF_FUNCTION;
      }
      matrixToFill[(i * nrCols) + j] = floatValue;
    }

    /* When we reach the end of the line, we need to read the newline character so that we move past it in our
    ** file input stream. */
    floatValue = fscanf(fp, "\n");
  }


  /* If we reach this point, then no error has occured and we can this set the exit status value to 0. */
  exitStatus = 0;


END_OF_FUNCTION:

  /* Close the file and return the function's return status. */
  if(fp != NULL) fclose(fp);
  return(exitStatus);
  
}
/*####################################################################################################################*/



















/*################################################################################################## 
### function writeMat() - Write a data matrix to file.
### *************************************************
### Note: This should eventually be merged with the above "mcReadMatrix" function, but we'll 
### keep it separate for now just to make sure the basic functionality works fine.
###
### Parameters:
###  -> fName: The name of the file to write to.
###  -> mat:   The data matrix to write.
###
### Returns: 0 if everything went fine, -1 otherwise.
##################################################################################################*/
int writeMat (char *fName, int **mat)
{
  int i, j, status;
  FILE *fp;

  status = 0;
  fp = NULL;
  
  /* Open the file for writing. */
  if((fp = fopen(fName, "w")) == NULL){
    Rprintf ("Can't open data file %s for writing.\n", fName);
    status = -1;
    goto END_OF_FUNCTION;
  }

  /* Write the 'meta data'. */
  fprintf(fp, "ncols %d\n", nrCols);
  fprintf(fp, "nrows %d\n", nrRows);
  fprintf(fp, "xllcorner %.9f\n", xllCorner);
  fprintf(fp, "yllcorner %.9f\n", yllCorner);
  fprintf(fp, "cellsize %.9f\n", cellSize);
  fprintf(fp, "NODATA_value %d\n", noData);
  
  /* Write the data to file. */
  for(i = 0; i < nrRows; i++){
    for(j = 0; j < nrCols; j++){
      fprintf(fp, "%d ", mat[i][j]);
    }
    fprintf(fp, "\n");
  }



END_OF_FUNCTION:
  
  /* Close the file and return the status. */
  if (fp != NULL) fclose(fp);
  return (status);

}
/*####################################################################################################################*/

