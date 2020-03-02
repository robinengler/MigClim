####################################################################################################
### R functions for the MigClim R package.
### Authors: Robin Engler and Wim Hordijk.
###
####################################################################################################

####################################################################################################
### MigClim.migrate.
### ***************
### Initializes the MigClim method by writing the parameter values passed by the user to a file 
### on disk, and then returns this file.
###
### Input parameters:
###  -> iniDist:     string giving the full name and path of the species initial distribution.
###  -> hsMap:
###  -> rcThreshold:
###  -> envChgSteps:
###  -> dispSteps:
###  -> dispKernel:
###  -> barrier:
###  -> barrierType:
###  -> iniMatAge:
###  -> propaguleProd:
###  -> lddFreq:
###  -> lddMinDist:
###  -> lddMaxDist:
###  -> simulName:
###  -> replicateNb:
###  -> overWrite:
###  -> testMode:
###  -> fullOutput:
###  -> keepTempFiles:
###
MigClim.migrate <- function(iniDist=NULL, hsMap=NULL, rcThreshold=0,
                            envChgSteps=1, dispSteps=1, 
                            dispKernel=c(1.0,1.0), dispKernelMapNb=NULL,
                            barrier=NULL, barrierType='strong',
                            iniMatAge=1, propaguleProd=c(1.0),
                            lddFreq=0.0, lddMinDist=0, lddMaxDist=0,
                            simulName='MigClimTest', replicateNb=1, overWrite=FALSE,
                            testMode=FALSE, fullOutput=FALSE, keepTempFiles=FALSE, 
                            randomGeneratorSeed=NULL){


    ### Input data check.
    ### ****************
    # Verify the reclassification threshold value:
    #  -> must be single, integer number.
    #  -> must be in the range [0:1000].
    if(!is.numeric(rcThreshold)) stop("INPUT ERROR: 'rcThreshold' must be an integer number in the range [0:1000].")
    if(rcThreshold < 0 | rcThreshold > 1000) stop("INPUT ERROR: 'rcThreshold' must be an integer number in the range [0:1000].")
    if(rcThreshold %% 1 != 0) stop("'rcThreshold' must be an integer number.")
    #
    # Verify environmental change steps values.
    if(!is.numeric(envChgSteps)) stop("'envChgSteps' must be an integer number in the range [1:295].")
    if(envChgSteps<1 | envChgSteps > 295) stop("'envChgSteps' must be an integer number in the range [1:295].")
    if(envChgSteps%%1!=0) stop("'envChgSteps' must be a number an integer number.")
    if(!is.numeric(dispSteps)) stop("'dispSteps' must be a number in the range [1:99].")
    if(dispSteps<1 | dispSteps > 99) stop("'dispSteps' must be a number in the range [1:99].")
    if(dispSteps%%1!=0) stop("'dispSteps' must be a number an integer number.")


    ### Dispersal kernel data (dispKernel).
    # dispKernel must be either a vector of values in the range [0-1], or a string giving the 
    # basename of the raster(s) to be used as dispersal kernels.
    #
    # Case 1: dispersal kernel is a vector of numeric values. In this case we check that all values are in the
    #         range [0-1].
    if(is.numeric(dispKernel)){
        if(any(dispKernel > 1) | any(dispKernel <= 0)) stop("Values of 'dispKernel' must be numbers > 0 and <= 1")
        dispKernelLength = length(dispKernel)

    # Case 2: dispersal kernel is a series of raster grids, where each grid contains the PDisp 
    #         value associated to pixels for a given distance class. E.g. the first raster gives 
    #         the PDisp values of pixels for distance class 1, the second raster gives the values 
    #         for distance class 2, etc. This method allows to have a custom kernel for each cell 
    #         in the landscape.
    } else if(is.character(dispKernel)){

        if(is.null(dispKernelMapNb)) stop("ERROR: when 'dispKernel' is given as a series of dispersal kernel ",
                                           "ascii files, then a value for 'dispKernelMapNb' must be provided.")
        dispKernelLength = dispKernelMapNb

        # TEMPORARY: DISABLE CHECK.
        # Verify that all dispersal kernel ascii grid files are available on disk.
        #for(i in 1:dispKernelLength){
        #    # Test that the file exists on disk.
        #    fileName = paste(dispKernel, i, ".asc", sep='')
        #    if(!file.exists(fileName)) stop("ERROR: cannot find input file [", fileName, "].", sep='')
        #    
        #    # Test that the file can be loaded as a raster object.
        #    rst = try(raster(fileName), silent=T)
        #    if(class(rst)[1] != "RasterLayer") stop("ERROR: input raster [", fileName, "] could not be read.", 
        #                                            " The raster file must be in ascii grid (.asc) format")
        #}

    } else{
        stop("ERROR: 'dispKernel' must either be a vector of numeric values in the range [0-1], ", 
             "or a string giving the basename of the dispersal kernel grids.")
    }


    # Verify the seed value for the random number generator. If the user provided a seed for the
    # random number generator, the seed must be an integer in the range [0,32767] - the range of 
    # 16 bit 'integer' numbers.
    if(!is.null(randomGeneratorSeed)){
        if(!is.numeric(randomGeneratorSeed) || 
           randomGeneratorSeed %% 1 != 0 || 
           randomGeneratorSeed < 1 || 
           randomGeneratorSeed > 65535){
            stop("ERROR: 'randomGeneratorSeed' must be an integer number in the range [1,65535].")
        }
        set.seed(randomGeneratorSeed)
    } else{
        randomGeneratorSeed = "NA"
    }



    if( !is.null(barrier) && !barrierType %in% c('weak','strong') ) stop("ERROR: 'barrierType' must be either 'weak' or 'strong'.")

    if(!is.numeric(iniMatAge)) stop("INPUT ERROR: 'iniMatAge' must be an integer number > 0.")
    if(iniMatAge<=0 | iniMatAge%%1!=0) stop("INPUT ERROR: 'iniMatAge' must be an integer number > 0.")
    if(!is.numeric(propaguleProd)) stop("INPUT ERROR: values of 'propaguleProd' must be numbers > 0 and < 1.")
    if(any(propaguleProd>1) | any(propaguleProd<=0)) stop("INPUT ERROR: values for 'propaguleProd' must be numbers > 0 and <= 1.")
    if(length(propaguleProd)>1) if(propaguleProd[length(propaguleProd)]==1) stop("If the length of the 'propaguleProd'", 
                                 " vector is > 1, then the last value cannot be 1. See the MigClim user guide", 
                                 " 'MigClim.userGuide()' for detailed explanations on this paramter.")

    if(!is.numeric(lddFreq)) stop("INPUT ERROR: 'lddFreq' must be a numeric value.")
    if(lddFreq<0 | lddFreq>1) stop("'lddFreq' must be a number >= 0 and <= 1.")
    if(lddFreq>0){
        if(!is.numeric(lddMinDist)) stop("INPUT ERROR: 'lddMinDist' must be a numeric value.")
        if(!is.numeric(lddMaxDist)) stop("INPUT ERROR: 'lddMaxDist' must be a numeric value.")
        if(lddMinDist%%1!=0 | lddMaxDist%%1!=0) stop("'lddMinDist' and 'lddMaxDist' must be integer numbers.")
        if(lddMinDist <= dispKernelLength) stop("INPUT ERROR: 'lddMinDist' must be larger than the length of the 'dispKernel'.")
        if(lddMaxDist < lddMinDist) stop("INPUT ERROR: 'lddMaxDist' must be >= 'lddMinDist'.")
    } else lddMinDist = lddMaxDist = 0

    
    # Verify that 'replicateNb' is an integer number >= 1.
    if(!is.numeric(replicateNb)) stop("INPUT ERROR: 'replicateNb' must be a numeric, integer, value.")
    if(replicateNb < 1 | replicateNb %% 1 != 0) stop("INPUT ERROR: 'replicateNb' must be an integer value >= 1.")

    
    # Verify that logical parameter have logical values (TRUE or FALSE).
    if(!is.logical(overWrite))     stop("INPUT ERROR: 'overWrite' must be either TRUE or FALSE.")
    if(!is.logical(testMode))      stop("INPUT ERROR: 'testMode' must be either TRUE or FALSE.")
    if(!is.logical(fullOutput))    stop("INPUT ERROR: 'fullOutput' must be either TRUE or FALSE.")
    if(!is.logical(keepTempFiles)) stop("INPUT ERROR: 'keepTempFiles' must be either TRUE or FALSE.")


    # Verify that initial distribution (iniDist), habitat suitability maps (hsMap) and barriers are provided either
    # as a string (name and path of ascii raster files), or as a matrix/data frame. Note that these 3 parameters must
    # all be of the same type: either all be strings, or all be matrices/data frames.
    if(!is.character(iniDist) & !is.matrix(iniDist) & !is.data.frame(iniDist)) stop("ERROR: 'iniDist' must be either a string, a data frame or a matrix.")
    if(!is.character(hsMap) & !is.matrix(hsMap) & !is.data.frame(hsMap)) stop("ERROR: 'hsMap' must be either a string, a data frame or a matrix.")
    if(!is.null(barrier)) if(!is.character(barrier) & !is.matrix(barrier) & !is.data.frame(barrier) & !is.vector(barrier)) stop("ERROR: 'barrier' must be either a string, a data frame or a matrix.")

    #
    # Verify that 'iniDist', 'hsMap' and 'barrier' are of the same type: all strings, or all matrices/dataframes.
    if(is.character(iniDist)){
        # If 'iniDist' was passed as a string, then we also check that all these 3 variables do have a length <= 1, 
        # just in case the user passes a vector of strings instead of a single string.
        if(length(iniDist) > 1 | length(hsMap) > 1 | length(barrier) > 1) stop("INPUT ERROR: iniDist', 'hsMap' and ", 
                                                  "'barrier' values must be a single string, not a vector of strings.")
        if(!is.character(hsMap)) stop("INPUT ERROR: 'iniDist' and 'hsMap' must have the same format: either both ", 
                                                                               "strings or both data frames/matrices.")
        if(!is.null(barrier) & !is.character(barrier)) stop("INPUT ERROR: 'iniDist' and 'barrier' must have the same ", 
                                                           "format: either both strings or both data frames/matrices.")
       

        # If the user has entered a file name (as opposed to a dataframe or matrix) then we remove any ".asc" or ".tif"
        # extension that the user may have specified in his/her filename.
        if(grepl('\\.(asc|tif)$', x=iniDist)) iniDist = strtrim(iniDist, nchar(iniDist)-4)
        if(grepl('\\.(asc|tif)$', x=hsMap))   hsMap   = strtrim(hsMap, nchar(hsMap)-4)
        if(!is.null(barrier)) if(grepl('\\.(asc|tif)$', x=barrier)) barrier = strtrim(barrier, nchar(barrier)-4)

    } else{
        if( is.character(hsMap) | is.character(barrier) ) stop("INPUT ERROR: 'iniDist', 'hsMap' ", 
                            "and 'barrier' must have the same format: either all strings or all data frames/matrices.")
    }



    ### Detect raster input type given by user.
    ### **************************************
    # The user can provide the input for "iniDist", "hsMap" and "barrier" either as:
    #  -> dataframe or matrix.
    #  -> raster in format: ascii grid (.asc), geo-tiff (.tif), ESRI raster (no extension), or R raster (no extension).
    #
    if(is.data.frame(iniDist) | is.matrix(iniDist)){
        # Case 1: user provided iniDist as matrix or dataframe. Note that if the user input is a matrix, we 
        # convert it to a data frame.
        if(is.matrix(iniDist)) iniDist = as.data.frame(iniDist)  
        RExt = ".DataFrame"

    } else{
        # Case 2: user provided iniDist as a string (file name and path). In this case we test whether we can find
        #         a file with the specified name. We try different extensions (.asc, .tif, nothing) and depending on 
        #         the extension we deduce the type of data.
        RExt = NA
        if(file.exists(iniDist)){
            rst = try(raster(iniDist), silent=T)
            if(class(rst)[1] == "RasterLayer") RExt = ""
            rm(rst)
        } else if(file.exists(paste(iniDist,".tif",sep=""))){
            rst = try(raster(paste(iniDist,".tif",sep="")), silent=T)
            if(class(rst)[1] == "RasterLayer") RExt = ".tif"
            rm(rst)
        } else if(file.exists(paste(iniDist,".asc",sep=""))){
            rst = try(raster(paste(iniDist,".asc",sep="")), silent=T)
            if(class(rst)[1] == "RasterLayer") RExt = ".asc"
            rm(rst)
        } else stop("ERROR: the 'iniDist' raster could not be found. Make sure the file ", getwd(), "/", iniDist, 
                                                                               " exists and that its path is correct.")
        if(is.na(RExt)) stop("ERROR: input data not recognized. The 'iniDist' input raster data must be in one of the", 
                        " following formats: ascii grid (.asc), geoTiff (.tif), ESRI grid or R raster (no extension).")        
    }
    


    # If the user chose to not allow overwriting of existing files (overWrite==F)
    # then we check that no future output file already exists.
    if(overWrite==F){

        # Check if output directory exists
        if(file.exists(simulName)) stop("The output directory '", getwd(), "/", simulName, "' already exists. Delete this directory or set 'overWrite=TRUE' in the function's parameters.\n")

        ### Check if any output ".asc" files already exist.
        if(RExt!=".asc"){
            if(file.exists(paste(basename(iniDist),".asc",sep=""))) stop("The output file '", getwd(), "/", paste(basename(iniDist),".asc",sep=""), "' already exists. \n Delete this file or set 'overWrite=TRUE' in the function's parameters.\n")
            for(J in 1:envChgSteps) if(file.exists(paste(basename(hsMap), J,".asc",sep=""))) stop("The output file '", getwd(), "/", paste(basename(hsMap), J,".asc",sep=""), "' already exists. \n Delete this file or set 'overWrite=TRUE' in the function's parameters.\n")
            if(!is.null(barrier)) if(file.exists(paste(basename(barrier),".asc",sep=""))) stop("The output file '", getwd(), "/", paste(basename(barrier),".asc",sep=""), "' already exists. \n Delete this file or set 'overWrite=TRUE' in the function's parameters.\n")
        }
        if(RExt==".asc"){
            if(iniDist!=basename(iniDist)) if(file.exists(paste(basename(iniDist),".asc",sep=""))) stop("The output file '", getwd(), "/", paste(basename(iniDist),".asc",sep=""), "' already exists. \n Delete this file or set 'overWrite=TRUE' in the function's parameters.\n")
            if(hsMap!=basename(hsMap)) for(J in 1:envChgSteps) if(file.exists(paste(basename(hsMap), J,".asc",sep=""))) stop("The output file '", getwd(), "/", paste(basename(hsMap), J,".asc",sep=""), "' already exists. \n Delete this file or set 'overWrite=TRUE' in the function's parameters.\n")
            if(!is.null(barrier)) if(barrier!=basename(barrier)) if(file.exists(paste(basename(barrier),".asc",sep=""))) stop("The output file '", getwd(), "/", paste(basename(barrier),".asc",sep=""), "' already exists. \n Delete this file or set 'overWrite=TRUE' in the function's parameters.\n")
        }
        if(RExt==".DataFrame"){
            if(file.exists(paste(simulName, ".InitialDist.asc", sep=""))) stop("The output file '", getwd(), "/", paste(simulName, ".InitialDist.asc", sep=""), "' already exists. \n Delete this file or set 'overWrite=TRUE' in the function's parameters.\n")
            for(J in 1:envChgSteps) if(file.exists(paste(simulName, ".HSmap", J, ".asc", sep=""))) stop("The output file '", getwd(), "/", paste(simulName, ".HSmap", J, ".asc", sep=""), "' already exists. \n Delete this file or set 'overWrite=TRUE' in the function's parameters.\n")
            if(!is.null(barrier)) if(file.exists(paste(simulName, ".Barrier.asc", sep=""))) stop("The output file '", getwd(), "/", paste(simulName, ".Barrier.asc", sep=""), "' already exists. \n Delete this file or set 'overWrite=TRUE' in the function's parameters.\n")
        }
    }



    # If the user has given the input as a matrix/dataframe, then we verify that
    # the data has the correct format. The correct format is as follows:
    # 'iniDist' = a data frame with 3 columns: X coordinate, Y coordinate and the species' initial distribution (0 or 1 values only).
    # 'hsMap' = a dataframe with ncol = envChgSteps. values must be in the range [0:1000]
    # 'barrier' = a dataframe or vector containing only values of 0 or 1.
    #
    if(RExt==".DataFrame"){
        cat("Converting data to ascii grid format... \n")

        # Convert all input data to data frame objects.
        if(is.matrix(hsMap)) hsMap = as.data.frame(hsMap)
        if(is.vector(hsMap)) hsMap = as.data.frame(hsMap)
        if(is.null(barrier) | is.character(barrier)) useBarrier = FALSE else useBarrier = TRUE
        if(is.matrix(barrier)) barrier = as.data.frame(barrier)
        if(is.vector(barrier)) barrier = as.data.frame(barrier)

        ### Verify all inputs are of the data frame type (note: iniDist has already been checked earlier).
        if(!is.data.frame(hsMap)) stop("INPUT ERROR: the 'hsMap' data could not be converted to a dataframe. All inputs must be of the same type. \n")
        if(useBarrier) if(!is.data.frame(barrier)) stop("INPUT ERROR: the 'barrier' data could not be converted to a dataframe. all inputs must be of the same type. \n")

        ### Verify all data frames have the correct number of rows and columns.
        if(ncol(iniDist)!=3) stop("Data input error. When entering 'iniDist' as a data frame or matrix, the data frame must have exactly 3 columns (in this order): X and Y coordinates, Initial distribution of the species. \n")
        if(ncol(hsMap)!=envChgSteps) stop("Data input error. When entering 'hsMap' as a data frame or matrix, the data frame must have a number of columns equal to envChgSteps. \n")
        if(nrow(hsMap)!=nrow(iniDist))  stop("Data input error. 'iniDist' and 'hsMap' must have the same number of rows.\n")
        if(useBarrier){
            if(ncol(barrier)!=1) stop("Data input error. When entering 'barrier' as a data frame, matrix or vector, the data must have a excatly 1 column. \n")
            if(nrow(barrier)!=nrow(iniDist))  stop("Data input error. 'iniDist' and 'barrier' must have the same number of rows.\n")
        }

        ### Verify all data frames contain meaningful values.
        if(any(is.na(match(unique(iniDist[,3]), c(0,1))))) stop("INPUT ERROR: the 3rd column of 'iniDist' should contain only values of 0 or 1. \n")
        if(any(hsMap<0) | any(hsMap>1000)) stop("INPUT ERROR: all values in 'hsMap' must be in the range [0:1000]. \n")
        if(useBarrier) if(any(is.na(match(unique(barrier[,1]), c(0,1))))) stop("INPUT ERROR: 'barrier' should contain only values of 0 or 1. \n")

        ### Convert data frames to ascii grid files.
        CreatedASCII = paste(simulName, c("InitialDist.asc", paste("HSmap", 1:envChgSteps, ".asc", sep="")), sep=".")
        dataframe2asc(cbind(iniDist[,c(2,1,3)], hsMap), outdir=getwd(), filenames=CreatedASCII, gz=FALSE)
        if(useBarrier){
            dataframe2asc(cbind(iniDist[,c(2,1)], barrier), outdir=getwd(), filenames=paste(simulName, ".Barrier", sep=""), gz=FALSE)
            CreatedASCII = c(CreatedASCII, paste(simulName, ".Barrier.asc", sep=""))
            barrier = paste(simulName, ".Barrier", sep="")
        }
        iniDist = paste(simulName, ".InitialDist", sep="")
        hsMap = paste(simulName, ".HSmap", sep="")
        RExt = ".asc"
    }


    # Verify that all the input raster files exist.
    if(!file.exists(paste(iniDist,RExt,sep=""))) stop('missing input file: ', iniDist, RExt)
    for(i in 1:envChgSteps){
    
        # TEMPORARY: DISABLE CHECK.
        #if(!file.exists(paste(hsMap,i,RExt,sep=""))) stop(paste("The 'hsMap' file '", hsMap, i, RExt, "' could not be found.\n",
        #                                                        "The naming convention for hsMap files is 'hsMap basename + 1', 'hsMap basename + 2', etc...\n",
        #                                                        "e.g. if you set 'hsMap='habitatSuitMap'' then your first hsMap file must be named 'habitatSuitMap1'.\n",
        #                                                        "the following hsMap file must be named 'habitatSuitMap2', 'habitatSuitMap3' and so on.\n", sep=""))
    }
    if(!is.null(barrier)) if(!file.exists(paste(barrier,RExt,sep=""))) stop(paste("The 'barrier' file '", barrier, RExt, "' could not be found.\n", sep=""))


    # If the input format is not ascii grid, then we convert the files to ascii grid format.
    # Note that we store the names of the created ascii files in the "CreatedASCII" object.
    if(RExt!=".asc"){
        cat("Converting data to ascii grid format... \n")
        rst = raster(paste(iniDist,RExt,sep=""))
        iniDist = basename(iniDist)
        rst2 = writeRaster(rst, filename=paste(iniDist,".asc",sep=""), format="ascii", overwrite=TRUE, datatype="INT2S", NAflag=-9999)
        CreatedASCII = paste(iniDist,".asc",sep="")
        for (J in 1:envChgSteps){
          rst = raster(paste(hsMap,J,RExt,sep=""))
          rst2 = writeRaster(rst, filename=paste(basename(hsMap),J,".asc",sep=""), format="ascii", overwrite=TRUE, datatype="INT2S", NAflag=-9999)
          CreatedASCII = c(paste(basename(hsMap),J,".asc",sep=""), CreatedASCII)
        }
        hsMap = basename(hsMap)
        if(!is.null(barrier)){
            rst = raster(paste(barrier,RExt,sep=""))
            barrier = basename(barrier)
            rst2 = writeRaster(rst, filename=paste(barrier,".asc",sep=""), format="ascii", overwrite=TRUE, datatype="INT2S", NAflag=-9999)
            CreatedASCII = c(paste(barrier,".asc",sep=""), CreatedASCII)
        }
        rm(rst,rst2)
    }

    
    
    # Raster data (ascii grid)
    # ***********************
    # Verify that all input raster files can be found in current working directory.
    # If not, we copy the files to the working directory.
    if (RExt==".asc"){
        if(iniDist!=basename(iniDist)){
            file.copy(from=paste(iniDist,".asc",sep=""), to=paste(basename(iniDist),".asc",sep=""), overwrite=T)
            iniDist = basename(iniDist)
            if(exists("CreatedASCII")) CreatedASCII = c(paste(iniDist,".asc",sep=""), CreatedASCII) else CreatedASCII = paste(iniDist,".asc",sep="")
        }
        if(hsMap!=basename(hsMap)){
            for(i in 1:envChgSteps){
                file.copy(from=paste(hsMap,i,".asc",sep=""), to=paste(basename(hsMap),i,".asc",sep=""), overwrite=T)
                if(exists("CreatedASCII")){
                    CreatedASCII = c(paste(basename(hsMap),i,".asc",sep=""), CreatedASCII)
                } else CreatedASCII = paste(basename(hsMap),i,".asc",sep="")
    		}
    		hsMap = basename(hsMap)
    	}
        if(!is.null(barrier)){
            if(barrier!=basename(barrier)){
               file.copy(from=paste(barrier,".asc",sep=""), to=paste(basename(barrier),".asc",sep=""), overwrite=T)
               barrier = basename(barrier)
               if(exists("CreatedASCII")) CreatedASCII = c(paste(barrier,".asc",sep=""), CreatedASCII) else CreatedASCII = paste(barrier,".asc",sep="")
            }
        }
    }

    # Verify that all ascii grid files have a correct structure and that their NoData value, 
    # if any, is set to -9999.
    noDataVal = get_nodata_value(paste(iniDist,".asc",sep=""))
    if(!is.na(noDataVal)){
        if(noDataVal == "ErrorInFile") stop("INPUT ERROR: 'iniDist' ascii grid has incorrect structure.")
        if(noDataVal >= 0)             stop("INPUT ERROR: 'iniDist' ascii grid must have 'NoData' values < 0.")
    }
    for(i in 1:envChgSteps){
        # TEMPORARY: DISABLE CHECK.
    	#noDataVal = get_nodata_value(paste0(hsMap, i, '.asc'))
        #if(!is.na(noDataVal)){
        #    if(noDataVal == "ErrorInFile") stop("INPUT ERROR: one or more 'hsMap' ascii grid files have incorrect structure.")
        #    if(noDataVal >= 0)             stop("INPUT ERROR: all 'hsMap' ascii grid files must have 'NoData' values set to a number < 0.")
        #}
    }
    if(!is.null(barrier)){
        noDataVal = get_nodata_value(paste(barrier,".asc",sep=""))
        if(!is.na(noDataVal)){
            if(noDataVal == "ErrorInFile") stop("INPUT ERROR: 'barrier' ascii grid has incorrect structure.")
            if(noDataVal >= 0)             stop("INPUT ERROR: 'barrier' ascii grid must have 'NoData' values < 0.")
        }
    }
    rm(noDataVal)


    # Verify that all raster have exactly the same dimensions and that they contain appropriate 
    # values: 'iniDist' and 'barrier' should contain only values of 0 or 1. 'hsMap' should contain 
    # only values in the range [0:1000].
    rst = raster(paste(iniDist,'.asc',sep=''))
    nrRows = nrow(rst)
    nrCols = ncol(rst)
    if(any(is.na(match(raster::unique(rst), c(0,1))))) stop("INPUT ERROR: the 'iniDist' raster should contain only values of 0 or 1.")
    for(i in 1:envChgSteps){
        # TEMPORARY: DISABLE CHECK.
        #rst = raster(paste(hsMap,i,'.asc',sep=''))
        #if(nrow(rst) != nrRows | ncol(rst) != nrCols) stop("INPUT ERROR: not all your rasters input data have the same dimensions.")
        #if(cellStats(rst,'min') < 0 | cellStats(rst,'max') > 1000) stop("INPUT ERROR: all habitat suitability rasters must have values in the range [0:1000].")
        #rm(rst)
    }
    if(!is.null(barrier)){
        rst = raster(paste(barrier,".asc",sep=""))
        if(nrow(rst)!=nrRows | ncol(rst)!=nrCols) stop("INPUT ERROR: not all your rasters input data have the same dimensions.")
        if(any(is.na(match(raster::unique(rst), c(0,1))))) stop("INPUT ERROR: the 'barrier' raster should contain only values of 0 or 1.")
        rm(rst)
    }



    # Create output directory.
    if(file.exists(simulName)) unlink(simulName, recursive=T)
    if(!dir.create(simulName)) stop("unable to create a '", simulName,"'subdirectory in the current workspace. ", 
                                     "Make sure the '", simulName,"'subdirectory does not already exists and that ", 
                                     "you have write permission in the current workspace.")

    # Write the "simulName_params.txt" file to disk.
    fileName = paste(simulName, "/", simulName, "_params.txt", sep="")
    write(paste("nrRows", nrRows), file=fileName, append=F)
    write(paste("nrCols", nrCols), file=fileName, append=T)
    write(paste("iniDist", iniDist), file=fileName, append=T)
    write(paste("hsMap", hsMap), file=fileName, append=T)
    write(paste("rcThreshold", rcThreshold), file=fileName, append=T)
    write(paste("envChgSteps", envChgSteps), file=fileName, append=T)
    write(paste("dispSteps", dispSteps), file=fileName, append=T)
    write(paste("dispDist", dispKernelLength), file=fileName, append=T)
    if(is.numeric(dispKernel)) write(c("dispKernel", dispKernel), file=fileName, append=T, ncolumns=dispKernelLength+1)
    if(is.character(dispKernel)) write(c("dispKernelFile", dispKernel), file=fileName, append=T)
    if(!is.null(barrier)) write(paste("barrier", barrier), file=fileName, append=T)
    if(is.null(barrier)) write("barrier NA", file=fileName, append=T)
    write(paste("barrierType", barrierType), file=fileName, append=T)
    write(paste("iniMatAge", iniMatAge), file=fileName, append=T)
    write(paste("fullMatAge", iniMatAge + length(propaguleProd)), file=fileName, append=T)
    write(c("propaguleProd", propaguleProd), file=fileName, append=T, ncolumns=length(propaguleProd)+1)
    write(paste("lddFreq", lddFreq), file=fileName, append=T)
    write(paste("lddMinDist", lddMinDist), file=fileName, append=T)
    write(paste("lddMaxDist", lddMaxDist), file=fileName, append=T)
    write(paste("fullOutput", fullOutput), file=fileName, append=T)
    write(paste("replicateNb", replicateNb), file=fileName, append=T)
    write(paste("simulName", simulName), file=fileName, append=T)
    write(paste("randomGeneratorSeed", randomGeneratorSeed), file=fileName, append=T)

    
    # End of test mode.
    # ****************
    # If the command is running in test mode, exit here.
    if(testMode){
       cat("### Test for", simulName, "completed sucessfully.\n")
       # Delete the created output directory and exit.
       unlink(simulName, recursive=T)
       return(envChgSteps)
    }

    
    # Call the main migclim C function.
    # ********************************
    cat("Starting simulation for ", simulName, "...\n")
    input_prefix = file.path(simulName, simulName)
    migrate = .C('mcMigrate', paste0(input_prefix, '_params.txt'), nr=integer(1))

    # If ascii grids were created by the MigClim.init() function, then we delete them here, unless 
    # the user has set "keepTempFiles = TRUE".
    if(keepTempFiles) rm(CreatedASCII)
    if(exists("CreatedASCII")){
        for (i in 1:length(CreatedASCII)) unlink(CreatedASCII[i])
        rm(CreatedASCII)
    }

    
    # Create average output across all replicates.
    # *******************************************
    # If the user has set replicateNb > 1 then we generate a final, averaged, output.
    # The individual outputs are conserved, though.
    if(replicateNb > 1){

        # Average '_stats.txt' files.
        merged_table = merge_files(prefix       = paste0(input_prefix, '_'), 
                                   suffix       = '_stats.txt', 
                                   replicate_nb = replicateNb, 
                                   merge_mode   = 'sum')
        write.table(round(merged_table/replicateNb, 2), 
                    file=paste0(input_prefix, '_stats.txt'), quote=F, row.names=F, sep='\t')

        # Average '_summary.txt' files
        merged_table = merge_files(prefix       = paste0(input_prefix, '_'), 
                                   suffix       = '_summary.txt', 
                                   replicate_nb = replicateNb, 
                                   merge_mode   = 'append')
        merged_table[replicateNb + 1, 1] = simulName
        merged_table[replicateNb + 1, 2:ncol(merged_table)] = round(
            apply(merged_table[1:replicateNb,2:ncol(merged_table)], 2, mean), 2)
        write.table(merged_table, 
                    file=paste0(input_prefix,'_summary.txt'), quote=F, row.names=F, sep='\t')
    }


    # Return the number of output files created.
    if(migrate$nr==envChgSteps) cat('### Simulation ', simulName, " completed successfully. Outputs stored in ", 
                                                                                 getwd(), "/", simulName, "\n", sep="")
    return(migrate$nr)
}

########################################################################################################################



####################################################################################################
get_nodata_value <- function(raster_file){
    ### Checks the structure of an ascii grid file, and if the structure is correct, returns its
    ### "NoData" value. If the structure of the file is not correct, the function returns the 
    ### string "ErrorInFile". If the "NoData" value is missing, which is possible, since this field
    ### is optional, the function returns NA.
    ###
    ### Input parameters:
    ###  -> raster_file: full name and path of the ascii grid file to check.
    ###

    # Verify the ascii file has the correct structure.
    # ***********************************************
    noDataVal = "ErrorInFile"
    fileStruct = c("ncols","nrows","xllcorner","yllcorner","cellsize")
    fileStruct2 = c(5,5,9,9,8)
    for(J in 0:4){
        lineVal = scan(file=raster_file, what="character", nlines=1, skip=J, quiet=TRUE)
        if(length(lineVal)!=2) return(noDataVal)
        if(nchar(lineVal[1])!=fileStruct2[J+1]) return(noDataVal)
        if(length(grep(fileStruct[J+1], lineVal[1], ignore.case=T))!=1) return(noDataVal)
    }


    # Get "NoData" value.
    # ******************
    # This line is optional in the file. If the line is missing we return NA.
    lineVal = scan(file=raster_file, what="character", nlines=1, skip=5, quiet=TRUE)
    if(length(lineVal)<2) return(noDataVal)
    noDataVal = NA
    if(length(grep("NODATA_value", lineVal[1], ignore.case=T)) == 1 & 
       nchar(lineVal[1])==12) noDataVal = as.numeric(lineVal[2])

    return(noDataVal)
}
####################################################################################################



####################################################################################################
merge_files <- function(prefix=NULL, suffix=NULL, replicate_nb=1, merge_mode='append'){
    # Merge the content of files
    
    
    # Verify input.
    stopifnot(is.numeric(replicate_nb))
    stopifnot(replicate_nb >= 1)
    stopifnot(merge_mode %in% c('append','sum'))
    
    # Loop through all replicates and merge data.
    for(i in 1:replicate_nb){
        # Load input table.
        input_file = paste0(prefix, i, suffix)
        if(!file.exists(input_file)) stop('missing input file:', input_file)
        input_table = read.table(input_file, header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)
        
        # Merge table with other replicates.
        if(i == 1){
            merged_table = input_table
        } else{
            if(merge_mode == 'append') merged_table = rbind(merged_table, input_table)
            if(merge_mode == 'sum')    merged_table = merged_table + input_table
        }
    }
    return(merged_table)
}
####################################################################################################
