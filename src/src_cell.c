
#include "migclim.h"

bool mcSrcCell (int i, int j, int **curState, int **pxlAge,
                int loopID, int habSuit, int **barriers){
    /* ************************************************************************
    ** For a given sink cell, search for a suitable source cell that can
    ** colonize it.
    **
    ** Parameters:
    **   i, j    : coordinate (row, column) of the sink cell.
    **   curState: matrix of the current state of the cellular automaton.
    **   pxlAge  : matrix of the 'age' of each colonized cell.
    **   loopID  : coded value of current environmental change step and
    **             dispersal step loop.
    **   habSuit : habitat suitability of the sink cell.
    **   barriers: pointer to the barriers matrix.
    **
    ** Return value:
    **   'true' if a source cell was found, 'false' otherwise.
    ** ********************************************************************* */
    int    k, l, realDist;
    double colonization_probability, rnd, kernel_value;

    /* Search for a potential source cell
    ** **********************************
    ** i, j are the coordinates of the sink cell.
    ** k, l are the coordinates of the potential source cell. */
    for(k = i - dispDist; k <= i + dispDist; k++){
    for(l = j - dispDist; l <= j + dispDist; l++){

        /* 1. Test of basic conditions to see if cell is a potential source
        **    cell:
        **     -> cell must be within the limits of the matrix's extent.
        **     -> cell must be colonized, but not during the current loop.
        **     -> cell must have reached "initial maturity", otherwise it
        **        cannot produce seeds.
        */
        if((k >= 0) && (k < nrRows) && (l >= 0) && (l < nrCols)){
        if((curState[k][l] > 0) && (curState[k][l] != loopID)){
        if (pxlAge[k][l] >= iniMatAge){

            /* 2. Compute the distance between sink and source cell in cell
            **    units. The distance must be <= maximum dispersal distance. */
            realDist = (int)round(sqrt((k-i)*(k-i) + (l-j)*(l-j)));
            if((realDist > 0) && (realDist <= dispDist)){

                /* 3. Compute the probability of colonization of the sink cell
                **    based on:
                **     - distance between source and sink cells (dispersal
                **       kernel value).
                **     - habitat suitability of the sink cell.
                **     - age of the source cell. */
                if(singleKernelMode){
                    kernel_value = dispKernel[realDist-1];
                } else{
                    kernel_value = dispKernel[(realDist-1) * nrRows * nrCols + k * nrCols + l];
                }
                colonization_probability = kernel_value * (habSuit / 1000.0);

                /* If the source cell has not reached its full maturity,
                ** rescale the probability by the source cell's current
                ** propagule production potential. */
                if(pxlAge[k][l] < fullMatAge){
                    colonization_probability = colonization_probability *
                                               propaguleProd[pxlAge[k][l] - iniMatAge];
                }

                /* Generate a random number */
                rnd = UNIF01;
                if(rnd < colonization_probability || colonization_probability == 1.0){
                    /* If this point is reached, the last thing to check is
                    ** whether there is a 'barrier' cell between the source
                    ** and sink cells.
                    ** This check is done last because it requires significant
                    ** compute time. */
                    if(!useBarrier) return(true);
                    if(!mcIntersectsBarrier (i, j, k, l, barriers)) return(true);
                }
            }
        }
        }
        }
    }
    }

    /* No source cell was found: return false. */
    return(false);
}

