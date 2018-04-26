/*******************************************************************************
 ** ArrayOperations.c
 ** Part of the mutual information toolbox
 **
 ** Contains functions to floor arrays, and to merge arrays into a joint
 ** state.
 ** 
 ** Author: Adam Pocock
 ** Created 17/2/2010
 ** Updated - 22/02/2014 - Added checking on calloc.
 **
 ** Copyright 2010-2017 Adam Pocock, The University Of Manchester
 ** www.cs.manchester.ac.uk
 **
 ** This file is part of MIToolbox, licensed under the 3-clause BSD license.
 *******************************************************************************/

#include <errno.h>
#include "MIToolbox/MIToolbox.h"
#include "MIToolbox/ArrayOperations.h"

void* checkedCalloc(size_t vectorLength, size_t sizeOfType) {
    void *allocated = CALLOC_FUNC(vectorLength, sizeOfType);
    if(allocated == NULL) {
#ifdef MEX_IMPLEMENTATION
        /* This call returns control to Matlab, with the associated error message */
        mexErrMsgTxt("Failed to allocate memory\n");
#elif defined(C_IMPLEMENTATION)
        fprintf(stderr, "Error: %s\nAttempted to allocate %lu length of size %lu\n", strerror(errno), vectorLength, sizeOfType);
        exit(EXIT_FAILURE);
#endif
    }
    return allocated;
}

void incrementVector(double* vector, int vectorLength) {
    /*This is used to map from C indices to MATLAB indices*/
    int i = 0;
    for (i = 0; i < vectorLength; i++) {
        vector[i]++;
    }/*for length of array */
}/* incrementVector(double*,int) */

void printDoubleVector(double *vector, int vectorLength) {
    int i;
    for (i = 0; i < vectorLength; i++) {
        if (vector[i] > 0) {
            printf("Value at i=%d, is %f\n",i,vector[i]);
        }
    }/*for number of items in vector*/
}/*printDoubleVector(double*,int)*/

void printIntVector(int *vector, int vectorLength) {
    int i;
    for (i = 0; i < vectorLength; i++) {
        printf("Value at i=%d, is %d\n",i,vector[i]);
    }/*for number of items in vector*/
}/*printIntVector(int*,int)*/

void printUintVector(uint *vector, int vectorLength) {
    int i;
    for (i = 0; i < vectorLength; i++) {
        printf("Value at i=%d, is %d\n",i,vector[i]);
    }/*for number of items in vector*/
}/*printUintVector(int*,int)*/

uint **generateIntIndices(uint *featureMatrix, uint noOfSamples, uint noOfFeatures) {
    int j;

    uint **feature2D = (uint **) checkedCalloc(noOfFeatures,sizeof(uint *));
    
    for (j = 0; j < noOfFeatures; j++) {
        feature2D[j] = featureMatrix + j*noOfSamples;
    }

    return feature2D;
}

double **generateDoubleIndices(double *featureMatrix, uint noOfSamples, uint noOfFeatures) {
    int j;

    double **feature2D = (double **) checkedCalloc(noOfFeatures,sizeof(double *));
    
    for (j = 0; j < noOfFeatures; j++) {
        feature2D[j] = featureMatrix + j*noOfSamples;
    }

    return feature2D;
}

int maxState(uint *vector, int vectorLength) {
    int i, max;
    max = 0;
    for (i = 0; i < vectorLength; i++) {
        if (vector[i] > max) {
            max = vector[i];
        }
    }
    return max + 1;
}

int numberOfUniqueValues(double *featureVector, int vectorLength) {
    int uniqueValues = 0;
    double *valuesArray = (double *) checkedCalloc(vectorLength,sizeof(double));

    int found = 0;
    int j = 0;
    int i;

    for (i = 0; i < vectorLength; i++) {
        found = 0;
        j = 0;
        while ((j < uniqueValues) && (found == 0)) {
            if (valuesArray[j] == featureVector[i]) {
                found = 1;
                featureVector[i] = (double) (j+1);
            }
            j++;
        }
        if (!found) {
            valuesArray[uniqueValues] = featureVector[i];
            uniqueValues++;
            featureVector[i] = (double) uniqueValues;
        }
    }/*for vectorlength*/

    FREE_FUNC(valuesArray);
    valuesArray = NULL;

    return uniqueValues;
}/*numberOfUniqueValues(double*,int)*/

/******************************************************************************* 
 ** normaliseArray takes an input vector and writes an output vector
 ** which is a normalised version of the input, and returns the number of states
 ** A normalised array has min value = 0, max value = number of states
 ** and all values are integers
 **
 ** length(inputVector) == length(outputVector) == vectorLength otherwise there
 ** is a memory leak
 *******************************************************************************/
int normaliseArray(double *inputVector, uint *outputVector, int vectorLength) {
    int minVal = 0;
    int maxVal = 0;
    int currentValue;
    int i;

    if (vectorLength > 0) {
        int* tempVector = (int*) checkedCalloc(vectorLength,sizeof(int));
        minVal = (int) floor(inputVector[0]);
        maxVal = (int) floor(inputVector[0]);

        for (i = 0; i < vectorLength; i++) {
            currentValue = (int) floor(inputVector[i]);
            tempVector[i] = currentValue;

            if (currentValue < minVal) {
                minVal = currentValue;
            } else if (currentValue > maxVal) {
                maxVal = currentValue;
            }
        }/*for loop over vector*/

        for (i = 0; i < vectorLength; i++) {
            outputVector[i] = tempVector[i] - minVal;
        }

        maxVal = (maxVal - minVal) + 1;

        FREE_FUNC(tempVector);
        tempVector = NULL;
    }

    return maxVal;
}/*normaliseArray(double*,double*,int)*/


/*******************************************************************************
 ** mergeArrays takes in two arrays and writes the joint state of those arrays
 ** to the output vector, returning the number of joint states
 **
 ** the length of the vectors must be the same and equal to vectorLength
 ** outputVector must be malloc'd before calling this function
 *******************************************************************************/
int mergeArrays(uint *firstVector, uint *secondVector, uint *outputVector, int vectorLength) {
    int firstNumStates = maxState(firstVector,vectorLength);
    int secondNumStates = maxState(secondVector,vectorLength);
    uint *stateMap = (uint *) checkedCalloc(firstNumStates*secondNumStates,sizeof(uint));
    int stateCount = 1;
    int i, curIndex;

    for (i = 0; i < vectorLength; i++) {
        curIndex = firstVector[i] + (secondVector[i] * firstNumStates);
        if (stateMap[curIndex] == 0) {
            stateMap[curIndex] = stateCount;
            stateCount++;
        }
        outputVector[i] = stateMap[curIndex];
    }

    FREE_FUNC(stateMap);
    stateMap = NULL;

    return stateCount;
}/*mergeArrays(double *,double *,int *, int)*/

/*******************************************************************************
 ** discAndMergeArrays takes in two arrays, discretises them and writes the joint
 ** state of those arrays to the output vector, returning the number of joint 
 ** states.
 **
 ** the length of the vectors must be the same and equal to vectorLength
 ** outputVector must be malloc'd before calling this function
 *******************************************************************************/
int discAndMergeArrays(double *firstVector, double *secondVector, uint *outputVector, int vectorLength) {
    uint *firstNormalisedVector;
    uint *secondNormalisedVector;
    int stateCount;

    firstNormalisedVector = (uint *) checkedCalloc(vectorLength,sizeof(uint));
    secondNormalisedVector = (uint *) checkedCalloc(vectorLength,sizeof(uint));

    normaliseArray(firstVector,firstNormalisedVector,vectorLength);
    normaliseArray(secondVector,secondNormalisedVector,vectorLength);

    stateCount = mergeArrays(firstNormalisedVector,secondNormalisedVector,outputVector,vectorLength);

    FREE_FUNC(firstNormalisedVector);
    FREE_FUNC(secondNormalisedVector);

    firstNormalisedVector = NULL;
    secondNormalisedVector = NULL;

    return stateCount;
}/*discAndMergeArrays(double *,double *,int *, int)*/

int mergeArraysArities(uint *firstVector, int numFirstStates, uint *secondVector, int numSecondStates, uint *outputVector, int vectorLength) {
    int i;
    int totalStates;
    int firstStateCheck, secondStateCheck;

    firstStateCheck = maxState(firstVector,vectorLength);
    secondStateCheck = maxState(secondVector,vectorLength);

    if ((firstStateCheck <= numFirstStates) && (secondStateCheck <= numSecondStates)) {
        for (i = 0; i < vectorLength; i++) {
            outputVector[i] = firstVector[i] + (secondVector[i] * numFirstStates) + 1;
        }
        totalStates = numFirstStates * numSecondStates;
    } else {
        totalStates = -1;
    }

    return totalStates;
}/*mergeArraysArities(int *,int,int *,int,int *,int)*/

int discAndMergeArraysArities(double *firstVector, int numFirstStates, double *secondVector, int numSecondStates, uint *outputVector, int vectorLength) {
    uint *firstNormalisedVector;
    uint *secondNormalisedVector;
    int i;
    int totalStates;
    int firstStateCheck, secondStateCheck;

    firstNormalisedVector = (uint *) checkedCalloc(vectorLength,sizeof(uint));
    secondNormalisedVector = (uint *) checkedCalloc(vectorLength,sizeof(uint));

    firstStateCheck = normaliseArray(firstVector,firstNormalisedVector,vectorLength);
    secondStateCheck = normaliseArray(secondVector,secondNormalisedVector,vectorLength);

    if ((firstStateCheck <= numFirstStates) && (secondStateCheck <= numSecondStates)) {
        for (i = 0; i < vectorLength; i++) {
            outputVector[i] = firstNormalisedVector[i] + (secondNormalisedVector[i] * numFirstStates) + 1;
        }
        totalStates = numFirstStates * numSecondStates;
    } else {
        totalStates = -1;
    }

    FREE_FUNC(firstNormalisedVector);
    FREE_FUNC(secondNormalisedVector);

    firstNormalisedVector = NULL;
    secondNormalisedVector = NULL;

    return totalStates;
}/*mergeArraysArities(double *,int,double *,int,double *,int)*/

int mergeMultipleArrays(double *inputMatrix, uint *outputVector, int matrixWidth, int vectorLength) {
    int i = 0;
    int currentIndex;
    int currentNumStates;
    uint *normalisedVector = (uint *) checkedCalloc(vectorLength, sizeof(uint));

    if (matrixWidth > 1) {
        currentNumStates = discAndMergeArrays(inputMatrix, (inputMatrix + vectorLength), outputVector, vectorLength);
        for (i = 2; i < matrixWidth; i++) {
            currentIndex = i * vectorLength;
            normaliseArray(inputMatrix+currentIndex, normalisedVector, vectorLength);
            currentNumStates = mergeArrays(outputVector, normalisedVector, outputVector, vectorLength);
        }
    } else {
        currentNumStates = normaliseArray(inputMatrix, normalisedVector, vectorLength);
        for (i = 0; i < vectorLength; i++) {
            outputVector[i] = normalisedVector[i];
        }
    }

    FREE_FUNC(normalisedVector);
    normalisedVector = NULL;

    return currentNumStates;
}/*mergeMultipleArrays(double *, double *, int, int)*/

int mergeMultipleArraysArities(double *inputMatrix, uint *outputVector, int matrixWidth, int *arities, int vectorLength) {
    int i = 0;
    int currentIndex;
    int currentNumStates;
    uint *normalisedVector = (uint *) checkedCalloc(vectorLength, sizeof(uint));

    if (matrixWidth > 1) {
        currentNumStates = discAndMergeArraysArities(inputMatrix, arities[0], (inputMatrix + vectorLength), arities[1], outputVector, vectorLength);
        for (i = 2; i < matrixWidth; i++) {
            currentIndex = i * vectorLength;
            normaliseArray(inputMatrix+currentIndex, normalisedVector, vectorLength);
            currentNumStates = mergeArraysArities(outputVector, currentNumStates, normalisedVector, arities[i], outputVector, vectorLength);
            if (currentNumStates == -1) {
                break;
            }
        }
    } else {
        currentNumStates = normaliseArray(inputMatrix, normalisedVector, vectorLength);
        for (i = 0; i < vectorLength; i++) {
            outputVector[i] = normalisedVector[i];
        }
    }

    FREE_FUNC(normalisedVector);
    normalisedVector = NULL;

    return currentNumStates;
}/*mergeMultipleArraysArities(double *, double *, int, int)*/

