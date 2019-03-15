/*******************************************************************************
** ArrayOperations.h
** Part of the mutual information toolbox
**
** Contains functions to floor arrays, and to merge arrays into a joint
** state.
** 
** Author: Adam Pocock
** Created 17/2/2010
** Updated - 22/02/2014 - Added checking on calloc, and an increment array function.
**
**  Copyright 2010-2017 Adam Pocock, The University Of Manchester
**  www.cs.manchester.ac.uk
**
**  This file is part of MIToolbox, licensed under the 3-clause BSD license.
*******************************************************************************/

#ifndef __ArrayOperations_H
#define __ArrayOperations_H

#include "MIToolbox/MIToolbox.h"

#ifdef __cplusplus
extern "C" {
#endif 

/*******************************************************************************
** A version of calloc which checks to see if memory was allocated.
*******************************************************************************/
void* checkedCalloc(size_t vectorLength, size_t sizeOfType);

/*******************************************************************************
** Increments each value in a double array
*******************************************************************************/
void incrementVector(double* vector, int vectorLength);

/*******************************************************************************
** Simple print functions for debugging
*******************************************************************************/
void printDoubleVector(double *vector, int vectorlength);
void printIntVector(int *vector, int vectorLength);
void printUintVector(uint *vector, int vectorLength);

/*******************************************************************************
** Generates a 2D representation of a contiguous fortran style matrix
*******************************************************************************/
uint **generateIntIndices(uint *featureMatrix, uint noOfSamples, uint noOfFeatures);
double **generateDoubleIndices(double *featureMatrix, uint noOfSamples, uint noOfFeatures);

/*******************************************************************************
** Finds the maximum state of an int array.
*******************************************************************************/
int maxState(uint* vector, int vectorLength);

/*******************************************************************************
** numberOfUniqueValues finds the number of unique values in an array by 
** repeatedly iterating through the array and checking if a value has been 
** seen previously
*******************************************************************************/
int numberOfUniqueValues(double *featureVector, int vectorLength);

/******************************************************************************* 
** normaliseArray takes an input vector and writes an output vector
** which is a normalised version of the input, and returns the number of states
** A normalised array has min value = 0, max value = number of states
** and all values are integers
**
** length(inputVector) == length(outputVector) == vectorLength otherwise there
** is a memory leak
*******************************************************************************/
int normaliseArray(double *inputVector, uint *outputVector, int vectorLength);

/*******************************************************************************
** mergeArrays takes in two arrays and writes the joint state of those arrays
** to the output vector
**
** the length of the vectors must be the same and equal to vectorLength
*******************************************************************************/
int mergeArrays(uint *firstVector, uint *secondVector, uint *outputVector, int vectorLength);
int discAndMergeArrays(double *firstVector, double *secondVector, uint *outputVector, int vectorLength);
int mergeArraysArities(uint *firstVector, int numFirstStates, uint *secondVector, int numSecondStates, uint *outputVector, int vectorLength);
int discAndMergeArraysArities(double *firstVector, int numFirstStates, double *secondVector, int numSecondStates, uint *outputVector, int vectorLength);

/*******************************************************************************
** mergeMultipleArrays takes in a matrix and repeatedly merges the matrix using 
** merge arrays and writes the joint state of that matrix
** to the output vector
**
** the length of the vectors must be the same and equal to vectorLength
** matrixWidth = the number of columns in the matrix
*******************************************************************************/
int mergeMultipleArrays(double *inputMatrix, uint *outputVector, int matrixWidth, int vectorLength);
int mergeMultipleArraysArities(double *inputMatrix, uint *outputVector, int matrixWidth, int *arities, int vectorLength);

#ifdef __cplusplus
}
#endif

#endif

