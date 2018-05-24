/*******************************************************************************
** CalculateProbability.h
** Part of the mutual information toolbox
**
** Contains functions to calculate the probability of each state in the array
** and to calculate the probability of the joint state of two arrays
** 
** Author: Adam Pocock
** Created 17/2/2010
**
**  Copyright 2010-2017 Adam Pocock, The University Of Manchester
**  www.cs.manchester.ac.uk
**
**  This file is part of MIToolbox, licensed under the 3-clause BSD license.
*******************************************************************************/

#ifndef __CalculateProbability_H
#define __CalculateProbability_H

#include "MIToolbox/MIToolbox.h"

#ifdef __cplusplus
extern "C" {
#endif 

typedef struct jpState
{
  double *jointProbabilityVector;
  int numJointStates;
  double *firstProbabilityVector;
  int numFirstStates;
  double *secondProbabilityVector;
  int numSecondStates;
} JointProbabilityState;

typedef struct pState
{
  double *probabilityVector;
  int numStates;
} ProbabilityState;

typedef struct wjpState
{
  double *jointProbabilityVector;
  double *jointWeightVector;
  int numJointStates;
  double *firstProbabilityVector;
  double *firstWeightVector;
  int numFirstStates;
  double *secondProbabilityVector;
  double *secondWeightVector;
  int numSecondStates;
} WeightedJointProbState;

typedef struct wpState
{
  double *probabilityVector;
  double *stateWeightVector;
  int numStates;
} WeightedProbState;

/*******************************************************************************
** calculateJointProbability returns the joint probability vector of two vectors
** and the marginal probability vectors in a struct.
** It is the base operation for all information theory calculations involving 
** two or more variables.
**
** length(firstVector) == length(secondVector) == vectorLength
** otherwise it will crash
*******************************************************************************/
JointProbabilityState calculateJointProbability(uint *firstVector, uint *secondVector, int vectorLength);

/*******************************************************************************
** discAndCalcJointProbability discretises the double vectors into int vectors,
** then generates a JointProbabilityState.
**
** length(firstVector) == length(secondVector) == vectorLength otherwise there
** will be a segmentation fault
*******************************************************************************/
JointProbabilityState discAndCalcJointProbability(double *firstVector, double *secondVector, int vectorLength);

/*******************************************************************************
** calculateWeightedJointProbability returns the same joint state object
** as calculateJointProbability, with additional weightVectors giving the 
** weight assigned to each state.
*******************************************************************************/
WeightedJointProbState calculateWeightedJointProbability(uint *firstVector, uint *secondVector, double *weightVector, int vectorLength);

/*******************************************************************************
** discAndCalcWeightedJointProbability returns the same joint state object
** as discAndCalcJointProbability, with additional weightVectors giving the 
** weight assigned to each state.
*******************************************************************************/
WeightedJointProbState discAndCalcWeightedJointProbability(double *firstVector, double *secondVector, double *weightVector, int vectorLength);

/*******************************************************************************
** calculateProbability returns the probability vector from one vector.
** It is the base operation for all information theory calculations involving 
** one variable
**
** length(dataVector) == vectorLength otherwise there
** will be a segmentation fault
*******************************************************************************/
ProbabilityState calculateProbability(uint *dataVector, int vectorLength);

/*******************************************************************************
** discAndCalcProbability discretises the double vector into an int vector,
** then generates a ProbabilityState.
**
** length(dataVector) == vectorLength otherwise there
** will be a segmentation fault
*******************************************************************************/
ProbabilityState discAndCalcProbability(double *dataVector, int vectorLength);

/*******************************************************************************
** calculateWeightedProbability returns the same state object
** as calculateProbability, with an additional weightVector giving the 
** weight assigned to each state.
*******************************************************************************/
WeightedProbState calculateWeightedProbability(uint *dataVector, double *weightVector, int vectorLength);

/*******************************************************************************
** discAndCalcWeightedProbability returns the same state object
** as discAndCalcProbability, with an additional weightVector giving the 
** weight assigned to each state.
*******************************************************************************/
WeightedProbState discAndCalcWeightedProbability(double *dataVector, double *weightVector, int vectorLength);

/*******************************************************************************
** Frees the struct members and sets all pointers to NULL.
*******************************************************************************/
void freeProbabilityState(ProbabilityState state);
void freeJointProbabilityState(JointProbabilityState state);
void freeWeightedProbState(WeightedProbState state);
void freeWeightedJointProbState(WeightedJointProbState state);

#ifdef __cplusplus
}
#endif

#endif

