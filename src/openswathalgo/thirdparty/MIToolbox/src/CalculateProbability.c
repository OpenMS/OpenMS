/*******************************************************************************
** CalculateProbability.c
** Part of the mutual information toolbox
**
** Contains functions to calculate the probability of each state in the array
** and to calculate the probability of the joint state of two arrays
** 
** Author: Adam Pocock
** Created: 17/02/2010
** Modified - 04/07/2011 - added weighted probability functions
** Updated - 22/02/2014 - Added checking on calloc.
**
** Copyright 2010-2017 Adam Pocock, The University Of Manchester
** www.cs.manchester.ac.uk
**
** This file is part of MIToolbox, licensed under the 3-clause BSD license.
*******************************************************************************/

#include "MIToolbox/MIToolbox.h"
#include "MIToolbox/ArrayOperations.h"
#include "MIToolbox/CalculateProbability.h"

JointProbabilityState calculateJointProbability(uint *firstVector, uint *secondVector, int vectorLength) {
  int *firstStateCounts;
  int *secondStateCounts;
  int *jointStateCounts;
  double *firstStateProbs;
  double *secondStateProbs;
  double *jointStateProbs;
  int firstNumStates;
  int secondNumStates;
  int jointNumStates;
  int i;
  double length = vectorLength;
  JointProbabilityState state;

  firstNumStates = maxState(firstVector,vectorLength);
  secondNumStates = maxState(secondVector,vectorLength);
  jointNumStates = firstNumStates * secondNumStates;
  
  firstStateCounts = (int *) checkedCalloc(firstNumStates,sizeof(int));
  secondStateCounts = (int *) checkedCalloc(secondNumStates,sizeof(int));
  jointStateCounts = (int *) checkedCalloc(jointNumStates,sizeof(int));
  
  firstStateProbs = (double *) checkedCalloc(firstNumStates,sizeof(double));
  secondStateProbs = (double *) checkedCalloc(secondNumStates,sizeof(double));
  jointStateProbs = (double *) checkedCalloc(jointNumStates,sizeof(double));
    
  /* Optimised for number of FP operations now O(states) instead of O(vectorLength) */
  for (i = 0; i < vectorLength; i++) {
    firstStateCounts[firstVector[i]] += 1;
    secondStateCounts[secondVector[i]] += 1;
    jointStateCounts[secondVector[i] * firstNumStates + firstVector[i]] += 1;
  }
  
  for (i = 0; i < firstNumStates; i++) {
    firstStateProbs[i] = firstStateCounts[i] / length;
  }
  
  for (i = 0; i < secondNumStates; i++) {
    secondStateProbs[i] = secondStateCounts[i] / length;
  }
  
  for (i = 0; i < jointNumStates; i++) {
    jointStateProbs[i] = jointStateCounts[i] / length;
  }

  FREE_FUNC(firstStateCounts);
  FREE_FUNC(secondStateCounts);
  FREE_FUNC(jointStateCounts);

  firstStateCounts = NULL;
  secondStateCounts = NULL;
  jointStateCounts = NULL;
  
  /*
  **typedef struct 
  **{
  **  double *jointProbabilityVector;
  **  int numJointStates;
  **  double *firstProbabilityVector;
  **  int numFirstStates;
  **  double *secondProbabilityVector;
  **  int numSecondStates;
  **} JointProbabilityState;
  */
  
  state.jointProbabilityVector = jointStateProbs;
  state.numJointStates = jointNumStates;
  state.firstProbabilityVector = firstStateProbs;
  state.numFirstStates = firstNumStates;
  state.secondProbabilityVector = secondStateProbs;
  state.numSecondStates = secondNumStates;

  return state;
}/*calcJointProbability(uint *,uint *, int)*/

JointProbabilityState discAndCalcJointProbability(double *firstVector, double *secondVector, int vectorLength) {
  uint *firstNormalisedVector = (uint *) checkedCalloc(vectorLength,sizeof(uint));
  uint *secondNormalisedVector = (uint *) checkedCalloc(vectorLength,sizeof(uint));
  JointProbabilityState state;
  
  normaliseArray(firstVector,firstNormalisedVector,vectorLength);
  normaliseArray(secondVector,secondNormalisedVector,vectorLength);

  state = calculateJointProbability(firstNormalisedVector,secondNormalisedVector,vectorLength);

  FREE_FUNC(firstNormalisedVector);
  FREE_FUNC(secondNormalisedVector);
    
  firstNormalisedVector = NULL;
  secondNormalisedVector = NULL;

  return state;
}/*discAndCalcJointProbability(double *,double *, int)*/

WeightedJointProbState calculateWeightedJointProbability(uint *firstVector, uint *secondVector, double *weightVector, int vectorLength) {
  int *firstStateCounts;
  int *secondStateCounts;
  int *jointStateCounts;
  double *firstStateProbs;
  double *secondStateProbs;
  double *jointStateProbs;
  double *firstWeightVec;
  double *secondWeightVec;
  double *jointWeightVec;
  int firstNumStates;
  int secondNumStates;
  int jointNumStates;
  int i;
  double length = vectorLength;
  WeightedJointProbState state;

  firstNumStates = maxState(firstVector,vectorLength);
  secondNumStates = maxState(secondVector,vectorLength);
  jointNumStates = firstNumStates * secondNumStates;
  
  firstStateCounts = (int *) checkedCalloc(firstNumStates,sizeof(int));
  secondStateCounts = (int *) checkedCalloc(secondNumStates,sizeof(int));
  jointStateCounts = (int *) checkedCalloc(jointNumStates,sizeof(int));
  
  firstStateProbs = (double *) checkedCalloc(firstNumStates,sizeof(double));
  secondStateProbs = (double *) checkedCalloc(secondNumStates,sizeof(double));
  jointStateProbs = (double *) checkedCalloc(jointNumStates,sizeof(double));
    
  firstWeightVec = (double *) checkedCalloc(firstNumStates,sizeof(double));
  secondWeightVec = (double *) checkedCalloc(secondNumStates,sizeof(double));
  jointWeightVec = (double *) checkedCalloc(jointNumStates,sizeof(double));
    
  for (i = 0; i < vectorLength; i++) {
    firstStateCounts[firstVector[i]] += 1;
    secondStateCounts[secondVector[i]] += 1;
    jointStateCounts[secondVector[i] * firstNumStates + firstVector[i]] += 1;

    firstWeightVec[firstVector[i]] += weightVector[i];
    secondWeightVec[secondVector[i]] += weightVector[i];
    jointWeightVec[secondVector[i] * firstNumStates + firstVector[i]] += weightVector[i];
  }
  
  for (i = 0; i < firstNumStates; i++) {
    if (firstStateCounts[i]) {
      firstStateProbs[i] = firstStateCounts[i] / length;
      firstWeightVec[i] /= firstStateCounts[i];
    }
  }
  
  for (i = 0; i < secondNumStates; i++) {
    if (secondStateCounts[i]) {
      secondStateProbs[i] = secondStateCounts[i] / length;
      secondWeightVec[i] /= secondStateCounts[i];
    }
  }
  
  for (i = 0; i < jointNumStates; i++) {
    if (jointStateCounts[i]) {
      jointStateProbs[i] = jointStateCounts[i] / length;
      jointWeightVec[i] /= jointStateCounts[i];
    }
  }

  FREE_FUNC(firstStateCounts);
  FREE_FUNC(secondStateCounts);
  FREE_FUNC(jointStateCounts);
    
  firstStateCounts = NULL;
  secondStateCounts = NULL;
  jointStateCounts = NULL;
  
  /*
  **typedef struct 
  **{
  **  double *jointProbabilityVector;
  **  double *jointWeightVector;
  **  int numJointStates;
  **  double *firstProbabilityVector;
  **  double *firstWeightVector;
  **  int numFirstStates;
  **  double *secondProbabilityVector;
  **  double *secondWeightVector;
  **  int numSecondStates;
  **} WeightedJointProbState;
  */
  
  state.jointProbabilityVector = jointStateProbs;
  state.jointWeightVector = jointWeightVec;
  state.numJointStates = jointNumStates;
  state.firstProbabilityVector = firstStateProbs;
  state.firstWeightVector = firstWeightVec;
  state.numFirstStates = firstNumStates;
  state.secondProbabilityVector = secondStateProbs;
  state.secondWeightVector = secondWeightVec;
  state.numSecondStates = secondNumStates;

  return state;
}

WeightedJointProbState discAndCalcWeightedJointProbability(double *firstVector, double *secondVector, double *weightVector, int vectorLength) {
  uint *firstNormalisedVector = (uint *) checkedCalloc(vectorLength,sizeof(uint));
  uint *secondNormalisedVector = (uint *) checkedCalloc(vectorLength,sizeof(uint));
  WeightedJointProbState state;
  
  normaliseArray(firstVector,firstNormalisedVector,vectorLength);
  normaliseArray(secondVector,secondNormalisedVector,vectorLength);
  
  state = calculateWeightedJointProbability(firstNormalisedVector,secondNormalisedVector,weightVector,vectorLength);
  FREE_FUNC(firstNormalisedVector);
  FREE_FUNC(secondNormalisedVector);
    
  firstNormalisedVector = NULL;
  secondNormalisedVector = NULL;

  return state;
}/*discAndCalcWeightedJointProbability(double *,double *,double *,int)*/

ProbabilityState calculateProbability(uint* dataVector, int vectorLength) {
  int numStates;
  int *stateCounts;
  double *stateProbs;
  ProbabilityState state;
  int i;
  double length = vectorLength;

  numStates = maxState(dataVector,vectorLength);
  
  stateCounts = (int *) checkedCalloc(numStates,sizeof(int));
  stateProbs = (double *) checkedCalloc(numStates,sizeof(double));
  
  /* Optimised for number of FP operations now O(states) instead of O(vectorLength) */
  for (i = 0; i < vectorLength; i++) {
    stateCounts[dataVector[i]] += 1;
  }
  
  for (i = 0; i < numStates; i++) {
    stateProbs[i] = stateCounts[i] / length;
  }
  
  FREE_FUNC(stateCounts);
  stateCounts = NULL;
  
  state.probabilityVector = stateProbs;
  state.numStates = numStates;

  return state;
}

ProbabilityState discAndCalcProbability(double *dataVector, int vectorLength) {
  ProbabilityState state;
  uint *normalisedVector = (uint *) checkedCalloc(vectorLength,sizeof(uint));
  
  normaliseArray(dataVector,normalisedVector,vectorLength);
  
  state = calculateProbability(normalisedVector,vectorLength);

  FREE_FUNC(normalisedVector);
  normalisedVector = NULL;

  return state;
}/*discAndCalcProbability(double *,int)*/


WeightedProbState calculateWeightedProbability(uint *dataVector, double *weightVector, int vectorLength) {
  int *stateCounts;
  double *stateProbs;
  double *stateWeights;
  int numStates, i;
  WeightedProbState state;
  double length = vectorLength;
  
  numStates = maxState(dataVector,vectorLength);
  
  stateCounts = (int *) checkedCalloc(numStates,sizeof(int));
  stateProbs = (double *) checkedCalloc(numStates,sizeof(double));
  stateWeights = (double *) checkedCalloc(numStates,sizeof(double));

  for (i = 0; i < vectorLength; i++) {
    stateCounts[dataVector[i]] += 1;
    stateWeights[dataVector[i]] += weightVector[i];
  }
  
  for (i = 0; i < numStates; i++) {
    stateProbs[i] = stateCounts[i] / length;
    stateWeights[i] /= stateCounts[i];
  }
  
  FREE_FUNC(stateCounts);
  stateCounts = NULL;
  
  state.probabilityVector = stateProbs;
  state.stateWeightVector = stateWeights;
  state.numStates = numStates;

  return state;
}/*calculateWeightedProbability(int *, double *, int)*/

WeightedProbState discAndCalcWeightedProbability(double *dataVector, double *weightVector, int vectorLength) {
  WeightedProbState state;
  uint *normalisedVector = (uint *) checkedCalloc(vectorLength,sizeof(uint));
  normaliseArray(dataVector,normalisedVector,vectorLength);
  
  state = calculateWeightedProbability(normalisedVector,weightVector,vectorLength);
  
  FREE_FUNC(normalisedVector);
  normalisedVector = NULL;

  return state;
}/*discAndCalcWeightedProbability(double *, double *, int)*/

void freeProbabilityState(ProbabilityState state) {
    FREE_FUNC(state.probabilityVector);
    state.probabilityVector = NULL;
}

void freeJointProbabilityState(JointProbabilityState state) {
    FREE_FUNC(state.firstProbabilityVector);
    state.firstProbabilityVector = NULL;
    FREE_FUNC(state.secondProbabilityVector);
    state.secondProbabilityVector = NULL;
    FREE_FUNC(state.jointProbabilityVector);
    state.jointProbabilityVector = NULL;
}

void freeWeightedProbState(WeightedProbState state) {
  FREE_FUNC(state.probabilityVector);
  state.probabilityVector = NULL;
  FREE_FUNC(state.stateWeightVector);
  state.stateWeightVector = NULL;
}

void freeWeightedJointProbState(WeightedJointProbState state) {
  FREE_FUNC(state.firstProbabilityVector);
  state.firstProbabilityVector = NULL;
  FREE_FUNC(state.secondProbabilityVector);
  state.secondProbabilityVector = NULL;
  FREE_FUNC(state.jointProbabilityVector);
  state.jointProbabilityVector = NULL;
  FREE_FUNC(state.firstWeightVector);
  state.firstWeightVector = NULL;
  FREE_FUNC(state.secondWeightVector);
  state.secondWeightVector = NULL;
  FREE_FUNC(state.jointWeightVector);
  state.jointWeightVector = NULL;
}
