/*******************************************************************************
** MutualInformation.c
** Part of the mutual information toolbox
**
** Contains functions to calculate the mutual information of 
** two variables X and Y, I(X;Y), to calculate the joint mutual information
** of two variables X & Z on the variable Y, I(XZ;Y), and the conditional
** mutual information I(x;Y|Z)
** 
** Author: Adam Pocock
** Created 19/2/2010
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
#include "MIToolbox/Entropy.h"
#include "MIToolbox/MutualInformation.h"

double mi(JointProbabilityState state) {
  double mutualInformation = 0.0;
  int firstIndex,secondIndex;
  int i;
    
  /*
  ** I(X;Y) = \sum_x \sum_y p(x,y) * \log (p(x,y)/p(x)p(y))
  */
  for (i = 0; i < state.numJointStates; i++) {
    firstIndex = i % state.numFirstStates;
    secondIndex = i / state.numFirstStates;
    
    if ((state.jointProbabilityVector[i] > 0) && (state.firstProbabilityVector[firstIndex] > 0) && (state.secondProbabilityVector[secondIndex] > 0)) {
      /*double division is probably more stable than multiplying two small numbers together
      ** mutualInformation += state.jointProbabilityVector[i] * log(state.jointProbabilityVector[i] / (state.firstProbabilityVector[firstIndex] * state.secondProbabilityVector[secondIndex]));
      */
      mutualInformation += state.jointProbabilityVector[i] * log(state.jointProbabilityVector[i] / state.firstProbabilityVector[firstIndex] / state.secondProbabilityVector[secondIndex]);
    }
  }
  
  mutualInformation /= log(LOG_BASE);
 
  return mutualInformation;
}/*mi(JointProbabilityState)*/

double calcMutualInformation(uint *dataVector, uint *targetVector, int vectorLength) {
  JointProbabilityState state = calculateJointProbability(dataVector,targetVector,vectorLength);
    
  double mutualInformation = mi(state);
  
  freeJointProbabilityState(state);
  
  return mutualInformation;
}/*calculateMutualInformation(uint *,uint *,int)*/

double discAndCalcMutualInformation(double *dataVector, double *targetVector, int vectorLength) {
  JointProbabilityState state = discAndCalcJointProbability(dataVector,targetVector,vectorLength);
    
  double mutualInformation = mi(state);
  
  freeJointProbabilityState(state);
  
  return mutualInformation;
}/*discAndCalcMutualInformation(double *,double *,int)*/

double calcConditionalMutualInformation(uint *dataVector, uint *targetVector, uint *conditionVector, int vectorLength) {
  double mutualInformation = 0.0;
  double firstCondition, secondCondition;
  uint *mergedVector = (uint *) checkedCalloc(vectorLength,sizeof(uint));
  
  mergeArrays(targetVector,conditionVector,mergedVector,vectorLength);
  
  /* I(X;Y|Z) = H(X|Z) - H(X|YZ) */
  /* double calculateConditionalEntropy(double *dataVector, double *conditionVector, int vectorLength); */
  firstCondition = calcConditionalEntropy(dataVector,conditionVector,vectorLength);
  secondCondition = calcConditionalEntropy(dataVector,mergedVector,vectorLength);
  
  mutualInformation = firstCondition - secondCondition;
  
  FREE_FUNC(mergedVector);
  mergedVector = NULL;
  
  return mutualInformation;
}/*calculateConditionalMutualInformation(double *,double *,double *,int)*/

double discAndCalcConditionalMutualInformation(double *dataVector, double *targetVector, double *conditionVector, int vectorLength) {
  double mutualInformation = 0.0;
  double firstCondition, secondCondition;
  uint *dataNormVector = (uint *) checkedCalloc(vectorLength,sizeof(uint));
  uint *targetNormVector = (uint *) checkedCalloc(vectorLength,sizeof(uint));
  uint *conditionNormVector = (uint *) checkedCalloc(vectorLength,sizeof(uint));
  uint *mergedVector = (uint *) checkedCalloc(vectorLength,sizeof(uint));
  
  normaliseArray(dataVector,dataNormVector,vectorLength);
  normaliseArray(targetVector,targetNormVector,vectorLength);
  normaliseArray(conditionVector,conditionNormVector,vectorLength);
  mergeArrays(targetNormVector,conditionNormVector,mergedVector,vectorLength);
  
  /* I(X;Y|Z) = H(X|Z) - H(X|YZ) */
  /* double calculateConditionalEntropy(double *dataVector, double *conditionVector, int vectorLength); */
  firstCondition = calcConditionalEntropy(dataNormVector,conditionNormVector,vectorLength);
  secondCondition = calcConditionalEntropy(dataNormVector,mergedVector,vectorLength);
  
  mutualInformation = firstCondition - secondCondition;
  
  FREE_FUNC(dataNormVector);
  FREE_FUNC(targetNormVector);
  FREE_FUNC(conditionNormVector);
  FREE_FUNC(mergedVector);
  dataNormVector = NULL;
  targetNormVector = NULL;
  conditionNormVector = NULL;
  mergedVector = NULL;
  
  return mutualInformation;
}/*calculateConditionalMutualInformation(double *,double *,double *,int)*/
