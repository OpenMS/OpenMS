/*******************************************************************************
** WeightedMutualInformation.c
** Part of the mutual information toolbox
**
** Contains functions to calculate the mutual information of 
** two variables X and Y, I(X;Y), to calculate the joint mutual information
** of two variables X & Z on the variable Y, I(XZ;Y), and the conditional
** mutual information I(X;Y|Z), while using a weight vector to modify the
** calculation.
** 
** Author: Adam Pocock
** Created: 20/06/2011
**
**  Copyright 2010-2017 Adam Pocock, The University Of Manchester
**  www.cs.manchester.ac.uk
**
**  This file is part of MIToolbox, licensed under the 3-clause BSD license.
*******************************************************************************/

#include "MIToolbox/MIToolbox.h"
#include "MIToolbox/ArrayOperations.h"
#include "MIToolbox/CalculateProbability.h"
#include "MIToolbox/WeightedEntropy.h"
#include "MIToolbox/WeightedMutualInformation.h"

double wmi(WeightedJointProbState state) {
  double mutualInformation = 0.0;
  int firstIndex,secondIndex;
  int i;
  
  /*
  ** I_w(X;Y) = \sum_x \sum_y w(x,y)p(x,y) * \log (p(x,y)/p(x)p(y))
  */
  for (i = 0; i < state.numJointStates; i++) {
    firstIndex = i % state.numFirstStates;
    secondIndex = i / state.numFirstStates;
    
    if ((state.jointProbabilityVector[i] > 0) && (state.firstProbabilityVector[firstIndex] > 0) && (state.secondProbabilityVector[secondIndex] > 0)) {
      mutualInformation += state.jointWeightVector[i] * state.jointProbabilityVector[i] * log(state.jointProbabilityVector[i] / state.firstProbabilityVector[firstIndex] / state.secondProbabilityVector[secondIndex]);
    }
  }
  
  mutualInformation /= log(LOG_BASE);

  return mutualInformation;
}

double calcWeightedMutualInformation(uint *dataVector, uint *targetVector, double *weightVector, int vectorLength) {
  WeightedJointProbState state = calculateWeightedJointProbability(dataVector,targetVector,weightVector,vectorLength);
  double mutualInformation = wmi(state);
    
  freeWeightedJointProbState(state);
  
  return mutualInformation;
}/*calcWeightedMutualInformation(uint *,uint *,double *,int)*/

double discAndCalcWeightedMutualInformation(double *dataVector, double *targetVector, double *weightVector, int vectorLength) {
  WeightedJointProbState state = discAndCalcWeightedJointProbability(dataVector,targetVector,weightVector,vectorLength);
  double mutualInformation = wmi(state);
    
  freeWeightedJointProbState(state);
  
  return mutualInformation;
}/*discAndCalcWeightedMutualInformation(double *,double *,double *,int)*/

double calcWeightedConditionalMutualInformation(uint *dataVector, uint *targetVector, uint *conditionVector, double *weightVector, int vectorLength) {
  double mutualInformation = 0.0;
  double firstCondition, secondCondition;
  uint *mergedVector = (uint *) checkedCalloc(vectorLength,sizeof(uint));
  
  mergeArrays(targetVector,conditionVector,mergedVector,vectorLength);
  
  /* I(X;Y|Z) = H(X|Z) - H(X|YZ) */
  /* double calculateWeightedConditionalEntropy(double *dataVector, double *conditionVector, double *weightVector, int vectorLength); */
  firstCondition = calcWeightedConditionalEntropy(dataVector,conditionVector,weightVector,vectorLength);
  secondCondition = calcWeightedConditionalEntropy(dataVector,mergedVector,weightVector,vectorLength);
  
  mutualInformation = firstCondition - secondCondition;
  
  FREE_FUNC(mergedVector);
  mergedVector = NULL;
  
  return mutualInformation;
}/*calcWeightedConditionalMutualInformation(double *,double *,double *,double *,int)*/

double discAndCalcWeightedConditionalMutualInformation(double *dataVector, double *targetVector, double *conditionVector, double *weightVector, int vectorLength) {
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
  
  /* I_w(X;Y|Z) = H_w(X|Z) - H_w(X|YZ) */
  /* double calculateWeightedConditionalEntropy(double *dataVector, double *conditionVector, double *weightVector, int vectorLength); */
  firstCondition = calcWeightedConditionalEntropy(dataNormVector,conditionNormVector,weightVector,vectorLength);
  secondCondition = calcWeightedConditionalEntropy(dataNormVector,mergedVector,weightVector,vectorLength);
  
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
}/*discAndCalcWeightedConditionalMutualInformation(double *,double *,double *,double *,int)*/
