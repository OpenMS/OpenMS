/*******************************************************************************
** WeightedEntropy.c
** Part of the mutual information toolbox
**
** Contains functions to calculate the entropy of a single variable H(X), 
** the joint entropy of two variables H(X,Y), and the conditional entropy
** H(X|Y), while using a weight vector to modify the calculation.
** 
** Author: Adam Pocock
** Created: 20/06/2011
**
** Copyright 2010-2017 Adam Pocock, The University Of Manchester
** www.cs.manchester.ac.uk
**
** This file is part of MIToolbox, licensed under the 3-clause BSD license.
*******************************************************************************/

#include "MIToolbox/MIToolbox.h"
#include "MIToolbox/CalculateProbability.h"
#include "MIToolbox/WeightedEntropy.h"

double wEntropy(WeightedProbState state) {
  double entropy = 0.0;
  double tempValue = 0.0;
  int i;
  
  /*H_w(X) = - \sum_x w(x)p(x) \log p(x)*/
  for (i = 0; i < state.numStates; i++) {
    tempValue = state.probabilityVector[i];
    if (tempValue > 0) {
      entropy -= state.stateWeightVector[i] * tempValue * log(tempValue);
    }
  }
  
  entropy /= log(LOG_BASE);

  return entropy;
}

double calcWeightedEntropy(uint *dataVector, double *weightVector, int vectorLength) {
  WeightedProbState state = calculateWeightedProbability(dataVector,weightVector,vectorLength);
  double entropy = wEntropy(state);
  
  freeWeightedProbState(state);
  
  return entropy;
}/*calcWeightedEntropy(uint *,double *,int)*/

double discAndCalcWeightedEntropy(double *dataVector, double *weightVector, int vectorLength) {
  WeightedProbState state = discAndCalcWeightedProbability(dataVector,weightVector,vectorLength);
  double entropy = wEntropy(state);
  
  freeWeightedProbState(state);
  
  return entropy;
}/*discAndCalcWeightedEntropy(double *,double *,int)*/

double wJointEntropy(WeightedJointProbState state) {
  double jointEntropy = 0.0;  
  double tempValue = 0.0;
  int i;

  /*H_w(X,Y) = - \sum_x \sum_y w(x,y) p(x,y) \log p(x,y)*/
  for (i = 0; i < state.numJointStates; i++) {
    tempValue = state.jointProbabilityVector[i];
    if (tempValue > 0) {
      jointEntropy -= state.jointWeightVector[i] * tempValue * log(tempValue);
    }
  }
  
  jointEntropy /= log(LOG_BASE);

  return jointEntropy;
}

double calcWeightedJointEntropy(uint *firstVector, uint *secondVector, double *weightVector, int vectorLength) {
  WeightedJointProbState state = calculateWeightedJointProbability(firstVector,secondVector,weightVector,vectorLength);
  double jointEntropy = wJointEntropy(state);
  
  freeWeightedJointProbState(state);
  
  return jointEntropy;
}/*calcWeightedJointEntropy(uint *,uint *,double *,int)*/

double discAndCalcWeightedJointEntropy(double *firstVector, double *secondVector, double *weightVector, int vectorLength) {
  WeightedJointProbState state = discAndCalcWeightedJointProbability(firstVector,secondVector,weightVector,vectorLength);
  double jointEntropy = wJointEntropy(state);
  
  freeWeightedJointProbState(state);
  
  return jointEntropy;
}/*discAndCalcWeightedJointEntropy(double *,double *,double *,int)*/

double wCondEntropy(WeightedJointProbState state) {
  double condEntropy = 0.0;  
  double jointValue = 0.0;
  double condValue = 0.0;
  int i;
  
  /*H_w(X|Y) = - \sum_x \sum_y w(x,y)p(x,y) \log (p(x,y)/p(y))*/
  /* to index by numFirstStates use modulus of i
  ** to index by numSecondStates use integer division of i by numFirstStates
  */
  for (i = 0; i < state.numJointStates; i++) {
    jointValue = state.jointProbabilityVector[i];
    condValue = state.secondProbabilityVector[i / state.numFirstStates];
    if ((jointValue > 0) && (condValue > 0)) {
      condEntropy -= state.jointWeightVector[i] * jointValue * log(jointValue / condValue);
    }
  }
  
  condEntropy /= log(LOG_BASE);
  
  return condEntropy;
}

double calcWeightedConditionalEntropy(uint *dataVector, uint *conditionVector, double *weightVector, int vectorLength) {
  WeightedJointProbState state = calculateWeightedJointProbability(dataVector,conditionVector,weightVector,vectorLength);
  double condEntropy = wCondEntropy(state);
  
  freeWeightedJointProbState(state);
  
  return condEntropy;
}/*calcWeightedConditionalEntropy(uint *,uint *,double *,int)*/

double discAndCalcWeightedConditionalEntropy(double *dataVector, double *conditionVector, double *weightVector, int vectorLength) {
  WeightedJointProbState state = discAndCalcWeightedJointProbability(dataVector,conditionVector,weightVector,vectorLength);
  double condEntropy = wCondEntropy(state);
  
  freeWeightedJointProbState(state);
  
  return condEntropy;
}/*discAndCalcWeightedConditionalEntropy(double *,double *,double *,int)*/
