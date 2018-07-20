/*******************************************************************************
** RenyiEntropy.c
** Part of the mutual information toolbox
**
** Contains functions to calculate the Renyi alpha entropy of a single variable 
** H_\alpha(X), the Renyi joint entropy of two variables H_\alpha(X,Y), and the 
** conditional Renyi entropy H_\alpha(X|Y)
** 
** Author: Adam Pocock
** Created 26/3/2010
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

double renyiEntropy(ProbabilityState state, double alpha) {
  double entropy = 0.0;
  double tempValue = 0.0;
  int i;
  
  /*H_\alpha(X) = 1/(1-alpha) * \log(\sum_x p(x)^alpha)*/
  for (i = 0; i < state.numStates; i++) {
    tempValue = state.probabilityVector[i];
    if (tempValue > 0) {
      entropy += pow(tempValue,alpha);
    }
  }
  
  entropy = log(entropy);
  entropy /= log(LOG_BASE);
  entropy /= (1.0-alpha);
  
  return entropy;
}

double calcRenyiEntropy(double alpha, uint *dataVector, int vectorLength) {
  ProbabilityState state = calculateProbability(dataVector,vectorLength);
  double h = renyiEntropy(state,alpha);
  
  freeProbabilityState(state);
  
  return h;
}/*calcRenyiEntropy(double,uint*,int)*/

double discAndCalcRenyiEntropy(double alpha, double *dataVector, int vectorLength) {
  ProbabilityState state = discAndCalcProbability(dataVector,vectorLength);
  double h = renyiEntropy(state,alpha);
  
  freeProbabilityState(state);
  
  return h;
}/*discAndCalcRenyiEntropy(double,double*,int)*/

double jointRenyiEntropy(JointProbabilityState state, double alpha) {
  double jointEntropy = 0.0;  
  double tempValue = 0.0;
  int i;
  
  /*H_\alpha(XY) = 1/(1-alpha) * log(2)(sum p(xy)^alpha)*/
  for (i = 0; i < state.numJointStates; i++) {
    tempValue = state.jointProbabilityVector[i];
    if (tempValue > 0) {
      jointEntropy += pow(tempValue,alpha);
    }
  }
  
  jointEntropy = log(jointEntropy);
  jointEntropy /= log(LOG_BASE);
  jointEntropy /= (1.0-alpha);

  return jointEntropy;
}

double calcJointRenyiEntropy(double alpha, uint *firstVector, uint *secondVector, int vectorLength) {
  JointProbabilityState state = calculateJointProbability(firstVector,secondVector,vectorLength);
  double h = jointRenyiEntropy(state,alpha);
  
  freeJointProbabilityState(state);
  
  return h;
}/*calcJointRenyiEntropy(double,uint*,uint*,int)*/

double discAndCalcJointRenyiEntropy(double alpha, double *firstVector, double *secondVector, int vectorLength) {
  JointProbabilityState state = discAndCalcJointProbability(firstVector,secondVector,vectorLength);
  double h = jointRenyiEntropy(state,alpha);
  
  freeJointProbabilityState(state);
  
  return h;
}/*discAndCalcJointRenyiEntropy(double,double*,double*,int)*/
