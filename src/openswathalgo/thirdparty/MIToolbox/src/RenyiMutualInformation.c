/*******************************************************************************
** RenyiMutualInformation.c
** Part of the mutual information toolbox
**
** Contains functions to calculate the Renyi mutual information of 
** two variables X and Y, I_\alpha(X;Y), using the Renyi alpha divergence and 
** the joint entropy difference
** 
** Author: Adam Pocock
** Created 26/3/2010
**
** Copyright 2010-2017 Adam Pocock, The University Of Manchester
** www.cs.manchester.ac.uk
**
** This file is part of MIToolbox, licensed under the 3-clause BSD license.
*******************************************************************************/

#include "MIToolbox/MIToolbox.h"
#include "MIToolbox/ArrayOperations.h"
#include "MIToolbox/CalculateProbability.h"
#include "MIToolbox/RenyiEntropy.h"
#include "MIToolbox/RenyiMutualInformation.h"

double renyiMI(JointProbabilityState state, double alpha) {
  int firstIndex,secondIndex;
  int i;
  double jointTemp, marginalTemp;
  double invAlpha = 1.0 - alpha;
  double mutualInformation = 0.0;

  /* standard MI is D_KL(p(x,y)||p(x)p(y))
  ** which expands to
  ** D_KL(p(x,y)||p(x)p(y)) = sum(p(x,y) * log(p(x,y)/(p(x)p(y))))
  **
  ** Renyi alpha divergence D_alpha(p(x,y)||p(x)p(y))
  ** expands to
  ** D_alpha(p(x,y)||p(x)p(y)) = 1/(alpha-1) * log(sum((p(x,y)^alpha)*((p(x)p(y))^(1-alpha))))
  */
  
  for (i = 0; i < state.numJointStates; i++) {
    firstIndex = i % state.numFirstStates;
    secondIndex = i / state.numFirstStates;
    
    if ((state.jointProbabilityVector[i] > 0) && (state.firstProbabilityVector[firstIndex] > 0) && (state.secondProbabilityVector[secondIndex] > 0)) {
      jointTemp = pow(state.jointProbabilityVector[i],alpha);
      marginalTemp = state.firstProbabilityVector[firstIndex] * state.secondProbabilityVector[secondIndex];
      marginalTemp = pow(marginalTemp,invAlpha);
      mutualInformation += (jointTemp * marginalTemp);
    }
  }

  mutualInformation = log(mutualInformation);
  mutualInformation /= log(LOG_BASE);  
  mutualInformation /= (alpha-1.0); 

  return mutualInformation;
}

double calcRenyiMIDivergence(double alpha, uint *dataVector, uint *targetVector, int vectorLength) {
  JointProbabilityState state = calculateJointProbability(dataVector,targetVector,vectorLength);
  double mutualInformation = renyiMI(state,alpha);
  
  freeJointProbabilityState(state);
  
  return mutualInformation;
}/*calcRenyiMIDivergence(double, uint *, uint *, int)*/

double discAndCalcRenyiMIDivergence(double alpha, double *dataVector, double *targetVector, int vectorLength) {
  JointProbabilityState state = discAndCalcJointProbability(dataVector,targetVector,vectorLength);
  double mutualInformation = renyiMI(state,alpha);
  
  freeJointProbabilityState(state);
  
  return mutualInformation;
}/*discAndCalcRenyiMIDivergence(double, double *, double *, int)*/

double calcRenyiMIJoint(double alpha, uint *dataVector, uint *targetVector, int vectorLength) {
  double hY = calcRenyiEntropy(alpha, targetVector, vectorLength);
  double hX = calcRenyiEntropy(alpha, dataVector, vectorLength);
  
  double hXY = calcJointRenyiEntropy(alpha, dataVector, targetVector, vectorLength);
  
  double answer = hX + hY - hXY;
  
  return answer;
}/*calcRenyiMIJoint(double, uint*, uint*, int)*/

double discAndCalcRenyiMIJoint(double alpha, double *dataVector, double *targetVector, int vectorLength) {
  double mi;
  uint *dataNormVector = (uint *) checkedCalloc(vectorLength, sizeof(uint));
  uint *targetNormVector = (uint *) checkedCalloc(vectorLength, sizeof(uint));

  normaliseArray(dataVector,dataNormVector,vectorLength);
  normaliseArray(targetVector,targetNormVector,vectorLength);
  
  mi = calcRenyiMIJoint(alpha,dataNormVector,targetNormVector,vectorLength);

  FREE_FUNC(dataNormVector);
  FREE_FUNC(targetNormVector);
  dataNormVector = NULL;
  targetNormVector = NULL;

  return mi;
}/*discAndCalcRenyiMIJoint(double, double*, double*, int)*/
