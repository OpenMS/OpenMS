/*******************************************************************************
** Entropy.h
** Part of the mutual information toolbox
**
** Contains functions to calculate the entropy of a single variable H(X), 
** the joint entropy of two variables H(X,Y), and the conditional entropy
** H(X|Y)
** 
** Author: Adam Pocock
** Created 19/2/2010
**
**  Copyright 2010-2017 Adam Pocock, The University Of Manchester
**  www.cs.manchester.ac.uk
**
**  This file is part of MIToolbox, licensed under the 3-clause BSD license.
*******************************************************************************/

#ifndef __Entropy_H
#define __Entropy_H

#include "MIToolbox/MIToolbox.h"
#include "MIToolbox/CalculateProbability.h"

#ifdef __cplusplus
extern "C" {
#endif 

/*******************************************************************************
** calculateEntropy returns the entropy in log base LOG_BASE of dataVector
** H(X)
**
** length(dataVector) == vectorLength otherwise there
** will be a segmentation fault
*******************************************************************************/
double discAndCalcEntropy(double *dataVector, int vectorLength);
double calcEntropy(uint *dataVector, int vectorLength);

/*******************************************************************************
** calculateJointEntropy returns the entropy in log base LOG_BASE of the joint 
** variable of firstVector and secondVector H(X,Y)
**
** length(firstVector) == length(secondVector) == vectorLength otherwise there
** will be a segmentation fault
*******************************************************************************/
double discAndCalcJointEntropy(double *firstVector, double *secondVector, int vectorLength);
double calcJointEntropy(uint *firstVector, uint *secondVector, int vectorLength);

/*******************************************************************************
** calculateConditionalEntropy returns the entropy in log base LOG_BASE of dataVector
** conditioned on conditionVector, H(X|Y)
**
** length(dataVector) == length(conditionVector) == vectorLength otherwise there
** will be a segmentation fault
*******************************************************************************/
double discAndCalcConditionalEntropy(double *dataVector, double *conditionVector, int vectorLength);
double calcConditionalEntropy(uint *dataVector, uint *conditionVector, int vectorLength);

/*******************************************************************************
** Inner functions which operate on state structs.
*******************************************************************************/
double entropy(ProbabilityState state);
double jointEntropy(JointProbabilityState state);
double condEntropy(JointProbabilityState state);

#ifdef __cplusplus
}
#endif

#endif

