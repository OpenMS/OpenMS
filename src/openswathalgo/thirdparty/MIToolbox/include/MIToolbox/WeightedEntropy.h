/*******************************************************************************
** WeightedEntropy.h
** Part of the mutual information toolbox
**
** Contains functions to calculate the entropy of a single variable H(X), 
** the joint entropy of two variables H(X,Y), and the conditional entropy
** H(X|Y), while using a weight vector to modify the calculation.
** 
** Author: Adam Pocock
** Created: 20/6/2011 
**
**  Copyright 2010-2017 Adam Pocock, The University Of Manchester
**  www.cs.manchester.ac.uk
**
**  This file is part of MIToolbox, licensed under the 3-clause BSD license.
*******************************************************************************/

#ifndef __WeightedEntropy_H
#define __WeightedEntropy_H

#include "MIToolbox/MIToolbox.h"
#include "MIToolbox/CalculateProbability.h"

#ifdef __cplusplus
extern "C" {
#endif 

/*******************************************************************************
** calculateWeightedEntropy returns the entropy in log base LOG_BASE of dataVector
** H_w(X), weighted by weightVector
**
** length(vectors) == vectorLength otherwise it will segmentation fault
*******************************************************************************/
double calcWeightedEntropy(uint *dataVector, double *weightVector, int vectorLength);
double discAndCalcWeightedEntropy(double *dataVector, double *weightVector, int vectorLength);

/*******************************************************************************
** calculateWeightedJointEntropy returns the entropy in log base LOG_BASE of the joint 
** variable of firstVector and secondVector, H_w(XY), weighted by weightVector
**
** length(vectors) == vectorLength otherwise it will segmentation fault
*******************************************************************************/
double calcWeightedJointEntropy(uint *firstVector, uint *secondVector, double *weightVector, int vectorLength);
double discAndCalcWeightedJointEntropy(double *firstVector, double *secondVector, double *weightVector, int vectorLength);

/*******************************************************************************
** calculateWeightedConditionalEntropy returns the entropy in log base LOG_BASE of 
** dataVector conditioned on conditionVector, H_w(X|Y), weighted by weightVector
**
** length(vectors) == vectorLength otherwise it will segmentation fault
*******************************************************************************/
double calcWeightedConditionalEntropy(uint *dataVector, uint *conditionVector, double *weightVector, int vectorLength);
double discAndCalcWeightedConditionalEntropy(double *dataVector, double *conditionVector, double *weightVector, int vectorLength);

/*******************************************************************************
** Inner functions which operate on state structs.
*******************************************************************************/
double wEntropy(WeightedProbState state);
double wJointEntropy(WeightedJointProbState state);
double wCondEntropy(WeightedJointProbState state);

#ifdef __cplusplus
}
#endif

#endif

