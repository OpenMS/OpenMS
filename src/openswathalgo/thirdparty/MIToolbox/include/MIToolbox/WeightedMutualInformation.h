/*******************************************************************************
** WeightedMutualInformation.h
** Part of the mutual information toolbox
**
** Contains functions to calculate the mutual information of 
** two variables X and Y, I(X;Y), to calculate the joint mutual information
** of two variables X & Z on the variable Y, I(XZ;Y), and the conditional
** mutual information I(x;Y|Z), while using a weight vector to modify
** the calculation.
** 
** Author: Adam Pocock
** Created: 20/6/2011
**
**  Copyright 2010-2017 Adam Pocock, The University Of Manchester
**  www.cs.manchester.ac.uk
**
**  This file is part of MIToolbox, licensed under the 3-clause BSD license.
*******************************************************************************/

#ifndef __WeightedMutualInformation_H
#define __WeightedMutualInformation_H

#include "MIToolbox/MIToolbox.h"
#include "MIToolbox/CalculateProbability.h"

#ifdef __cplusplus
extern "C" {
#endif 

/*******************************************************************************
** calculateWeightedMutualInformation returns the LOG_BASE mutual information 
** between dataVector and targetVector, I_w(X;Y), weighted by weightVector
**
** length(vectors) == vectorLength otherwise it will segmentation fault
*******************************************************************************/
double calcWeightedMutualInformation(uint *dataVector, uint *targetVector, double *weightVector, int vectorLength);
double discAndCalcWeightedMutualInformation(double *dataVector, double *targetVector, double *weightVector, int vectorLength);

/*******************************************************************************
** calculateWeightedConditionalMutualInformation returns the LOG_BASE 
** mutual information between dataVector and targetVector, conditioned on 
** conditionVector, I(X;Y|Z), weighted by weightVector
**
** length(vectors) == vectorLength otherwise it will segmentation fault
*******************************************************************************/
double calcWeightedConditionalMutualInformation(uint *dataVector, uint *targetVector, uint *conditionVector, double *weightVector, int vectorLength);
double discAndCalcWeightedConditionalMutualInformation(double *dataVector, double *targetVector, double *conditionVector, double *weightVector, int vectorLength);

/*******************************************************************************
** Inner functions which operate on state structs.
*******************************************************************************/
double wmi(WeightedJointProbState state);

#ifdef __cplusplus
}
#endif

#endif

