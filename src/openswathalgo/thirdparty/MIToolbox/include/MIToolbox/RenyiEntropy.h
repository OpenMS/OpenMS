/*******************************************************************************
** RenyiEntropy.h
** Part of the mutual information toolbox
**
** Contains functions to calculate the Renyi alpha entropy of a single variable 
** H_\alpha(X), the Renyi joint entropy of two variables H_\alpha(X,Y), and the 
** conditional Renyi entropy H_\alpha(X|Y)
** 
** Author: Adam Pocock
** Created 26/3/2010
**
**  Copyright 2010-2017 Adam Pocock, The University Of Manchester
**  www.cs.manchester.ac.uk
**
**  This file is part of MIToolbox, licensed under the 3-clause BSD license.
*******************************************************************************/

#ifndef __Renyi_Entropy_H
#define __Renyi_Entropy_H

#include "MIToolbox/MIToolbox.h"
#include "MIToolbox/CalculateProbability.h"

#ifdef __cplusplus
extern "C" {
#endif 

/*******************************************************************************
** calculateRenyiEntropy returns the Renyi entropy in log base LOG_BASE of dataVector
** H_{\alpha}(X), for \alpha != 1
**
** length(dataVector) == vectorLength otherwise there
** will be a segmentation fault
*******************************************************************************/
double calcRenyiEntropy(double alpha, uint *dataVector, int vectorLength);
double discAndCalcRenyiEntropy(double alpha, double *dataVector, int vectorLength);

/*******************************************************************************
** calculateJointRenyiEntropy returns the Renyi entropy in log base LOG_BASE of the 
** joint variable of firstVector and secondVector H_{\alpha}(XY), 
** for \alpha != 1
**
** length(firstVector) == length(secondVector) == vectorLength otherwise there
** will be a segmentation fault
*******************************************************************************/
double calcJointRenyiEntropy(double alpha, uint *firstVector, uint *secondVector, int vectorLength);
double discAndCalcJointRenyiEntropy(double alpha, double *firstVector, double *secondVector, int vectorLength);

/*******************************************************************************
** Inner functions which operate on state structs.
*******************************************************************************/
double renyiEntropy(ProbabilityState state, double alpha);
double jointRenyiEntropy(JointProbabilityState state, double alpha);

#ifdef __cplusplus
}
#endif

#endif

