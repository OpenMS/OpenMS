/*******************************************************************************
** RenyiMutualInformation.h
** Part of the mutual information toolbox
**
** Contains functions to calculate the Renyi mutual information of 
** two variables X and Y, I_\alpha(X;Y), using the Renyi alpha divergence and 
** the joint entropy difference
** 
** Author: Adam Pocock
** Created 26/3/2010
**
**  Copyright 2010-2017 Adam Pocock, The University Of Manchester
**  www.cs.manchester.ac.uk
**
**  This file is part of MIToolbox, licensed under the 3-clause BSD license.
*******************************************************************************/

#ifndef __Renyi_MutualInformation_H
#define __Renyi_MutualInformation_H

#include "MIToolbox/MIToolbox.h"
#include "MIToolbox/CalculateProbability.h"

#ifdef __cplusplus
extern "C" {
#endif 

/*******************************************************************************
** calculateRenyiMIDivergence returns the log base LOG_BASE Renyi mutual information
** between dataVector and targetVector, I_{\alpha}(X;Y), for \alpha != 1
** This uses Renyi's generalised alpha divergence as the difference measure
** instead of the KL-divergence as in Shannon's Mutual Information
**
** length(dataVector) == length(targetVector) == vectorLength otherwise there
** will be a segmentation fault
*******************************************************************************/
double calcRenyiMIDivergence(double alpha, uint *dataVector, uint *targetVector, int vectorLength);
double discAndCalcRenyiMIDivergence(double alpha, double *dataVector, double *targetVector, int vectorLength);

/****************************************************************************** 
** This function returns a different value to the alpha divergence mutual 
** information, and thus is not a correct mutual information.
** It is maintained to show how different Renyi's Information Theory is.
******************************************************************************/
double calcRenyiMIJoint(double alpha, uint *dataVector, uint *targetVector, int vectorLength);
double discAndCalcRenyiMIJoint(double alpha, double *dataVector, double *targetVector, int vectorLength);

/*******************************************************************************
** Inner functions which operate on state structs.
*******************************************************************************/
double renyiMI(JointProbabilityState state, double alpha);

#ifdef __cplusplus
}
#endif 

#endif

