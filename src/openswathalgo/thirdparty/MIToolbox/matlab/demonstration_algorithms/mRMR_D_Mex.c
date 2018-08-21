/*******************************************************************************
** Demonstration feature selection algorithm - MATLAB r2009a
**
** Initial Version - 13/06/2008
** Updated - 07/07/2010
** based on mRMR_D.m
**
** Minimum Relevance Maximum Redundancy
** in
** "Feature Selection Based on Mutual Information: Criteria of Max-Dependency, Max-Relevance, and Min-Redundancy"
** H. Peng et al. (2005)
**
** Author - Adam Pocock
** Demonstration code for MIToolbox
*******************************************************************************/

#include "mex.h"
#include "MIToolbox/MIToolbox.h"
#include "MIToolbox/MutualInformation.h"

void mRMRCalculation(int k, int noOfSamples, int noOfFeatures,double *featureMatrix, double *classColumn, double *outputFeatures)
{
  double **feature2D = (double**) mxCalloc(noOfFeatures,sizeof(double*));
  /*holds the class MI values*/
  double *classMI = (double *) mxCalloc(noOfFeatures,sizeof(double));
  int *selectedFeatures = (int *) mxCalloc(noOfFeatures,sizeof(int));
  /*holds the intra feature MI values*/
  int sizeOfMatrix = k*noOfFeatures;
  double *featureMIMatrix = (double *)mxCalloc(sizeOfMatrix,sizeof(double));
  
  double maxMI = 0.0;
  int maxMICounter = -1;
  
  /*init variables*/
  
  double score, currentScore, totalFeatureMI;
  int currentHighestFeature;
  
  int arrayPosition, i, j, x;
  
  for(j = 0; j < noOfFeatures; j++)
  {
    feature2D[j] = featureMatrix + (int)j*noOfSamples;
  }
  
  for (i = 0; i < sizeOfMatrix;i++)
  {
    featureMIMatrix[i] = -1;
  }/*for featureMIMatrix - blank to -1*/
  

  for (i = 0; i < noOfFeatures;i++)
  {
    classMI[i] = discAndCalcMutualInformation(feature2D[i], classColumn, noOfSamples);
    if (classMI[i] > maxMI)
    {
      maxMI = classMI[i];
      maxMICounter = i;
    }/*if bigger than current maximum*/
  }/*for noOfFeatures - filling classMI*/
  
  selectedFeatures[maxMICounter] = 1;
  outputFeatures[0] = maxMICounter;
  
  /*************
  ** Now we have populated the classMI array, and selected the highest
  ** MI feature as the first output feature
  ** Now we move into the mRMR-D algorithm
  *************/
  
  for (i = 1; i < k; i++)
  {
    /*to ensure it selects some features
    **if this is zero then it will not pick features where the redundancy is greater than the 
    **relevance
    */
    score = -1000.0;
    currentHighestFeature = 0;
    currentScore = 0.0;
    totalFeatureMI = 0.0;
    
    for (j = 0; j < noOfFeatures; j++)
    {
      /*if we haven't selected j*/
      if (selectedFeatures[j] == 0)
      {
        currentScore = classMI[j];
        totalFeatureMI = 0.0;
        
        for (x = 0; x < i; x++)
        {
          arrayPosition = x*noOfFeatures + j;
          if (featureMIMatrix[arrayPosition] == -1)
          {
            /*work out intra MI*/
            
            /*double discAndCalcMutualInformation(double *firstVector, double *secondVector, int vectorLength);*/
            featureMIMatrix[arrayPosition] = discAndCalcMutualInformation(feature2D[(int) outputFeatures[x]], feature2D[j], noOfSamples);
          }
          
          totalFeatureMI += featureMIMatrix[arrayPosition];
        }/*for the number of already selected features*/
        
        currentScore -= (totalFeatureMI/i);
			  if (currentScore > score)
			  {
				  score = currentScore;
				  currentHighestFeature = j;
			  }
			}/*if j is unselected*/
	  }/*for number of features*/
  
    selectedFeatures[currentHighestFeature] = 1;
    outputFeatures[i] = currentHighestFeature;
  
  }/*for the number of features to select*/
  
  for (i = 0; i < k; i++)
  {
    outputFeatures[i] += 1; /*C indexes from 0 not 1*/
  }/*for number of selected features*/
  
}/*mRMRCalculation(double[][],double[])*/

/*entry point for the mex call
**nlhs - number of outputs
**plhs - pointer to array of outputs
**nrhs - number of inputs
**prhs - pointer to array of inputs
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*************************************************************
  ** this function takes 3 arguments:
  ** k = number of features to select,
  ** featureMatrix[][] = matrix of features
  ** classColumn[] = targets
  ** the arguments should all be discrete integers.
  ** and has one output:
  ** selectedFeatures[] of size k
  *************************************************************/
  
  int k, numberOfFeatures, numberOfSamples, numberOfTargets;
  double *featureMatrix, *targets, *output;
  
  
  if (nlhs != 1)
  {
    printf("Incorrect number of output arguments");
  }/*if not 1 output*/
  if (nrhs != 3)
  {
    printf("Incorrect number of input arguments");
  }/*if not 3 inputs*/
  
  /*get the number of features to select, cast out as it is a double*/
  k = (int) mxGetScalar(prhs[0]);

  numberOfFeatures = mxGetN(prhs[1]);
  numberOfSamples = mxGetM(prhs[1]);
  
  numberOfTargets = mxGetM(prhs[2]);
  
  if (numberOfTargets != numberOfSamples)
  {
    printf("Number of targets must match number of samples\n");
    printf("Number of targets = %d, Number of Samples = %d, Number of Features = %d\n",numberOfTargets,numberOfSamples,numberOfFeatures);
    
    plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
  }/*if size mismatch*/
  else
  {
    
    featureMatrix = mxGetPr(prhs[1]);
    targets = mxGetPr(prhs[2]);
    
    plhs[0] = mxCreateDoubleMatrix(k,1,mxREAL);
    output = (double *)mxGetPr(plhs[0]);
    
    /*void mRMRCalculation(int k, int noOfSamples, int noOfFeatures,double *featureMatrix, double *classColumn, double *outputFeatures)*/
    mRMRCalculation(k,numberOfSamples,numberOfFeatures,featureMatrix,targets,output);
  }
  
  return;
}/*mexFunction()*/
