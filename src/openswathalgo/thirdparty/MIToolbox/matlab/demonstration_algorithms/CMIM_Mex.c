/*******************************************************************************
** Demonstration feature selection algorithm - MATLAB r2009a
**
** Initial Version - 13/06/2008
** Updated - 07/07/2010
** based on CMIM.m
**
** Conditional Mutual Information Maximisation
** in
** "Fast Binary Feature Selection using Conditional Mutual Information Maximisation
** F. Fleuret (2004)
**
** Author - Adam Pocock
** Demonstration code for MIToolbox
*******************************************************************************/

#include "mex.h"
#include "MIToolbox/MIToolbox.h"
#include "MIToolbox/MutualInformation.h"
  
void CMIMCalculation(int k, int noOfSamples, int noOfFeatures,double *featureMatrix, double *classColumn, double *outputFeatures)
{
  /*holds the class MI values
  **the class MI doubles as the partial score from the CMIM paper
  */
  double *classMI = (double *)mxCalloc(noOfFeatures,sizeof(double));
  /*in the CMIM paper, m = lastUsedFeature*/
  int *lastUsedFeature = (int *)mxCalloc(noOfFeatures,sizeof(int));
  
  double score, conditionalInfo;
  int iMinus, currentFeature;
  
  double maxMI = 0.0;
  int maxMICounter = -1;
  
  int j,i;

  double **feature2D = (double**) mxCalloc(noOfFeatures,sizeof(double*));

  for(j = 0; j < noOfFeatures; j++)
  {
    feature2D[j] = featureMatrix + (int)j*noOfSamples;
  }
  
  for (i = 0; i < noOfFeatures;i++)
  {
    classMI[i] = discAndCalcMutualInformation(feature2D[i], classColumn, noOfSamples);
    
    if (classMI[i] > maxMI)
    {
      maxMI = classMI[i];
      maxMICounter = i;
    }/*if bigger than current maximum*/
  }/*for noOfFeatures - filling classMI*/
  
  outputFeatures[0] = maxMICounter;
  
  /*****************************************************************************
  ** We have populated the classMI array, and selected the highest
  ** MI feature as the first output feature
  ** Now we move into the CMIM algorithm
  *****************************************************************************/
  
  for (i = 1; i < k; i++)
  {
    score = 0.0;
    iMinus = i-1;
    
    for (j = 0; j < noOfFeatures; j++)
    {
      while ((classMI[j] > score) && (lastUsedFeature[j] < i))
      {
        /*double discAndCalcConditionalMutualInformation(double *firstVector, double *targetVector, double *conditionVector, int vectorLength);*/
        currentFeature = (int) outputFeatures[lastUsedFeature[j]];
        conditionalInfo = discAndCalcConditionalMutualInformation(feature2D[j],classColumn,feature2D[currentFeature],noOfSamples);
        if (classMI[j] > conditionalInfo)
        {
          classMI[j] = conditionalInfo;
        }/*reset classMI*/
        /*moved due to C indexing from 0 rather than 1*/
        lastUsedFeature[j] += 1;
      }/*while partial score greater than score & not reached last feature*/
      if (classMI[j] > score)
      {
        score = classMI[j];
        outputFeatures[i] = j;
      }/*if partial score still greater than score*/
	  }/*for number of features*/
  }/*for the number of features to select*/
  
  
  for (i = 0; i < k; i++)
  {
    outputFeatures[i] += 1; /*C indexes from 0 not 1*/
  }/*for number of selected features*/
  
}/*CMIMCalculation*/

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
    
    /*void CMIMCalculation(int k, int noOfSamples, int noOfFeatures,double *featureMatrix, double *classColumn, double *outputFeatures)*/
    CMIMCalculation(k,numberOfSamples,numberOfFeatures,featureMatrix,targets,output);
  }  
  
  return;
}/*mexFunction()*/
