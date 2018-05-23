/*******************************************************************************
** MIToolboxMex.c
** is the MATLAB entry point for the MIToolbox functions when called from
** a MATLAB/OCTAVE script.
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
#include "MIToolbox/MutualInformation.h"

/*******************************************************************************
** entry point for the mex call
** nlhs - number of outputs
** plhs - pointer to array of outputs
** nrhs - number of inputs
** prhs - pointer to array of inputs
*******************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*****************************************************************************
  ** this function takes a flag and a variable number of arguments
  ** depending on the value of the flag and returns either a construct 
  ** containing probability estimates, a merged vector or a double value 
  ** representing an entropy or mutual information
  *****************************************************************************/
  
  int flag, i, numberOfSamples, checkSamples, thirdCheckSamples, numberOfFeatures, checkFeatures, thirdCheckFeatures;
  int  numArities, errorTest;
  double *dataVector, *condVector, *targetVector, *firstVector, *secondVector, *output, *numStates;
  double *matrix, *arities;
  uint *mergedVector;
  uint *outputIntVector, *intArities;
  
  double *jointOutput, *numJointStates, *firstOutput, *numFirstStates, *secondOutput, *numSecondStates;
  
  ProbabilityState state;
  JointProbabilityState jointState;
  
  /*if (nlhs != 1)
  {
    printf("Incorrect number of output arguments\n");
  }//if not 1 output
  */
  switch (nrhs)
  {
    case 2:
    {
        /*printf("Must be H(X), discAndCalcProbability(X), merge(X), normaliseArray(X)\n");*/
        break;
    }
    case 3:
    {
        /*printf("Must be H(XY), H(X|Y), discAndCalcJointProbability(XY), I(X;Y)\n");*/
        break;
    }
    case 4:
    {
        /*printf("Must be I(X;Y|Z)\n");*/
        break;
    }
    default:
    {
        printf("Incorrect number of arguments, format is MIToolbox(\"FLAG\",varargin)\n");
        break;
    }
  }
  
  /* number to function map
  ** 1 = discAndCalcProbability
  ** 2 = discAndCalcJointProbability
  ** 3 = mergeArrays
  ** 4 = H(X)
  ** 5 = H(XY)
  ** 6 = H(X|Y)
  ** 7 = I(X;Y)
  ** 8 = I(X;Y|Z)
  ** 9 = normaliseArray
  */
  
  flag = *mxGetPr(prhs[0]);
  
  switch (flag)
  {
    case 1:
    {
      /*
      **discAndCalcProbability
      */
      numberOfSamples = mxGetM(prhs[1]);
      dataVector = (double *) mxGetPr(prhs[1]);

      /*ProbabilityState discAndCalcProbability(double *dataVector, int vectorLength);*/
      state = discAndCalcProbability(dataVector,numberOfSamples);
      
      plhs[0] = mxCreateDoubleMatrix(state.numStates,1,mxREAL);
      plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
      output = (double *)mxGetPr(plhs[0]);
      numStates = (double *) mxGetPr(plhs[1]);
      
      *numStates = state.numStates;
      
      for (i = 0; i < state.numStates; i++)
      {
        output[i] = state.probabilityVector[i];
      }

      freeProbabilityState(state);
      
      break;
    }/*case 1 - discAndCalcProbability*/
    case 2:
    {
      /*
      **discAndCalcJointProbability
      */
      numberOfSamples = mxGetM(prhs[1]);
      firstVector = (double *) mxGetPr(prhs[1]);
      secondVector = (double *) mxGetPr(prhs[2]);

      /*JointProbabilityState discAndCalcJointProbability(double *firstVector, double *secondVector int vectorLength);*/
      jointState = discAndCalcJointProbability(firstVector,secondVector,numberOfSamples);
      
      plhs[0] = mxCreateDoubleMatrix(jointState.numJointStates,1,mxREAL);
      plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
      plhs[2] = mxCreateDoubleMatrix(jointState.numFirstStates,1,mxREAL);
      plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
      plhs[4] = mxCreateDoubleMatrix(jointState.numSecondStates,1,mxREAL);
      plhs[5] = mxCreateDoubleMatrix(1,1,mxREAL);
      jointOutput = (double *)mxGetPr(plhs[0]);
      numJointStates = (double *) mxGetPr(plhs[1]);
      firstOutput = (double *)mxGetPr(plhs[2]);
      numFirstStates = (double *) mxGetPr(plhs[3]);
      secondOutput = (double *)mxGetPr(plhs[4]);
      numSecondStates = (double *) mxGetPr(plhs[5]);
      
      *numJointStates = jointState.numJointStates;
      *numFirstStates = jointState.numFirstStates;
      *numSecondStates = jointState.numSecondStates;
      
      for (i = 0; i < jointState.numJointStates; i++)
      {
        jointOutput[i] = jointState.jointProbabilityVector[i];
      }
      for (i = 0; i < jointState.numFirstStates; i++)
      {
        firstOutput[i] = jointState.firstProbabilityVector[i];
      }
      for (i = 0; i < jointState.numSecondStates; i++)
      {
        secondOutput[i] = jointState.secondProbabilityVector[i];
      }

      freeJointProbabilityState(jointState);
      
      break;
    }/*case 2 - discAndCalcJointProbability */
    case 3:
    {
      /*
      **mergeArrays
      */
      numberOfSamples = mxGetM(prhs[1]);
      numberOfFeatures = mxGetN(prhs[1]);
            
      numArities = 0;
      if (nrhs > 2)
      {
        numArities = mxGetN(prhs[2]);
        /*printf("arities = %d, features = %d, samples = %d\n",numArities,numberOfFeatures,numberOfSamples);*/
      }
      
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
      
      if (numArities == 0)
      {
        /*
        **no arities therefore compress output
        */
        if ((numberOfFeatures > 0) && (numberOfSamples > 0))
        { 
          matrix = (double *) mxGetPr(prhs[1]);
          mergedVector = (uint *) mxCalloc(numberOfSamples,sizeof(uint));
            
          plhs[0] = mxCreateDoubleMatrix(numberOfSamples,1,mxREAL);
          output = (double *)mxGetPr(plhs[0]);
          
          /*int mergeMultipleArrays(double *inputMatrix, double *outputVector, int matrixWidth, int vectorLength)*/
          mergeMultipleArrays(matrix, mergedVector, numberOfFeatures, numberOfSamples);
          for (i = 0; i < numberOfSamples; i++)
          {
            output[i] = mergedVector[i];
          }
          
          mxFree(mergedVector);
          mergedVector = NULL;
        }
      }
      else if (numArities == numberOfFeatures)
      {
        if ((numberOfFeatures > 0) && (numberOfSamples > 0))
        { 
          
          matrix = (double *) mxGetPr(prhs[1]);
          mergedVector = (uint *) mxCalloc(numberOfSamples,sizeof(uint));
          
          arities = (double *) mxGetPr(prhs[2]);
          intArities = (uint *) mxCalloc(numberOfFeatures,sizeof(uint));
          for (i = 0; i < numArities; i++)
          {
            intArities[i] = (uint) floor(arities[i]);
          }
          
          /*int mergeMultipleArrays(double *inputMatrix, double *outputVector, int matrixWidth, int *arities, int vectorLength);*/
          errorTest = mergeMultipleArraysArities(matrix, mergedVector, numberOfFeatures, intArities, numberOfSamples);
           
          if (errorTest != -1)
          {
            plhs[0] = mxCreateDoubleMatrix(numberOfSamples,1,mxREAL);
            output = (double *)mxGetPr(plhs[0]);
            for (i = 0; i < numberOfSamples; i++)
            {
              output[i] = mergedVector[i];
            }
          }
          else
          {
            printf("Incorrect arities supplied. More states in data than specified\n");
          }
          
          mxFree(mergedVector);
          mergedVector = NULL;
        }
      }
      else
      {
        printf("Number of arities does not match number of features, arities should be a row vector\n");
      }
      
      break;
    }/*case 3 - mergeArrays*/
    case 4:
    {
      /*
      **H(X)
      */
      numberOfSamples = mxGetM(prhs[1]);
      numberOfFeatures = mxGetN(prhs[1]);
      
      dataVector = (double *) mxGetPr(prhs[1]);

      plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
      output = (double *)mxGetPr(plhs[0]);

      if (numberOfFeatures == 1)
      {
        /*double discAndCalcEntropy(double *dataVector, int vectorLength);*/
        *output = discAndCalcEntropy(dataVector,numberOfSamples);
      }
      else
      {
        printf("No columns in input\n");
        *output = -1.0;
      }
      
      break;
    }/*case 4 - H(X)*/
    case 5:
    {
      /*
      **H(XY)
      */
      numberOfSamples = mxGetM(prhs[1]);
      checkSamples = mxGetM(prhs[2]);
      
      numberOfFeatures = mxGetN(prhs[1]);
      checkFeatures = mxGetN(prhs[2]);

      firstVector = mxGetPr(prhs[1]);
      secondVector = mxGetPr(prhs[2]);

      plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
      output = (double *)mxGetPr(plhs[0]);


      if ((numberOfFeatures == 1) && (checkFeatures == 1))
      {
        if ((numberOfSamples == 0) && (checkSamples == 0))
        {
          *output = 0.0;
        }
        else if (numberOfSamples == 0)
        {
          *output = discAndCalcEntropy(secondVector,numberOfSamples);
        }
        else if (checkSamples == 0)
        {
          *output = discAndCalcEntropy(firstVector,numberOfSamples);
        }
        else if (numberOfSamples == checkSamples)
        {
          /*double discAndCalcJointEntropy(double *firstVector, double *secondVector, int vectorLength);*/
          *output = discAndCalcJointEntropy(firstVector,secondVector,numberOfSamples);
        }
        else
        {
          printf("Vector lengths do not match, they must be the same length\n");
          *output = -1.0;
        }
      }
      else
      {
        printf("No columns in input\n");
        *output = -1.0;
      }
      
      break;
    }/*case 5 - H(XY)*/
    case 6:
    {
      /*
      **H(X|Y)
      */
      numberOfSamples = mxGetM(prhs[1]);
      checkSamples = mxGetM(prhs[2]);
      
      numberOfFeatures = mxGetN(prhs[1]);
      checkFeatures = mxGetN(prhs[2]);

      dataVector = mxGetPr(prhs[1]);
      condVector = mxGetPr(prhs[2]);

      plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
      output = (double *)mxGetPr(plhs[0]);

      if ((numberOfFeatures == 1) && (checkFeatures == 1))
      {
        if (numberOfSamples == 0)
        {
          *output = 0.0;
        }
        else if (checkSamples == 0)
        {
          *output = discAndCalcEntropy(dataVector,numberOfSamples);
        }
        else if (numberOfSamples == checkSamples)
        {
          /*double discAndCalcConditionalEntropy(double *dataVector, double *condVector, int vectorLength);*/
          *output = discAndCalcConditionalEntropy(dataVector,condVector,numberOfSamples);
        }
        else
        {
          printf("Vector lengths do not match, they must be the same length\n");
          *output = -1.0;
        }
      }
      else
      {
        printf("No columns in input\n");
        *output = -1.0;
      }
      break;
    }/*case 6 - H(X|Y)*/
    case 7:
    {
      /*
      **I(X;Y)
      */
      numberOfSamples = mxGetM(prhs[1]);
      checkSamples = mxGetM(prhs[2]);

      numberOfFeatures = mxGetN(prhs[1]);
      checkFeatures = mxGetN(prhs[2]);
      
      firstVector = mxGetPr(prhs[1]);
      secondVector = mxGetPr(prhs[2]);

      plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
      output = (double *)mxGetPr(plhs[0]);

      if ((numberOfFeatures == 1) && (checkFeatures == 1))
      {
        if ((numberOfSamples == 0) || (checkSamples == 0))
        {
          *output = 0.0;
        }
        else if (numberOfSamples == checkSamples)
        {
          /*double discAndCalcMutualInformation(double *firstVector, double *secondVector, int vectorLength);*/
          *output = discAndCalcMutualInformation(firstVector,secondVector,numberOfSamples);
        }
        else
        {
          printf("Vector lengths do not match, they must be the same length\n");
          *output = -1.0;
        }
      }
      else
      {
        printf("No columns in input\n");
        *output = -1.0;
      }
      break;
    }/*case 7 - I(X;Y)*/
    case 8:
    {
      /*
      **I(X;Y|Z)
      */
      numberOfSamples = mxGetM(prhs[1]);
      checkSamples = mxGetM(prhs[2]);
      thirdCheckSamples = mxGetM(prhs[3]);
      
      numberOfFeatures = mxGetN(prhs[1]);
      checkFeatures = mxGetN(prhs[2]);
      thirdCheckFeatures = mxGetN(prhs[3]);

      firstVector = mxGetPr(prhs[1]);
      targetVector = mxGetPr(prhs[2]);
      condVector = mxGetPr(prhs[3]);

      plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
      output = (double *)mxGetPr(plhs[0]);
      
      if ((numberOfFeatures == 1) && (checkFeatures == 1))
      {
        if ((numberOfSamples == 0) || (checkSamples == 0))
        {
          *output = 0.0;
        }
        else if ((thirdCheckSamples == 0) || (thirdCheckFeatures != 1))
        {
          *output = discAndCalcMutualInformation(firstVector,targetVector,numberOfSamples);
        }
        else if ((numberOfSamples == checkSamples) && (numberOfSamples == thirdCheckSamples))
        {
          /*double discAndCalcConditionalMutualInformation(double *firstVector, double *targetVector, double *condVector, int vectorLength);*/
          *output = discAndCalcConditionalMutualInformation(firstVector,targetVector,condVector,numberOfSamples);
        }
        else
        {
          printf("Vector lengths do not match, they must be the same length\n");
          *output = -1.0;
        }
      }
      else
      {
        printf("No columns in input\n");
        *output = -1.0;
      }
      break;
    }/*case 8 - I(X;Y|Z)*/
    case 9:
    {
      /*
      **normaliseArray
      */
      numberOfSamples = mxGetM(prhs[1]);
      dataVector = (double *) mxGetPr(prhs[1]);
      
      outputIntVector = (uint *) mxCalloc(numberOfSamples,sizeof(uint));

      plhs[0] = mxCreateDoubleMatrix(numberOfSamples,1,mxREAL);
      plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
      output = (double *)mxGetPr(plhs[0]);
      numStates = (double *) mxGetPr(plhs[1]);
      
      /*int normaliseArray(double *inputVector, int *outputVector, int vectorLength);*/
      *numStates = normaliseArray(dataVector, outputIntVector, numberOfSamples);
      
      for (i = 0; i < numberOfSamples; i++)
      {
        output[i] = outputIntVector[i];
      }

      mxFree(outputIntVector);
      outputIntVector = NULL;
      
      break;
    }/*case 9 - normaliseArray*/
    default:
    {
      printf("Unrecognised flag\n");
      break;
    }/*default*/
  }/*switch(flag)*/
  
  return;
}/*mexFunction()*/

