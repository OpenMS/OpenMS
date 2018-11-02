#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "MIToolbox/ArrayOperations.h"
#include "MIToolbox/Entropy.h"
#include "MIToolbox/MutualInformation.h"

int main(int argc, char *argv[])
{
  int i;
  double length, miTarget, entropyTarget, cmiTarget;
  double firstEntropy, secondEntropy, thirdEntropy, targetEntropy;
  double firstMItarget, secondMItarget, thirdMItarget, targetMItarget;
  int *testFirstVector, *testSecondVector, *testThirdVector, *testMergedVector;
  struct timeval start,end;

  int *firstVector = (int *) calloc(4,sizeof(int));
  int *secondVector = (int *) calloc(4,sizeof(int));
  int *thirdVector = (int *) calloc(4,sizeof(int));
  int *targetVector = (int *) calloc(4,sizeof(int));
  
  firstVector[0] = 0;
  firstVector[1] = 0;
  firstVector[2] = 1;
  firstVector[3] = 1;
  
  secondVector[0] = 0;
  secondVector[1] = 1;
  secondVector[2] = 0;
  secondVector[3] = 1;
  
  thirdVector[0] = 0;
  thirdVector[1] = 1;
  thirdVector[2] = 1;
  thirdVector[3] = 1;
  
  targetVector[0] = 0;
  targetVector[1] = 1;
  targetVector[2] = 1;
  targetVector[3] = 0;
  
  firstEntropy = calcEntropy(firstVector,4);
  secondEntropy = calcEntropy(secondVector,4);
  thirdEntropy = calcEntropy(thirdVector,4);
  targetEntropy = calcEntropy(targetVector,4);
  
  printf("Entropies - first: %f, second: %f, third: %f, target %f\n",firstEntropy,secondEntropy,thirdEntropy,targetEntropy);
    
  firstMItarget = calcMutualInformation(firstVector,targetVector,4);
  secondMItarget = calcMutualInformation(secondVector,targetVector,4);
  thirdMItarget = calcMutualInformation(thirdVector,targetVector,4);
  targetMItarget = calcMutualInformation(targetVector,targetVector,4);
  
  printf("MIs - first: %f, second: %f, third: %f, target %f\n",firstMItarget,secondMItarget,thirdMItarget,targetMItarget);
  
  testFirstVector = (int *) calloc(10000,sizeof(int));
  testSecondVector = (int *) calloc(10000,sizeof(int));
  testThirdVector = (int *) calloc(10000,sizeof(int));
  testMergedVector = (int *) calloc(10000,sizeof(int));
  
  for (i = 0; i < 10000; i++)
  {
    testFirstVector[i] = i % 2;
    testSecondVector[i] = i % 4;
    testThirdVector[i] = i % 3;
  }
  /* struct timeval
   * {
   *    time_t         tv_sec      seconds
   *    suseconds_t    tv_usec     microseconds
   * }
   */

  gettimeofday(&start, NULL);
  for (i = 0; i < 1000; i++)
  {
    miTarget = calcMutualInformation(testFirstVector,testSecondVector,10000);
    entropyTarget = calcEntropy(testFirstVector,10000);
    cmiTarget = calcConditionalMutualInformation(testFirstVector,testSecondVector,testThirdVector,10000);
    mergeArrays(testFirstVector,testSecondVector,testMergedVector,10000);
  }
  gettimeofday(&end, NULL);
  printf("I(X;Y) = %f, H(X) = %f, I(X;Y|Z) = %f\n",miTarget,entropyTarget,cmiTarget);
  
  length = end.tv_sec - start.tv_sec;
  length = length + (end.tv_usec - start.tv_usec) / 1000000.0;
  
  printf("Time taken for a thousand I(X;Y), H(X), I(X;Y|Z), merge(X,Y) is %lf seconds\n",length);
}/*main(int, char **)*/
