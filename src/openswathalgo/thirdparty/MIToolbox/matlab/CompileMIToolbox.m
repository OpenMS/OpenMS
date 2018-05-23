% Compiles the MIToolbox functions

mex -I../include MIToolboxMex.c ../src/MutualInformation.c ../src/Entropy.c ../src/CalculateProbability.c ../src/ArrayOperations.c
mex -I../include RenyiMIToolboxMex.c ../src/RenyiMutualInformation.c ../src/RenyiEntropy.c ../src/CalculateProbability.c ../src/ArrayOperations.c
mex -I../include WeightedMIToolboxMex.c ../src/WeightedMutualInformation.c ../src/WeightedEntropy.c ../src/CalculateProbability.c ../src/ArrayOperations.c
