function selectedFeatures = DISR(k, featureMatrix, classColumn)
%function selectedFeatures = DISR(k, featureMatrix, classColumn)
%
%Computers optimal features according to DISR algorithm from 
%On the Use of variable "complementarity for feature selection" 
%by P Meyer, G Bontempi (2006)
%
%Computes the top k features from
%a dataset featureMatrix with n training examples and m features
%with the classes held in classColumn.
%
%DISR - arg(Xi) max(sum(Xj mem XS)(SimRel(Xij,Y)))
%where SimRel = MI(Xij,Y) / H(Xij,Y)

totalFeatures = size(featureMatrix,2);
classMI = zeros(totalFeatures,1);
unselectedFeatures = ones(totalFeatures,1);
score = 0;
currentScore = 0;
innerScore = 0;
iMinus = 0;
answerFeatures = zeros(k,1);
highestMI = 0;
highestMICounter = 0;
currentHighestFeature = 0;

%create a matrix to hold the SRs of a feature pair. 
%initialised to -1 as you can't get a negative SR.
featureSRMatrix = -(ones(k,totalFeatures));

for n = 1 : totalFeatures
	classMI(n) = mi(featureMatrix(:,n),classColumn);
	if classMI(n) > highestMI
		highestMI = classMI(n);
		highestMICounter = n;
	end
end

answerFeatures(1) = highestMICounter;
unselectedFeatures(highestMICounter) = 0;

for i = 2 : k
	score = 0;
	currentHighestFeature = 0;
	iMinus = i-1;
  for j = 1 : totalFeatures
		if unselectedFeatures(j) == 1
			%DISR - arg(Xi) max(sum(Xj mem XS)(SimRel(Xij,Y)))
			%where SimRel = MI(Xij,Y) / H(Xij,Y)
			currentScore = 0;
			for m = 1 : iMinus
        if featureSRMatrix(m,j) == -1
          unionedFeatures = joint([featureMatrix(:,answerFeatures(m)),featureMatrix(:,j)]);
				  tempUnionMI = mi(unionedFeatures,classColumn);
				  tempTripEntropy = h([unionedFeatures,classColumn]);
				  featureSRMatrix(m,j) = tempUnionMI/tempTripEntropy;
        end
        
				currentScore =  currentScore + featureSRMatrix(m,j);
			end
			if (currentScore > score)
				score = currentScore;
				currentHighestFeature = j;
			end
		end
	end
	%now highest feature is selected in currentHighestFeature
	%store it
	unselectedFeatures(currentHighestFeature) = 0;
	answerFeatures(i) = currentHighestFeature;
end

selectedFeatures = answerFeatures;
