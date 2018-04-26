function selectedFeatures = mRMR_D(k, featureMatrix, classColumn)
%function selectedFeatures = mRMR_D(k, featureMatrix, classColumn)
%
%Selects optimal features according to the mRMR-D algorithm from 
%"Feature Selection Based on Mutual Information: Criteria of Max-Dependency, Max-Relevance, and Min-Redundancy"
%by H. Peng et al. (2005)
%
%Calculates the top k features
%a dataset featureMatrix with n training examples and m features
%with the classes held in classColumn (an n x 1 vector)

noOfTraining = size(classColumn,1);
noOfFeatures = size(featureMatrix,2);
unselectedFeatures = ones(noOfFeatures,1);

classMI = zeros(noOfFeatures,1);
answerFeatures = zeros(k,1);
highestMI = 0;
highestMICounter = 0;
currentHighestFeature = 0;

featureMIMatrix = -(ones(k,noOfFeatures));

%setup the mi against the class
for n = 1 : noOfFeatures
  classMI(n) = mi(featureMatrix(:,n),classColumn);
	if classMI(n) > highestMI
		highestMI = classMI(n);
		highestMICounter = n;
	end
end

answerFeatures(1) = highestMICounter;
unselectedFeatures(highestMICounter) = 0;

%iterate over the number of features to select
for i = 2:k
  score = -100;
	currentHighestFeature = 0;
	iMinus = i-1;
  for j = 1 : noOfFeatures
		if unselectedFeatures(j) == 1
			currentMIScore = 0;
			for m = 1 : iMinus
        if featureMIMatrix(m,j) == -1
          featureMIMatrix(m,j) = mi(featureMatrix(:,j),featureMatrix(:,answerFeatures(m)));
        end
				currentMIScore =  currentMIScore + featureMIMatrix(m,j);
			end
			currentScore = classMI(j) - (currentMIScore/iMinus);
			
			if (currentScore > score)
				score = currentScore;
				currentHighestFeature = j;
			end
		end
	end
	
  if score < 0 
    disp(['at selection ' int2str(j) ' mRMRD is negative with value ' num2str(score)]);
  end

	%now highest feature is selected in currentHighestFeature
	%store it
	unselectedFeatures(currentHighestFeature) = 0;
	answerFeatures(i) = currentHighestFeature;
end

selectedFeatures = answerFeatures;
