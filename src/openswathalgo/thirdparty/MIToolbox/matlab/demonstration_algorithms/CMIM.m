function selectedFeatures = CMIM(k, featureMatrix, classColumn)
%function selectedFeatures = CMIM(k, featureMatrix, classColumn)
%Computes conditional mutual information maximisation algorithm from
%"Fast Binary Feature Selection with Conditional Mutual Information"
%by F. Fleuret (2004)

%Computes the top k features from
%a dataset featureMatrix with n training examples and m features
%with the classes held in classColumn.

noOfTraining = size(classColumn,1);
noOfFeatures = size(featureMatrix,2);

partialScore = zeros(noOfFeatures,1);
m = zeros(noOfFeatures,1);
score = 0;
answerFeatures = zeros(k,1);
highestMI = 0;
highestMICounter = 0;

for n = 1 : noOfFeatures
  partialScore(n) = mi(featureMatrix(:,n),classColumn);
  if partialScore(n) > highestMI
    highestMI = partialScore(n);
    highestMICounter = n;
  end
end

answerFeatures(1) = highestMICounter;

for i = 2 : k
  score = 0;
  limitI = i - 1;
  for n = 1 : noOfFeatures
    while ((partialScore(n) >= score) && (m(n) < limitI))
      m(n) = m(n) + 1;
      conditionalInfo = cmi(featureMatrix(:,n),classColumn,featureMatrix(:,answerFeatures(m(n))));
      if partialScore(n) > conditionalInfo
        partialScore(n) = conditionalInfo;
      end
    end
    if partialScore(n) >= score
      score = partialScore(n);
      answerFeatures(i) = n;
    end
  end
end

selectedFeatures = answerFeatures;
