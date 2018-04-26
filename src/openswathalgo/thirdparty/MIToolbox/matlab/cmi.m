function output = cmi(X,Y,Z)
%function output = cmi(X,Y,Z)
%X, Y & Z can be matrices which are converted into a joint variable
%before computation
%
%expects variables to be column-wise
%
%returns the mutual information between X and Y conditioned on Z, I(X;Y|Z)

if nargin == 3
  if (~isa(X,'double') || ~isa(Y,'double') || ~isa(Z,'double'))
    error('Error, inputs must be double vectors or matrices')
  end
  if (size(X,2)>1)
	mergedFirst = MIToolboxMex(3,X);
  else
	mergedFirst = X;
  end
  if (size(Y,2)>1)
	mergedSecond = MIToolboxMex(3,Y);
  else
	mergedSecond = Y;
  end
  if (size(Z,2)>1)
	mergedThird = MIToolboxMex(3,Z);
  else
	mergedThird = Z;
  end
  [output] = MIToolboxMex(8,mergedFirst,mergedSecond,mergedThird);
elseif nargin == 2
  if (~isa(X,'double') || ~isa(Y,'double'))
    error('Error, inputs must be double vectors or matrices')
  end
  output = mi(X,Y);
else
  output = 0;
end
