function output = condh(X,Y)
%function output = condh(X,Y)
%X & Y can be matrices which are converted into a joint variable
%before computation
%
%expects variables to be column-wise
%
%returns the conditional entropy of X given Y, H(X|Y)

if nargin == 2
  if (~isa(X,'double') || ~isa(Y,'double'))
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
  [output] = MIToolboxMex(6,mergedFirst,mergedSecond);
elseif nargin == 1
  if (~isa(X,'double'))
    error('Error, inputs must be double vectors or matrices')
  end
  output = h(X);
else
  output = 0;
end
