function output = mi(X,Y)
%function output = mi(X,Y)
%X & Y can be matrices which are converted into a joint variable
%before computation
%
%expects variables to be column-wise
%
%returns the mutual information between X and Y, I(X;Y)

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
[output] = MIToolboxMex(7,mergedFirst,mergedSecond);
