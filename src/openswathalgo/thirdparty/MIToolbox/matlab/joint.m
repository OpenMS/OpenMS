function output = joint(X,arities)
%function output = joint(X,arities)
%returns the joint random variable of the matrix X
%assuming the variables are in columns
%
%if passed a vector of the arities then it produces a correct
%joint variable, otherwise it may not include all states
%
%if the joint variable is only compared with variables using the same samples,
%then arity information is not required

if (nargin == 2)
  if (~isa(X,'double') || ~isa(arities,'double'))
    error('Error, inputs must be double vectors or matrices')
  end
  [output] = MIToolboxMex(3,X,arities);
else
  if (~isa(X,'double'))
    error('Error, input must be a double vector or matrix')
  end
  [output] = MIToolboxMex(3,X);
end
