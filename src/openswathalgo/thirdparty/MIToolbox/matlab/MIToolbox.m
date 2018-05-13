function [varargout] = MIToolbox(functionName, varargin)
%function [varargout] = MIToolbox(functionName, varargin)
%
%Provides access to the functions in MIToolboxMex
%
%Expects column vectors, will not work with row vectors
%
%Function list
%"joint" = joint variable of the matrix
%"entropy" or "h" = H(X)
%"ConditionalEntropy" or "condh" = H(X|Y)
%"mi" = I(X;Y)
%"ConditionalMI" or "cmi" = I(X;Y|Z)
%
%Arguments and returned values
%[jointVariable] = joint(matrix)
%[entropy] = H(X) = H(vector)
%[entropy] = H(X|Y) = H(vector,condition)
%[mi] = I(X;Y) = I(vector,target)
%[mi] = I(X;Y|Z) = I(vector,target,condition)
%
%Internal MIToolbox function number
%Joint                = 3
%Entropy              = 4
%Conditional Entropy  = 6
%Mutual Information   = 7
%Conditional MI       = 8

for i = 1:length(varargin)
    if (~isa(varargin{i},'double'))
        error('Error, MIToolbox requires inputs to be double vector or matrices')
    end
end

if (strcmpi(functionName,'Joint') || strcmpi(functionName,'Merge'))
    [varargout{1}] = MIToolboxMex(3,varargin{1});
elseif (strcmpi(functionName,'Entropy') || strcmpi(functionName,'h'))
    %disp('Calculating Entropy');
    if (size(varargin{1},2)>1)
			mergedVector = MIToolboxMex(3,varargin{1});
		else
			mergedVector = varargin{1};
    end
    [varargout{1}] = MIToolboxMex(4,mergedVector);
elseif ((strcmpi(functionName,'ConditionalEntropy')) || strcmpi(functionName,'condh'))
    if (size(varargin{1},2)>1)
			mergedFirst = MIToolboxMex(3,varargin{1});
		else
			mergedFirst = varargin{1};
    end
    if (size(varargin{2},2)>1)
			mergedSecond = MIToolboxMex(3,varargin{2});
		else
			mergedSecond = varargin{2};
    end
    [varargout{1}] = MIToolboxMex(6,mergedFirst,mergedSecond);
elseif (strcmpi(functionName,'mi'))
    if (size(varargin{1},2)>1)
			mergedFirst = MIToolboxMex(3,varargin{1});
		else
			mergedFirst = varargin{1};
    end
    if (size(varargin{2},2)>1)
			mergedSecond = MIToolboxMex(3,varargin{2});
		else
			mergedSecond = varargin{2};
    end
    [varargout{1}] = MIToolboxMex(7,mergedFirst,mergedSecond);
elseif (strcmpi(functionName,'ConditionalMI') || strcmpi(functionName,'cmi'))
    if (size(varargin{1},2)>1)
			mergedFirst = MIToolboxMex(3,varargin{1});
		else
			mergedFirst = varargin{1};
    end
    if (size(varargin{2},2)>1)
			mergedSecond = MIToolboxMex(3,varargin{2});
		else
			mergedSecond = varargin{2};
    end
    if (size(varargin{3},2)>1)
			mergedThird = MIToolboxMex(3,varargin{3});
		else
			mergedThird = varargin{3};
    end
    [varargout{1}] = MIToolboxMex(8,mergedFirst,mergedSecond,mergedThird);
else
    varargout{1} = 0;
    disp(['Unrecognised functionName ' functionName]);
end
