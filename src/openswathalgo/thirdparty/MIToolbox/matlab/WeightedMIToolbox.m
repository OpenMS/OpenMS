function [varargout] = WeightedMIToolbox(functionName, weightVector, varargin)
%function [varargout] = WeightedMIToolbox(functionName, weightVector, varargin)
%
%Provides access to the functions in WeightedMIToolboxMex
%
%Expects column vectors, will not work with row vectors
%
%Function list
%"entropy" or "h" = H(X)
%"ConditionalEntropy" or "condh" = H(X|Y)
%"mi" = I(X;Y)
%"ConditionalMI" or "cmi" = I(X;Y|Z)
%
%Arguments and returned values
%[entropy] = H_w(X) = H_w(vector)
%[entropy] = H_w(X|Y) = H_w(vector,condition)
%[mi] = I_w(X;Y) = I_w(vector,target)
%[mi] = I_w(X;Y|Z) = I_w(vector,target,condition)
%
%Internal MIToolbox function number
%Entropy              = 1
%Conditional Entropy  = 3
%Mutual Information   = 4
%Conditional MI       = 5

for i = 1:length(varargin)
    if (~isa(varargin{i},'double'))
        error('Error, MIToolbox requires inputs to be double vector or matrices')
    end
end

if (strcmpi(functionName,'Entropy') || strcmpi(functionName,'h'))
    %disp('Calculating Entropy');
    if (size(varargin{1},2)>1)
			mergedVector = MIToolboxMex(3,varargin{1});
		else
			mergedVector = varargin{1};
    end
    [varargout{1}] = WeightedMIToolboxMex(1,weightVector,mergedVector);
elseif ((strcmpi(functionName,'JointEntropy')) || strcmpi(functionName,'jointh'))
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
    [varargout{1}] = WeightedMIToolboxMex(2,weightVector,mergedFirst,mergedSecond);
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
    [varargout{1}] = WeightedMIToolboxMex(3,weightVector,mergedFirst,mergedSecond);
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
    [varargout{1}] = WeightedMIToolboxMex(4,weightVector,mergedFirst,mergedSecond);
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
    [varargout{1}] = WeightedMIToolboxMex(5,weightVector,mergedFirst,mergedSecond,mergedThird);
else
    varargout{1} = 0;
    disp(['Unrecognised functionName ' functionName]);
end
