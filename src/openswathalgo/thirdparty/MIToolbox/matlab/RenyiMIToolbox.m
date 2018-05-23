function [varargout] = RenyiMIToolbox(functionName, alpha, varargin)
%function [varargout] = RenyiMIToolbox(functionName, alpha, varargin)
%
%Provides access to the functions in RenyiMIToolboxMex
%
%Expects column vectors, will not work with row vectors
%
%Function list
%"Entropy" = H_{\alpha}(X) = 1
%"MI" = I_{\alpha}(X;Y) = 3
%
%Arguments and returned values
%[entropy] = H_\alpha(X) = H(alpha,vector)
%[mi] = I_\alpha(X;Y) = I(alpha,vector,target)
%
%Internal RenyiMIToolbox function number
%Renyi Entropy = 1;
%Renyi MI = 3;

for i = 1:length(varargin)
    if (~isa(varargin{i},'double'))
        error('Error, MIToolbox requires inputs to be double vector or matrices')
    end
end

if (alpha ~= 1)
  if (strcmpi(functionName,'Entropy') || strcmpi(functionName,'h'))
      %disp('Calculating Entropy');
      if (size(varargin{1},2)>1)
  			mergedVector = MIToolboxMex(3,varargin{1});
  		else
  			mergedVector = varargin{1};
      end
      [varargout{1}] = RenyiMIToolboxMex(1,alpha,mergedVector);
  elseif (strcmpi(functionName,'MI'))
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
      [varargout{1}] = RenyiMIToolboxMex(3,alpha,mergedFirst,mergedSecond);
  else
      varargout{1} = 0;
      disp(['Unrecognised functionName ' functionName]);
  end
else
  disp('For alpha = 1 use functions in MIToolbox.m');
  disp('as those functions are the implementation of Shannon''s Information Theory');
end
