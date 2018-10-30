function [cmb association] = IAMB( data, targetindex, THRESHOLD)
%function [cmb association] = IAMB( data, targetindex, THRESHOLD)
%
%Performs the IAMB algorithm of Tsmardinos et al. (2003)
%from "Towards principled feature selection: Relevancy, filters and wrappers"

if (nargin == 2)
  THRESHOLD = 0.02;
end

numf = size(data,2);
targets = data(:,targetindex);
data(:,targetindex) = -10;

cmb = [];

finished = false;
while ~finished
    for n = 1:numf
        cmbVector = joint(data(:,cmb));
        if isempty(cmb)
            association(n) = mi( data(:,n), targets );
        end
        
        if ismember(n,cmb)
            association(n) = -10; %arbtirary large negative constant
        else
            association(n) = cmi( data(:,n), targets, cmbVector);
        end
    end

    [maxval maxidx] = max(association);
    if maxval < THRESHOLD
        finished = true;
    else
        cmb = [ cmb maxidx ];
    end
end

finished = false;
while ~finished && ~isempty(cmb)
    association = [];
    for n = 1:length(cmb)
        cmbwithoutn = cmb;
        cmbwithoutn(n)=[];
        association(n) = cmi( data(:,cmb(n)), targets, data(:,cmbwithoutn) );
    end
    
    [minval minidx] = min(association);
    if minval > THRESHOLD
        finished = true;
    else        
        cmb(minidx) = [];
    end
end

