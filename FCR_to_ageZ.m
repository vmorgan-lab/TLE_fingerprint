function [fc_corr,fc_res] = FCR_to_ageZ(FCr,age, ts_length)
% This program will convert the FC matrix from r values to z values
% corrected for age matched control. The indices in the correction matrix
% can be found in the node_info.csv under connROI and can be used as is if
% the input matrix has order of ROIS set by ID
%
%  
% INPUTS - 
% FCr = matrix of FC r correlation values. ROIS need to match the order of ID in node_info.csv or file must be edited
% age = age of the patient in years
% ts_length = time series length
%
% The program will load in node_info.csv for indicies info and MA117_FC_conage_fit_85.mat
% with correction parameters derived from 85 healthy controls

load('MA117_FC_conage_fit.mat');
node_all=readtable('node_info.csv');
numROIS=length(node_all.ID);

% need age and length of time series for conversion
if nargin < 3
    ts_length = input('Enter length of time series for Z conversion: ');
    if nargin < 2
        age = input('Enter age of the patient in years: ');
    end
end


% load in connectivity matrix in r values and convert to Z values
pat_Z=(0.5.* log((1+FCr(:,:))./(1-FCr(:,:))))*sqrt(ts_length-3); % convert to Z
fc_corr=NaN(numROIS,numROIS);
fc_res=NaN(numROIS,numROIS);

for ic=1:numROIS
    i=node_all.connROI(ic);
    for jc=1:numROIS
        j=node_all.connROI(jc);
        if i~=j 
            if con_age_fit(i,j,3)~=0 % if fit is not available then don't do, make corr NaN
                co=pat_Z(ic,jc);
                fc_corr(ic,jc)=[co-(age.*con_age_fit(i,j,1)+con_age_fit(i,j,2))]./con_age_fit(i,j,3);
                fc_res(ic,jc)=[co-(age.*con_age_fit(i,j,1)+con_age_fit(i,j,2))];
            end
        end
    end
end

end

