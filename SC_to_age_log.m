function [sc_corr,sc_res] = SC_to_age_log(SC,age)
% This program will convert the SC matrix from raw SC to log SC
% corrected for age matched control. The indices in the correction matrix
% can be found in the node_info.csv under connROI and can be used as is if
% the input matrix has order of ROIS set by ID
%
%  
% INPUTS - 
% SC = matrix of SC streamline count matrix from mrtrix3. ROIS need to match the order of ID in node_info.csv or file must be edited
% age = age of the patient in years
%
% The program will load in node_info.csv for indicies info and MA117_SC_conage_fit_85.mat
% with correction parameters derived from 85 healthy controls

load('MA117_SC_conage_fit.mat');
node_all=readtable('node_info.csv');
numROIS=length(node_all.ID);

% need age and length of time series for conversion
if nargin < 2
    age = input('Enter age of the patient in years: ');
end


sc_corr=NaN(numROIS,numROIS);
sc_res=NaN(numROIS,numROIS);

for ic=1:numROIS
    i=node_all.connROI(ic);
    for jc=1:numROIS
        j=node_all.connROI(jc);
        if i~=j 
            if con_age_fit(i,j,3)~=0 % if fit is not available then don't do, make corr NaN
                co=log(SC(ic,jc));
                sc_corr(ic,jc)=[co-(age.*con_age_fit(i,j,1)+con_age_fit(i,j,2))]./con_age_fit(i,j,3);
                sc_res(ic,jc)=[co-(age.*con_age_fit(i,j,1)+con_age_fit(i,j,2))];
            end
        end
    end
end

end


