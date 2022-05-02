function [MEFC_pat, MESC_pat, MEFCSC_pat, pat_profile] = TLE_fingerprint_distance(model,pat_data, side_1)
% This function will plot patient profiles and the TLE fingerprint on a
% radar plot as shown in Morgan et al. Brain Communications 2022 Figure 5
% subplots
%
% INPUTs - 
% model = matrix of the TLE fingerprint in 40 x 6 x 2
%   40 are 40 permutations of seizure free patients
%   6 are lobes = prefrontal, parietal, occipital, temporal,
%   motor/somatosensory and subcortical
%   2 are FC and SC
% pat_data = matrix of the invidual patient data in n xn n x 2 = nodes x
% nodes x FC and SC
%
% side_1 = side of focus - left = 0, right = 1
%
% this program loads in the excel spreadsheet - node_info.xlsx which
% contains info about the nodes, their order, their side and their lobe.
% Interest column should be one for the nodes of interest.
%
% pat_data(:,:,1)=all_subj_FC(36,:,:);
% pat_data(:,:,2)=all_subj_SC(36,:,:);
% model(:,:,1)=w_perm_FC_lobe(:,1:6,1);
% model(:,:,2)=w_perm_FC_lobe(:,1:6,2);
%
% OUTPUTS - 
% MEFC = functional connectivity distance to model profile
% MESC = structural connectivity distance to model profile
% total = total distance to model profile
% pat_profile = 6 lobes x 2 (FC and SC) that can be plotted with plot_profile.m 

% created by VL Morgan April 2022

% check to see if side of focus is entered
if nargin < 3
    side_1 = input('Enter side of focus (0 = left, 1 = right): ');
    while (side_1 ~= 0) & (side_1 ~= 1) 
        side_1 = input('Enter side of focus (0 = left, 1 = right): ');
    end
end
 

% load in node data - edit as needed - to compute the weights matrices
% for now nodes of interest for FC and SC are the same
node_all=dataset('XLSFile','node_info.xlsx','Sheet','BC_2022');
num_rois=length(node_all.ID);

side = node_all.side_right;
if side_1 == 0 % swap if left sided
    side = node_all.side_left;
end

node_intFC=find(node_all.interest==1); %regular 14
node_intSC=node_intFC;

all_nodesFC=zeros(num_rois,1);
all_nodesFC(node_intFC)=1;
skip_nodesFC=find(all_nodesFC==0);

all_nodesSC=zeros(num_rois,1);
all_nodesSC(node_intSC)=1;
skip_nodesSC=find(all_nodesSC==0);

FCweight_mat=ones(num_rois, num_rois);
SCweight_mat=ones(num_rois, num_rois);

ip_ip=2; sfc(1,1)=ip_ip;
ip_con=1; sfc(1,2)=ip_con; sfc(2,1)=ip_con;
con_con=1; sfc(2,2)=con_con;

ip_ip2=2; ssc(1,1)=ip_ip2;
ip_con2=1; ssc(1,2)=ip_con2; ssc(2,1)=ip_con2;
con_con2=1; ssc(2,2)=con_con2;

for i=1:num_rois
    for j=1:num_rois
        FCweight_mat(i,j)=FCweight_mat(i,j)*sfc(side(i),side(j));
        SCweight_mat(i,j)=SCweight_mat(i,j)*ssc(side(i),side(j));
    end
end

for i=1:num_rois
    FCweight_mat(i,skip_nodesFC)=0;
    SCweight_mat(i,skip_nodesSC)=0;
end


% model data is input as 40 "subjects" x 6 lobes X 2 types of connectivity
% (FC then SC)

% set parameters for Malahanobis distance to model
clear muFC muSC CFC iCFC CSC iSCS
muFC=nanmean(model(:,:,1),1)';
muSC=nanmean(model(:,:,2),1)';
[CFC,rho] = shrinkage_cov(model(:,:,1),'rblw');
iCFC=inv(CFC); 
[CSC,rho] = shrinkage_cov(model(:,:,2),'rblw');
iCSC=inv(CSC); 


% load in patient data as full patient matrix - all nodes to all nodes x 2
% types of connectivity (FC then SC)
clear pat_FC pat_SC
pat_FC(:,:)=pat_data(:,:,1);
pat_SC(:,:)=pat_data(:,:,2);

% multiply by weighting - 2 x nodes of interest ipsi and 1 x nodes of
% interest contra
clear w_pat_FC w_pat_SC
w_pat_FC(:,:)=squeeze(pat_FC(:,:)).*(FCweight_mat);
w_pat_SC(:,:)=squeeze(pat_SC(:,:)).*(SCweight_mat);

% compute nodes degree -
clear w_patFC_node w_pat_SC_node 
w_pat_FC_node=squeeze(sum((w_pat_FC(:,:)),2,'omitnan'));
w_pat_SC_node=squeeze(sum((w_pat_SC(:,:)),2,'omitnan')); 
 
% compute pat average per ipsilateral lobe
clear w_pat_data_FC_lobe w_pat_data_SC_lobe 
for lo=1:6
        ind_side_lobe=find(side==1 & node_all.lobe==lo);
        w_pat_FC_lobe(lo)=nanmean(w_pat_FC_node(ind_side_lobe));
        w_pat_SC_lobe(lo)=nanmean(w_pat_SC_node(ind_side_lobe));
end

% compute pat distances to model
clear MFC_pat MSC_pat EFC_pat ESC_pat MEFC_pat MESC_pat MEFCSC_pat 
dd = 6; % number of lobes used
MFC_pat=sqrt((w_pat_FC_lobe'-muFC)'*iCFC*(w_pat_FC_lobe'-muFC));  % Mahalnobis functional
MSC_pat=sqrt((w_pat_SC_lobe'-muSC)'*iCSC*(w_pat_SC_lobe'-muSC));  % Mahalnobis structural
EFC_pat=pdist2(w_pat_FC_lobe,muFC','euclidean');                                              % Euclidean functional 
ESC_pat=pdist2(w_pat_SC_lobe,muSC','euclidean');                                              % Euclidean structural
MEFC_pat=MFC_pat+EFC_pat./dd;                                                                 % total functional  
MESC_pat=MSC_pat+ESC_pat./dd;                                                                 % total structural
MEFCSC_pat=sqrt((MEFC_pat.^2)+(MESC_pat.^2));                                                 % total functional and structural distance  

% save profile = 6 lobes x 2 types of conn
% lobes = pref par occ temp mot/som sub
% connectivity = FC then SC

clear pat_profile
pat_profile=[w_pat_FC_lobe; w_pat_SC_lobe]';
end

